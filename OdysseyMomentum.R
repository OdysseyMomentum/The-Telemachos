##### ODYSSEY MOMENTUM #####
#####  THE TELEMACHOS  #####

# We will model the deterioration of a patient by formulating it as a Markov chain.

# We have the following possible actions:
# - No additional attention needed (DN) 
# - For now, the upcoming actions will be indicated by performing the action "inspection" (IN), 
#    during the pilot we will differentiate between these actions.
#     - Patients can fill in a questionnaire (A1) with corresponding cost c1
#     - Patients can be seen by a "wijkverpleegkundige" (A2) with corresponding cost c2
#     - Patients needs to be seen by the GP (A3) with corresponding cost c3
#     - Patients needs to be seen by a specialist (A4) with corresponding cost c4
# - Patients need to go the emergency departement (ED) with corresponding cost c5 

# We do not know the actual health status of the patient with certainty 
#  until the patient fills in the questionnaire or is seen by a health professional
# Due to this, the best approach is to model the problem as a partially observable Markov decision process (POMDP).

# We let theta denote the knowledge state. So theta_i,j denote the knowledge state that represents that the 
#  last observed health status is i, and that it is j time periods ago that this health status has been observed.

# Input
S <- c("1", "2", "3", "4", "5")       # Five possible health states
m <- 4                                # The number of health stats before the ED
A <- c("DN", "IN", "ED")              # The three possible actions

# Cost parameters / Negative rewards
c_inspect <- 20
c_GP      <- 50
c_ED      <- 10000
N  <- 50   # Maximum number of time periods between actions. (This is needed as the state space has to be bounded) 

# For now we will assume patients with "Disease A" deteriorate according to a gamma process with parameters a and b.
#  This will be adjusted when we get insight into the data. 
a <- 0.02
b <- 0.05
L <- 1     # We normalize the deterioration level for which you enter state 5 to one.
dt <- 1    # Each time period in our model consists of one day.

# To find the thresholds corresponding to the other health statuses we discretize the gamma deterioration process. 
#  The discrete-time Markov chain is now described by the transition probability matrix P. 
TPM_gammaprocess <- function(a, b, L, m, dt){
  P <- matrix(0, nrow=m+1, ncol = m+1)
  dx <-  L/m
  
  prob <- function(i, dx, dt, a, b){
    (1/dx)*integrate(pgamma, lower = (i-1)*dx, upper =i*dx, shape = a*dt, scale = b)$value
  }
  
  prob = Vectorize(prob, "i")
  
  q = prob(1:m, dx, dt, a, b)
  p = diff(c(0,q))
  
  for(i in 1:m){
    P[i,i:m] = p[1:(m+1-i)]
    P[i, m+1] = 1-sum(P[i,1:m])
  }
  P[m+1, m+1] = 1
  
  P[P<0] = 0
  return(P)
}

P <- TPM_gammaprocess(a, b, L, m, dt)
# Note that P is upper-diagonal and thus we assume patients cannot improve their health status without interference of a health professional.
#  If data shows this assumption is wrong, this can easily be changed. 

# We have to take into account that new patients join our system (and automatically enter class 1) and people leave the system (due to mortality). 
#  We make the assumption that the amount of patients joining and leaving the system is equal, such that our system size stays the same. 
#  Furthermore, for now, we assume that the patients that leave the system are uniformly distributed among the states and the total percentage of patients leaving is 1e-7. 
P.test[2,2:5] <- P[2,2:5] - 1e-7/(3*4)
P.test[2,1]   <- 1e-7/3 
P.test[3,3:5] <- P[3,3:5] - 1e-7/(3*3)
P.test[3,1]   <- 1e-7/3 
P.test[4,4:5] <- P[4,4:5] - 1e-7/(3*2)
P.test[4,1]   <- 1e-7/3 

P <- P.test

# Creating the knowledge states theta
theta <- array(NA, dim = c(m, N+1, length(S)))
q     <- array(NA, dim = c(m, N+1))  # Probability of going to the ED when choosing DN in state theta_{i,j}

for(i in 1:m){
  theta[i,1,] <- c(rep(0,i-1), 1, rep(0,m-i+1)) # Calculating the knowledge state of the patients that had their last inspection today. 
                                                #  We know their health status with certainty.
  for(j in 1:N){
    theta[i,j+1,] <- c(theta[i,j,]%*%P[,1:m]/sum(theta[i,j,]%*%P[,1:m]),0) # Calculating the knowledge states for all patients that did not have 
                                                                           #  an inspection today. Hence, answer the question: 
                                                                           #  what is likely to be their current health status?
    
    q[i,j]        <- 1 - sum(theta[i,j,]%*%P[,1:m])                        # Probability of going to the ED when choosing DN in state theta_{i,j}
  }
  q[i, N+1] <- 1 - sum(theta[i,N+1,]%*%P[,1:m])
}

# Now we move to definining the different transaction probabilities per action. 
P.new <- array(0, dim = c(length(A), m*(N+1), m*(N+1)))

# Transition matrix when choosing DN
# Note that we have got m*(N+1) possible states
# We use the following order of the states: c(theta^{1,0}, theta^{1,1}, ..., theta^{2,0}, theta^{2,1}, ...)
# DN is possible in all states except at theta^{i,N}
# You go back to theta^{1,0} if failure occurs, otherwise you continue to the next knowledge state
for(i in 1:(m*(N+1))){
  P.new[1,i,3] <- as.vector(t(q))[i]          # For now, we assume a patient coming out of the ED will be 
                                              #  a class 3 patient the rest of the time period until the next inspection.
  for(j in 2:(m*(N+1))){
    if(j==i+1)
      P.new[1,i,j] <- 1 - as.vector(t(q))[i] # If you do not inspect the patient and he did not come to the ED themselves, 
                                             #  the next knowledge state is theta^{i, j + 1}
  }
  if((i%%(N+1) == 0) == 1){
    P.new[1,i,] <- 0
  }
}

# Transition matrix when choosing IN
# A1 is possible in all states except theta^{i,0}
# You move to one of the theta^{i,0} as after IN you know the precise health status of the patient 
# To simplify the calculations we first put our knowledge states in a matrix
theta.matrix <- c()
for(i in 1:m){
  for(j in 1:(N+1)){
    theta.matrix <- rbind(theta.matrix,theta[i,j,]) 
  }
}

for(j in seq(1,m*(N+1), by = N+1)){
  P.new[2,,j]  <- theta.matrix[,which(seq(1,m*(N+1), by = N+1) == j)]
  P.new[2,j,]  <- 0
}

# Furthermore, there is a probability that this action improves the current health state of the patient.
#  This will be incoorporated during the pilot.

# Transition matrix when choosing ED
# The new health status is random:
P.new[3,,2] <- 0.2
P.new[3,,3] <- 0.6
P.new[3,,4] <- 0.2

# The Value Iteration Algorithm
epsilon <- 0.001  
v0 <- rep(0, m*(N+1)) # Initial value vector

val.iter.algo <- function(epsilon, A, v0, c_inspect, c_GP, c_ED, P.new){
  n <- 0
  v <- matrix(v0, ncol = m*(N+1), nrow = 1)
  
  span <- 2*epsilon
  while(span>epsilon){
    n = n+1
    v_all = matrix(nrow=m*(N+1), ncol=length(A))
    # Cost-to-go
    v_all[,1] = P.new[1,,] %*% v[n,] + as.vector(t(q*c_ED))
    v_all[,2] = P.new[2,,] %*% v[n,] + c_inspect
    v_all[,3] = P.new[3,,] %*% v[n,] + c_GP
    
    # Setting rewards of impossible actions to Inf
    for(j in seq(N+1,m*(N+1), by = N+1)){
      v_all[j,1] <- Inf
    }
    for(j in seq(1,m*(N+1), by = N+1)){
      v_all[j,2] <- Inf
    }
    
    v <- rbind(v, apply(v_all, 1, min))
    span <- max(v[n+1,]-v[n,]) - min(v[n+1,]-v[n,])
  }
  
  d <- A[max.col(-v_all, "first")]
  rownames(v) <- 0:n
  
  return <- list(v=v, d=d, span = span)
  return
}

result <- val.iter.algo(epsilon, A, v0, c_inspect, c_GP, c_ED, P.new)
# The d array in result indicates the best possible action in each state.
result



