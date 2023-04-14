library(ggplot2)
#-------------------------------------------------------------------------------
# Global Parameters
b_0 = 1.15
d = 1
#-------------------------------------------------------------------------------
# Cell specs
mu <- 0.055 # Driver mutation rate
nu <- 0.3 # Passenger mutation rate
s_p <- 0.002 # Passenger strength
L <- 2
#-------------------------------------------------------------------------------
# Birth rate as a function of type i
b <- function(i) {
  return(b_0*( (1+s_p)^i ))
}

# Probability generating functions
pgf <- function(i,x,y,z) {
  denom = b(i)+mu+nu+d
  return((d/denom) + (nu*x/denom) + (b(i)*(y^2)/denom) + (mu*z/denom))
}
#-------------------------------------------------------------------------------
# Simulation specs
n = 100 # >= t+nL+1 or <= t-nL-1 are taboo types
iter = 1000 # Number of iterations for Hautphenne's algorithm
N = 50*L # Decide the number of initial types
ini_types = seq(from = (-N), to = N, by = 1) # Initial cell types
ini_types_len = 2*N + 1

# Initialize extinction probabilities
ext_prob = rep(0,ini_types_len)

for (t in 1:ini_types_len) {
  w = rep(0,(2*n*L)+1) # Reset w when initial type changes
  current_type = ini_types[t]
  non_taboo_types = current_type+((-(n*L)):(n*L))
  
  for (i in 1:iter){
    w_holder = vector() # For updating purpose
    
    for (j in non_taboo_types) {
      index = j + (n*L) - current_type + 1
      
      if (j == (current_type-n*L)) {
        w_holder[index] = pgf(j,0,w[index],w[index+L])
      } else if (!((j+L) %in% non_taboo_types)) {
        w_holder[index] = pgf(j,w[index-1],w[index],0)
      } else {
        w_holder[index] = pgf(j,w[index-1],w[index],w[index+L])
      }
    }
    
    w = w_holder # Update w and it approaches ext prob
  }
  
  ext_prob[t] = w[(n*L)+1] # Record ext prob when initialize by current_type
}

type_vs_ext_P1 <- as.data.frame(cbind(ini_types, ext_prob))

ggplot(type_vs_ext_P1, aes(x=ini_types, y=ext_prob)) + 
  geom_point() + xlab("Initial Types")+ylab("Extinction Probability")
ggsave("ext_prob_P1.png")

write.csv(type_vs_ext_P1,"type_vs_ext_P1.csv")






























