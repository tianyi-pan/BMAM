set.seed(123)
p <- 0.5 # proportion of treatment group
beta <- c(-0.8, 1.2, 1.5, 0.4)
N_list <- c(300)
n_list <- c(10)
mean.formula <- ~time*Xe
lv.formula <- ~1+time
Sigma_myfun <- matrix(c(1^2, 0, 0, 0.5^2), nrow = 2)

N <- 100
n <- 20

nclust <- rep(n,N) # number of units in each cluster
id <- rep(seq(N), nclust) # id of each unit
Xe <- rep(rep(c(1, 0), times = c(N*p, N-N*p)), times = nclust) # treatment: first half cluster 1, second half 0
time <- as.numeric(sapply( nclust, function(ZZ) {rnorm(ZZ)} )) # draw time from N(0,1)
data  <- data.frame(id, time, Xe) # generate data (without y)
data  <- data[order(data$id, data$time),] # order data by id and time

