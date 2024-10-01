# Source the main function file
# the main function file contains Wu & Lu functions
# as well as functions written by myself
source('01_functions/functions.R')

# Source the QA function file
# this file contains in-progress functions
# being tested/actively worked on
source('01_functions/QA_functions.R')


# ---- 0.1 Parameter setting ----
#===========================================#
# Parameter setting
#===========================================#

N <- 4000; n <- 100
mm <- 100; nsim <- 2

# ---- 0.2 Finite population setting ----

# ==========================================#
# Finite population generation
# ==========================================#

X1 <- rlnorm(N); X2 <- rchisq(N,2)
X3 <- X2*(X2<2); X4 <- sqrt(X1)
X <- cbind(X1, X2, X3, X4)

Y1 <- 1+X1+X2+X3+X4+2*rnorm(N); Y1 <- Y1-min(Y1)
Y2 <- 1+X1+X2+X3+X4+26*rnorm(N); Y2 <- Y2-min(Y2)
Y3 <- rexp(N)
Y <- cbind(Y1, Y2, Y3)

vpi <- (max(X2)-mm*min(X2))/(mm-1)
W <- X2+vpi; W <- n*W/sum(W)
d <- 1/W; qs <- rep(1,n)

pop <- cbind(X, Y, d, W)
Tx <- colSums(X)
My <- colMeans(Y)
Ty <- colSums(Y)

# ---- 0.3 Generate tables to store results ----

## Matrices to save estimates of mu_Y
my_HT <- matrix(0, nsim, ncol(Y))
my_gr <- my_HT
my_gd1 <- my_HT; my_gd2 <- my_HT; my_gd3 <- my_HT; my_gd4 <- my_HT

## Matrices to save distances from aux. totals
dist_x_gd1 <- matrix(0, nsim, ncol(X))
dist_x_gd2 <- dist_x_gd1; dist_x_gd3 <- dist_x_gd1; dist_x_gd4 <- dist_x_gd1
dist_x_gr <- dist_x_gd1; dist_rm <- dist_x_gd1

## Vectors to save distance from design weights
dist_w_gd1 <- rep(0,nsim); dist_w_gd2 <- rep(0,nsim); dist_w_gd3 <- rep(0,nsim); dist_w_gd4 <- rep(0,nsim); 
dist_w_gr <- rep(0, nsim)

## Vectors for non-convergence flags and no. of iterations
flags1 <- rep(0, nsim); flags2 <- flags1; flags3 <- flags1; flags4 <- flags1
iter1 <- flags1; iter2 <- flags1; iter3 <- flags1; iter4 <- flags1;
obj1 <- flags1; obj2 <- flags1; obj3 <- flags1; obj4 <- flags1

# ---- 1.0 Simulation Study ----

pb <- txtProgressBar(min = 0, max = nsim, style = 3, width = 50, char = "=")
for (i in 1:nsim){
  
  ## ---- 1.1 Sampling and calculate weights using GREG ----
  # Take the first sample using the Rao-Sampford method
  idx <- RSsample(W,n)
  xs <- X[idx,]
  ys <- Y[idx,]
  ds <- d[idx]
  
  # Calculate estimates using the GREG calibrator
  GR <- cal_est('GR', xs, ys, ds, qs, Tx, maxiter = 100000, initial = 0, c1, c2)
  
  ## ---- 1.2 Use iterative method ----
  # k1 = k2 = 1
  GD1 <- grad_desc(GR$w_cal, ds, xs, Tx,maxiter = 100000, 
                   tol = 1e-10, k2 = 1, k1 = 1, alpha = 0.01, decay = 1e-8)
  # k1 = 1, k2 = 8
  GD2 <- grad_desc(GR$w_cal, ds, xs, Tx,maxiter = 100000, 
                   tol = 1e-10, k2 = 8, k1 = 1, alpha = 0.01, decay = 1e-8)
  # k1 = 1, k2 = 2
  GD3 <- grad_desc(GR$w_cal, ds, xs, Tx,maxiter = 100000, 
                   tol = 1e-10, k2 = 2, k1 = 1, alpha = 0.01, decay = 1e-8)
  flags3[i] <- GD3$flag
  iter3[i] <- GD3$iter
  obj3[i] <- GD3$obj
  
  # k1 = 4, k2 = 1
  GD4 <- grad_desc(GR$w_cal, ds, xs, Tx,maxiter = 100000, 
                   tol = 1e-10, k2 = 1, k1 = 4, alpha = 0.01, decay = 1e-8)
  
  ## ---- 1.3 Store results ----
  
  # Flags (0 = converged, 1 = fail to converge)
  flags1[i] <- GD1$flag
  flags2[i] <- GD2$flag
  flags3[i] <- GD3$flag
  flags4[i] <- GD4$flag
  
  # Iterations
  iter1[i] <- GD1$iter
  iter2[i] <- GD2$iter 
  iter3[i] <- GD3$iter
  iter4[i] <- GD4$iter
  
  # Value of the objective function
  obj1[i] <- GD1$obj
  obj2[i] <- GD2$obj
  obj3[i] <- GD3$obj
  obj4[i] <- GD4$obj
  
  ## ---- 1.4 Calculate eval. metrics ----
  my_HT[i,] <- colSums(ds*ys)/N
  dist_rm[i,] <- colSums(xs*ds)-Tx
  
  my_gr[i,] <- colSums(GR$w_cal*ys)/N
  dist_w_gr[i] <- sum((GR$w_cal-ds)^2)
  dist_x_gr[i,] <- colSums(GR$w_cal*xs)-Tx
  
  # mean of Y
  my_gd1[i,] <- colSums(GD1$w*ys)/N
  my_gd2[i,] <- colSums(GD2$w*ys)/N
  my_gd3[i,] <- colSums(GD3$w*ys)/N
  my_gd4[i,] <- colSums(GD4$w*ys)/N
  
  # distance from aux. totals
  dist_x_gd1[i,] <- colSums(GD1$w*xs)-Tx
  dist_x_gd2[i,] <- colSums(GD2$w*xs)-Tx
  dist_x_gd3[i,] <- colSums(GD3$w*xs)-Tx
  dist_x_gd4[i,] <- colSums(GD4$w*xs)-Tx
  
  # distance from design weights
  dist_w_gd1[i] <- sum((GD1$w-ds)^2)
  dist_w_gd2[i] <- sum((GD2$w-ds)^2)
  dist_w_gd3[i] <- sum((GD3$w-ds)^2)
  dist_w_gd4[i] <- sum((GD4$w-ds)^2)
  
  setTxtProgressBar(pb, i)
}

close(pb)


# ---- 2.0 Results ---- 

## ---- 2.1 Means of Y ----
means <- rbind(My, colMeans(my_HT), colMeans(my_gr), colMeans(my_gd1[flags1==0,]), colMeans(my_gd2[flags2==0,]),
               colMeans(my_gd3[flags3==0,]), colMeans(my_gd4[flags4==0,]))
names <- rbind('True value', 'HT', 'GR', 'k2 = k1 = 1', 'k2 = 8, k1 = 1', 'k2 = 2, k1 = 1',
               'k2 = 1, k1 = 4')

dat <- data.frame(means)
rownames(dat) <- names

## ---- 2.2 Distance from design weights ----
dist_w <- rbind(mean(dist_w_gr), mean(dist_w_gd1[flags1==0]), mean(dist_w_gd2[flags2==0]), 
                mean(dist_w_gd3[flags3==0]),mean(dist_w_gd4[flags4==0]))

dat_w <- data.frame(dist_w)
names_w <- rbind('GR', 'k2 = k1 = 1', 'k2 = 8, k1 = 1/8', 'k2 = 2, k1 = 1/2',
                 'k2 = 1/4, k1 = 4')
rownames(dat_w) <- names_w
colnames(dat_w) <- ''

## ---- 2.3 Distance from aux. totals ----
dist_x <- rbind(colMeans(dist_x_gr), colMeans(dist_x_gd1[flags1==0,]), colMeans(dist_x_gd2[flags2==0,]), 
                colMeans(dist_x_gd3[flags3==0,]),colMeans(dist_x_gd4[flags4==0,]))


dat_x <- data.frame(dist_x)
rownames(dat_x) <- names_w

## ---- 2.4 Iterations ----
iter_mean <- rbind(mean(iter1), mean(iter2), mean(iter3), mean(iter4))
names_iter <- rbind('k2 = k1 = 1', 'k2 = 8, k1 = 1/8', 'k2 = 2, k1 = 1/2',
                    'k2 = 1/4, k1 = 4')
dat_it <- data.frame(iter_mean)
rownames(dat_it) <- names_iter
colnames(dat_it) <- ''

## ---- 2.5 Ratio of first moment ----
dat_rm <- data.frame(dat_x/colMeans(dist_rm))
rownames(dat_rm) <- names_w

## ---- 2.6 Final Summary (all metrics) ----
{
  cat('Mean of Y:\n')
  print(dat)
  cat('\n\n')
  cat('Distance from design weights:\n')
  print(dat_w)
  cat('\n\n')
  cat('Distance from aux totals:\n')
  print(dat_x)
  cat('\n\n')
  cat('Mean number of iterations:\n')
  print(dat_it)
  cat('\n\n')
  cat('First ratio moments:\n')
  print(dat_rm)
}


