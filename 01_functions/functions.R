# ---- GREG Lagrange Multiplier ----

#================================================#
# R Function for GREG Lagrange Multiplier        #
#                                                #
# Input: xs = (x_1, x_2,..., x_n)                #
#        ds = (d_1, d_2,..., d_n)                #
#        qs = (q_1, q_2,..., q_n)                #
#        Tx = (T_x1, T_x2,..., T_xn)             #
#                                                #
# Output: lambda = lambda                        #
#                                                #
#                                                #
# Written by: W. Sutherland                      #
#================================================#

lambda_GREG <- function(xs,ds,qs,Tx){
  eps <- 1e-8
  txHT <- colSums(ds*xs)
  d1 <- t(xs)%*%diag(ds*qs)%*%xs
  det <- abs(det(d1))
  if (det <= eps){
    lambda <- rep(9999, length(Tx))
  } else {
    lambda <- solve(d1,Tx-txHT) 
  }
  return(lambda = lambda)
}

# ---- ET Lagrange Multiplier ----

#================================================#
# R Function for ET Lagrange Multiplier          #
#                                                #
# Input: xs = (x_1, x_2,..., x_n)                #
#        ds = (d_1, d_2,..., d_n)                #
#        qs = (q_1, q_2,..., q_n)                #
#        Tx = (T_x1, T_x2,..., T_xn)             #
#        maxiter = maximum no. of iterations     #
#        initial = starting value of lambda      #
#                                                #
# Output: lambda = lambda                        #
#         iterations = no. of iterations         #
#                                                #
# lambda_ET is modified from lag4                #
#                                                #
# Modified by: W. Sutherland                     #
#================================================#

lambda_ET <- function(xs, ds, qs, Tx, maxiter = 50, initial = 0){
  lambda <- rep(initial, length(Tx))
  eps <- 1e-8
  eps_1 <- .Machine$double.eps
  run <- 0
  for (i in 1:maxiter){
    d1 <- colSums(ds*c(exp(xs%*%lambda*qs))*xs)-Tx
    d2 <- t(xs)%*%(c(qs*ds*exp((xs%*%lambda)*qs))*xs)
    det <- abs(det(d2))
    rec <- rcond(d2)
    if (is.na(det) || rec < eps_1){
      lambda <- lambda+9999
      break
    }
    if (det <= eps){
      lambda <- lambda + 9999
      break
    }
    if (i == maxiter){
      lambda <- lambda + 9999
      break
    }
    if (det > eps){
      dd <- solve(d2,d1,eps_1)
      dif <- max(abs(dd))
      run <- run+1
    }
    if (dif <= eps){
      break
    }
    lambda <- lambda-dd
  }
  return(list(lambda = as.vector(lambda), iterations = run))
}

# ---- EL Lagrange Multiplier ----

#================================================#
# R Function for EL Lagrange Multiplier          #
#                                                #
# Input: xs = (x_1, x_2,..., x_n)                #
#        ds = (d_1, d_2,..., d_n)                #
#        qs = (q_1, q_2,..., q_n)                #
#        Tx = (T_x1, T_x2,..., T_xn)             #
#        maxiter = maximum no. of iterations     #
#        initial = starting value of lambda      #
#                                                #
# Output: lambda = lambda                        #
#         iterations = no. of iterations         #
#                                                #
# lambda_EL is modified from lag3                #
#                                                #
# Modified by: W. Sutherland                     #
#================================================#

lambda_EL <- function(xs, ds, qs, Tx, maxiter = 50, initial = 0){
  lambda <- rep(initial, length(Tx))
  eps <- 1e-8
  eps_1 <- .Machine$double.eps
  run <- 0
  for (i in 1:maxiter){
    d1 <- t(xs)%*%((ds/(1+qs*xs%*%lambda))*qs)-Tx
    d2 <- -t(xs)%*%(c((qs*ds/(1+qs*(xs%*%lambda))^2))*xs)
    det <- abs(det(d2))
    rec <- rcond(d2)
    if (is.na(det) || rec < eps_1){
      lambda <- lambda+9999
      break
    }
    if (det <= eps){
      lambda <- lambda + 9999
      break
    }
    if (det > eps){
      dd <- solve(d2,d1,eps_1)
      dif <- max(abs(dd))
      run <- run+1
    }
    if (i == maxiter){
      lambda <- lambda + 9999
      break
    }
    if (min(1+qs*(xs%*%(lambda-dd)))<= 0){
      dd <- dd/2
    }
    if (dif <= eps){
      break
    }
    lambda <- lambda-dd
  }
  return(list(lambda = as.vector(lambda), iterations = run))
}

# ---- Trimmed weights ----

#================================================#
# R Function for trimmed weights                 #
#                                                #
# Input: w_cal = (w_1, w_2,..., w_n)             #
#        g_cal = (g_1, g_2,..., g_n)             #
#        ys = (y_1, y_2,..., y_n)                #
#        ds = (d_1, d_2,..., d_n)                #
#        c1 & c2 = the weight constraints        #
#                                                #
# Output: trim_w_cal = trim_w_cal                #
#         ys = ys                                #
#                                                #
# trim_wts is modified from TRIM                 #
#                                                #
# Modified by: W. Sutherland                     #
#================================================#

trim_wts <- function(w_cal,g_cal,ys,ds,c1,c2){
  n <- length(g_cal)
  min <- min(g_cal)
  max <- max(g_cal)
  
  id1 <- (1:n)[g_cal < c1]
  id2 <- (1:n)[g_cal >= c1 & g_cal <= c2]
  id3 <- (1:n)[g_cal > c2]
  idx_cal =c(id1,id2,id3)
  
  w_cal_t1 <- sum(w_cal[id1])
  w_cal_t2 <- sum(w_cal[id2])
  w_cal_t3 <- sum(w_cal[id3])
  
  ds_t1 <- sum(ds[id1])
  ds_t3 <- sum(ds[id3])
  if(max >= c1 & min <=c2){
    trim_w_cal  <- w_cal
    ys0 <-  ys
  }
  k <- 1-(c1*ds_t1+c2*ds_t3-w_cal_t1-w_cal_t3)/w_cal_t2
  trim_w_cal <- c(c1*ds[id1],k*w_cal[id2],c2*ds[id3])
  if (!is.null(ncol(ys))){
    ys0 <- ys[idx_cal,]
  }
  if (is.null(ncol(ys))){
    ys0 <- ys[idx_cal]
  }
  return(list(trim_w_cal = trim_w_cal ,ys0 = ys0, idx_cal = idx_cal))
}

# ---- Calibration of estimators ----

#==========================================================#
# R Function for calibration of our estimators             #
#                                                          #
# Input: method = 'GR' or 'ET' or 'EL'                     # 
#        xs = (x_1, x_2,..., x_n)                          # 
#        ys = (y_1, y_2,..., y_n)                          #  
#        ds = (d_1, d_2,..., d_n)                          #
#        qs = (q_1, q_2,..., q_n)                          #
#        Tx = (T_x1, T_x2,..., T_xn)                       #
#        maxiter = maximum no. of iterations               #
#        initial = starting value of lambda                # 
#        c1 & c2 = the weight constraints                  #
#                                                          #
# Output: ty_cal = calibrated total                        #
#         trim_ty_cal = constrained calibrated total       #
#                                                          #
# cal_est is modified from cal_est                         #
#                                                          #
# Modified by: W. Sutherland                               #
#==========================================================#

cal_est <- function(method, xs, ys, ds, qs, Tx, maxiter = 50, initial = 0, c1, c2){
  flag <- 0
  if (method == 'GR'){
    lambda <- lambda_GREG(xs,ds,qs,Tx)
    if (max(lambda)>9900){
      flag <- 1
      g_cal <- qs
    }
    g_cal <- 1+as.vector(xs%*%lambda)*qs
  }
  if (method == 'ET'){
    lambda <- lambda_ET(xs,ds,qs,Tx,maxiter,initial)$lambda
    if (max(lambda)>9900){
      flag <- 1
      g_cal <- qs
    }
    g_cal <- as.vector(exp(xs%*%lambda)*qs)
  }
  if (method == 'EL'){
    lambda <- lambda_EL(xs,ds,qs,Tx,maxiter,initial)$lambda
    if (max(lambda)>9900){
      flag <- 1
      g_cal <- qs
    }
    g_cal <- as.vector(1/(1+(xs%*%lambda)*qs))
  }
  w_cal <- ds*g_cal
  trimmed <- trim_wts(w_cal, g_cal, ys, ds, c1, c2)
  trim_w_cal <- trimmed$trim_w_cal
  ys0 <- trimmed$ys0
  idx_cal <- trimmed$idx_cal
  if (flag == 0 & !is.null(ncol(ys))){
    ty_cal <- colSums(w_cal*ys)
    trim_ty_cal <- colSums(trim_w_cal*ys0)
  }
  if (flag == 0 & is.null(ncol(ys))){
    ty_cal <- sum(w_cal*ys)
    trim_ty_cal <- sum(trim_w_cal*ys0)
  }
  if (flag == 1){
    ty_cal <- 0
    trim_ty_cal <- 0
  }
  return(list(flag = flag, ty_cal = ty_cal , trim_ty_cal = trim_ty_cal, 
              trim_w_cal = trim_w_cal, w_cal = w_cal, idx_cal = idx_cal))
}






















# ---- Rao-Sampford sampling method ----
#====================================#
## R function to draw a pps sample   ##
## Rao-Sampford pps sampling method  ##
##                                   ##
## Input: p is the Nx1 popu vector   ##
##        for the size variable      ##
##        n is the sample size       ##
## Output: the set of sampled units  ##
##                                   ##
## Written by Changbao Wu, July 2004 ##
#====================================#
RSsample<-function(p,n){
  N<-length(p)
  p=p/sum(p)
  lam<-rep(0,N+1)
  for(i in 1:N) lam[i+1]<-lam[i]+p[i]
  q<-rep(0,N)
  for(i in 1:N) q[i]<-p[i]/(1-n*p[i])
  q<-q/sum(q)
  lam2<-rep(0,N+1)
  for(i in 1:N) lam2[i+1]<-lam2[i]+q[i]
  ntot<-1
  while(ntot<1963){
    sam<-NULL
    rand<-runif(1)
    dif<-1
    L<-1
    U<-N+1
    while(dif>0){
      M<-floor((U-L)/2)
      if(lam[L+M]>rand) U<-U-M
      if(lam[L+M]<=rand) L<-L+M
      if(lam[U]>=rand & lam[U-1]<rand){
        si<-U-1
        dif<-0
      }
      if(lam[L]<rand & lam[L+1]>=rand){
        si<-L
        dif<-0
      }
    }
    sam<-c(sam,si)
    nn<-0
    while(nn<n-1){
      dif<-1
      rand<-runif(1)
      L<-1
      U<-N+1
      while(dif>0){
        M<-floor((U-L)/2)
        if(lam2[L+M]>rand) U<-U-M
        if(lam2[L+M]<=rand) L<-L+M
        if(lam2[U]>=rand & lam2[U-1]<rand){
          si<-U-1
          dif<-0
        }
        if(lam2[L]<rand & lam2[L+1]>=rand){
          si<-L
          dif<-0
        }
      }
      if(min(abs(si-sam))==0) nn<-n+1
      if(min(abs(si-sam))>0){
        sam<-c(sam,si)
        nn<-nn+1
      }
    }
    if(nn==n+1) ntot<-1
    if(nn==n-1) ntot<-1989
  }
  return(sam)
}