# ---- QA functions ----

normalize <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    col_means <- colMeans(x)
    col_stds <- apply(x, 2, sd)
    normalized_x <- scale(x, center = col_means, scale = col_stds)
  } else {
    mean_val <- mean(x)
    std_val <- sd(x)
    normalized_x <- (x - mean_val) / std_val
  }
  normalized_x <- replace(normalized_x, is.nan(normalized_x), 0)
  return(normalized_x)
}

objective <- function(w, xs, ds, Tx, k1, k2) {
  p1 <- sum((w-ds)^2)
  p2 <- sum((colSums(w * xs) - Tx)^2)
  
  k1*p1 + k2*p2
}

gradient <- function(w, ds, xs, Tx, k1, k2) {
  dw_p1 <- 2 * (w - ds)
  dw_p2 <- 2 * xs %*% (colSums(w * xs) - Tx)
  
  dw_p1 <- normalize(dw_p1)
  dw_p2 <- normalize(dw_p2)

  
  k1*dw_p1+k2*dw_p2
}

#######################################################
# R Function for gradient descent method              #
#                                                     #
# Function to minimize f(w)                           #
# f = k1 * (sum(w-ds)^2) + k2 * sum((sum(w*xs)-Tx)^2) #
#                                                     #
#                                                     #
# Input: w = (w_1, w_2,..., w_n)                      #
#        ds = (d_1, d_2,..., d_n)                     #
#        xs = (x_1, x_2,..., x_n)                     #
#        Tx = (T_x1, T_x2,..., T_xn)                  #
#        alpha = learning rate                        #
#        decay = decay rate                           #
#        k1 & k2 = scaling constants                  #
#                                                     #
#                                                     #
# Output: w = (w_1, w_2,..., w_n)                     #
#                                                     #
# Written by: W. Sutherland                           #
#######################################################


grad_desc <- function(w, ds, xs, Tx, alpha = 0.0001, decay = 0, maxiter = 10000, tol = 1e-8, k1, k2) {
  #pb <- txtProgressBar(min = 0,max = maxiter, style = 3,  width = 50,char = "=")
  
  iter <- 0
  flag <- 0
  prev_obj <- objective(w, xs, ds, Tx, k1, k2)
  
  while (iter < maxiter) {
    gradient <- gradient(w, ds, xs, Tx, k1, k2)
    
    w_new <- as.vector(w - alpha * gradient)
    
    
    obj <- objective(w_new, xs, ds, Tx, k1, k2)
    
    obj_change <- abs(prev_obj - obj)
    
    
    if (is.na(obj_change)){
      flag <- 1
      break
    }
    
    if (obj_change < tol) {
      break
    }
    
    w <- w_new
    prev_obj <- obj
    
    alpha <- ifelse(iter >= 10000, alpha / (1 + decay * iter), alpha )
    iter <- iter + 1
    
    if (iter == maxiter) {
      flag <- 1
    }
    #setTxtProgressBar(pb, iter)
  }
  #close(pb)
  return(list(w = w, flag = flag, iter = iter, obj = obj, obj_change = obj_change, alpha = alpha))
}
