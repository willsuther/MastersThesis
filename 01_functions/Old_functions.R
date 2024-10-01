##################################################
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
##################################################
dw_part1OLD <- function(w, ds) {
  2 * (w - ds)
}

dw_part2OLD <- function(w, xs, Tx) {
  2 * xs %*% (colSums(w * xs) - Tx)
}

gradOLD <- function(w, ds, xs, Tx, k1, k2) {
  dw_part1 <- k1*dw_part1OLD(w,ds)
  dw_part2 <- k2*dw_part2OLD(w,xs,Tx)
  
  return(list(dw_part1 = dw_part1, dw_part2 = dw_part2))
}

objectiveOLD <- function(w, xs, ds, Tx, k1, k2) {
  k1*sum((w-ds)^2) + k2*sum((colSums(w * xs) - Tx)^2)
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


grad_descOLD <- function(w, ds, xs, Tx, alpha = 0.0001, decay = 0, maxiter = 10000, tol = 1e-8, k1, k2) {
  iter <- 0
  flag <- 0
  prev_obj <- objectiveOLD(w, xs, ds, Tx, k1, k2)
  obj_out <- rep(0,maxiter)
  
  while (iter < maxiter) {
    gradient <- gradOLD(w, ds, xs, Tx, k1, k2)
    
    w_new <- as.vector(w - alpha * (gradient$dw_part1 + gradient$dw_part2))
    
    
    obj <- objectiveOLD(w_new, xs, ds, Tx, k1, k2)
    
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
    
    iter <- iter + 1
    alpha <- alpha / (1 + decay * iter)
    obj_out[iter] <- obj
    
    if (iter == maxiter) {
      flag <- 1
    }
  }
  
  return(list(w = w, flag = flag, iter = iter, obj_change = obj_change))
}