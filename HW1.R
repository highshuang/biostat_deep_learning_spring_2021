## PHS597 Deep learning HW1 implement PPR 
set.seed(1)
library("scatterplot3d") 

# Simulate data 
n <- 300
x1 <- rnorm(n, 0, 1)
x2 <- rnorm(n, 0, 1)
y <- 1/(1+exp(-5*((x1+x2)/sqrt(2)-0.5))) + rnorm(n, 0, 0.2)
#y <- 1/(1+exp(x1+x2))
data <- data.frame(y,x1,x2)

scatterplot3d(x=data$x1, y = data$x2, z = data$y)

# implement PPR iteritively with three knots 
PPR_imp2 <- function(data, n_iter=1000) {
  # 
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  X <- as.matrix(data[, -1])
  Y <- data[, 1]
  
  # randomly initialize omega 
  w_old <- rep(1, p)
  
  # normalize vector 
  w_old <- w_old / sqrt(sum(w_old^2))
  
  for (iter in 1:n_iter) {
    # Step 1: given omega get g() with natural cubic splines 
    V <- X %*% w_old
    
    # set the knots based on quantile 
    # get the knots based on the number of data points 
    if (length(V)<50) {
      knot <- sort(V) # type is vector or matrix check 
    } else {
      knot <- unname(quantile(V, probs = seq(0, 1, length.out = 50)))
    }
    
    n_knot <- length(knot)
    
    
    # calcualte indicator matrix for each knot 
    Indicator <- matrix(rep(0, n*n_knot), nrow = n)
    
    for (j in 1:n_knot) {
      Indicator[,j] <- ifelse(V>= knot[j], (V-knot[j])^3, 0)
    }
    
    # calcualte the basis
    M <- matrix(rep(0, n*n_knot), nrow = n)
    M[,1] <- rep(1, n)
    M[,2] <- V
    
    for (j in 1:(n_knot-2)) {
      M[,(j+2)] <- (Indicator[,j] - Indicator[,n_knot])/(knot[j]-knot[n_knot]) - 
        (Indicator[,(n_knot-1)] - Indicator[,n_knot])/(knot[(n_knot-1)]-knot[n_knot])
    }
    
    # least square estimates 
    beta_g <- ginv(t(M) %*% M) %*% t(M) %*% Y
    
    # Step 2: given g estimates omega 
    # get g'(omega_oldTxi) = g'(v)
    g_der <- sapply(X = V, FUN=g_derivative, knot = knot, beta_g=beta_g)
    
    # update omega 
    # calculate g(omega_oldTxi)
    scaler_vec <- g_der^2 * (V + (Y - M %*% beta_g) / g_der)
    num_sum <- rep(0, p)
    den_sum <- matrix(rep(0, p^2), nrow = p)
    
    # calculate sum for denominator and numerator 
    for (i in 1:n) {
      num_sum <- num_sum + scaler_vec[i] * X[i, ]
      xi <- as.vector(X[i, ])
      den_sum <- den_sum + xi %*% t(xi) * (g_der^2)[i]
    }
    
    # calculate new omega 
    w_new <- ginv(den_sum) %*% as.vector(num_sum)
    
    # normalize w_new 
    w_new <- w_new / sqrt(sum(w_new^2))
    
    tol <- max(abs(w_new - w_old)) 
    if(tol < 1e-6) {
      print("Method2: Convergence !!!")
      return(list("g"=beta_g, "omega"=w_new, "tol"=tol, "knot"=knot))
    }
    
    print("Not Convergence")
    w_old <- w_new 
    
    
  }#loop end 
  return(list("g_not"=beta_g, "omega"=w_new, "tol"=tol, "knot"=knot))
} # function end




g_derivative <- function(x, knot, beta_g) {
  derivative <- beta_g[2]#
  K <- length(knot)
  pos <- findInterval(x, knot)
  
  if (pos <= (K-2)) {
    for (j in 1:pos) {
      derivative <- derivative + before_K_1(x = x, knot_j = knot[j], knot_K = knot[K],
                                            beta_j = beta_g[(j+2)])
    }
  } else {
    derivative <- derivative + after_K_1(x = x, knot = knot, beta_g = beta_g)
  }
  return(derivative)
}

# two additional functions for calcualte derivative of g function 
# for x is between knot j and knot K-1
before_K_1 <- function(x, knot_j, knot_K, beta_j) {
  return(3*beta_j*(x-knot_j)^2/(knot_j-knot_K))
}

after_K_1 <- function(x, knot, beta_g) {
  K <- length(knot)
  sum <- 0 
  
  for (j in 1:(K-2)) {
    sum <- sum + 3*beta_g[(j+2)]*(x-knot[j])^2/(knot[j]-knot[K]) - 3*beta_g[(j+2)]*(x-knot[(K-1)])^2/(knot[(K-1)]-knot[K])
  }
  
  return(sum)
}


## test 
res <- PPR_imp2(data = data, n_iter = 1000)
res$omega

fit <- ppr(x = data[,-1], y=data[,1], nterms = 1, sm.method='spline', df = 3)

summary(fit)
fit$alpha

#### plotting data 
plotdata <- as.data.frame(rbind(data[,-1], data[,-1], data[,-1]))
plotdata$pred <- c(data[,1], as.matrix(data[,-1]) %*% res$omega, as.matrix(data[,-1]) %*% fit$alpha)
plotdata$group <- factor(c(rep("raw data", n), rep("self2", n), rep("ppr", n)))

colors <- c("#ea9df2", "#E69F00", "#56B4E9")
colors <- colors[as.numeric(plotdata$group)]


p1<-scatterplot3d(plotdata[,1:3],
              main="Comparison between Method2 and \n R functions for PPR",
              xlab = "x1",
              ylab = "x2",
              zlab = "y", pch=1, color = colors)

legend(p1$xyz.convert(2.5, 0, 0), legend = levels(plotdata$group),
       col =  c("#ea9df2", "#E69F00", "#56B4E9"), pch = 1)

levels(plotdata$group)[3]
