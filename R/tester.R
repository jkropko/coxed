library(coxed)
set.seed(22902)
data <- sim.survdata(beta = c(2,3,4))
model <- coxph(Surv(y, failed) ~ X1 + X2 + X3, data=data$data)
coef(model)


b <- data$betas
x <- as.matrix(data$xdata)
xb <- x %*% b
s <- data$baseline$survivor
survival <- t(sapply(xb, FUN=function(x){s^exp(x)}, simplify=TRUE))
y <- apply(survival, 1, FUN=function(x){
     z <- diff(x < runif(1))
     r <- ifelse(all(z==0), T, which.max(z))
     return(r)
})

table(xb == data$xb)


table(model$y[,1]==y) ## not the same???
table(data$data$y==y) ## not the same???

# problem is that data is not the same as y ....

set.seed(22902)
baseline <- baseline.build(T=100, knots=8, spline=TRUE)
xb <- generate.lm(baseline, X=NULL, beta=c(2,3,4), N=1000)

table(xb$data$y==y) ## not the same???
table(xb$data$y==model$y[,1]) ## not the same???
table(xb$data$y==data$data$y) ## not the same???
