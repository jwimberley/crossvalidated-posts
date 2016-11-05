Response <- function(fit) {
  return(fit$model[,1])
}

Pi <- function(fit) {
  return(fit$fitted.values)
}

DPi <- function(fit) {
  mylink <- make.link(fit$family$link)
  pi <- fit$fitted.values
  dpi <- model.matrix(fit)
  dpi <- dpi*mylink$mu.eta(mylink$linkfun(pi))
  return(dpi)
}

InformationMatrix <- function(pi,dpi) { 
  
  ppinv <- 1/(pi*(1-pi))
  mat <- apply(dpi,1,function(x){outer(x,x)})
  fisher <- mat %*% ppinv
  fisher <- matrix(fisher,ncol(dpi),ncol(dpi))
  return(fisher)
}

GenericScoreVector <- function(pi,dpi,metric) {
  
  t <- metric(TRUE,pi)
  f <- metric(FALSE,pi)
  val <- t-f
  score <- colSums(dpi*val)   
  
  return(score)
}

GenericMetric <- function(Y,pi,dpi,metric) {
  fisher <- InformationMatrix(pi,dpi)
  score <- GenericScoreVector(pi,dpi,metric)
  scorep <- solve(fisher,score)
  
  num <- length(pi)
  t <- metric(TRUE,pi)
  f <- metric(FALSE,pi)
  tf <- t-f
  obs <- mapply(function(y,t,f) { if (y) t else f },Y,t,f)
  exp <- pi*t + (1-pi)*f
  difference <- sum(obs-exp)
  var <- pi*(t-exp)^2 + (1-pi)*(f-exp)^2
  subvar <- tf %*% (dpi %*% scorep)
  variance <- sum(var)-subvar
  z <-  difference/sqrt(abs(variance))
  return(z)
}

FitGenericMetric <- function(fit,metric) {
  return(GenericMetric(Y=Response(fit),pi=Pi(fit),dpi=DPi(fit),metric))
}

deviance <- function(c, p) {
  if (c)
    return(-2*log(p))
  else
    return(-2*log(1-p))
}

pearson <- function(c, p) {
  den <- p*(1-p)
  if (c)
    return((1-p)^2/den)
  else
    return(p^2/den)
}

chch <- function(c, p) {
  if (c)
    return((1-p)*(1-p))
  else
    return(p*p)
}

RegressionFitCreator <- function(N,beta0,beta1,linkname) {
  mylink <- make.link(linkname)
  x <- NULL
  if (linkname != "identity")
    x <- rnorm(N)
  else {
    # let beta0 + beta1 x be in the range 0.2 to 0.8 to prevent numerical problems
    # thus x in the range (0.2-beta0,0.8-beta0)/beta1
    x <- runif(N,min=(0.2-beta0)/beta1,max=(0.8-beta0)/beta1)
  }
  probs <- mylink$linkinv(beta0+beta1*x)
  y <- rbinom(N,1,probs)
  fit <- glm(y~x,family=binomial(link=linkname),start=c(rnorm(1,beta0,0.1*beta0),rnorm(1,beta1,0.1*beta1)))
  #print(summary(fit))
  return(fit)
}

ToyStudy <- function(Ntoys=1000, Nevents=500, beta0=0.25,beta1=2,linkname="logit") {
  data <- data.frame("uG2"=numeric(Ntoys),"uX2"=numeric(Ntoys),"uS"=numeric(Ntoys))
  for (toy in 1:Ntoys) {
    if (toy %% 50 == 0) print(paste("On toy ",toy))
    fit <- RegressionFitCreator(Nevents,beta0,beta1,linkname)
    uG2 <- FitGenericMetric(fit,deviance)
    uX2 <- FitGenericMetric(fit,pearson)
    uS <- FitGenericMetric(fit,chch)
    data$uG2[toy] <- uG2
    data$uX2[toy] <- uX2
    data$uS[toy] <- uS
  }
  return(data)
}

BootstrappedToyStudy <- function(Ntoys=1000,Nevents=50,data,linkname="logit") {
  mylink = make.link(linkname)
  data <- data.frame("uG2"=numeric(Ntoys),"uX2"=numeric(Ntoys),"uS"=numeric(Ntoys))
  for (toy in 1:Ntoys) {
    if (toy %% 50 == 0) print(paste("On toy ",toy))
    ns <- sample(nrow(data),Nevents,replace=TRUE)
    subtable <- table[ns,]
    subtable$Eta <- mylink$linkfun(subtable$Eta)
    fit <- glm(Correct~Eta*I(Flavour/2),data=subtable,family=binomial(link=linkname))
    uG2 <- FitGenericMetric(fit,deviance)
    uX2 <- FitGenericMetric(fit,pearson)
    uS <- FitGenericMetric(fit,chch)
    data$uG2[toy] <- uG2
    data$uX2[toy] <- uX2
    data$uS[toy] <- uS
  }
  return(data)
}
