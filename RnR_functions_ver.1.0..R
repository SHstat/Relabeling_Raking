library(e1071)  # classifier: svm, naiveBayes
library(rpart)  # classifier: q tree
library(class)  # classifier: knn
library(RSNNS)  # classifier: MLP
library(ebmc)   # measure


### RnR algorithm for continuous variables
RnR <- function(train, N, N.pop=1){
  #### Preparation ####
  major=train[train$y==0,-ncol(train)]  # Drop the labeling after splitting the data.
  minor=train[train$y==1,-ncol(train)]
  major.id=which(train$y==0)
  minor.id=which(train$y==1)
  
  m <- ceiling(nrow(train)/2)
  n1 <- floor(m*0.5)
  n2 <- m - n1
  
  #### Calibration for minor class ####
  mx <- apply(minor,2,mean)
  w <- rep(1/nrow(major), nrow(major))
  
  for(i in 1:N){
    for(j in 1:length(mx)){
      w=Lagrange.HH(major[,j],mx[j],w)$weight
    }
  }

  rep.minor1.id = sample(minor.id, n1, replace=T)
  rep.minor2.id = sample(major.id, n2, replace=T, prob=w)
  rep.minor.id=c(rep.minor1.id,rep.minor2.id)
  
  #### Calibration for major class ####
  sloc=which(is.na(match(major.id, unique(rep.minor2.id)))==T)
  reduced.major.id=major.id[sloc]
  reduced.major=train[reduced.major.id,]
  
  mx=apply(major,2,mean)
  w=rep(1/nrow(reduced.major),nrow(reduced.major))
  
  for(i in 1:N){
    for(j in 1:length(mx)){
      w=Lagrange.HH(reduced.major[,j],mx[j],w)$weight
    }
  }
  
  rep.major.id <- sample(reduced.major.id,m,replace=T,prob=w)
  
  #### Balanced Data ####
  rnr.data <- train[c(rep.minor.id,rep.major.id),]
  rnr.data$y[1:ceiling(nrow(train)/2)]=1
  #table(rnr.data$y)
  
  return(rnr.data)
}


### RnR algorithm for categorical variables
RnR.cate <- function(train, N){
  #### Preparation ####
  major=train[train$y==0,-ncol(train)]  # Drop the labeling after splitting the data.
  minor=train[train$y==1,-ncol(train)]
  
  
  #### Calibration for minor class ####
  w <- rep(1, nrow(major))
  mx.li <- list()
  for (i in 1:ncol(minor)){mx.li[[names(mx[i])]] <- table(minor[,i])/nrow(minor)}
  
  for(i in 1:N){
    for(j in 1:ncol(minor)){
      w=sraking(w, major[,j], mx.li[[j]])
    }
  }
  
  n2 <- n1 <- floor(nrow(train)/4)
  major.id=which(train$y==0)
  minor.id=which(train$y==1)
  rep.minor1.id = sample(minor.id, n1, replace=T)
  rep.minor2.id = sample(major.id, n2, replace=T, prob=w)
  rep.minor.id=c(rep.minor1.id,rep.minor2.id)
  
  #### Calibration for major class ####
  sloc=which(is.na(match(major.id, unique(rep.minor2.id)))==T)
  reduced.major.id=major.id[sloc]
  reduced.major=train[reduced.major.id, -ncol(train)]
  
  w=rep(1,nrow(reduced.major))
  cri <- ceiling(nrow(train)/2)
  mx.li <- list()
  for (i in 1:ncol(major)){mx.li[[names(mx[i])]] <- table(major[,i])/nrow(major)}
  
  for(i in 1:N){
    for(j in 1:ncol(reduced.major)){
      w=sraking(w, reduced.major[,j], mx.li[[j]])
    }
  }
  
  
  rep.major.id <- sample(reduced.major.id,cri,replace=T,prob=w)
  
  #### Balanced Data ####
  rnr.data <- train[c(rep.minor.id,rep.major.id),]
  rnr.data$y[1:ceiling(nrow(train)/2)]=1
  #table(rnr.data$y)
  
  return(rnr.data)
}


### getting raking weights for categorical variables
sraking=function(W,dat,pdist){
  # W: input weight
  # dat: input (sample) data (should be categorical)
  # pdist: distribution of population
  # fw: final weight
  
  n=length(W); nu=length(pdist);
  frac=sum(W)/sum(pdist)
  ratio=pdist*frac                                  # (50, 50)
  subtotal=list()
  for(i in names(pdist)) subtotal[[i]]=sum(W[which(dat==i)])
  
  fw=vector(length=n)
  for(j in 1:n){
    v=dat[j]
    fw[j]=W[j]*ratio[v]/subtotal[[v]]
  }
  return(fw)
}


add <- function(x) Reduce("+", x)


Lagrange=function(x,mx,d,cri=0){
  U <- vector(); myList <- list()
  X <- cbind(x)
  lam <- rep(0, ncol(X)); W <- as.vector(d*exp(X%*%lam))
  
  if (length(which(sign(X)==sign(mx))) >  ceiling(cri)){
    repeat{
      lam0=lam; W0 <- W
      
      U <- colSums(W0*X)-c(mx)
      
      for (i in 1:dim(X)[1]){myList[[i]] <- W0[i] * X[i,] %*% t(X[i,])}
      S<-add(myList)
      
      # Newton's method
      lam <- lam0-MASS::ginv(as(S, "matrix"), tol=256*.Machine$double.eps)%*%U
      
      W <- as.vector(d*exp(X%*%lam))
      dif <- sum((colSums(W*X) - mx)^2)
      #print(dif)
      if(dif < 0.0000001){break}
      if(sum((lam0 - lam)^2)==0){break}
    }
  }
  else {
    print(paste(names(mx), "is hard to rake."))
    W <- d
  }
  return(list(lambda=lam,weight=W))
}


# getting gmean score and AUC score
Gmean_AUC <- function(tr_li, test, by=c("SVM", "tree", "logistic", "knn", "mlp", "NB")){  
  res1=c()
  res2=c()
  
  for (i in 1:length(tr_li)){
    train <- tr_li[[i]]
    
    if (by == "SVM"){
      fit <- svm(y~., data = train, probability=TRUE, type = "C")
      pred <- predict(fit, test, probability=TRUE)
      prob <- attr(pred, "probabilities")[,2]
    }
    
    if (by == "tree"){
      fit <- rpart(y~., data=train, method="class",
                   control = rpart.control(cp = 0.0001))
      bestcp <- fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]  # Step3: Prune the tree using the best cp.
      fit<- prune(fit, cp = bestcp)
      prob <- as.vector(predict(fit, test, type="prob")[,2])
      pred <- rep(0,length(prob))
      pred[prob>=0.5] <- 1
    }
    
    if (by == "logistic"){
      fit <- glm(y~., data=train, family = binomial)
      prob <- predict(fit, test, type = 'response')
      pred <- prob>0.5
    }
    
    res <- M(test$y, pred, prob)
    res1 <- c(res1, res[[1]])
    res2 <- c(res2, res[[2]])
  }
  return(list(Gmean = res1, AUC = res2))
}

# nested function for Gmean_AUC function.
M = function(real, pred, prob){
  con.mat <- table(real, pred)
  nc=ncol(con.mat)
  
  if(nc==2){
    TNrate=con.mat[1,1]/sum(con.mat[1,])
    TPrate=con.mat[2,2]/sum(con.mat[2,])
    Gmean=sqrt(TPrate*TNrate)
  }
  
  if(nc==1){    
    Gmean=0
  }
  
  auc<-measure(label = real, probability = prob, metric = "auc")
  return(list(Gmean=Gmean, Auc=auc))
}
