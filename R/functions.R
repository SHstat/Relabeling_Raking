######################################################################
# Functions: for index of training data set
######################################################################
generate_train_test <- function(ind_maj, ind_min, seed = 1){
  
  #' Make train and test dataset
  #'
  #' @param data dataset
  #' @param ind_maj index of major instances
  #' @param ind_min index of minor instances
  #' @param seed seed for sampling training dataset
  
  n_tr_maj = floor(0.5*length(ind_maj))      # the number of major instances in training set 
  n_tr_min = floor(0.5*length(ind_min))      # the number of minor instances in training set
  
  set.seed(seed)
  training_maj <- sample(ind_maj, n_tr_maj)
  training_min <- sample(ind_min, n_tr_min)
  tr.id <- c(training_maj, training_min)     # index of training dataset
  
  return(tr.id)
}

######################################################################
# Functions: for RnR Algorithm
######################################################################
add <- function(x) Reduce("+", x)

Lagrange <- function(x,mx,d,cri=0, type){
  
  # numeric
  if (type == 'numeric'){
    U <- vector(); myList <- list()
    X <- cbind(x)
    lam <- rep(0, ncol(X)); W <- as.vector(d*exp(X%*%lam))
    if (length(which(sign(X)==sign(mx[[1]]))) >  ceiling(cri)){
      repeat{
        lam0=lam; W0 <- W
        
        U <- colSums(W0*X)-c(mx[[1]])
        
        for (i in 1:dim(X)[1]){myList[[i]] <- W0[i] * X[i,] %*% t(X[i,])}
        S<-add(myList)
        
        # Newton's method
        lam <- lam0-MASS::ginv(as(S, "matrix"), tol=256*.Machine$double.eps)%*%U
        
        W <- as.vector(d*exp(X%*%lam))
        dif <- sum((colSums(W*X) - mx[[1]])^2)
        #print(dif)
        if(dif < 0.0000001){break}
        if(sum((lam0 - lam)^2)==0){break}
      }
    }else {
      print(paste(names(mx), "is hard to rake."))
      W <- d
    }
  }else {# sraking
    n=length(d); nu=length(mx);
    frac=sum(d)/sum(mx)
    ratio=mx*frac
    subtotal=list()
    for(i in names(mx)) subtotal[[i]]=sum(d[which(x==i)])
    
    W=vector(length=n)
    for(j in 1:n){
      v=x[j]
      if (!v %in% names(ratio)){
        ratio[[v]] <- 1
        subtotal[[v]] <- 1
      }
      W[j]=d[j]*ratio[v]/subtotal[[v]]
    }
  }
  return(W)
}

RnR <- function(train, N = 10, v_id = NULL){
  
  #' Run the RnR algorithm 
  #'
  #' @param train train dataset
  #' @param N the number of iteration in raking.
  #' @param v_id the column name of selected variables.
  #' default = 10
  
  #####################################
  #### Preparation ####
  major=train[train$y==0,-ncol(train)]  # major instances
  minor=train[train$y==1,-ncol(train)]  # minor instances
  major.id=which(train$y==0)            # location index of major instances (not index)
  minor.id=which(train$y==1)            # location index of minor instances (not index)
  
  #### Subset columns based on variable selection results
  tmp <- major %>% summary.default %>% as.data.frame %>% 
    dplyr::group_by(Var1) %>%  tidyr::spread(key = Var2, value = Freq)
  
  if (length(v_id) != 0){
    tmp <- filter(tmp, Var1 %in% v_id)
    ncol <- length(v_id)
  }else{ncol <- ncol(major)}
  
  #####################################
  #### Calibration for minor class ####
  w <- rep(1/nrow(major), nrow(major))  # default weights for major instances
  
  for(i in 1:N){
    for(j in 1:ncol){
      type <- subset(tmp, Var1 == tmp$Var1[j])$Mode  #'numeric' or 'categorical'
      if (type == "numeric"){
        pdist <- list()
        pdist[[as.character(tmp$Var1[j])]] <- mean(minor[[tmp$Var1[j]]])
      }else pdist <- table(minor[[tmp$Var1[j]]])/nrow(minor)
      w <- Lagrange(major[[tmp$Var1[j]]], pdist, w, type=type)
    }
  }
  
  set.seed(2)
  m <- ceiling(nrow(train)/2)           # the number of major/minor instances after RnR (50%/50%)
  n1 <- floor(m*0.5)                    # the number of minor instances from minor instances (25%)
  n2 <- m - n1                          # the number of minor instances from major instances (25%)
  
  rep.minor1.id <- sample(minor.id, n1, replace=TRUE)         # minor instances from minor instances
  rep.minor2.id <- sample(major.id, n2, replace=TRUE, prob=w) # minor instances from major instances
  rep.minor.id <- c(rep.minor1.id,rep.minor2.id)
  
  #####################################
  #### Calibration for major class ####
  sloc <- which(is.na(match(major.id, unique(rep.minor2.id)))==TRUE)
  reduced.major.id <- major.id[sloc]
  reduced.major <- train[reduced.major.id,]                   # major instances after relabeling
  
  w <- rep(1/nrow(reduced.major),nrow(reduced.major))
  
  for(i in 1:N){
    for(j in 1:ncol){
      type <- subset(tmp, Var1 == tmp$Var1[j])$Mode
      if (type == "numeric"){
        pdist <- list()
        pdist[[as.character(tmp$Var1[j])]] <- mean(reduced.major[[tmp$Var1[j]]])
      }else pdist <- table(major[,j])/nrow(major)
      w <- Lagrange(reduced.major[[tmp$Var1[j]]], pdist, w, type=type)
    }
  }
  
  rep.major.id <- sample(reduced.major.id,m,replace=TRUE,prob=w)
  
  #####################################
  #### Balanced Data ####
  rnr.data <- train[c(rep.minor.id,rep.major.id),]
  rnr.data$y[1:ceiling(nrow(train)/2)] <- 1   
  #table(rnr.data$y)
  
  return(rnr.data)
}

######################################################################
# Functions: for Alternative Algorithm
######################################################################
Adasyn=function(train, k=1){
  Adasyn_dt <- ADAS(subset(train, select=-c(y)), train$y, K = min(table(train$y)[2]-1, 5))[["data"]]
  names(Adasyn_dt)[length(names(Adasyn_dt))]<-"y"
  Adasyn_dt[,"y"] <- as.numeric(Adasyn_dt[,"y"])
  
  return(Adasyn_dt)
}

ROS=function(train, k=1){
  ### Split train into major and minor class, respectively. 
  # size: sample size for balanced sample    
  osample1=train[train$y==0,]      
  dat1=train[train$y==1,]   
  n0=nrow(osample1)     
  n1=nrow(dat1)   
  size=k*n0
  
  #############################################################
  #### Oversampling with size of size =n0
  rid2=sort(sample(1:n1,size,replace=TRUE))
  osample2=dat1[rid2,] 
  #############################################################
  #### Balanced Sample
  mat2=rbind(osample1, osample2)
  
  return(mat2)
}

Smote=function(train, k=1){
  Smote_dt <- SMOTE(subset(train, select=-c(y)), train$y, K = min(table(train$y)[2]-1, 5))[["data"]]
  names(Smote_dt)[length(names(Smote_dt))]<-"y"
  Smote_dt[,"y"] <- as.numeric(Smote_dt[,"y"])
  
  return(Smote_dt)
}

######################################################################
# Functions: for evaluation
######################################################################
M2 = function(test, pred){
  con.mat <- table(test$y, pred)
  nc=ncol(con.mat)
  
  # negative = majority classes
  # positive = minority class
  if(nc==2){
    TN <- con.mat[1,1]
    FP <- con.mat[1,2]
    FN <- con.mat[2,1]
    TP <- con.mat[2,2]
  }
  
  if(nc==1){    
    TN <- con.mat[1,1]
    FP <- 0
    FN <- con.mat[2,1]
    TP <- 0
  }
  
  specificity <- TN/(TN+FP)
  recall <- TP/(FN+TP)
  precision <- TP/(FP+TP)
  if (is.nan(precision)) {precision <- 0}
  
  Gmean = sqrt(recall*specificity)
  F1 = (2*precision*recall)/(precision+recall)
  if (is.nan(F1)) {F1 <- 0}
  auc<-roc.curve(test$y, pred, plotit = FALSE)$auc
  
  return(list(Gmean = Gmean, Fmeasure = F1, Auc = auc))
}

Eval <- function(train, test, by, alter, cate.col){
  
  if (by == "SVM"){
    # Since every alternative methods had already changed categorical variables into dummy variables.
    # Only RnR based methods need this procedure.
    if (alter == FALSE){
      for (col in cate.col){
        a = length(table(train[col]))
        b = length(table(test[col]))
        
        if (a>b) {large <- train; small <- test}
        if (a<b) {large <- test; small <- train}
        
        if (a!=b) {
          large[[col]] <- factor(large[[col]])
          small[[col]] <- factor(small[[col]], levels = levels(large[[col]]))
        }
        
        if (a>b) {train <- large; test <- small}
        if (a<b) {test <- large; train <- small}
      }
    }
    set.seed(1)
    fit <- svm(y~., data = train, type = "C")
    pred <- predict(fit, test)
  }
  
  if (by == "Rf"){
    rf <- randomForest(as.factor(y) ~ ., data=train)
    pred = predict(rf, newdata=test)
    #prob <- as.vector(predict(rf, test, type="prob")[,2])
  }
  
  if (by == "NB"){
    fit <- naiveBayes(y ~ ., data = train)
    pred <- predict(fit, test)
  }
  
  return(M2(test, pred))
}

######################################################################