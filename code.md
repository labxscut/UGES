# UGS-brca breast cancer intrinsic subtype classifier

Here, we complete the UGS-brca pipeline on a small sample dataset. Our starting point is based on pre-processed genomics/epigenomics/expression/clinical data that have been segmented by sample. We demonstrate how the data can be imported into R software for analysis.

## Getting ready

First we load the packages.

```{R}
library(glmnet)
library(pROC)
library(dplyr)
```

Then we import the example data, and divide it into four parts that denote combined data and three single-feature data, repectively.

```{R}
data <- read.csv("~/R/Breastcancer/newdata/cna_gene_methylation.csv")       # original combined data
standard <- read.csv("~/R/Breastcancer/newdata/standard_cna_mu_me.csv")     # standardize combined data
mutation <- data[,2:16772]
cna <- data[,c(2,16773:42366)]
methylation <- data[,c(2,42367:50833)]
```

## Lasso-Logistic regression

We evaluate three competitive hierarchical stepwise classification strategies: (B,H,LA,LB), (B,H,(LA,LB)) and ((B,H),(LA,LB)). Among these strategies, we use Lasso-Logistic regression method to construct different classifiers and compared their performance.

### (B,H,LA,LB) strategy

```{R}
s1 <- which(data$SUBTYPE == "1")
s2 <- which(data$SUBTYPE == "2")
s3 <- which(data$SUBTYPE == "3")
s4 <- which(data$SUBTYPE == "4")

# Lasso-Logistic regression
lasso_four <- function(x){        
  train_data <- x[train,]
  test_data <- x[test,]
  x.train <- as.matrix(sapply(data.frame(train_data[,c(2:ncol(train_data))]),as.numeric))
  y.train <- as.matrix(sapply(data.frame(train_data[,1]),as.numeric))
  x.test <- as.matrix(sapply(data.frame(test_data[,c(2:ncol(test_data))]),as.numeric))
  y.test <- as.matrix(sapply(data.frame(test_data[,1]),as.numeric))
  
  #glmnet()
  fit <- glmnet(x.train, y.train, family="multinomial", alpha=1)
  
  #cv.glmnet()
  set.seed(1)
  cv.fit <- cv.glmnet(x.train, y.train, family="multinomial", alpha=1, 
                      type.measure = "class", nfolds = 10)
  bestlam <- cv.fit$lambda.min
  
  #training result
  train_pred <- predict(fit, s = bestlam, newx = x.train, type="class")
  train_pred_prob <- predict(cv.fit, newx=x.train , s=bestlam, type="response")
  train_t <- table(y.train, train_pred, dnn=c("True","Prediction"))
  train_acc <- rep(0,4)
  if (dim(train_t)[2] == 4){
    for (i in 1:4) {
      train_acc[i] <- train_t[i,i]/sum(train_t[i,])
    }
  } else{
    for (i in 1:3) {
      train_acc[i] <- train_t[i,i]/sum(train_t[i,])
    }
  }
  train_auc <- auc(multiclass.roc(train_data$SUBTYPE,as.numeric(train_pred)))
  
  #testing result
  test_pred <- predict(fit, s = bestlam, newx = x.test, type="class")
  test_pred_prob <- predict(cv.fit, newx=x.test , s=bestlam, type="response")
  test_t <- table(y.test, test_pred, dnn=c("True","Prediction"))
  test_acc <- rep(0,4)
  if (dim(test_t)[2] == 4){
    for (i in 1:4) {
      test_acc[i] <- test_t[i,i]/sum(test_t[i,])
    }
  } else{
    for (i in 1:3) {
      test_acc[i] <- test_t[i,i]/sum(test_t[i,])
    }
  }
  test_auc <- auc(multiclass.roc(test_data$SUBTYPE,as.numeric(test_pred)))
  
  result <- list(fit = fit, bestlam = bestlam, train_pred = train_pred, 
                 train_pred_prob = train_pred_prob, train_t = train_t, 
                 train_acc =train_acc, train_auc = train_auc, 
                 test_pred = test_pred, test_pred_prob = test_pred_prob,
                 test_t = test_t, test_acc = test_acc, test_auc = test_auc)
  return(result)
}

# divide the training and testing set, ensuring sample balance
set.seed(900)
train <- c(sample(s1,200),sample(s2,200),sample(s3,200),sample(s4,200))   #sample balance
test <- c(1:nrow(data))[-train]

# generate ROC curves of the (B,H,LA,LB) strategy
IV_OVR_train <- function(x){
  train_data <- x[train,]
  test_data <- x[test,]
  IV <- lasso_four(x)
  IV$train_pred_prob <- as.matrix(data.frame(IV$train_pred_prob))
  colnames(IV$train_pred_prob) <- c(1,2,3,4)
  IV$train_label <- matrix(0,dim(train_data)[1],4)
  IV$train_label[,1:4] <- train_data$SUBTYPE
  IV$train_label[,1][which(IV$train_label[,1] != 1)] <- 0
  IV$train_label[,2][which(IV$train_label[,2] != 2)] <- 0
  IV$train_label[,2][which(IV$train_label[,2] == 2)] <- 1
  IV$train_label[,3][which(IV$train_label[,3] != 3)] <- 0
  IV$train_label[,3][which(IV$train_label[,3] == 3)] <- 1
  IV$train_label[,4][which(IV$train_label[,4] != 4)] <- 0
  IV$train_label[,4][which(IV$train_label[,4] == 4)] <- 1
  
  train_roc1 <- roc(IV$train_label[,1], IV$train_pred_prob[,1], direction = "<", )
  train_roc2 <- roc(IV$train_label[,2], IV$train_pred_prob[,2], direction = "<")
  train_roc3 <- roc(IV$train_label[,3], IV$train_pred_prob[,3], direction = "<")
  train_roc4 <- roc(IV$train_label[,4], IV$train_pred_prob[,4], direction = "<")
  train_roc <- train_roc1
  train_roc$sensitivities <- (train_roc1$sensitivities + train_roc2$sensitivities + train_roc3$sensitivities + train_roc4$sensitivities)/4
  train_roc$specificities <- (train_roc1$specificities + train_roc2$specificities + train_roc3$specificities + train_roc4$specificities)/4
  train_roc$auc <- (train_roc1$auc + train_roc2$auc + train_roc3$auc + train_roc4$auc)/4
  roc <- plot(train_roc1,col='1')          #select 0 as control, 1 as case
  plot(train_roc2,col='2',add = TRUE)
  plot(train_roc3,col='3',add = TRUE)
  plot(train_roc4,col='4',add = TRUE)
  plot(train_roc,col='6',add = TRUE)
  legend('bottomright',
         legend = c(paste("Basal: ", round(auc(train_roc1),3), sep =""), 
                    paste("Her2: ", round(auc(train_roc2),3), sep =""), 
                    paste("LumA: ", round(auc(train_roc3),3), sep =""), 
                    paste("LumB: ", round(auc(train_roc4),3), sep =""),
                    paste("Overall: ", round(train_roc$auc,3), sep ="")),
         col = c("1", "2", "3", "4","6"), lty = 1)
  
  return(roc)
}

IV_OVR_test <- function(x){
  train_data <- x[train,]
  test_data <- x[test,]
  IV <- lasso_four(x)
  IV$test_pred_prob <- as.matrix(data.frame(IV$test_pred_prob))
  colnames(IV$test_pred_prob) <- c(1,2,3,4)
  IV$test_label <- matrix(0,dim(test_data)[1],4)
  IV$test_label[,1:4] <- test_data$SUBTYPE
  IV$test_label[,1][which(IV$test_label[,1] != 1)] <- 0
  IV$test_label[,2][which(IV$test_label[,2] != 2)] <- 0
  IV$test_label[,2][which(IV$test_label[,2] == 2)] <- 1
  IV$test_label[,3][which(IV$test_label[,3] != 3)] <- 0
  IV$test_label[,3][which(IV$test_label[,3] == 3)] <- 1
  IV$test_label[,4][which(IV$test_label[,4] != 4)] <- 0
  IV$test_label[,4][which(IV$test_label[,4] == 4)] <- 1
  
  test_roc1 <- roc(IV$test_label[,1], IV$test_pred_prob[,1], direction = "<", )
  test_roc2 <- roc(IV$test_label[,2], IV$test_pred_prob[,2], direction = "<")
  test_roc3 <- roc(IV$test_label[,3], IV$test_pred_prob[,3], direction = "<")
  test_roc4 <- roc(IV$test_label[,4], IV$test_pred_prob[,4], direction = "<")
  test_roc <- test_roc1
  test_roc$sensitivities <- (test_roc1$sensitivities + test_roc2$sensitivities + test_roc3$sensitivities + test_roc4$sensitivities)/4
  test_roc$specificities <- (test_roc1$specificities + test_roc2$specificities + test_roc3$specificities + test_roc4$specificities)/4
  test_roc$auc <- (test_roc1$auc + test_roc2$auc + test_roc3$auc + test_roc4$auc)/4
  roc <- plot(test_roc1,col='1')          #select 0 as control, 1 as case
  plot(test_roc2,col='2',add = TRUE)
  plot(test_roc3,col='3',add = TRUE)
  plot(test_roc4,col='4',add = TRUE)
  plot(test_roc,col='6',add = TRUE)
  legend('bottomright',
         legend = c(paste("Basal: ", round(auc(test_roc1),3), sep =""), 
                    paste("Her2: ", round(auc(test_roc2),3), sep =""), 
                    paste("LumA: ", round(auc(test_roc3),3), sep =""), 
                    paste("LumB: ", round(auc(test_roc4),3), sep =""),
                    paste("Overall: ", round(test_roc$auc,3), sep ="")),
         col = c("1", "2", "3", "4","6"), lty = 1)
  
  return(roc)
}

par(mfrow = c(2,4))
IV_OVR_train(mutation)
IV_OVR_train(cna)
IV_OVR_train(methylation)
IV_OVR_train(standard[,-1])
IV_OVR_test(mutation)
IV_OVR_test(cna)
IV_OVR_test(methylation)
IV_OVR_test(standard[,-1])

IV_OVR <- function(x){
  train_data <- x[train,]
  test_data <- x[test,]
  IV <- lasso_four(x)
  
  data <- rbind(train_data,test_data)
  IV$pred_prob <- as.matrix(rbind(data.frame(IV$train_pred_prob),data.frame(IV$test_pred_prob)))
  colnames(IV$pred_prob) <- c(1,2,3,4)
  IV$label <- matrix(0,dim(data)[1],4)
  IV$label[,1:4] <- data$SUBTYPE
  IV$label[,1][which(IV$label[,1] != 1)] <- 0
  IV$label[,2][which(IV$label[,2] != 2)] <- 0
  IV$label[,2][which(IV$label[,2] == 2)] <- 1
  IV$label[,3][which(IV$label[,3] != 3)] <- 0
  IV$label[,3][which(IV$label[,3] == 3)] <- 1
  IV$label[,4][which(IV$label[,4] != 4)] <- 0
  IV$label[,4][which(IV$label[,4] == 4)] <- 1
  
  roc1 <- roc(IV$label[,1], IV$pred_prob[,1], direction = "<")
  roc2 <- roc(IV$label[,2], IV$pred_prob[,2], direction = "<")
  roc3 <- roc(IV$label[,3], IV$pred_prob[,3], direction = "<")
  roc4 <- roc(IV$label[,4], IV$pred_prob[,4], direction = "<")
  averroc <- roc1
  averroc$sensitivities <- (roc1$sensitivities + roc2$sensitivities + roc3$sensitivities + roc4$sensitivities)/4
  averroc$specificities <- (roc1$specificities + roc2$specificities + roc3$specificities + roc4$specificities)/4
  averroc$auc <- (roc1$auc + roc2$auc + roc3$auc + roc4$auc)/4
  roc <- plot(roc1,col='1')          #select 0 as control, 1 as case
  plot(roc2,col='2',add = TRUE)
  plot(roc3,col='3',add = TRUE)
  plot(roc4,col='4',add = TRUE)
  plot(averroc,col='6',add = TRUE)
  legend('bottomright',
         legend = c(paste("Basal: ", round(auc(roc1),3), sep =""), 
                    paste("Her2: ", round(auc(roc2),3), sep =""), 
                    paste("LumA: ", round(auc(roc3),3), sep =""), 
                    paste("LumB: ", round(auc(roc4),3), sep =""),
                    paste("Overall: ", round(averroc$auc,3), sep ="")),
         col = c("1", "2", "3", "4","6"), lty = 1)
  
  return(roc)
}

par(mfrow = c(2,4))
IV_OVR(mutation)
IV_OVR(cna)
IV_OVR(methylation)
IV_OVR(standard[,-1])
```

### (B,H,(LA,LB)) strategy

