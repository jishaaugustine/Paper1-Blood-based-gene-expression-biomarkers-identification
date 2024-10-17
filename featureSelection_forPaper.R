library(parallel)
library(doParallel)
library(data.table)
setwd("~/Desktop/imple/microarray_integration_new")
load("CombindZ_4datafull.RData")
load("GSE72267_validation.RData")
data1 <- fread(file="existing_signature.csv",header=T)
b = data1[,data1$Falchetti59]
#load("GSE63060.RData")
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
########Resampling into two group without replacement####
combat_edataF=as.data.frame(combat_edata)
combat_edataF$label=phenoDta$label

library(ROCR)
library(caret)
library(FSelector)
library(MASS)
library(FSinR)

feature_varImp=function(model,trData,label,threshold)
{
  set.seed(123,sample.kind = "default")
  model <- train(x=trData,y=label, method=model)
  gbmImp <- varImp(model, scale = FALSE)
  imp=as.data.frame(gbmImp$importance)
  subset=rownames(imp)[order(imp, decreasing=TRUE)][1:threshold]
  return(subset)
}

library(MLeval)
predictionModel=function(model,trData,tstData,tstLabel)
{
  set.seed(123,sample.kind = "default")
  #Classification with 5 fold cross validation
  trainControl <- trainControl(method="cv",number=10 ,savePredictions = TRUE, classProbs = TRUE,summaryFunction = twoClassSummary,allowParallel = TRUE)
  
  metric <- "ROC"
  
  #Fitting Model
  fit.model <- train(label~., data=trData, method=model, metric=metric,
                     trControl=trainControl)

  x <- evalm(fit.model,showplots = FALSE)
  x1=x$optres$`Group 1`["AUC-ROC",1]

  #Prediction with validation set
  predictions <-predict(fit.model, val,type="prob")
  #Calculate AUC for independent validation set
  pred=prediction(predictions$X1,label)
  perf_AUC=performance(pred,"auc") #Calculate the AUC value
  x2=perf_AUC@y.values[[1]]
  fit.model=NULL

  return(list(x1,x2))
  
}


set.seed(123,sample.kind = "default")
##############################################################################
result=list()
summary_auc=data.frame(row.names=1)
sig_genes=list()#List of significant genes from feature selection
sig_genes1=list()
sig_genes2=list()
sig_genes3=list()
sig_genes4=list()
dataset=list()
datasetT=list()
val=list()
auc_gbm=c()
#####################################################################
set.seed(5,sample.kind = "default")
trainIndex <- createDataPartition(combat_edataF$label, p = .8,list=FALSE)
# # select 10% of the data for validation
validation <- combat_edataF[-trainIndex,]

# # use the remaining 90% of data to training and testing the models
train <- combat_edataF[trainIndex,]

#########################################################################
train=combat_edataF
train$label=combat_edataF$label
validation=as.data.frame(validationset)
validation$label=validationset_label
##############################################################
library(gbm)
set.seed(123,sample.kind = "default")
sig1=list()
thresh=ncol(train)-1
sig_genes1=feature_varImp("gbm",train[,-ncol(train)],as.factor(train[,ncol(train)]),thresh)
########################################################
train_new=train[,colnames(train) %in% sig_genes1[1:60]]
train_new$label=train$label

#Get validation data with selected features
test=validation[,colnames(validation) %in% sig_genes1[1:60]]
test$label=validation[,ncol(validation)]
####################################################################################
x=train_new[,-ncol(train_new)]
y=as.factor(make.names(train_new$label))
# construct rfeControl object
rfe_control = rfeControl(functions = caretFuncs, #caretFuncs here
                         method="cv",
                         number=10,allowParallel=TRUE)

# construct trainControl object for your train method 
fit_control = trainControl(classProbs=T,
                           search="random")

# get results
rfe_fit = rfe(x,y,
              sizes = 1:60,
              rfeControl = rfe_control,
              method="glm",
              trControl=fit_control)
sig_genes2=rfe_fit$optVariables

#######################################################################

ctrl <- gafsControl(functions = caretGA,allowParallel=TRUE)
obj <- gafs(x = train_new[,-ncol(train_new)], 
            y = train_new$label,
            iters = 100,
            gafsControl = ctrl,
            method = "glm")
sig_genes3 <- obj$ga$final # Get features selected by GA
####################################################################################
#Stepwise logistic regression
set.seed(123,sample.kind = "default")
library(MASS)
library(dbplyr)
library(plyr)
# Fit the model
model <- glm(label ~., data = train_new, family = binomial) %>%
  stepAIC(trace = FALSE, direction="both")
coeff=as.data.frame(model$coefficients)
sig_genes4=row.names(as.data.frame(model$coefficients))[2:nrow(coeff)]
###############################################################
set.seed(123,sample.kind = "default")
library(MASS)
step.model <- train(label ~., data = train_new,
                    method = "glmStepAIC", 
                    trace = FALSE
)
# Model accuracy
step.model$results
# Final model coefficients
int=names(step.model$finalModel$coefficients)[2:nrow(coeff)]




#################################################################
#intersect
b=intersect(sig_genes2,intersect(sig_genes3,sig_genes4))
#############################################################
train1=train[,colnames(train) %in% b]
train1$label=train$label
library(DMwR)
library(ROSE)
set.seed(123,sample.kind = "default")
train_rose <- ovun.sample(label ~ ., data = train1, N=nrow(train1),  p=0.5, seed=123,method="both")$data
train_smote=SMOTE(label~ ., data=train1,perc.over = 100,perc.under = 200)
train=train_rose
#####################################################################
for (th in seq(10,150,10))
{
 
  sig_genes=sig_genes1[1:th]
  print(length(sig_genes))
  dataset=train[,colnames(train) %in% sig_genes]
  datasetT=as.data.frame(dataset)
  datasetT$label=make.names(train[,ncol(train)])
  
  #Get validation data with selected features
  val=validation[,colnames(validation) %in% sig_genes]
  
  
  label=make.names(validation[,ncol(validation)])
  
  auc_gbm=predictionModel("gbm",datasetT, val,label)
  
  auc_model<-as.data.frame(x=c(auc_gbm[[1]]))
  auc_validation<-as.data.frame(x=c(auc_gbm[[2]]))

  names(auc_model)=th
  names(auc_validation)=th
  auc_list=list(auc_model,auc_validation)
  summary_auc=cbind(summary_auc,auc_list)
}
l=list(sig_genes,summary_auc)
save(l,file="GBM10_150.RData")
##########################################################################
for(i in list(123,111,500,100))
{
  i=500

  library(caret)
  set.seed(i,sample.kind = "default")
  trainIndex <- createDataPartition(combat_edataF$label, p = .8,list=FALSE)
  # # select 10% of the data for validation
  validation <- combat_edataF[-trainIndex,]

  # # use the remaining 90% of data to training and testing the models
  train <- combat_edataF[trainIndex,]
  
  train1=train[,colnames(train) %in% b]
  train1$label=train$label
  colnames(train1)[24]="HLA"
  library(DMwR)
  library(ROSE)
  set.seed(i,sample.kind = "default")
  train_rose <- ovun.sample(label ~ ., data = train1, N=nrow(train1),  p=0.5, seed=500,method="both")$data
  train_smote=SMOTE(label~ ., data=train1,perc.over = 100,perc.under = 200)
  train=train_rose
  colnames(train1)[24]="HLA-A"
  # train=combat_edataF
  # train$label=combat_edataF$label
  # validation=as.data.frame(validationset)
  # validation$label=validationset_label  
  # 
dataset=train[,colnames(train) %in% b]
datasetT=as.data.frame(dataset)
datasetT$label=make.names(train[,ncol(train)])

#Get validation data with selected features
val=validation[,colnames(validation) %in% b]
label=make.names(validation[,ncol(validation)])
#######################################################
library(caretEnsemble)
library(gbm)
library("caTools")
trainControl <- trainControl(method="repeatedcv",repeats=3, number=10 ,savePredictions = TRUE, classProbs = TRUE,summaryFunction = twoClassSummary,allowParallel = TRUE)
metric <- "ROC"
#algorithmList <- c('xgbTree','knn','nb','rf','rpart','svmLinear','svmRadial','gbm')
#algorithmList <- c('rf','xgbTree','svmLinear','svmRadial','gbm')
algorithmList <- c('svmLinear','svmRadial')
set.seed(i,sample.kind = "default")
models <- caretList(label~., data=datasetT, trControl=trainControl, methodList=algorithmList)
model_preds <- lapply(models, predict, newdata=val, type="prob")
model_preds <- lapply(model_preds, function(x) x[,"X1"])
model_preds <- data.frame(model_preds)
print(i)
print(caTools::colAUC(model_preds, label,plotROC=FALSE,alg="ROC"))
}

#Get specificity and sensitivity
label=validation[,ncol(validation)]
threshold=0.575
pred1rpart <- predict(models$svmRadial, val,type = "prob")
pred1rpart1<-ifelse(pred1rpart$X1 > threshold,1,0)
confusionMatrix (as.factor(pred1rpart1), as.factor(label))


#Plot ROC
roc_pred <- prediction(predictions = model_preds$svmRadial  , labels = validationset_label)
roc_perf <- performance(roc_pred , "tpr" , "fpr")
plot(roc_perf)
####################################################################################
#Deep Neural network
library(keras)
library(recipes)
library(keras)
use_session_with_seed(123)
X_train2=train[,colnames(train) %in% b]
X_train2$label=train$label
X_test2=as.data.frame(validationset[,colnames(validationset) %in% b])
X_test2$label=validationset_label

rec_obj <- 
  recipe (label ~ ., data = X_train2) %>%
  step_center (all_predictors(), -all_outcomes()) %>%
  step_scale (all_predictors(), -all_outcomes()) %>%
prep (data = X_train2)
X_train=as.matrix(X_train2[,-ncol(X_train2)])

rec_obj1 <- 
  recipe (label ~ ., data = X_test2) %>%
  step_center (all_predictors(), -all_outcomes()) %>%
  step_scale (all_predictors(), -all_outcomes()) %>%
prep (data = X_test2)
X_test=as.matrix(X_test2[,-ncol(X_test2)])

y_train_vec <- ifelse (pull (X_train2, label) == "1", 1, 0)
y_train=y_train_vec

y_test_vec  <- ifelse (pull (X_test2, label) == "1", 1, 0)
y_test=y_test_vec

#X_train1=train[,colnames(train) %in% b]
#X_train <- X_train1 %>% 
#  scale()
#y_train=to_categorical(train$label)

#X_test1=validationset[,colnames(validationset) %in% b]
#X_test <- X_test1 %>% 
#  scale()
#y_test=to_categorical(validationset_label)


model <- keras_model_sequential()

# Add layers to the model
model %>% 
  layer_dense(units = 25, activation = 'relu',kernel_regularizer = regularizer_l2(0.01),kernel_initializer = 'glorot_normal',input_shape = ncol(X_train)) %>% 
  layer_dropout(0.5) %>%
  layer_dense(units = 10,kernel_regularizer = regularizer_l2(0.01),activation = 'relu') %>%
  layer_dropout(0.5) %>%
  layer_dense(units = 1, activation = 'sigmoid')

library(tensorflow)
library(tidyverse)
model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = 'Nadam',
  metrics = list(tf$keras$metrics$AUC(),tf$keras$metrics$TrueNegatives(),tf$keras$metrics$TruePositives(),tf$keras$metrics$FalseNegatives(),tf$keras$metrics$FalsePositives())

  #metrics = "accuracy"
)

es = callback_early_stopping(monitor='val_loss', mode='auto', verbose=1)

model %>% fit(
  X_train, 
  y_train, 
  epochs = 100, 
  batch_size = 32, 
  shuffle=TRUE,
  callbacks = es,
  validation_data = list(X_test,y_test)
  #validation_split = 0.2
)

model %>% save_model_hdf5("my_model.h5")

new_model <- load_model_hdf5("my_model.h5")

new_model %>% evaluate(X_test, y_test, verbose = 1)

pred_class <- predict_classes(object = model,
                              x = X_test) %>% as.vector()
pred_prob <- predict_proba(object = model,
                           x = X_test) %>% as.vector()

library(tibble)
library(tidyverse)
library(yardstick)
library(pROC)
library(ROCR)
predict_value <- tibble(
  truth = as.factor(validationset_label) %>% fct_recode(Yes = "1", No = "0"),
  estimate = as.factor(pred_class) %>% fct_recode(Yes = "1", No = "0"),
  pred_prob = pred_prob[1:59]
)

print(predict_value)
predict_value %>% metrics(truth, estimate)
estimate_tbl2 %>% roc_auc(truth, class_prob)

#######################################
#cross-validation in DNN
set.seed(123,sample.kind = "default")
library(keras)
library(recipes)
library(tensorflow)
library(tidyverse)
library(caret)
use_session_with_seed(123)
trainData=train[,colnames(train) %in% b]
trainData$label=train$label
testData=as.data.frame(validationset[,colnames(validationset) %in% b])
testData$label=validationset_label

rec_obj1 <- 
  recipe (label ~ ., data = testData) %>%
  step_center (all_predictors(), -all_outcomes()) %>%
  step_scale (all_predictors(), -all_outcomes()) %>%
  prep (data = testData)
X_test=as.matrix(testData[,-ncol(testData)])


y_test_vec  <- ifelse (pull (testData, label) == "1", 1, 0)
y_test=y_test_vec

# define 10-fold cross validation test harness
kfold <- createFolds(trainData$label, k = 10)

# running through each fold of the cross-validation
for (fold in kfold){
  
  X_train2=trainData[-fold,]
  X_val2=trainData[fold,]
  rec_obj <- 
    recipe (label ~ ., data = X_train2) %>%
    step_center (all_predictors(), -all_outcomes()) %>%
    step_scale (all_predictors(), -all_outcomes()) %>%
    prep (data = X_train2)
  X_train=as.matrix(X_train2[,-ncol(X_train2)])
  
  y_train_vec <- ifelse (pull (X_train2, label) == "1", 1, 0)
  y_train=y_train_vec
  
  rec_obj <- 
    recipe (label ~ ., data = X_val2) %>%
    step_center (all_predictors(), -all_outcomes()) %>%
    step_scale (all_predictors(), -all_outcomes()) %>%
    prep (data = X_val2)
  X_val=as.matrix(X_val2[,-ncol(X_val2)])
  
  y_val_vec <- ifelse (pull (X_val2, label) == "1", 1, 0)
  y_val=y_val_vec
  
  
  model <- keras_model_sequential()
  
  # Add layers to the model
  model %>% 
    layer_dense(units = 25, activation = 'relu',kernel_regularizer = regularizer_l2(0.01),kernel_initializer = 'glorot_normal',input_shape = ncol(X_train)) %>% 
    layer_dropout(0.5) %>%
    layer_dense(units = 10,kernel_regularizer = regularizer_l2(0.01),activation = 'relu') %>%
    layer_dropout(0.5) %>%
    layer_dense(units = 1, activation = 'sigmoid')
  
  model %>% compile(
    loss = 'binary_crossentropy',
    optimizer = 'Nadam',
    metrics = list(tf$keras$metrics$AUC(),tf$keras$metrics$TrueNegatives(),tf$keras$metrics$TruePositives(),tf$keras$metrics$FalseNegatives(),tf$keras$metrics$FalsePositives())
    #metrics = tf$keras$metrics$AUC()
    #metrics = "accuracy"
  )
  
  es = callback_early_stopping(monitor='val_loss', mode='auto', verbose=1)
  
  model %>% fit(
    X_train, 
    y_train, 
    epochs = 100, 
    batch_size = 32, 
    shuffle=TRUE,
    callbacks = es,
    validation_data = list(X_val,y_val),
    verbose = 0
  )
  
  model %>% evaluate(X_val, y_val, verbose = 1)
  # evaluating the performance of the model
  model %>% evaluate(X_test, y_test, verbose = 1)
  # model %>% predict_classes(x_test)
}

