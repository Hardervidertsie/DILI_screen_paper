# testing consistency
output.confus = alist()
output_testComps = alist()

output_tune = alist()

for( j in 1:10){
  caseInd <- createDataPartition(all.treatsnoDMSO_dili$class,p=.8,list=FALSE)
  
 
  trainData_long<- as.data.frame(ML_data_final_data_long[ treatment %in% all.treatsnoDMSO_dili[caseInd, treatment] ])
  testData_long <- as.data.frame(ML_data_final_data_long[ treatment %in% all.treatsnoDMSO_dili[-caseInd, treatment] ])
  
  trainData_long <- as.data.table(trainData_long)
  
  results.ks = alist()
  for( i in seq_along(all.ks.feats)){
    results.ks[[i]] <- ks.test(x = trainData_long[ variable == all.ks.feats[i] & class == "nonSevere" , value ],
                               y = trainData_long[ variable == all.ks.feats[i] & class == "severe" , value ],
                               alternative = "two.sided", exact = FALSE)
    
  }
  
  ks_results <- data.table(data.frame(feature = all.ks.feats, 
                                      statistic = sapply(results.ks, "[[", "statistic"),
                                      p.value = sapply(results.ks, "[[", "p.value")))
  
  
  ks_results <- ks_results[ order(ks_results$p.value), ]
  sign_ks_results <- ks_results[ p.value <= 0.05, ]
  
  
  ML_data_final_data$treatment[!unique(ML_data_final_data$treatment) %in% unique(BMC_data_add$treatment)] # dmso
  
  
  ind <- colnames(ML_data_final_data) %in% sign_ks_results$feature
  # sel_data is training data
  sel_data <- as.data.frame(
    ML_data_final_data[ treatment %in% all.treatsnoDMSO_dili[caseInd, treatment], ind, with = FALSE]
  )
 cor_data <- cor(sel_data)
 
  
 
  df2 = cor(sel_data)
  
  hc = findCorrelation(abs(df2), cutoff=0.8) # putt any value as a "cutoff" 
  hc = sort(hc)
  reduced_Data = sel_data[,-c(hc)]
  sel_data <- reduced_Data
 
 
  selFeats <- colnames(sel_data)
 
  indKeep <-  colnames(ML_data_final_data) %in% selFeats
  sum(indKeep)
  
  indKeep <- indKeep | colnames(ML_data_final_data) %in% c("diliClass")
  
  ML_data_final_data[, diliClass := ML_data_final_annot$diliClass ]
  ML_data_final_data[, diliClass := factor(diliClass) ]
  
  trainData <- as.data.frame(
    ML_data_final_data[ treatment %in% all.treatsnoDMSO_dili[caseInd, treatment] , indKeep, with = FALSE]
  )
  testData <- as.data.frame(
    ML_data_final_data[ treatment %in% all.treatsnoDMSO_dili[-caseInd, treatment] , indKeep, with = FALSE]
  )
  
  trainData$diliClass <- factor(trainData$diliClass, levels =c("severe","nonSevere"))
  testData$diliClass <- factor(testData$diliClass, levels =c("severe","nonSevere"))
  dim(testData)
  
  # feature selection on 26 features using ga selection with svm
  
  trainX <-trainData[, !colnames(trainData) %in% 'diliClass'] # Create training feature data frame
  testX <- testData[,!colnames(trainData) %in% 'diliClass'] # Create test feature data frame 
  y=trainData$diliClass # Target variable for training
  
  
  trainX2 <- as.data.frame(trainX) # training data: selected features
  testX2 <- as.data.frame(testX) # test data: selected features
  
  
  
  ## SUPPORT VECTOR MACHINE MODEL
  dim(testX2)
  #Note the default method of picking the best model is accuracy and Cohen's Kappa 
  # Set up training control
  ctrl <- trainControl(method="repeatedcv", # 10fold cross validation
                       number = 10,
                       repeats=10, # do 5 repititions of cv
                       summaryFunction=twoClassSummary, # Use AUC to pick the best model
                       classProbs=TRUE)
  
  #Use the expand.grid to specify the search space 
  #Note that the default search grid selects 3 values of each tuning parameter
  
  grid <- expand.grid(interaction.depth = c(1,2,3), #tree depths from 1 to 4
                      n.trees=seq(10,100,by=10), # let iterations go from 10 to 100
                      shrinkage=c(0.01,0.1), # Try 2 values fornlearning rate 
                      n.minobsinnode = 20)
  # 
  # Set up for parallel processing
  #set.seed(1951)
  
  registerDoParallel(7,cores=8)
  
  #support vector machine or random forest that deals well with a large number of predictors.
  #Train and Tune the SVM
  psvmTuneGrid <- expand.grid(C=seq(0.005,0.25, length.out=4), sigma=seq(0.001, 0.1, length.out=8))
  svm.tune <- train(x=trainX2,
                    y= trainData$diliClass,
                    method = "svmRadial",
                    degree = 2,
                    tuneLength = 3, # 9 values of the cost function
                    preProc = c("center","scale"),
                    metric="ROC",
                    trControl=ctrl
                    #tuneGrid=psvmTuneGrid
  ) 
  
  
  #closeAllConnections()
 
  #svm.tune$results[rev(order(svm.tune$results$ROC)), ]
  #Finally, assess the performance of the model using the test data set.
  
  #Make predictions on the test data with the SVM Model
  
  svm.pred <- predict(svm.tune, newdata = testX2, type = "raw")
  
  
  output_tune[[j]] <- svm.tune
  output_testComps[[j]] <- all.treatsnoDMSO_dili[-caseInd, treatment]
output.confus[[j]] <-   confusionMatrix(svm.pred,testData$diliClass)
print(j)
closeAllConnections()
}
output_testComps[[2]]
output.confus[[1]]$overall[1]
mean(sapply(lapply(output.confus, '[[', 'overall'),'[[',1))
max(sapply(lapply(output.confus, '[[', 'overall'),'[[',1))
min(sapply(lapply(output.confus, '[[', 'overall'),'[[',1))
myresult = alist()
for( i in seq_along(output_tune)){
myresult[[i]] <- (output_tune[[i]]$results[which(max(output_tune[[i]]$results$ROC) == output_tune[[i]]$results$ROC),][, 3:5])
}

colMeans(do.call('rbind', myresult))

