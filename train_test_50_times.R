# testing consistency
output.confus = alist()
output_testComps = alist()
output_feats = alist()
output_tune = alist()
output_test = alist()
output_predict = alist()
for( j in 1:200){
  caseInd <- createDataPartition(all.treatsnoDMSO_dili$class,p=0.8,list=FALSE)
  
 
  trainData_long<- as.data.frame(ML_data_final_data_long[ treatment %in% all.treatsnoDMSO_dili[caseInd, treatment] ])
  testData_long <- as.data.frame(ML_data_final_data_long[ treatment %in% all.treatsnoDMSO_dili[-caseInd, treatment] ])
  
  trainData_long <- as.data.table(trainData_long)
  
  # results.ks = alist()
  # for( i in seq_along(all.ks.feats)){
  #   results.ks[[i]] <- ks.test(x = trainData_long[ variable == all.ks.feats[i] & class == "nonSevere" , value ],
  #                              y = trainData_long[ variable == all.ks.feats[i] & class == "severe" , value ],
  #                              alternative = "two.sided", exact = FALSE)
  #   
  # }
  # 
  # ks_results <- data.table(data.frame(feature = all.ks.feats, 
  #                                     statistic = sapply(results.ks, "[[", "statistic"),
  #                                     p.value = sapply(results.ks, "[[", "p.value")))
  # 
  # 
  # ks_results <- ks_results[ order(ks_results$p.value), ]
  # sign_ks_results <- ks_results[ p.value <= 0.05, ]
  # 
  # 
  # ML_data_final_data$treatment[!unique(ML_data_final_data$treatment) %in% unique(BMC_data_add$treatment)] # dmso
  # 
  # 
  #sel_feats_200run<-selFeats filtered on correlation using all data
  ind <- colnames(ML_data_final_data) %in% sel_feats_200run
  # sel_data is training data
  sel_data <- as.data.frame(
    ML_data_final_data[ treatment %in% all.treatsnoDMSO_dili[caseInd, treatment], ind, with = FALSE]
  )
 cor_data <- cor(sel_data)
 
  
 
  df2 = cor(sel_data)
  
 #  hc = findCorrelation(abs(df2), cutoff=0.8) # putt any value as a "cutoff" 
 #  hc = sort(hc)
 #  reduced_Data = sel_data[,-c(hc)]
 #  sel_data <- reduced_Data
 # dim(sel_data)
 # 
 #  selFeats <- colnames(sel_data)
  selFeats <- sel_feats_200run
  #write.table(sel_feats_200run, file ="../tmp/sel_feats_200run.txt", sep ="\t", col.names = NA)
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
  psvmTuneGrid <- expand.grid(C=seq(0.005,0.25, length.out=3), sigma=seq(0.001, 0.1, length.out=8))
  svm.tune <- train(x=trainX2,
                    y= trainData$diliClass,
                    method = "svmRadial",
                    degree = 2,
                    tuneLength = 3, # 9 values of the cost function
                    preProc = c("center","scale"),
                    metric="ROC",
                    trControl=ctrl,
                    tuneGrid=psvmTuneGrid
  ) 
  
  
  #closeAllConnections()
 
  #svm.tune$results[rev(order(svm.tune$results$ROC)), ]
  #Finally, assess the performance of the model using the test data set.
  
  #Make predictions on the test data with the SVM Model
  
  svm.pred <- predict(svm.tune, newdata = testX2, type = "raw")
  output_predict[[j]] <- svm.pred
  output_test[[j]] <- testData
  output_tune[[j]] <- svm.tune
  output_testComps[[j]] <- all.treatsnoDMSO_dili[-caseInd, treatment]
output.confus[[j]] <-   confusionMatrix(svm.pred,testData$diliClass)
#output_feats[[j]] <- sign_ks_results$feature
output_feats[[j]] <- selFeats
print(j)
closeAllConnections()
}

compound_prediction <- data.frame(unlist(output_testComps), 
           unlist(output_predict), 
           unlist(lapply(output_test, "[[", 'diliClass')))# to count how often compounds are correctly tested
colnames(compound_prediction) <- c("treatment", "predicted", "diliClass")
compound_prediction$count <- 1
compound_prediction$correct <- compound_prediction$predicted == compound_prediction$diliClass


count_comps_correct <- aggregate(data = compound_prediction, .~treatment, sum)
count_comps_correct$frac_correct <- count_comps_correct$correct / count_comps_correct$count
# check out what dili class they are
count_comps_correct$diliClass <- NULL


count_comps_correct <- merge(count_comps_correct, annotationData, by = 'treatment', all.x = TRUE)
colnames(count_comps_correct)
count_comps_correct<- count_comps_correct[, c("treatment", "count", "frac_correct", 
                                              "SeverityLabel", "DILIConcern", "dissolveInfo")]
write.table(count_comps_correct,file = '../generated/results/SVM/compounds_pred_correct.txt', sep ="\t", col.names=NA)
head(count_comps_correct)
names(output.confus[[1]])
output.confus[[1]]$byClass[1]

median(sapply(lapply(output.confus, '[[', 'overall'),'[[','Accuracy')) # 0.6818182
median(sapply(lapply(output.confus, '[[', 'byClass'),'[[','Sensitivity')) # 0.6
median(sapply(lapply(output.confus, '[[', 'byClass'),'[[','Specificity')) # 0.75

median(sapply(lapply(output.confus, '[[', 'byClass'),'[[','Balanced Accuracy')) #  0.67

max(sapply(lapply(output.confus, '[[', 'overall'),'[[',1))
min(sapply(lapply(output.confus, '[[', 'overall'),'[[',1))
myresult = alist()
for( i in seq_along(output_tune)){
myresult[[i]] <- (output_tune[[i]]$results[which(max(output_tune[[i]]$results$ROC) == output_tune[[i]]$results$ROC),][, 3:5])
}

colMeans(do.call('rbind', myresult))

output_tune[[1]]["trainingData"]
svm.probs = alist()
svm.ROC = alist()
for( i in seq_along(output_tune)){
svm.probs[[i]] <- predict(output_tune[[i]], newdata = 
                       output_test[[i]][, colnames(output_test[[i]]) != "diliClass"],type="prob") # Gen probs for ROC

svm.ROC[[i]] <- roc(predictor=svm.probs[[i]]$severe,
               response=output_test[[i]]$diliClass,
               levels=rev(levels(output_test[[i]]$diliClass)))
}
svm.ROC

plot(svm.ROC[[1]],main="ROC for SVM", ylim = c(0,1), lt = 2)
names(svm.ROC[[1]])

sensdata <- as.data.frame(sapply(svm.ROC, "[[", 'sensitivities' ))
specdata <- as.data.frame(sapply(svm.ROC, "[[", 'specificities' ))
colnames(sensdata) <- paste0("run", 1:ncol(sensdata))
colnames(specdata) <- paste0("run", 1:ncol(specdata))
head(specdata)

sensdata <- melt(sensdata)
specdata <- melt(specdata)
colnames(sensdata)[2] <- "Sensitivity"
colnames(specdata)[2] <- "Specificity"

ROCdata <- cbind(sensdata, specdata)
head(ROCdata)
colnames(ROCdata)[3] <- "run"
lapply(svm.ROC, lines, pch = '.')
pdf(file = '../generated/results/SVM/final/roc200final.pdf', width = 8, height = 8)
ggplot(ROCdata,aes(1-Specificity, Sensitivity))+geom_line(aes(group = run), stat="smooth",method = "loess",
                                                          size = 1.5,
                                                          linetype ="dotted",
                                                          alpha = 0.4)+
  labs(title= "ROC curves 200 test-set runs: median ROC = 0.73", 
       x = "False Positive Rate (1-Specificity)", 
       y = "True Positive Rate (Sensitivity)") + theme_minimal() + coord_cartesian(ylim = c(0,1.2))
dev.off()



# figure for how often predicted

comps_correct_hm <- data.frame(count_comps_correct$frac_correct)
rownames(comps_correct_hm) <- count_comps_correct$treatment
correct_comps_Colors <- myColors
correct_comps_Colors$type <- NULL
correct_comps_Colors$diliClass <- NULL
correct_comps_Colors$SeverityLabel <- SeverityLabel


pdf("../generated/results/SVM/final/correct_predictions.pdf", width= 6, height = 18)
aheatmap(comps_correct_hm,  Colv=NA, 
         color = my.colors, #breaks = colBreaks,
         fontsize = 10, width = 10, height=10, legend = TRUE,
         annRow = count_comps_correct[, c("SeverityLabel", "DILIConcern")],annColors = correct_comps_Colors   #,txt =hm_data_data[ indOrder ,colIndex]
)
dev.off()






names(svm.ROC[[1]])
mean(sapply(svm.ROC, "[[", 'auc'))
median(sapply(svm.ROC, "[[", 'auc'))

which(sapply(lapply(output.confus, '[[', 'overall'),'[[','Accuracy') == 0.4)
output_feats[[158]]
length(table(as.character(unlist(output_feats)))[table(as.character(unlist(output_feats)))>150])
# small test
# select the >150 feats from selection of 200 runs. filter with corr filter
sel_feats_200run <- names(table(as.character(unlist(output_feats)))[table(as.character(unlist(output_feats)))>150])

