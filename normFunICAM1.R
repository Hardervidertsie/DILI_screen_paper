normFunIcam <- function(dataIn, normVar) {
  
  buffer <- dataIn[ fingerprints%in% normVar] # select 1 variable
  #subtract DMSO per plate-time followed by min max single plate within replicates
  meanDMSO <-  buffer[ treatment %in% "dmso"  ,  
                       mean(value, na.rm=TRUE),
                       by = c("plateID", "timeID")
                       ]
  
  setkeyv(meanDMSO, c("plateID", "timeID"))
  setkeyv(buffer, c("plateID" ,"timeID"))
  
  buffer<- meanDMSO[buffer]
  buffer[, DMSO_sub.value:= value - V1 ]
  buffer[, value:=NULL]
  buffer[, fingerprints:= paste("mmnDMSO.sub", normVar, sep ="_")]
  buffer[, V1:=NULL]
  
  # min and max per replicate (because plates contain compound set for specific cmax's)
  
  minValue <- buffer[ , min(DMSO_sub.value, na.rm = TRUE), by = replID] 
  maxValue <- buffer[ , max(DMSO_sub.value, na.rm = TRUE), by = replID]
  
  
  setkey(buffer, "replID")
  setkey(maxValue, "replID")
  buffer <- buffer[maxValue]
  setnames(buffer, "V1", "maxValue")
  
  setkey(buffer, "replID")
  setkey(minValue, "replID")
  buffer <- buffer[minValue]
  setnames(buffer, "V1", "minValue")
  
  buffer[, mmnvalue:= 2* ((DMSO_sub.value - minValue)/ 
           (maxValue-minValue)) - 1 ]
  
  setnames(buffer, "mmnvalue", "value")
  buffer[, DMSO_sub.value:=NULL]
  buffer[, minValue:=NULL]
  buffer[, maxValue:=NULL]
  output <- buffer
  return(output)  
}

