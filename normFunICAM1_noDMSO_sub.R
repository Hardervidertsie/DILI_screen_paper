normFunIcamnodmso <- function(dataIn, normVar) {
  
  buffer <- dataIn[ fingerprints%in% normVar] # select 1 variable
  
  # min and max per replicate (because plates contain compound set for specific cmax's)
  
  minValue <- buffer[ , min(value, na.rm = TRUE), by = replID] 
  maxValue <- buffer[ , max(value, na.rm = TRUE), by = replID]
  
  setkey(buffer, "replID")
  setkey(maxValue, "replID")
  buffer <- buffer[maxValue]
  setnames(buffer, "V1", "maxValue")
  
  setkey(buffer, "replID")
  setkey(minValue, "replID")
  buffer <- buffer[minValue]
  setnames(buffer, "V1", "minValue")
  
  buffer[ , mmnvalue := (value - minValue)/ 
           (maxValue-minValue)]
  buffer[, value:=NULL]
  setnames(buffer, "mmnvalue", "value")
  buffer[, fingerprints := paste( "mmn" ,fingerprints, sep = "_")]
  
  buffer[, minValue:=NULL]
  buffer[, maxValue:=NULL]
  output <- buffer
  return(output)  
}
