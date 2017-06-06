
require(data.table)
require("fda.usc")

 load("../tmp/model.results.GFP_data.df.Rdata")
model.results.GFP_data.df <- as.data.table(model.results.GFP_data.df)

all.fingerprints <- unique(model.results.GFP_data.df[, fingerprints])

model.results.GFP_data.df[ , TD := paste0(treatment, "(", cmax, ")")]
model.results.GFP_data.df[ treatment == "dmso", TD := "dmso" ]

metaTime <- unique(model.results.GFP_data.df[ , timeAfterExposure])

all.fingerprints.out = list() # stores all the fingerprints lists

for ( j in seq_along( all.fingerprints ) ) {

  anova.perm.out = list() # new list for each fingerprint

fingerprintSubset <-  model.results.GFP_data.df[ fingerprints %in% all.fingerprints[j] , ]

all.TD <- unique(fingerprintSubset[, TD])

all.TD <- all.TD[ all.TD != "dmso"]

#p21 50 cmax only 1 replicate. skip for now. Combine 2 doses in a following analysis

# if( grepl("p21", all.fingerprints[j])  ){ # if a p21 bases fingerprint: skip all 50 cmax related TD's 
#   ind50cmax <- grepl("50cmax", all.TD)
# all.TD <- all.TD[!ind50cmax]
#   }

for ( i in seq_along(all.TD)) {
  
  TDSubset <- fingerprintSubset[TD == all.TD[i] | TD == "dmso", ]
  

  TDSubset_wide <- dcast(data = TDSubset, formula = TD + replID + cmax ~ timeAfterExposure, value.var = "mod")


Group <- grepl("dmso", TDSubset_wide[ , "TD"] )
ind <- Group
Group[ ind ] <- "dmso"
Group[! ind ] <- unique(TDSubset_wide$TD[ !grepl("dmso", TDSubset_wide$TD)])
Group <- factor(Group)


# check if dcast worked ( could be wrong cmax values for dmso .. )
if( !length(Group) ==  ( nrow(TDSubset) / (length( metaTime ) ) ) ) {
  stop( paste( "dcast failed at i: ", i, 'and j: ', j ) )
}

  mdata <- as.matrix(TDSubset_wide[,-c(1, 2, 3) ])
  argvals <-metaTime
  fdata.imagedata <- fdata( mdata = mdata, argvals = argvals, rangeval = range(argvals))

  
  if( length( Group[!ind]) > 1 ) {
    anova.perm.out[[i]] <- anova.onefactor( object = fdata.imagedata, group = Group, nboot = 100, plot = FALSE)  
  } else {
    anova.perm.out[[i]] <- list(pvalue = NA, stat = NA) # no replicate...
  }
    anova.perm.out[[i]]$wm <- NULL
    
  names(anova.perm.out)[i] <- paste0( unique(as.character(Group)), collapse = "_")


    }
  all.fingerprints.out[[j]] <- anova.perm.out
  names(all.fingerprints.out)[j] <-  as.character( all.fingerprints[j] )
print(j)
  


}


result[[1]]
result <- sapply(all.fingerprints.out, "[" ) 
result.2 <- unlist(result)

write.table( result.2, file = "data/permutation_functionaldataanalysis_dili.txt", sep ="\t")


