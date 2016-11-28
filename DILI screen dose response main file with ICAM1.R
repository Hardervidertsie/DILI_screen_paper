



#analyze plan:
#1) single cell data van ICAM1 dataset gebruiken om cell boven een bepaalde threshhold te tellen
  #i) Eerst per tijdpunt gemiddelde DMSO (alle onder 0.2) van waarden afhalen
  #ii) Dan gemiddelde van "achtergrond" waarden bepalen adv eerste tijdpunt DMSO waarden
  #iii) deze achtergrond waarden gebruiken om single cell data te tellen dat boven en onder deze threshholds zit

#2) gemiddelden van ICAM1 berekenen, dan combineren met andere reporters
#3) time courses modelen (smoothen van data en gelijke tijdpunten sampelen)
#4) plaat normalizaties (tov DMSO per tijdpunt - DMSO is dan dus altijd plat): dit is niet nodig voor threshhold-getelde cellen
#5) alleen van de responsen een time-course heatmap (oranje-tinten positief/ blauw-tinten negatief)
# - geen cell dood time courses want deze is er alleen voor ICAM1 dataset live

#6) Dan tijd dimensie verwijderen adv absolute maximum tijd-responsen (want icam kan ook naar beneden)
#7) dose responsen modelen + cell dood + viability markers (aantal cell, speed, etc)
#8) cell dood en viability features normalizeren
#9) dose response heatmap




# dose response script for DILI data
# first figure to generate: per compound dose response fingerprints. normalize all features per plate relative to DMSO/DMEM t0
# then summarize the time dimension to meaningfull features

# Ha Steven,
# hierbij de data. Alles staat in het mapje:
#   \DILI Exposure experiments\Nikon3 - ICAM1 DILI screen\data Steven Wink
# Op de harde schijf die nu in je computer geplugd gaat worden 
# 
# 20150608 − 10cMax_wells heeft langer geimaged --> stukje afhalen
# 
# Dit zijn de correcte files:
#   
#   20160114	100	ICAM1
# 20150830	100	ICAM1
# 20140701	50	ICAM1
# 20150828	50	ICAM1
# 20150604	10	ICAM1
# 20150608	10	ICAM1
# 20150601	5	ICAM1
# 20151231	5	ICAM1
# 20150829	1	ICAM1
# 20150914	1	ICAM1
# 
# Het wijst zich redelijk vanzelf. 
# Mocht je nog andere data/iets anders nodig hebben, laat het even weten.
# Succes!
#   Suus


# analyze plan per 31-1-2016
#1) van een set compounds is alle (dus srxn1, chop en p21 en ook ICAM1) data veranderd. inladen en oude compounds weg
#2) overal moet ICAM1 bij


# 3) ik heb ICAM1 alle data op "H:\suz icam1" (mijn 3TB interne HD) gezet. rdata files zijn de single cell response data. summary data CD zijn de CD data laatste tijdpunt
# 4) data ICAM1 data van 3)  moeten deels overschreven worden door: H:\DILI screen\meta analyse DILI screen\nieuwe cmax data
# 5) de summary data files van srxn1/chop en p21 staan op H:/DILI screen/meta analyse DILI screen/Summary data files
# 6) ook 5) deels overschreven door H:\DILI screen\meta analyse DILI screen\nieuwe cmax data
# 7) celdood data van screen: H:/DILI screen/all cytotoxicity output/cytotox Summary Data.txt"
# 8) ook 7) moet deels overschreven worden door H:\DILI screen\meta analyse DILI screen\nieuwe cmax data

# deze moeten deels overschreven worden door


rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("E:/suz icam1") # laad we workspace


.libPaths( c( .libPaths(), "D:/R-3.2.2/library") )

options(stringsAsFactors = FALSE)
library(data.table)
library(ggplot2)
library(plyr)
require(grid)
require(splines)
require(stats); require(graphics)
library(reshape2)
library(NMF)


rm(list=ls())
source("E:/R_WORK/R scipts/theme_sharp.R")

path.to.nonICAM <- ("E:/DILI screen/meta analyse DILI screen/Summary data files")
dir()


raw.data =list()
dir.files <- dir(path.to.nonICAM)
dir.files <- dir.files[grepl("(Summary Data)", dir.files)]

for (i in seq_along(dir.files)){
  raw.data[[i]] <- read.table(file = paste(path.to.nonICAM, dir.files[i], sep ="/"), header = TRUE, sep ="\t")
}
head(raw.data[[1]])
# first 3 no timeAfterExposure
# 2013-07-10 1.2  0.8
# 2013-07-24 1.076 tr  40min del
# 2013-07-29  59 min  20 min


#load ICAM single cell data 
ICAM.files <- dir("Rdata files/")[grepl(".Rdata", dir("Rdata files/") )]

ICAM.raw =list()

for(i in seq_along(ICAM.files)) 
{
print(i)
    load(paste("Rdata files/",ICAM.files[i], sep ="/"))
  ICAM.raw[[i]] <- outputList
rm("outputList")
  }







#fix column name inconsistencies 

colnames(ICAM.raw[[4]]$myDT)
#img -> image
colnames(ICAM.raw[[1]]$myDT) <- gsub("img", "image", colnames(ICAM.raw[[1]]$myDT))
colnames(ICAM.raw[[2]]$myDT) <- gsub("img", "image", colnames(ICAM.raw[[2]]$myDT))
colnames(ICAM.raw[[3]]$myDT) <- gsub("img", "image", colnames(ICAM.raw[[3]]$myDT))

#obj_cytoplasm_  -> obj_cyto_only_  
colnames(ICAM.raw[[1]]$myDT) <- gsub("obj_cytoplasm_", "obj_cyto_only_", colnames(ICAM.raw[[1]]$myDT))
colnames(ICAM.raw[[2]]$myDT) <- gsub("obj_cytoplasm_", "obj_cyto_only_", colnames(ICAM.raw[[2]]$myDT))
colnames(ICAM.raw[[3]]$myDT) <- gsub("obj_cytoplasm_", "obj_cyto_only_", colnames(ICAM.raw[[3]]$myDT))

for(i in seq_along(ICAM.raw)) {
  colnames(ICAM.raw[[i]]$myDT) <- gsub("DistanceTraveled_[0-9]{2}", "DistanceTraveled", colnames(ICAM.raw[[i]]$myDT))
}

#column names of variables to be counted
TH.cols <- colnames(ICAM.raw[[1]]$myDT)[c(16, 17, 18, 19 )]


#calculate mean of DMSO for <= 0.2 concentration for each time point
mean.dmso.colnames <- paste("meanDMSO_", TH.cols, sep="")

unique(as.numeric(ICAM.raw[[1]]$myDT[, timeID]))

#ICAM.raw[[i]]$myDT[, mean.dmso.colnames]
meanDMSO = list()
for(i in seq_along(ICAM.raw)) {
  print(paste("i: ", i ))
  ICAM.raw[[i]]$myDT[, dose_uM:= as.numeric(as.character(ICAM.raw[[i]]$myDT[, dose_uM]))]# change factor to numeric of dose 
    meanDMSO[[i]] <-  ICAM.raw[[i]]$myDT[ treatment %in% "DMSO" & dose_uM <= 0.2 ,  
                                                               lapply(.SD, function(x) { mean(x, na.rm=TRUE) }),
                                                               by = c("plateID",   "timeID"),
                                                               .SDcols =
                                                                 TH.cols
                                                               ]
    setkeyv(ICAM.raw[[i]]$myDT, c("plateID", "timeID"))
    setkeyv(meanDMSO[[i]], c("plateID", "timeID"))
    
    ICAM.raw[[i]]$myDT <- ICAM.raw[[i]]$myDT[meanDMSO[[i]]]
    
          }

# subtract these values
normCols <- paste("norm", TH.cols, sep ="_")

for(i in seq_along(ICAM.raw)) {
  print(paste("i: ", i ))
 
  for(j in seq_along(normCols)) {
    print(j)
    print(normCols[j])
    print(TH.cols[j])
    print(paste("i", TH.cols,sep=".")[j])
    
       ICAM.raw[[i]]$myDT[, normCols[j]:= (get(TH.cols[j]) - get(paste("i", TH.cols,sep=".")[j] )) ] # the"i"  columns are the DMSO averages
  }
}

TH.cols.mean <- paste("backgroundValue.",TH.cols ,sep = "")
mean1above <- paste( "mean1above", TH.cols, sep ="_")
mean2above <- paste( "mean2above", TH.cols, sep ="_")
mean3above <- paste( "mean3above", TH.cols, sep ="_")
mean1below <- paste( "mean1below", TH.cols, sep ="_")
mean2below <- paste( "mean2below", TH.cols, sep ="_")
mean3below <- paste( "mean3below", TH.cols, sep ="_")
# count cells above and below background (of DMSO time point 1)
# check if timepoint 1 has indeed average of around 0 for all the dmso' s <=0.2
#hier gebelven, kan code opnieuw gebruiken\]-

#remove redundent data to save some memory
colnames(ICAM.raw[[5]]$myDT)

rmCols <- c("imageNumber", "plateWellID", "groupNumber","imageID","plateWellID","obj_cells_Number_Object_Number",
            "i.obj_icam1_Intensity_IntegratedIntensity_image_GFP","i.obj_icam1_Intensity_MeanIntensity_image_GFP",
            "i.obj_cyto_only_Intensity_IntegratedIntensity_image_GFP","i.obj_cyto_only_Intensity_MeanIntensity_image_GFP")
for(i in seq_along(ICAM.raw)) {
  print(paste("i: ", i ))
  for(j in seq_along(rmCols)){
  ICAM.raw[[i]]$myDT[, rmCols[j]:= NULL ]
  }
}
#make backup of all single cell data (had wat corrupte Rdata files van ICAM1 op image data schijf....)
#save.image("bkupICAM1_1_2_2016.RData")
#hier Dinsdag avond verder. Probeer al een flink stuk achter de kiezen te hebben voor donderdag.
#load("bkupICAM1_1_2_2016.RData")

# count norm cols above and bellow thresholds.
for(i in seq_along(ICAM.raw)) {
 print(paste("i: ", i ))

  
#ICAM.raw[[i]][["myDT"]][, get("TH.cols.mean"):=NULL]

 ICAM.raw[[i]][["myDT"]][, get(
    "TH.cols.mean"):= list(mean(get(TH.cols[1])[treatment %in% "DMSO" & dose_uM <= 0.2 & timeID==1], na.rm = TRUE), 
                           mean(get(TH.cols[2])[treatment %in% "DMSO" & dose_uM <= 0.2 & timeID==1], na.rm = TRUE),
                           mean(get(TH.cols[3])[treatment %in% "DMSO" & dose_uM <= 0.2 & timeID==1], na.rm = TRUE),
                           mean(get(TH.cols[4])[treatment %in% "DMSO" & dose_uM <= 0.2 & timeID==1], na.rm = TRUE))]

  
  
  for(j in seq_along(normCols)) { # dus aantal cellen die boven achtergrond uitsteken nadat achtergrond er al af is gehaald
    print(j) 
  # normcol is 1 achtergrond waarde reeds afgehaald. dus -1 aan rechter zijde gedaan
    ICAM.raw[[i]]$myDT[, mean1above[j]:= ((get(normCols[j])) >= (0.5*get(TH.cols.mean[j]))  )] 
    ICAM.raw[[i]]$myDT[, mean2above[j]:= ((get(normCols[j])) >= (1*get(TH.cols.mean[j]))  )] 
    ICAM.raw[[i]]$myDT[, mean3above[j]:= ((get(normCols[j])) >= (2*get(TH.cols.mean[j])) )] 
    ICAM.raw[[i]]$myDT[, mean1below[j]:= ((get(normCols[j])) <= (-0.5*get(TH.cols.mean[j]) ) )] 
    ICAM.raw[[i]]$myDT[, mean2below[j]:= ((get(normCols[j])) <= (-1*get(TH.cols.mean[j]))  )] 
    ICAM.raw[[i]]$myDT[, mean3below[j]:= ((get(normCols[j])) <= (-2*get(TH.cols.mean[j]) ) )] 
    
  }
}

# calculate mean

all.icam.feats<-c("imageCountTracked", "obj_icam1_AreaShape_Area",
             "obj_icam1_Intensity_IntegratedIntensity_image_GFP","obj_icam1_Intensity_MeanIntensity_image_GFP",
             "obj_cyto_only_Intensity_IntegratedIntensity_image_GFP","obj_cyto_only_Intensity_MeanIntensity_image_GFP",
             "obj_nuclei_AreaShape_Area","obj_nuclei_Intensity_MeanIntensity_image_hoechst",
             "obj_nuclei_TrackObjects_DistanceTraveled", "obj_cells_Children_obj_icam1_Count",
             normCols,
             TH.cols.mean, mean1above, mean2above, mean3above, mean1below, mean2below, mean3below
             
)

mean.data = list()

all.icam.feats[ !all.icam.feats%in% colnames(ICAM.raw[[2]]$myDT)] # check if column names correct now


#calculate mean, however what about NA values? no problem I think: mean(c(NA,NA,NA,1,0),na.rm=TRUE) = 0.5
for(i in seq_along(ICAM.raw)){
print(i)
    mean.data[[i]] <-  ICAM.raw[[i]]$myDT[ , lapply(.SD, function(x) { mean(x, na.rm=TRUE) }),
                                            by = c("treatment", "dose_uM", "cell_line", "plateID", "timeID", "timeAfterExposure"),
                                  .SDcols = all.icam.feats
                            ]
  }

mean.data

all.mean.d <- do.call('rbind', mean.data)



#25.8 27.2 28.6 30.0 31.4 32.8 34.2 35.6 37.0 38.4 moet weg in  20150608 - 10cMax_wells
dim(all.mean.d)
all.mean.d <- all.mean.d[ !( plateID %in% "20150608 - 10cMax_wells" &
                            timeAfterExposure %in% c(25.8, 27.2, 28.6, 30.0, 31.4, 32.8, 34.2, 35.6, 37.0, 38.4)) ,]

all.mean.d.long<- melt(all.mean.d,  id.vars = c("treatment", "dose_uM", "cell_line", "plateID", "timeID", "timeAfterExposure"))
tail(all.mean.d.long)



head(all.mean.d)
head(all.mean.d.long)

all.treats<-unique(all.mean.d.long[, treatment])
unique(all.mean.d.long[, plateID])

#20150608 ??? 10cMax_wells timeAfterExposure niet goed exposureDelay <- c("00:35" ) # hh:mm  timeBetweenFrames <- c("01:27:00" )  verkeerd om?

#for plotting remove put DMSO as zero

all.mean.d.long.dmso0<- all.mean.d.long
all.mean.d.long.dmso0<-as.data.table(all.mean.d.long.dmso0)
all.mean.d.long.dmso0[ treatment %in% "DMSO", dose_uM:=0]
all.icam.feats<- unique(all.mean.d.long.dmso0[, variable])
for(i in seq_along(all.icam.feats)){

    pdf(paste( all.icam.feats[i],".pdf" ,sep=""), width = 100, height =80)
  p<-ggplot(all.mean.d.long.dmso0[variable%in% all.icam.feats[i] & !treatment %in% c("g.f. DMEM" )],
            aes(x=timeAfterExposure, y = value, color = as.character(plateID))) + 
    geom_line(aes(color=as.character(plateID), group = plateID)) + geom_point(aes(size = as.factor(dose_uM)))+  facet_wrap(~treatment)
  print(p)
  dev.off()
}



all.mean.d.long.dmso0[variable%in% all.icam.feats[21],]




# hier de nieuwe data inladen van aangepaste Cmax waarden. zelfde procedure als boven maar andere technieken ivm dose ranges op 
# ivm geheugen summary data wegschrijven en geheugen leeg maken

write.table(all.mean.d, file = "ICAM1_origall.mean.data.wide.txt", sep ="\t", col.names = TRUE)
write.table(all.mean.d.long, file = "ICAM1_orig.all.mean.data.long.txt", sep ="\t", col.names = TRUE)

rm(list=ls())
gc()
icam.orig.all.mean.d <- read.table("ICAM1_origall.mean.data.wide.txt", header= TRUE, sep ="\t")
icam.orig.all.mean.d.long <- read.table("ICAM1_orig.all.mean.data.long.txt", header = TRUE, sep ="\t" )

##### begin plak bovenstaande data voor de cmax aangepaste nieuwe ICAM RData

library(data.table)
library(ggplot2)
library(plyr)
require(grid)
require(splines)
require(stats); require(graphics)
library(reshape2)




#source("D:/R_HOME/R_WORK/R scipts/theme_sharp.R")

path.to.nonICAM <- ("E:/DILI screen/meta analyse DILI screen/Summary data files")
dir(path.to.nonICAM)






raw.data =list()
dir.files <- dir(path.to.nonICAM)
dir.files <- dir.files[grepl("(Summary Data)", dir.files)]

for (i in seq_along(dir.files)){
  raw.data[[i]] <- read.table(file = paste(path.to.nonICAM, dir.files[i], sep ="/"), header = TRUE, sep ="\t")
}
head(raw.data[[1]])
# first 3 no timeAfterExposure
# 2013-07-10 1.2  0.8
# 2013-07-24 1.076 tr  40min del
# 2013-07-29  59 min  20 min


#load new cmax compound ICAM single cell data 
# 930 en 1001 zijn rdata files van 1 plaat dat halverwege gecrasht was
ICAM.files <- dir("H:/DILI screen/meta analyse DILI screen/nieuwe cmax data/ICAM1/Rdata files")[grepl(".Rdata", 
dir("H:/DILI screen/meta analyse DILI screen/nieuwe cmax data/ICAM1/Rdata files/") )]

ICAM.raw =list()

for(i in seq_along(ICAM.files)) 
{
  load(paste("H:/DILI screen/meta analyse DILI screen/nieuwe cmax data/ICAM1/Rdata files",ICAM.files[i], sep ="/"))
  ICAM.raw[[i]] <- outputList
  rm("outputList")
}

head(ICAM.raw[[1]])


#fix column name inconsistencies 

colnames(ICAM.raw[[1]]$myDT)


for(i in seq_along(ICAM.raw)) {
  colnames(ICAM.raw[[i]]$myDT) <- gsub("DistanceTraveled_[0-9]{2}", "DistanceTraveled", colnames(ICAM.raw[[i]]$myDT))
}

#column names of variables to be counted
TH.cols <- colnames(ICAM.raw[[1]]$myDT)[c(16, 17, 18, 19 )]
TH.cols

head(ICAM.raw[[i]]$myDT)
unique(ICAM.raw[[3]]$myDT[, dose_uM])
head(all.mean.d)

#calculate mean of DMSO for <= 0.2 concentration for each time point
mean.dmso.colnames <- paste("meanDMSO_", TH.cols, sep="")




# De plaat met SRXN1 en ICAM1 van 0930 op 1001 is gecrashed na 5 uur en 45 minuten. 
# Het eerste experiment (dus alles in map 0930)  was 's middags gestart om 16:55
# Het experiment is herstart om 9:02
# Het tweede deel (dus alles in map 1001) duurt nog 7 uur. Tijdpunten zijn elke 1 uur en 35 minuten. 

# dus 5:45 --> gecrasht om 16:55 + 5u45 min = 22:40 gecrasht. Start 10:20 later

unique((ICAM.raw[[3]]$myDT[, timeID])) 
#change plateID
(ICAM.raw[[3]]$myDT[, plateID:="20150930_ICAM1"])
#modify time



#ICAM.raw[[i]]$myDT[, mean.dmso.colnames]

timeBetweenFrames <- "01:35:00"
exposureDelay <- "16:05"


timeBetweenFrames <- round(as.integer(strftime(strptime(timeBetweenFrames, format = "%H:%M:%S"), "%H")) + 
                             1/60 * as.integer(strftime(strptime(timeBetweenFrames, 
                                                                 format = "%H:%M:%S"), "%M")) +
                             1/3600 * as.integer(strftime(strptime(timeBetweenFrames, 
                                                                   format = "%H:%M:%S"), "%S"))
                           , digit =2 )

exposureDelay <- round(as.integer(strftime(strptime(exposureDelay, format = "%H:%M"), "%H")) + 
                         1/60 * as.integer(strftime(strptime(exposureDelay, 
                                                             format = "%H:%M"), "%M")), digit =1 )
ICAM.raw[[3]]$myDT[, timeAfterExposure:=round(as.numeric(timeID)*timeBetweenFrames+exposureDelay - timeBetweenFrames, digits = 1) ]
unique((ICAM.raw[[3]]$myDT[, timeAfterExposure])) 
unique((ICAM.raw[[3]]$myDT[, timeID])) 
((ICAM.raw[[3]]$myDT[, timeID:=as.numeric(timeID)+4])) 

#combine:
object.size(ICAM.raw[[2]]$myDT) #121289408 bytes
object.size(ICAM.raw[[3]]$myDT) #159178016 bytes
nrow(ICAM.raw[[2]]$myDT) #673454
nrow(ICAM.raw[[3]]$myDT) #864733

ICAM.raw[[2]]$myDT <- rbind(ICAM.raw[[2]]$myDT, ICAM.raw[[3]]$myDT)
object.size(ICAM.raw[[2]]$myDT)  #276963688 bytes
nrow(ICAM.raw[[2]]$myDT)  #1538187
#673454+864733 = 1538187
ICAM.raw[[3]] <- NULL
unique(ICAM.raw[[2]]$myDT[, timeAfterExposure])

ICAM.raw[[3]]$myDT[ treatment%in%"DMSO" & as.numeric(as.character(dose_uM)) <= 0.2]




meanDMSO = list()
for(i in seq_along(ICAM.raw)) {
  print(paste("i: ", i ))
  ICAM.raw[[i]]$myDT[, dose_uM:= as.numeric(as.character(ICAM.raw[[i]]$myDT[, dose_uM]))]# change factor to numeric of dose 
  meanDMSO[[i]] <-  ICAM.raw[[i]]$myDT[ treatment %in% "DMSO" & dose_uM <= 0.2 ,  
                                        lapply(.SD, function(x) { mean(x, na.rm=TRUE) }),
                                        by = c("plateID", "timeID"),
                                        .SDcols =
                                          TH.cols
                                        ]
  setkeyv(ICAM.raw[[i]]$myDT, c("plateID", "timeID"))
  setkeyv(meanDMSO[[i]], c("plateID", "timeID"))
  
  ICAM.raw[[i]]$myDT <- ICAM.raw[[i]]$myDT[meanDMSO[[i]]]
  
}
meanDMSO[[1]]
ICAM.raw[[1]]$myDT

# subtract these values
normCols <- paste("norm", TH.cols, sep ="_")

for(i in seq_along(ICAM.raw)) {
  print(paste("i: ", i ))
  
  for(j in seq_along(normCols)) {
    print(j)
    print(normCols[j])
    print(TH.cols[j])
    print(paste("i", TH.cols,sep=".")[j])
    
    ICAM.raw[[i]]$myDT[, normCols[j]:= (get(TH.cols[j]) - get(paste("i", TH.cols,sep=".")[j] )) ] # the"i"  columns are the DMSO averages
  }
}

TH.cols.mean <- paste("backgroundValue.",TH.cols ,sep = "")
mean1above <- paste( "mean1above", TH.cols, sep ="_")
mean2above <- paste( "mean2above", TH.cols, sep ="_")
mean3above <- paste( "mean3above", TH.cols, sep ="_")
mean1below <- paste( "mean1below", TH.cols, sep ="_")
mean2below <- paste( "mean2below", TH.cols, sep ="_")
mean3below <- paste( "mean3below", TH.cols, sep ="_")
# count cells above and below background (of DMSO time point 1)
# check if timepoint 1 has indeed average of around 0 for all the dmso' s <=0.2
#hier gebelven, kan code opnieuw gebruiken\]-

#remove redundent data to save some memory
colnames(ICAM.raw[[1]]$myDT)

rmCols <- c("imageNumber", "plateWellID", "groupNumber","imageID","plateWellID","obj_cells_Number_Object_Number",
            "i.obj_icam1_Intensity_IntegratedIntensity_image_GFP","i.obj_icam1_Intensity_MeanIntensity_image_GFP",
            "i.obj_cyto_only_Intensity_IntegratedIntensity_image_GFP","i.obj_cyto_only_Intensity_MeanIntensity_image_GFP")
for(i in seq_along(ICAM.raw)) {
  print(paste("i: ", i ))
  for(j in seq_along(rmCols)){
    ICAM.raw[[i]]$myDT[, rmCols[j]:= NULL ]
  }
}
#hiero
# count norm cols above and bellow thresholds. 
for(i in seq_along(ICAM.raw)) {
  print(paste("i: ", i ))

  #ICAM.raw[[i]][["myDT"]][, get("TH.cols.mean"):=NULL]
  
  ICAM.raw[[i]][["myDT"]][, get(
    "TH.cols.mean"):= list(mean(get(TH.cols[1])[treatment %in% "DMSO" & dose_uM <= 0.2 & timeID==1], na.rm = TRUE), 
                           mean(get(TH.cols[2])[treatment %in% "DMSO" & dose_uM <= 0.2 & timeID==1], na.rm = TRUE),
                           mean(get(TH.cols[3])[treatment %in% "DMSO" & dose_uM <= 0.2 & timeID==1], na.rm = TRUE),
                           mean(get(TH.cols[4])[treatment %in% "DMSO" & dose_uM <= 0.2 & timeID==1], na.rm = TRUE))]

  for(j in seq_along(normCols)) { # dus aantal cellen die boven achtergrond uitsteken nadat achtergrond er al af is gehaald
    print(j) 
    # normcol is 1 achtergrond waarde reeds afgehaald. dus -1 aan rechter zijde gedaan
    ICAM.raw[[i]]$myDT[, mean1above[j]:= ((get(normCols[j])) >= (0.5*get(TH.cols.mean[j]))  )] 
    ICAM.raw[[i]]$myDT[, mean2above[j]:= ((get(normCols[j])) >= (1*get(TH.cols.mean[j]))  )] 
    ICAM.raw[[i]]$myDT[, mean3above[j]:= ((get(normCols[j])) >= (2*get(TH.cols.mean[j])) )] 
    ICAM.raw[[i]]$myDT[, mean1below[j]:= ((get(normCols[j])) <= (-0.5*get(TH.cols.mean[j]) ) )] 
    ICAM.raw[[i]]$myDT[, mean2below[j]:= ((get(normCols[j])) <= (-1*get(TH.cols.mean[j]))  )] 
    ICAM.raw[[i]]$myDT[, mean3below[j]:= ((get(normCols[j])) <= (-2*get(TH.cols.mean[j]) ) )] 
    
  }
}

# calculate mean
all.icam.feats<-c("imageCountTracked", "obj_icam1_AreaShape_Area",
                  "obj_icam1_Intensity_IntegratedIntensity_image_GFP","obj_icam1_Intensity_MeanIntensity_image_GFP",
                  "obj_cyto_only_Intensity_IntegratedIntensity_image_GFP","obj_cyto_only_Intensity_MeanIntensity_image_GFP",
                  "obj_nuclei_AreaShape_Area","obj_nuclei_Intensity_MeanIntensity_image_hoechst",
                  "obj_nuclei_TrackObjects_DistanceTraveled", 
                  normCols,
                  TH.cols.mean, mean1above, mean2above, mean3above, mean1below, mean2below, mean3below
)

all( all.icam.feats %in% colnames(ICAM.raw[[i]]$myDT) )

mean.data = list()

#all.icam.feats[ !all.icam.feats%in% colnames(ICAM.raw[[2]]$myDT)] # check if column names correct now
#save.image("bkup31012016.RData")

#calculate mean, however what about NA values? no problem I think: mean(c(NA,NA,NA,1,0),na.rm=TRUE) = 0.5
for(i in seq_along(ICAM.raw)){
  print(i)
  mean.data[[i]] <-  ICAM.raw[[i]]$myDT[ , lapply(.SD, function(x) { mean(x, na.rm=TRUE) }),
                                         by = c("treatment", "dose_uM", "cell_line", "plateID", "timeID", "timeAfterExposure"),
                                         .SDcols = all.icam.feats
                                         ]
}



new.cmaxICAM1 <- do.call('rbind', mean.data)

##### einde plak code ICAM RDAta aangepaste compound

#heb nu dus:

icam.orig.all.mean.d <-as.data.table(icam.orig.all.mean.d )
icam.orig.all.mean.d.long <- as.data.table(icam.orig.all.mean.d.long)
raw.data # orig srxn1 chop & p21 summarized data

#Carmustine 0.037 ipv 0.0356
new.cmaxICAM1[ treatment%in% "Carmustine" & dose_uM %in% 0.356, dose_uM:= 0.037]
new.cmaxICAM1[ treatment%in% "Carmustine" & dose_uM %in% 1.78, dose_uM:= 5*0.037]
new.cmaxICAM1[ treatment%in% "Carmustine" & dose_uM %in% 3.56, dose_uM:= 10*0.037]
new.cmaxICAM1[ treatment%in% "Carmustine" & dose_uM %in% 8.90, dose_uM:= 25*0.037]
new.cmaxICAM1[ treatment%in% "Carmustine" & dose_uM %in% 17.80, dose_uM:= 50*0.037]
new.cmaxICAM1[ treatment%in% "Carmustine" & dose_uM %in% 35.60, dose_uM:= 100*0.037]

write.table(new.cmaxICAM1, file ="new.cmaxICAM1.txt", col.names = TRUE, sep ="\t")

# pdfjes om nieuwe cmax ICAM1 data te verifieren


#new cmax have changed and 1 extra
#1
#5
#10
#25 <-- new
#50
#100

newCmax <- read.table("newCmax.txt", sep ="\t")

newCmax$cmax <- 1
colnames(newCmax) <- c("treatment", "dose_uM", "cmax")
cmax5 <- newCmax
cmax10 <- newCmax
cmax25 <- newCmax
cmax50 <- newCmax
cmax100 <- newCmax

cmax5$dose_uM<-newCmax$dose_uM*5
cmax5$cmax <- 5
cmax10$dose_uM<- newCmax$dose_uM*10
cmax10$cmax<-10
cmax25$dose_uM<- newCmax$dose_uM*25
cmax25$cmax<-25
cmax50$dose_uM <- newCmax$dose_uM*50
cmax50$cmax<-50
cmax100$dose_uM <- newCmax$dose_uM*100
cmax100$cmax<-100


newCmax<- rbind(newCmax,cmax5,cmax10, cmax25, cmax50, cmax100)
newCmax<- as.data.table(newCmax)
setkeyv(new.cmaxICAM1, c("treatment", "dose_uM") )
setkeyv(newCmax, c("treatment", "dose_uM"))

new.cmaxICAM1<-newCmax[new.cmaxICAM1]

new.cmaxICAM1.long<- melt(new.cmaxICAM1,  id.vars = c("treatment", "dose_uM","cmax", "cell_line", "plateID", "timeID", "timeAfterExposure"))
head(new.cmaxICAM1.long)

sum(is.na(new.cmaxICAM1.long[, cmax]))

unique(new.cmaxICAM1.long[is.na(new.cmaxICAM1.long[, cmax]), cmax:=0.2])



all.treats<-unique(new.cmaxICAM1.long[, treatment])


new.cmaxICAM1.long.dmso0<- new.cmaxICAM1.long
new.cmaxICAM1.long.dmso0<-as.data.table(new.cmaxICAM1.long.dmso0)
new.cmaxICAM1.long.dmso0[ treatment %in% "DMSO", dose_uM:=0]
all.icam.feats<- unique(new.cmaxICAM1.long.dmso0[, variable])
new.cmaxICAM1.long.dmso0[, replID:= plateID]
new.cmaxICAM1.long.dmso0[, plateID:= paste(plateID, cmax, sep="_")]
for(i in seq_along(all.icam.feats)){
  
  pdf(paste( all.icam.feats[i],"_newCMAXICAM.pdf" ,sep=""), width = 30, height =20)
  p<-ggplot(new.cmaxICAM1.long.dmso0[variable%in% all.icam.feats[i] ],
            aes(x=timeAfterExposure, y = value, color = as.character(cmax))) + 
    geom_line(aes( group = plateID)) + geom_point(aes(shape=replID) ) + facet_wrap(~treatment)
  print(p)
  dev.off()
}


write.table(new.cmaxICAM1.long, file ="new.cmaxICAM1.long.txt", sep ="\t", col.names=TRUE)

# hier verder: eventueel normalizatie srxn1 chop p21 aanpassen. Moet nog eea min max normalizeren (de niet tel ICAM1 dingen)
# en niet vergeten dat nieuwe cmax oude moet vervangen
timeBetweenFrames <- 1.2
exposureDelay <- 0.8
raw.data[[1]]$timeAfterExposure<- as.integer(raw.data[[1]]$timeID) * timeBetweenFrames + exposureDelay - timeBetweenFrames
head(raw.data[[1]])

timeBetweenFrames <- 1.076
exposureDelay <- 0.667
raw.data[[2]]$timeAfterExposure<- as.integer(raw.data[[2]]$timeID) * timeBetweenFrames + exposureDelay - timeBetweenFrames
head(raw.data[[2]])

timeBetweenFrames <- 0.983
exposureDelay <- 0.333
raw.data[[3]]$timeAfterExposure<- as.integer(raw.data[[3]]$timeID) * timeBetweenFrames + exposureDelay - timeBetweenFrames
head(raw.data[[3]])

lapply(raw.data, ncol)

raw.data[[26]]$doseLevel <- NULL

raw.data <- do.call("rbind" ,raw.data)




raw.data$control <-NULL
sum(is.na(raw.data))


head(raw.data)

raw.data$variable <- gsub("_>_", "_larger_", raw.data$variable)
raw.data$variable <- gsub("_<_", "_smaller_", raw.data$variable)

# fix cell line names
unique(raw.data$cell_line)
raw.data$cell_line <- gsub("Chop", "DDIT3", raw.data$cell_line)

unique(raw.data$variable)

# add replicate ID
unique(raw.data$plateID)
raw.data$replID <-NA
raw.data$cmax <- NA


# first the 9 icam1 plates
 #cmax
 raw.data[ raw.data$plateID == "20140306 - 100cMax_wells", "cmax"] <-   "100cmax"  #icam1 100cmax repl1
 raw.data[ raw.data$plateID == "20150830 - 100cMax_wells", "cmax"] <-   "100cmax"  #icam1 100cmax repl2
 
 raw.data[ raw.data$plateID == "20140701 - 50cMax_wells", "cmax"] <-   "50cmax"  #icam1 50cmax repl1
 raw.data[ raw.data$plateID == "20150828 - 50cMax_wells", "cmax"] <-   "50cmax"  #icam1 50cmax repl2
 
 raw.data[ raw.data$plateID == "20150604 - 10cMax_wells", "cmax"] <-   "10cmax"  #icam1 10cmax repl1
 raw.data[ raw.data$plateID == "20150608 - 10cMax_wells", "cmax"] <-   "10cmax"  #icam1 10cmax repl2
 
 raw.data[ raw.data$plateID == "20150601 - 5cMax_wells", "cmax"] <-   "5cmax"  #icam1 5cmax repl1
 raw.data[ raw.data$plateID == "20140606 - 5cMax_wells", "cmax"] <-   "5cmax"  #icam1 5cmax repl2
 
 raw.data[ raw.data$plateID == "20150829 - 1cMax_wells", "cmax"] <-   "1cmax"  #icam1 1cmax repl1

 
 #replID
 raw.data[ raw.data$plateID == "20140306 - 100cMax_wells", "replID"] <-   "replicate_1"  #icam1 100cmax repl1
 raw.data[ raw.data$plateID == "20150830 - 100cMax_wells", "replID"] <-   "replicate_2"  #icam1 100cmax repl2
 
 raw.data[ raw.data$plateID == "20140701 - 50cMax_wells", "replID"] <-   "replicate_1"  #icam1 50cmax repl1
 raw.data[ raw.data$plateID == "20150828 - 50cMax_wells", "replID"] <-   "replicate_2"  #icam1 50cmax repl2
 
 raw.data[ raw.data$plateID == "20150604 - 10cMax_wells", "replID"] <-   "replicate_1"  #icam1 10cmax repl1
 raw.data[ raw.data$plateID == "20150608 - 10cMax_wells", "replID"] <-   "replicate_2"  #icam1 10cmax repl2
 
 raw.data[ raw.data$plateID == "20150601 - 5cMax_wells", "replID"] <-   "replicate_1"  #icam1 5cmax repl1
 raw.data[ raw.data$plateID == "20140606 - 5cMax_wells", "replID"] <-   "replicate_2"  #icam1 5cmax repl2
 
 raw.data[ raw.data$plateID == "20150829 - 1cMax_wells", "replID"] <-   "replicate_1"  #icam1 1cmax repl1
 
 
 

raw.data[ raw.data$plateID == "2013_07_10_Srxn1_z1", "replID"] <-   "replicate_1"  #Srxn1 100cmax
raw.data[ raw.data$plateID == "2013_07_10_Srxn1_z1", "cmax"] <-   "100cmax"  #Srxn1 100cmax

raw.data[ raw.data$plateID == "2013_07_24_Srxn1", "replID"] <-      "replicate_2"  #Srxn1 100cmax
raw.data[ raw.data$plateID == "2013_07_24_Srxn1", "cmax"] <-      "100cmax"  #Srxn1 100cmax

raw.data[ raw.data$plateID == "2013_07_29_Srxn1", "replID"] <-      "replicate_3"  #Srxn1 100cmax
raw.data[ raw.data$plateID == "2013_07_29_Srxn1", "cmax"] <-      "100cmax"  #Srxn1 100cmax

raw.data[ raw.data$plateID == "2013_08_01_Srxn1", "replID"] <-      "replicate_4"  #Srxn1 100cmax
raw.data[ raw.data$plateID == "2013_08_01_Srxn1", "cmax"] <-      "100cmax"  #Srxn1 100cmax

raw.data[ raw.data$plateID == "2013_08_14_Srxn1", "replID"] <-      "replicate_5"  #Srxn1 100cmax
raw.data[ raw.data$plateID == "2013_08_14_Srxn1", "cmax"] <-      "100cmax"  #Srxn1 100cmax

raw.data[ raw.data$plateID == "2013_08_15_Srxn1", "replID"] <-      "replicate_1"  #Srxn1 50cmax
raw.data[ raw.data$plateID == "2013_08_15_Srxn1", "cmax"] <-      "50cmax"  #Srxn1 50cmax

raw.data[ raw.data$plateID == "2013_08_21_Srxn1", "replID"] <-      "replicate_2"  #Srxn1 50cmax
raw.data[ raw.data$plateID == "2013_08_21_Srxn1", "cmax"] <-      "50cmax"  #Srxn1 50cmax

raw.data[ raw.data$plateID == "2013_08_22_Srxn1", "replID"] <-      "replicate_3"  #Srxn1 50cmax     
raw.data[ raw.data$plateID == "2013_08_22_Srxn1", "cmax"] <-      "50cmax"  #Srxn1 50cmax

raw.data[ raw.data$plateID == "2013_08_28_nii_Srxn1", "replID"] <-  "replicate_1"  #Srxn1 10cmax
raw.data[ raw.data$plateID == "2013_08_28_nii_Srxn1", "cmax"] <-  "10cmax"  #Srxn1 10cmax

raw.data[ raw.data$plateID == "2013_08_28_niii_Srxn1", "replID"] <- "replicate_2"  #Srxn1 10cmax  
raw.data[ raw.data$plateID == "2013_08_28_niii_Srxn1", "cmax"] <- "10cmax"  #Srxn1 10cmax

raw.data[ raw.data$plateID == "2013_09_04_Srxn1", "replID"] <-      "replicate_3"  #Srxn1 10cmax
raw.data[ raw.data$plateID == "2013_09_04_Srxn1", "cmax"] <-      "10cmax"  #Srxn1 10cmax

raw.data[ raw.data$plateID == "2013_09_05_Srxn1", "replID"] <-      "replicate_1"  #Srxn1 5cmax
raw.data[ raw.data$plateID == "2013_09_05_Srxn1", "cmax"] <-      "5cmax"  #Srxn1 5cmax

raw.data[ raw.data$plateID == "2013_09_09_Srxn1", "replID"] <-      "replicate_2"  #Srxn1 5cmax
raw.data[ raw.data$plateID == "2013_09_09_Srxn1", "cmax"] <-      "5cmax"  #Srxn1 5cmax

raw.data[ raw.data$plateID == "2013_09_16_Srxn1", "replID"] <-      "replicate_1"  #Srxn1 1cmax
raw.data[ raw.data$plateID == "2013_09_16_Srxn1", "cmax"] <-      "1cmax"  #Srxn1 1cmax


raw.data$plateID <- gsub("25-9-2013_Srxn1", "2013_09_25_Srxn1", raw.data$plateID)

raw.data[ raw.data$plateID == "2013_09_25_Srxn1", "replID"] <-       "replicate_2"  #Srxn1 1cmax
raw.data[ raw.data$plateID == "2013_09_25_Srxn1", "cmax"] <-       "1cmax"  #Srxn1 1cmax

raw.data[ raw.data$plateID == "2013_09_26_DDIT3", "replID"] <-      "replicate_1"  #DDIT3 100cmax
raw.data[ raw.data$plateID == "2013_09_26_DDIT3", "cmax"] <-      "100cmax"  #DDIT3 100cmax

raw.data[ raw.data$plateID == "2013_10_03_DDIT3", "replID"] <-      "replicate_2"  #DDIT3 100cmax
raw.data[ raw.data$plateID == "2013_10_03_DDIT3", "cmax"] <-      "100cmax"  #DDIT3 100cmax

raw.data[ raw.data$plateID == "2013_10_04_DDIT3", "replID"] <-      "replicate_1"  #DDIT3 50cmax
raw.data[ raw.data$plateID == "2013_10_04_DDIT3", "cmax"] <-      "50cmax"  #DDIT3 50cmax

raw.data[ raw.data$plateID == "2013_10_09_niii_p21", "replID"] <-   "replicate_1"  #p21 100cmax
raw.data[ raw.data$plateID == "2013_10_09_niii_p21", "cmax"] <-   "100cmax"  #p21 100cmax

raw.data[ raw.data$plateID == "2013_10_09_nii_p21", "replID"] <-    "replicate_2"  #p21 100cmax
raw.data[ raw.data$plateID == "2013_10_09_nii_p21", "cmax"] <-    "100cmax"  #p21 100cmax

raw.data[ raw.data$plateID == "2013_10_17_nii_DDIT3", "replID"] <-  "replicate_2"  #DDIT3 50cmax
raw.data[ raw.data$plateID == "2013_10_17_nii_DDIT3", "cmax"] <-  "50cmax"  #DDIT3 50cma

raw.data[ raw.data$plateID == "2013_10_17_niii", "replID"] <-       "replicate_3"  #DDIT3 50 cmax
raw.data[ raw.data$plateID == "2013_10_17_niii", "cmax"] <-       "50cmax"  #DDIT3 50 cm

raw.data[ raw.data$plateID == "2013_10_24_DDIT3", "replID"] <-      "replicate_1"  #DDIT3 10 cmax
raw.data[ raw.data$plateID == "2013_10_24_DDIT3", "cmax"] <-      "10cmax"  #DDIT3 10 cm

raw.data[ raw.data$plateID == "2013_11_20_DDIT3", "replID"] <-      "replicate_2"  #DDIT3 10 cmax
raw.data[ raw.data$plateID == "2013_11_20_DDIT3", "cmax"] <-      "10cmax"  #DDIT3 10 cm

raw.data[ raw.data$plateID == "2013_11_27_DDIT3", "replID"] <-      "replicate_1"  #DDIT3 5 cmax
raw.data[ raw.data$plateID == "2013_11_27_DDIT3", "cmax"] <-      "5cmax"  #DDIT3 5 cma

raw.data[ raw.data$plateID == "2013_12_04_DDIT3_niii", "replID"] <-      "replicate_1"  #DDIT3 1 cmax
raw.data[ raw.data$plateID == "2013_12_04_DDIT3_niii", "cmax"] <-      "1cmax"  #DDIT3 

raw.data[ raw.data$plateID == "2013_12_04_DDIT3_nii", "replID"] <-  "replicate_2"  #DDIT3 1 cmax
raw.data[ raw.data$plateID == "2013_12_04_DDIT3_nii", "cmax"] <-  "1cmax"  #DDIT3 1 cma

raw.data[ raw.data$plateID == "2013_12_11_DDIT3", "replID"] <-  "replicate_2"  #DDIT3 5 cmax
raw.data[ raw.data$plateID == "2013_12_11_DDIT3", "cmax"] <-  "5cmax"  #DDIT3 5 cmax

raw.data[ raw.data$plateID == "2014_02_05_nii_p21", "replID"] <-    "replicate_1" #p21 10cmax
raw.data[ raw.data$plateID == "2014_02_05_nii_p21", "cmax"] <-    "10cmax" #p21 10cma

raw.data[ raw.data$plateID == "2014_02_05_p21_niii", "replID"] <-   "replicate_2" #p21 10cmax
raw.data[ raw.data$plateID == "2014_02_05_p21_niii", "cmax"] <-   "10cmax" #p21 10cma

raw.data[ raw.data$plateID == "2014_02_11_p21", "replID"] <-        "replicate_1" #p21 50cmax
raw.data[ raw.data$plateID == "2014_02_11_p21", "cmax"] <-        "50cmax" #p21 50cma

raw.data[ raw.data$plateID == "2014_02_18_p21", "replID"] <-        "replicate_1" #p21 5cmax 
raw.data[ raw.data$plateID == "2014_02_18_p21", "cmax"] <-        "5cmax" #p21 5cmax

raw.data[ raw.data$plateID == "2014_02_19_p21", "replID"] <-        "replicate_2" #p21 5cmax
raw.data[ raw.data$plateID == "2014_02_19_p21", "cmax"] <-        "5cmax" #p21 5cmax

raw.data[ raw.data$plateID == "2014_02_25_p21", "replID"] <-        "replicate_1" #p21 1cmax
raw.data[ raw.data$plateID == "2014_02_25_p21", "cmax"] <-        "1cmax" #p21 1cmax

raw.data[ raw.data$plateID == "2014_03_05_p21", "replID"] <-        "replicate_2" #p21 1cmax
raw.data[ raw.data$plateID == "2014_03_05_p21", "cmax"] <-        "1cmax" #p21 1cmax

sum(is.na(raw.data$replID)) #0
sum(is.na(raw.data$cmax)) #0
unique(raw.data$dose_uM)


#fix dosages: should be exactly the same digits if same concentration
raw.data<- as.data.table(raw.data)
raw.data$comp_dose <- paste(raw.data$treatment, raw.data$dose_uM, sep="_")

raw.data$dose_uM <- round(raw.data$dose_uM, digits = 4)

ind.more10 <- raw.data$dose_uM >= 10
raw.data$dose_uM[ind.more10] <- round(raw.data$dose_uM[ind.more10], digits = 1)
ind.more1 <- raw.data$dose_uM >= 1
raw.data$dose_uM[ind.more1] <- round(raw.data$dose_uM[ind.more1], digits = 2)
head(raw.data)
unique(raw.data[, dose_uM])



# tabel met compound namen en concentraties maken
raw.data[, comp_dose:= paste(treatment, dose_uM, sep ="_")]
compDosetabel <- unique(raw.data[!cell_line %in% "ICAM1", list(treatment, cmax, dose_uM, comp_dose)])

compDosetabel <- as.data.frame(compDosetabel)
compDosetabel<- compDosetabel[ order(compDosetabel$treatment, compDosetabel$dose_uM),]
# count how many doses each compound has and fix if incorrect

compDosetabel$comp_dose <- paste(compDosetabel$treatment, compDosetabel$dose_uM, sep="_")
comps <- unique(compDosetabel[ , c("treatment", "dose_uM")])[ order(unique(compDosetabel$comp_dose)),  ]
comps <- comps[ order(comps$treatment, comps$dose_uM),]
tt <- ddply(comps, .(treatment), summarize, length(treatment))
tt[ tt[, "..1"] >5,]

unique(raw.data[, cell_line])

#fix compDosetabel, daarna merge operatie ervan uitgaande dat de cmax'n kloppen
getwd()

incDose <- comps[ comps$treatment %in%  c("digoxin","epinephrine", "haloperidol", "moxisylyte",
                               "paroxetine","phenacetin", "propranolol","ribavirin",
                               "simvastatin","sulindac", "tacrine", "ticlopidine"), ]
incDose
compDosetabel[ compDosetabel$treatment == "digoxin" & compDosetabel$dose_uM == 0.2830, "dose_uM" ] <-0.2800
compDosetabel[ compDosetabel$treatment == "epinephrine" & compDosetabel$dose_uM == 0.1763, "dose_uM" ] <-0.1800
compDosetabel[ compDosetabel$treatment == "haloperidol" & compDosetabel$dose_uM == 0.5321, "dose_uM" ] <-0.53
compDosetabel[ compDosetabel$treatment == "moxisylyte" & compDosetabel$dose_uM == 15.7000, "dose_uM" ] <-15.8
compDosetabel[ compDosetabel$treatment == "paroxetine" & compDosetabel$dose_uM == 3.0400, "dose_uM" ] <-3.05
compDosetabel[ compDosetabel$treatment == "phenacetin" & compDosetabel$dose_uM == 30.8000, "dose_uM" ] <-31.0
compDosetabel[ compDosetabel$treatment == "propranolol" & compDosetabel$dose_uM == 10.1000, "dose_uM" ] <-10.0
compDosetabel[ compDosetabel$treatment == "ribavirin" & compDosetabel$dose_uM == 130.6000, "dose_uM" ] <-130.7
compDosetabel[ compDosetabel$treatment == "simvastatin" & compDosetabel$dose_uM == 4.1200, "dose_uM" ] <-4.1
compDosetabel[ compDosetabel$treatment == "sulindac" & compDosetabel$dose_uM == 1599.2000, "dose_uM" ] <-1599.3
compDosetabel[ compDosetabel$treatment == "tacrine" & compDosetabel$dose_uM == 3.9800, "dose_uM" ] <-4
compDosetabel[ compDosetabel$treatment == "ticlopidine" & compDosetabel$dose_uM ==  403.7000, "dose_uM" ] <- 403.8
# wat niet klopt: 
#digoxin 
#epinephrine
#haloperidol
#moxisylyte
#paroxetine
#phenacetin
#propranolol
#ribavirin
#simvastatin
#sulindac
#tacrine
#ticlopidine

raw.data$comp_dose <- NULL
# normalization considerations
  # normalize over entire set per cell line?
  # or use the above DMSO background?
  # definitely don't use the per-plate min max normalized...
  # min max normalization over multiple plates will not work.. need per plate normalization

# --> first try the above DMSO background counts
# --> works fine, but cell number and cell speed normaliztion is not. Normalize per plate relative to DMSO wells per time point 

# collapse time to make dose response curves
# 1) max based
# 2) final time point based


head(compDosetabel)
compDosetabel<-as.data.table(compDosetabel)
compDosetabel[, comp_dose:=NULL]
compDosetabel<-unique(compDosetabel)

setkeyv(compDosetabel, c("treatment", "cmax") )
setkeyv(raw.data, c("treatment", "cmax"))

dim(raw.data)
compDosetabel<- unique(compDosetabel)
raw.data.test <- raw.data[compDosetabel]
dim(raw.data.test)


raw.data<-raw.data.test
# lets have a look at the cmax dose_uM levels with respect to each other
tail(raw.data)
raw.data[, dose_uM:=NULL]
setnames(raw.data, old="i.dose_uM", new = "dose_uM")
all.t<- unique(raw.data$treatment)
for(i in seq_along(all.t)){
print(dim(unique(raw.data[ raw.data$treatment == all.t[i], list(dose_uM,cmax) ])))
} # 5 2 for all?




unique(raw.data$variable)
require(stringr)
########lol
raw.data[ , variable:= gsub("count_Nuclei_Intensity_IntegratedIntensity_Image_GFP_larger_7.5.1",
                            "count_Nuclei_Intensity_IntegratedIntensity_Image_GFP_larger_7.51", variable)]

unique(raw.data[, variable])
raw.data[, ThreshN:=NULL]
raw.data[, ThreshN:= str_match(variable,"_larger_([0-9 .]{2,5})$")[,2]]

unique(raw.data[, ThreshN])

raw.data[, ThreshN:= as.numeric(ThreshN)]
unique(raw.data$ThreshN)
head(test)
test<- ddply(raw.data, .(plateID), transform, isMin= ThreshN == min(ThreshN, na.rm=TRUE))

test<- ddply(test, .(plateID), transform, isMax= ThreshN == max(ThreshN, na.rm=TRUE))

test<- ddply(test, .(plateID), transform, isMiddle= ThreshN != max(ThreshN, na.rm=TRUE) & ThreshN!=min(ThreshN, na.rm=TRUE))

####
#2013-08-28 nii eruit vanwege duidelijke fouten, nog uitzoeken wat er aan de hand was
test <- test[ test$plateID!="2013_08_28_nii_Srxn1",]

test[1:300,]


#change variables to consistent names

test<-as.data.table(test)
test[, variable:= gsub("obj_nc_", "Nuclei_", variable) ]
test[, variable:= gsub("_img_gfp", "_Image_GFP", variable) ]
  


unique(test$variable)

head(test)

as.data.frame(unique(test[, list(plateID,variable, ThreshN, isMin,isMax, isMiddle)]))


test[ isMin ==TRUE  & !grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_[0-9 .]{3,5}", 
                        "count_Cyto.2m",variable)]
test[ isMin ==TRUE  & grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Nuclei_Intensity_IntegratedIntensity_Image_GFP_larger_[0-9 .]{3,5}", 
                                                                                  "count_Nuclei.2m",variable)]


test[ isMiddle ==TRUE  & !grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_[0-9 .]{3,5}", 
                                                                                  "count_Cyto.3m",variable)]
test[ isMiddle ==TRUE  & grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Nuclei_Intensity_IntegratedIntensity_Image_GFP_larger_[0-9 .]{3,5}", 
                                                                                 "count_Nuclei.3m",variable)]

test[ isMax ==TRUE  & !grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_[0-9 .]{3,5}", 
                                                                                     "count_Cyto.m3sd",variable)]
test[ isMax ==TRUE  & grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Nuclei_Intensity_IntegratedIntensity_Image_GFP_larger_[0-9 .]{3,5}", 
                                                                                    "count_Nuclei.m3sd",variable)]
is.data.table(test)
unique(test$variable)

# nog aanpassen
test[, variable:= gsub("count_Cyto", "GFP_pos", variable)]
test[, variable:= gsub("count_Nuclei", "GFP_pos", variable)]
# maak grafiek, 3 cellijnen op 1 plot, eerst 3X bckground proberen


unique(test$variable)


test[, variable:= gsub("Nuclei_TrackObjects_DistanceTraveled_[0-9]{2}",
                       "CellSpeed", variable)]

test[, variable:= gsub("imageCountTracked",
                       "CellNumber", variable)]

unique(test[, variable])

# make wide
# calculate norm
# make long

test[, ThreshN:=NULL]
test[, isMin:=NULL]
test[, isMax:=NULL]
test[, isMiddle:=NULL]
test[, valueControl:=NULL]

unique(test$variable)


orig srxn1, chop & p21 summarized data qc data

# te maken plotjes
# een time course plot voor 1 per cell line gfp feature
# een heatmap time course plot voor 1 per cell lijn gfp features
# een heatmap voor time-viability features

# een dose response voor gfp features + cell dood (PI/Annexin fracties)
# een dose response heatmap voor viability featues
#####################
### === Nu eerst alle time course data plaat-normalizeren === ###
# normalizatie hetzelfde als ICAM1 data: per tijdpunt plaat-specifieke DMSO waarden eraf halen
# voor de time course plot is alleen de 2X boven belangrijk en nog 1 kiezen voor ICAM1
unique(test$fingerprints)


test[, fingerprints:=variable]
  
 
  test[variable %in% "GFP_pos.m3sd", fingerprints:= paste(variable, cell_line, sep = "_")]
  test[variable %in% "GFP_pos.2m", fingerprints:= paste(variable, cell_line, sep = "_")]
  test[variable %in% "GFP_pos.3m", fingerprints:= paste(variable, cell_line, sep = "_")]
  
  
  test[variable %in% "Cytoplasm_Intensity_IntegratedIntensity_Image_GFP", fingerprints:= paste(variable, cell_line, sep = "_")]
  test[variable %in% "Cytoplasm_Intensity_MeanIntensity_Image_GFP", fingerprints:= paste(variable, cell_line, sep = "_")]
  test[variable %in% "Nuclei_Intensity_IntegratedIntensity_Image_GFP", fingerprints:= paste(variable, cell_line, sep = "_")]
  test[variable %in% "Nuclei_Intensity_MeanIntensity_Image_GFP", fingerprints:= paste(variable, cell_line, sep = "_")]
  
  unique(test[,fingerprints])
  # remove in GUI normalized features
  test <- test[!variable %in%  c("min_maxNorm_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP",
                                 "min_maxNorm_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP",
                                 "min_maxNorm_CellNumber",
                                 "min_maxNorm_CellSpeed",
                                 "min_maxNorm_Nuclei_Intensity_IntegratedIntensity_Image_GFP")]
  
 
  
  
  test[, variable:=NULL]
  test[, fingerprints:= gsub("Nuclei_Intensity_MeanIntensity_img_nc", 
                             "Nuclei_Intensity_MeanIntensity_Image_Hoechst", fingerprints)]
  
  varsToNorm <- c( "Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1",
                  "Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1",
                  "Nuclei_Intensity_MeanIntensity_Image_Hoechst", 
                  "Nuclei_Intensity_IntegratedIntensity_Image_GFP_DDIT3", "Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3",
                  "Nuclei_Intensity_IntegratedIntensity_Image_GFP_p21", "Nuclei_Intensity_MeanIntensity_Image_GFP_p21"
  )

# ja wel DMSO eraftrekken ipv fold change. dan min-max
# normalization function: normalize based on internal control then scale between 0-1
# oude norm functie gebruikien? of anders proberen zonder fold change? Of DMSO waarden eraf?

  test[, dose_uM:= as.numeric(as.character(test[, dose_uM]))]# change factor to numeric of dose 

  
 
 
  options(width = 120)
 checkIt <- test[, max(value), by = c( "fingerprints", "cmax")]
 checkIt[order(fingerprints, cmax),]
  
  
  write.table(test, file = "bkupTest24_04_2016.txt", col.names=TRUE, sep ="\t")

  test <- read.delim(file = "C:/Users/steve_000/Documents/work/DILI paper/DILIpaper/data/bkupTest24_04_2016.txt", sep ="\t", header = TRUE)
  test<-data.table(test)
  
  
  # first add the nieuw cmax data - then do the normalization per replicate (new plan 23/07/2016)
  
  
  # dataIn <- test
  # normVar <- "Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1" 
  oldreporterdata
  #normVar<- "Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1"
  unique(test$treatment)
  unique(test[ treatment =="dmso",list(plateID, replID)])
  #dataIn <- test
  normFunAll <- function(dataIn, normVar) {
  
  buffer <- dataIn[ fingerprints%in% normVar] # select 1 variable
  #subtract DMSO per plate followed by min max multiple plates within replicates
  meanDMSO <-  buffer[ treatment %in% "dmso"  ,  
                     mean(value, na.rm=TRUE),
                     by = c("plateID", "timeID")
                     ]
 # buffer[is.na(value)]
#  meanDMSO[is.na(V1)]
  
  setkeyv(meanDMSO, c("plateID", "timeID"))
  setkeyv(buffer, c("plateID" ,"timeID"))
  
  # iets heeft geen dmso
 # a<-unique(buffer[ , list(plateID, replID, timeID) ])
#  b<-unique(buffer[ treatment == "dmso" , list(plateID, replID, timeID) ])
#  a[, poep:= paste(plateID,replID,timeID)]
#  b[, poep:= paste(plateID,replID,timeID)]
#  setdiff(a$poep,b$poep)
  buffer<- meanDMSO[buffer]
 # meanDMSO[buffer][is.na(V1)]
  
  buffer[, DMSO_sub.value:= value - V1 ]
  
  buffer[, value:=NULL]
  buffer[, fingerprints:= paste("mmnDMSO.sub", normVar, sep ="_")]
  
  
  buffer[, minValue:= min(DMSO_sub.value, na.rm= TRUE), by = replID]
  buffer[, maxValue:= max(DMSO_sub.value, na.rm= TRUE), by = replID]
  
  buffer[, mmnvalue:= (DMSO_sub.value - minValue)/ 
           (maxValue-minValue)]
  
  setnames(buffer, "mmnvalue", "value")
  buffer[, DMSO_sub.value:=NULL]
  buffer[, V1:=NULL]
  buffer[, maxValue:=NULL]
  buffer[, minValue:=NULL]
  output <- buffer
    return(output)  
}



# #adv  cell count kan helemaal weg uit dataset:
# doxycycline 1125.7  voor srxn1 en ddit3 en p21
# doxycycline 562.9 adv ddit3
# doxycycline 112.6 adv ddit3 (minder erg)
# doxycyline helemaal eruit wellicht?
# 
# tetracycline 2096.3 voor srxn1
# tetracycline 1048.2 voor ddit3
# aflatoxin B1 100 voor Srxn1
# dim(test)


#doxycycline eruit
test<- test[!(treatment %in% "doxycycline" ) ]
# tetracycline 2096.3 en 1048.2 eruit
test<- test[!(treatment %in% "tetracycline" & dose_uM %in% 2096.3) ]
test<- test[!(treatment %in% "tetracycline" & dose_uM %in% 1048.2) ]

test<- test[!(treatment %in% "aflatoxin B1" & dose_uM %in% 100) ]

test[ fingerprints %in% "Nuclei_AreaShape_Area" & value > 500, list(value,treatment, cmax, timeID, cell_line)]


unique(test[, fingerprints])

test<- test[!(treatment %in% "aflatoxin B1" & dose_uM %in% 100) ]
write.table(test, file =  "data/notNorm_oldreporterdata.txt", sep = "\t", col.names = TRUE)
# nu opnieuw normalizeren

# heb ook de niet genormaliseerde nieuwe cmax ( verder op in script...), geanalyseerd en weggeschreven
niewcmaxdata <- read.delim( file = "C:/Users/steve_000/Documents/work/DILI paper/DILIpaper/data/notNormalizedNieuweCmax.txt", sep = "\t", header = TRUE)
niewcmaxdata<- data.table(niewcmaxdata)
oldreporterdata <- read.delim( file =  "C:/Users/steve_000/Documents/work/DILI paper/DILIpaper/data/notNorm_oldreporterdata.txt", sep = "\t", header = TRUE)
oldreporterdata<- data.table(oldreporterdata)



colnames(oldreporterdata)
colnames(niewcmaxdata)
niewcmaxdata$variable <- NULL

# normalize, then pull appart again, to be able to keep old code
nieuwPlates <- unique(niewcmaxdata$plateID)

# nu iets verzinnen met die 5 replicates van srxn1, welke zijn goed? welke combineren? en chop"





test <- rbind(oldreporterdata, niewcmaxdata)
test <- test[ !is.na(value)]


# reduceer aantal replicaten
2013_07_29_Srxn1 replicate_3 # 100 cmax 182 comps  
2013_08_01_Srxn1 replicate_4 # 100 cmax  182 comps keep (make repl 1 )
2013_08_14_Srxn1 replicate_5 # 100 cmax # 176 comps keep (make repl 2)
2013_08_22_Srxn1 replicate_3 # 50 cmax 183 compounds
2013_09_04_Srxn1 replicate_3 # 10 cmax 184 compounds only 2 reps change 3 to 1 ( repl 1 was removed due to erronous looking data)

2013_10_17_niii replicate_3 # CHOP 183 50 cmax (make repl 1)

unique( test[plateID == "2013_10_17_niii" & replID == "replicate_3", list(treatment, cmax, replID)])
unique( test[cell_line == "DDIT3" & cmax == "50cmax", list( cmax, replID, plateID)])

# remove: 2013_07_10_Srxn1_z1 replicate_1 (86 comps 100 cmax)
# remove 2013_07_24_Srxn1 replicate_2 100cmax 123 compounds
# remove 2013_07_29_Srxn1 replicate_3 100 cmax 
# remove 2013_08_22_Srxn1 replicate_3  50 cmax
# remove 2013_10_04_DDIT3 replicate_1 50 cmax 

test <- test[ plateID != "2013_07_10_Srxn1_z1",]
test <- test[ plateID != "2013_07_24_Srxn1",]
test <- test[ plateID != "2013_07_29_Srxn1",]
test <- test[ plateID != "2013_08_22_Srxn1",]
test <- test[ plateID != "2013_10_04_DDIT3",]

test[plateID == "2013_08_01_Srxn1", replID:= "replicate_1" ]
test[plateID == "2013_08_14_Srxn1", replID:= "replicate_2" ]
test[plateID == "2013_09_04_Srxn1", replID:= "replicate_1" ]
test[plateID == "2013_10_17_niii", replID:= "replicate_1" ]


unique(test[,  list(plateID, replID)])

#pull appart:

niewcmaxdata <- test[plateID %in% nieuwPlates, ]
oldreporterdata <- test[ !plateID %in% nieuwPlates, ]

# check if each cell line and dose has two replicates

unique(niewcmaxdata[,  list(plateID, cmax, replID)])
unique(oldreporterdata[,  list(plateID, cmax, replID)])

# where is repl2 of p21 50 cmax?? : lost when external HD
unique(test$fingerprints)

test <- test[!grepl("mmnDMSO", fingerprints)]



mmn.list = list()

unique(test[, fingerprints])
test[, treatment:= tolower(treatment)]
for( i in seq_along(varsToNorm)){
  mmn.list[[i]] <- normFunAll(dataIn = test, normVar = varsToNorm[i] )
}

mmn.df <- do.call("rbind", mmn.list)
mmn.df <- as.data.table(mmn.df)
mmn.df <- mmn.df[!is.na(value)]

test <- rbind(test, mmn.df)

#pull appart:


niewcmaxdata <- test[plateID %in% nieuwPlates, ]
oldreporterdata <- test[ !plateID %in% nieuwPlates, ]



write.table(test, "data/norm_origandnew_reporters.txt", col.names= TRUE, sep ="\t")
dir()
# eens kijken hoe dit gegaan is:
fingerprintVars <- unique(test[, fingerprints])
# om ook even nuclei size en cellnumber per cell line te zien even dummy variabele aanmaken: omdat outliers verwijderen adv deze plotjes
test[, fingerprintDummy:= paste(fingerprints, cell_line, sep ="_")]
fingerprintVarsD <- unique(test[, fingerprintDummy])
writepdfs <- "H:/DILI screen/meta analyse DILI screen/graphs/graphs febr2016/time/"
# plot each variable seperately:
for( i in seq_along(fingerprintVarsD)){
  
  plotData <- test[ fingerprintDummy %in% fingerprintVarsD[i]]
  plotData[, treat.dose:= paste(treatment, dose_uM)]
  
  
  p <- ggplot(data =plotData, aes(x= timeAfterExposure, y=value)) + facet_wrap(~treat.dose, scales = "free_x")
  p<-p + geom_point(aes(color=replID, shape = replID)) +geom_line(aes(group=replID, color=replID)) + theme_sharp()
  pdf(file = paste(writepdfs, fingerprintVarsD[i], ".pdf", sep =""), width =60, height = 60)
  print(p)
  dev.off()
}
#save.image("03082016.Rdata")
# aan de hand van deze grafieken outliers verwijderen ( zie bijv cellnumber -- enorme outliers)
# vergeet niet ook modeled time course data op te slaan voor machine learning

 load("20160208.Rdata.RData")

# waar staan we nu? :
 
 test # orig reporter response data (weggeschreven. test wordt overschreven door nieuwe cmax later in script)
 new.cmaxICAM1.long # upgedate cmax ICAM1 data
 icam.orig.all.mean.d.long # originele ICAM1 data
 
 
 #nieuwe Srxn1, Chop en DDIT3 moet nog
 
 # celdood data moet nog
  # daarna modelen en illustraties
# nieuwe Srxn1 Chop en DDIT3 data. Zelfde als oude maar even oppasen met min-max normalizatie ivm meerdere concentraties op 1 plaat (dus per plaat doen)
 
 ###======||=====
 # hier
 #rm(list=ls())
 rm("test")
 path.to.nonICAM <- ("E:/DILI screen/meta analyse DILI screen/nieuwe cmax data/summaryFilesNieuwCmax")
 
 setwd(path.to.nonICAM)
 
 raw.data =list()
 dir.files <- dir(path.to.nonICAM)
 dir.files <- dir.files[grepl("(Summary Data)", dir.files)]
 
 for (i in seq_along(dir.files)){
   raw.data[[i]] <- as.data.table(read.table(file = paste(path.to.nonICAM, dir.files[i], sep ="/"), header = TRUE, sep ="\t"))
 }
 head(raw.data[[1]])
  
 
 
 unique(raw.data[[7]]$plateID)
 raw.data[[7]][, plateID:= "20150930_SRXN1" ]
 #modify time
 
  timeBetweenFrames <- "01:35:00"
 exposureDelay <- "16:05"
 
 timeBetweenFrames <- round(as.integer(strftime(strptime(timeBetweenFrames, format = "%H:%M:%S"), "%H")) + 
                              1/60 * as.integer(strftime(strptime(timeBetweenFrames, 
                                                                  format = "%H:%M:%S"), "%M")) +
                              1/3600 * as.integer(strftime(strptime(timeBetweenFrames, 
                                                                    format = "%H:%M:%S"), "%S"))
                            , digit =2 )
 
 exposureDelay <- round(as.integer(strftime(strptime(exposureDelay, format = "%H:%M"), "%H")) + 
                          1/60 * as.integer(strftime(strptime(exposureDelay, 
                                                              format = "%H:%M"), "%M")), digit =1 )
 
 raw.data[[7]][, timeAfterExposure:=round(as.numeric(timeID)*timeBetweenFrames+exposureDelay - timeBetweenFrames, digits = 1) ]
 
 unique((raw.data[[7]][, timeAfterExposure])) 
 unique((raw.data[[7]][, timeID])) 
 
 ((raw.data[[7]][, timeID:=as.numeric(timeID)+4])) 
 
 unique(raw.data[[6]][, plateID])
 raw.data[[1]]$replID <- "replicate_1"
 raw.data[[2]]$replID <- "replicate_2"
 raw.data[[3]]$replID <- "replicate_1"
 raw.data[[4]]$replID <- "replicate_2"
 raw.data[[5]]$replID <- "replicate_1"
 raw.data[[6]]$replID <- "replicate_2"
 raw.data[[7]]$replID <- "replicate_2"
 
 raw.data <- do.call("rbind" ,raw.data)
 head(raw.data)
 
 raw.data$control <-NULL
 head(raw.data)
 unique(raw.data$variable)
 raw.data$variable <- gsub("_>_", "_larger_", raw.data$variable)
 raw.data$variable <- gsub("_<_", "_smaller_", raw.data$variable)
 
 # fix cell line names
 unique(raw.data$cell_line)
 raw.data$cell_line <- gsub("CHOP", "DDIT3", raw.data$cell_line)
 raw.data$cell_line <- gsub("SRXN1", "Srxn1", raw.data$cell_line)
 raw.data$cell_line <- gsub("P21", "p21", raw.data$cell_line)
 unique(raw.data$variable)
 
 
 
 # moet worden  _image_hoechst  --> _Image_Hoechst  _image_GFP  -->  _Image_GFP
 raw.data[ , variable:= gsub("_image_hoechst", "_Image_Hoechst", variable  )]
 raw.data[ , variable:= gsub("_image_GFP", "_Image_GFP", variable  )]
  
 
 #fix dosages: should be exactly the same digits if same concentration
 raw.data<- as.data.table(raw.data)
 raw.data$comp_dose <- paste(raw.data$treatment, raw.data$dose_uM, sep="_")
 
 raw.data$dose_uM <- round(raw.data$dose_uM, digits = 4)
 
 ind.more10 <- raw.data$dose_uM >= 10
 raw.data$dose_uM[ind.more10] <- round(raw.data$dose_uM[ind.more10], digits = 1)
 ind.more1 <- raw.data$dose_uM >= 1
 raw.data$dose_uM[ind.more1] <- round(raw.data$dose_uM[ind.more1], digits = 2)
 head(raw.data)
 count(as.data.frame(unique(raw.data[, list(treatment, dose_uM) ]))[  order(as.data.frame(unique(raw.data[, list(treatment, dose_uM) ]))$treatment)     , ]$treatment)
 
 # create cmax' s
 
 compdose.table<-as.data.frame(unique(raw.data[, list(treatment, dose_uM) ]))
 compdose.table<-compdose.table[ order(compdose.table$treatment, compdose.table$dose_uM), ]
 compdose.table <- compdose.table[ compdose.table$treatment !="DMSO",]
 
 compdose.table$cmax <- rep( c("1cmax", "5cmax", "10cmax", "25cmax", "50cmax", "100cmax"),length(unique(compdose.table$treatment)) )
 setkeyv(raw.data, c("treatment", "dose_uM"))
 compdose.table<- as.data.table(compdose.table)
 setkeyv(compdose.table, c("treatment", "dose_uM"))
 
 raw.data <-compdose.table[raw.data]
 unique(raw.data[is.na(cmax), treatment])
 # excellent =)(=
 raw.data$comp_dose <- NULL
 
 
 
 
 
 require(stringr)
 unique(raw.data[, variable])
 raw.data[, ThreshN:=NULL]
 raw.data[, ThreshN:= str_match(variable,"_larger_([0-9 .]{2,6})$")[,2]]
 
 unique(raw.data[, ThreshN])
 # die plaat crash molt weer met script ff snelle oplossing zodat t weer werkt.
 
 raw.data[ , variable:= gsub("count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_8.074", "count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_9.582", variable)]
 raw.data[ , variable:= gsub("count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_12.111", "count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_14.373", variable)]
 raw.data[ , variable:= gsub("count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_38.816", "count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_34.524", variable)]
 raw.data[, ThreshN:= as.numeric(ThreshN)]
 unique(raw.data$ThreshN)
 
 test<- ddply(raw.data, .(plateID), transform, isMin= ThreshN == min(ThreshN, na.rm=TRUE))
 test<- ddply(test, .(plateID), transform, isMax= ThreshN == max(ThreshN, na.rm=TRUE))
 test<- ddply(test, .(plateID), transform, isMiddle= ThreshN != max(ThreshN, na.rm=TRUE) & ThreshN!=min(ThreshN, na.rm=TRUE))
 
  test[1:300,]
 
 test<-as.data.table(test)
 
 as.data.frame(unique(test[, list(plateID,variable, ThreshN, isMin,isMax, isMiddle)]))
 
 
 test[ isMin ==TRUE  & !grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_[0-9 .]{3,6}", 
                                                                                   "count_Cyto.2m",variable)]
 test[ isMin ==TRUE  & grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Nuclei_Intensity_MeanIntensity_Image_GFP_larger_[0-9 .]{3,6}", 
                                                                                  "count_Nuclei.2m",variable)]
 
 
 test[ isMiddle ==TRUE  & !grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_[0-9 .]{3,6}", 
                                                                                      "count_Cyto.3m",variable)]
 test[ isMiddle ==TRUE  & grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Nuclei_Intensity_MeanIntensity_Image_GFP_larger_[0-9 .]{3,6}", 
                                                                                     "count_Nuclei.3m",variable)]
 
 test[ isMax ==TRUE  & !grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_larger_[0-9 .]{3,6}", 
                                                                                   "count_Cyto.m3sd",variable)]
 test[ isMax ==TRUE  & grepl("count_Nuclei_Intensity",variable) , variable:= gsub("count_Nuclei_Intensity_MeanIntensity_Image_GFP_larger_[0-9 .]{3,6}", 
                                                                                  "count_Nuclei.m3sd",variable)]
 is.data.table(test)
 unique(test$variable)
 
 
 test[, variable:= gsub("count_Cyto", "GFP_pos", variable)]
 test[, variable:= gsub("count_Nuclei", "GFP_pos", variable)]
 # maak grafiek, 3 cellijnen op 1 plot, eerst 3X bckground proberen
 
 
 unique(test$variable)
 
 
 
 unique(test[, variable])
 
 # make wide
 # calculate norm
 # make long
 
 test[, ThreshN:=NULL]
 test[, isMin:=NULL]
 test[, isMax:=NULL]
 test[, isMiddle:=NULL]
 
 
 unique(test$variable)
 
 
 
 test[ , fingerprints:= variable ]
 test[variable %in% "GFP_pos.m3sd", fingerprints:= paste(variable, cell_line, sep = "_")]
 test[variable %in% "GFP_pos.2m", fingerprints:= paste(variable, cell_line, sep = "_")]
 test[variable %in% "GFP_pos.3m", fingerprints:= paste(variable, cell_line, sep = "_")]
 
 
 test[variable %in% "Cytoplasm_Intensity_IntegratedIntensity_Image_GFP", fingerprints:= paste(variable, cell_line, sep = "_")]
 test[variable %in% "Cytoplasm_Intensity_MeanIntensity_Image_GFP", fingerprints:= paste(variable, cell_line, sep = "_")]
 test[variable %in% "Nuclei_Intensity_IntegratedIntensity_Image_GFP", fingerprints:= paste(variable, cell_line, sep = "_")]
 test[variable %in% "Nuclei_Intensity_MeanIntensity_Image_GFP", fingerprints:= paste(variable, cell_line, sep = "_")]
 
 
 unique(test[,fingerprints])
 
 
 
 #test[, variable:=NULL]
 
 varsToNorm <- c("Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3", "Nuclei_Intensity_MeanIntensity_Image_Hoechst",
                 "Nuclei_Intensity_MeanIntensity_Image_GFP_p21",
                 "Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1", "Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1"
 )
 
 # ja wel DMSO eraftrekken ipv fold change. dan min-max
 # normalization function: normalize based on internal control then scale between 0-1
 # oude norm functie gebruikien? of anders proberen zonder fold change? Of DMSO waarden eraf?
 
 test[, dose_uM:= as.numeric(as.character(test[, dose_uM]))]# change factor to numeric of dose 
 
 
 #14-02-2016
 #save.image("21072016.RData")
 #getwd()
 
 #"I:/DILI screen/meta analyse DILI screen/nieuwe cmax data/summaryFilesNieuwCmax"
 # min max schaling naar per-plaat aanpassen is next
 
 
 
 unique(test[, fingerprints])
 # dataIn <- test
 #normVar <- "Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1"
 
 
 test[treatment =="DMSO", cmax:= "100cmax"]
 
 
 normFun <- function(dataIn, normVar) {
   
   buffer <- dataIn[ fingerprints%in% normVar] # select 1 variable
   #subtract DMSO per plate-time followed by min max single plate within replicates
   meanDMSO <-  buffer[ treatment %in% "DMSO"  ,  
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
   # min and max per plate (because multiple dosages per plate in this case)
   minValue <- buffer[ , min(DMSO_sub.value, na.rm = TRUE), by = plateID] 
   maxValue <- buffer[ , max(DMSO_sub.value, na.rm = TRUE), by = plateID]
   maxValue[, V1:=V1+0.2] # add constant term, to normalize for the difference between original dose-per-plate and all dose on 1 plate.
   # likely the minmax-normalization assumption of equal max value is not valid, beacause only 12 new cmax compounds.
   #0.1 is the difference in mean of these plates (see a few lines bellow the normFun). 2* 0.1 for a 0.1 in mean
   
   setkey(buffer, "plateID")
   setkey(minValue, "plateID")
   setkey(maxValue, "plateID")
   
   buffer <- buffer[minValue]
   buffer <- buffer[maxValue]
   
   setnames(buffer, "V1", "minValue")
   setnames(buffer, "i.V1", "maxValue")
   
   buffer[, mmnvalue:= (DMSO_sub.value - minValue)/ 
            (maxValue-minValue)]
   
   setnames(buffer, "mmnvalue", "value")
   buffer[, DMSO_sub.value:=NULL]
   buffer[, minValue:=NULL]
   buffer[, maxValue:=NULL]
   output <- buffer
   return(output)  
 } # heb script aangepast voor nieuwe cmax data
 
 
 
 # #adv  cell count kan helemaal weg uit dataset:
 # doxycycline 1125.7  voor srxn1 en ddit3 en p21
 # doxycycline 562.9 adv ddit3
 # doxycycline 112.6 adv ddit3 (minder erg)
 # doxycyline helemaal eruit wellicht?
 # 
 # tetracycline 2096.3 voor srxn1
 # tetracycline 1048.2 voor ddit3
 # aflatoxin B1 100 voor Srxn1
 # dim(test)
 
 
 #doxycycline eruit
 test<- test[!(treatment %in% "doxycycline" ) ]
 # tetracycline 2096.3 en 1048.2 eruit
 test<- test[!(treatment %in% "tetracycline" & dose_uM %in% 2096.3) ]
 test<- test[!(treatment %in% "tetracycline" & dose_uM %in% 1048.2) ]
 
 test<- test[!(treatment %in% "aflatoxin B1" & dose_uM %in% 100) ]
 
 test[ fingerprints %in% "Nuclei_AreaShape_Area" & value > 500, list(value,treatment, cmax, timeID, cell_line)]
 
 
 unique(test[, fingerprints])
 
 test<- test[!(treatment %in% "aflatoxin B1" & dose_uM %in% 100) ]
 # write on 26/07/2016
 write.table(test, file = "C:/Users/steve_000/Documents/work/DILI paper/DILIpaper/data/yesNormalizedOldandNieuweCmax.txt", sep = "\t", col.names = TRUE)
 # nu opnieuw normalizeren
  # status 21-07-2016: normalizatie al stukken beter. maar srxn1 int int nog te hoog. chop te laag. PI is goed genoeg.
 # nieuwe status: min max per replicate, inclusief de nieuwe cmax.
 mmn.list = list()
 for( i in seq_along(varsToNorm)){
   
   mmn.list[[i]] <- normFun(dataIn = test, normVar = varsToNorm[i] )
 }
 
 mmn.df <- do.call("rbind", mmn.list)
 
 
 head(mmn.df)
 
 
 
 # check distribution of normalized new-cmax data in comparison with rest of data
 mmn.df <- data.table(mmn.df)
 # Srxn int int nog niet goed. doe appart per cell line ook
 mean(test[ !grepl("Image_Hoechst", fingerprints), mean(value, na.rm = TRUE), by = list(plateID, fingerprints)]$V1)
 0.04544379 0.05407546 0.14137718 0.24612652 0.16187806 0.18455411 0.15512363 0.20427306
 0.1491065
 ##==##==##==##
 unique(mmn.df[, treatment])
 
 
 
 origrep.data <- read.delim( file ="E:/suz icam1/norm_orig_reporters.txt", sep = "\t", header=TRUE)
 origrep.data<- data.table( origrep.data )
 
 mean(origrep.data[ grepl("mmnDMSO", fingerprints) &
                 !grepl("Image_Hoechst", fingerprints),
                  mean(value, na.rm = TRUE), by = list(replID, fingerprints)]$V1)
 [1] 0.033642297 0.015468747 0.016354061 0.133233207 0.125315131 0.082739815 0.058331765 0.045914619 0.113831399 0.127299113 0.028474811 0.012734851
 [13] 0.026343123 0.013973739 0.009071702 0.016034666 0.062119731 0.049938942 0.052094196 0.063628847
 
 # grand mean: 0.05432724
 
 
 test <- rbind(test, mmn.df)
 
 test[, cmax:= factor(cmax, levels = c("1cmax", "5cmax", "10cmax","50cmax","100cmax"))]
 #remove 25cmax
 test <- test[!is.na(cmax),]
 
 pdf(file = "testNewCmaxNorn_.pdf", width = 20, height = 30 )
 ggplot(data = test[ fingerprints %in% c("Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1"
                                         )], aes(x= timeAfterExposure, y = value, color = fingerprints)) + geom_point() +
   facet_grid(treatment~cmax)
 dev.off()
 
 
 dir("E:/DILI screen/meta analyse DILI screen/nieuwe cmax data")
 write.table(test, "E:/DILI screen/meta analyse DILI screen/nieuwe cmax data/norm_nieuweCMAX_reporters21072016.txt", col.names= TRUE, sep ="\t")
 
 
 
 
 # eens kijken hoe dit gegaan is:
 fingerprintVars <- unique(test[, fingerprints])
 # om ook even nuclei size en cellnumber per cell line te zien even dummy variabele aanmaken: omdat outliers verwijderen adv deze plotjes
 test[, fingerprintDummy:= paste(fingerprints, cell_line, sep ="_")]
 fingerprintVarsD <- unique(test[, fingerprintDummy])
 
 writepdfs <- "I:/DILI screen/meta analyse DILI screen/graphs/graphs april2016/time_newcmax/"
 
 # plot each variable seperately:
 
 
 
 for( i in seq_along(fingerprintVarsD)){
   
   plotData <- test[ fingerprintDummy %in% fingerprintVarsD[i]]
   plotData[, treat.dose:= paste(treatment, dose_uM)]
   
   
   p <- ggplot(data =plotData, aes(x= timeAfterExposure, y=value)) + facet_wrap(~treat.dose, scales = "free_x")
   p<-p + geom_point(aes(color=replID, shape = replID)) +geom_line(aes(group=replID, color=replID)) + theme_sharp()
   pdf(file = paste(writepdfs, fingerprintVarsD[i], ".pdf", sep =""), width =60, height = 60)
   print(p)
   dev.off()
 }
 
 ###======||=====
 

 # nieuw plan mbt celdood data.
 # -> alleen en enkel PI positive fraction + number of cells. (this is missing in the new-cmax data)
 # nieuwe cmax data check (eerste set reporters en ICAM1)
 # oude data check
 # eerste ICAM check
 
 # # eerste oude data doorlopen
 # # daarna nieuwe cmax en eerste-ICAM1 laden, opschonen en samenvoegen.
 
 #ICAM1-eerst: I:\suz icam1\Summary files CD
 #niewe cmax: I:\DILI screen\meta analyse DILI screen\nieuwe cmax data\CDsummary
 # origninele reporters: I:/DILI screen/all cytotoxicity output/cytotox Summary Data.txt"
 
 # toxiciteit data orig reporters: E:\DILI screen\meta analyse DILI screen\nieuwe cmax data\summaryFilesNieuwCmax (cytotox_data_originalReporters.txt)


load("data_07april2016")
 # load cytotox data + verify with a plot
 cytotox.d <- read.delim(file ="E:/DILI screen/all cytotoxicity output/cytotox Summary Data.txt",
                         sep ="\t", header = T)
 
 #fix 2013_10_09_niii_p21 concentration
 cytotox.d<- as.data.table(cytotox.d)
 
 
 
 cytotox.d[ plateID%in%"2013_10_09_niii_p21" , dose_uM:=2*dose_uM ]
 
 cytotox.d <- cytotox.d[ !variable%in% "min_maxNorm_imageCountParentObj", ] # this normalization makes no sense
 cytotox.d <- cytotox.d[ !variable%in% "min_maxNorm_Nuclei_obj_AreaShape_Area", ] # this normalization makes no sense
 cytotox.d <- cytotox.d[ !variable%in% "min_maxNorm_Nuclei_obj_Intensity_MeanIntensity_Hoechst_Image", ] # this normalization makes no sense
 #change names
 cytotox.d[, variable:=gsub("imageCountParentObj", "cellCountCytotox", variable)]
 cytotox.d[, variable:=gsub("count_Annex_obj_AreaShape_Area.DIV.Nuclei_obj_AreaShape_Area_larger_0.02_", "apoptosis positive", variable)]
 cytotox.d[, variable:=gsub("count_Annex_obj_masked_primaryID_AreaShape_Area.DIV.Nuclei_obj_AreaShape_Area_larger_0.02_", "apoptosis positive pID", variable)]
 cytotox.d[, variable:=gsub("count_PI_obj_AreaShape_Area.DIV.Nuclei_obj_AreaShape_Area_larger_0.02_", "necrosis positive", variable)]
 cytotox.d[, variable:=gsub("count_PI_obj_masked_primaryID_AreaShape_Area.DIV.Nuclei_obj_AreaShape_Area_larger_0.02_", "necrosis positive pID", variable)]
 cytotox.d[, variable:=gsub("Nuclei_obj_Intensity_MeanIntensity_Hoechst_Image", "NucleiHoechstIntCytotox", variable)]
 cytotox.d[, variable:=gsub("Nuclei_obj_AreaShape_Area", "NucleiAreaCytotox", variable)]
 
 unique(cytotox.d[, variable])
 # keep: "cellCountCytotox"   "necrosis positive"
 
 cytotox.d <- cytotox.d[ variable %in% c("cellCountCytotox", "necrosis positive", "NucleiAreaCytotox","NucleiHoechstIntCytotox")]
 
 tox.vars <- unique(cytotox.d[, variable])
 setnames(cytotox.d, "variable", "fingerprints")
 
tox.vars



# normalization of cell number(to noisy elsewhise) requires alternative..
# not per plate due to dose differences per plate (this is the case for original cytotox data)


normFunCS <- # function for normalizing single time point data with first a fold change wrt plate specific DMSO. also includes outlier removal
  function(dataIn, normVar) {
    
    buffer <- dataIn[ fingerprints%in% normVar] # select 1 variable
    Checkcontrol <- buffer[ treatment %in% "DMSO"] # select the DMSO control
    controlV<-Checkcontrol[ ,mean(value, na.rm=TRUE), by = plateID]   # only single time point exist for this data
    setkey(controlV, "plateID")
    setkey(buffer, "plateID")
    buffer<- buffer[controlV]
    # remove outliers that are < 1.5 X DMSO values
    buffer <- buffer[ value < ( 1.5 * V1 ) ]
    
    buffer[, FC_value:= value/V1 ]
    buffer[, value:=NULL]
    buffer[, fingerprints:= paste("mmnFC", normVar, sep ="_")]
    minValue <- min(buffer[, FC_value], na.rm=TRUE)
    maxValue <- max(buffer[ , FC_value], na.rm=TRUE)
    buffer[, mmnvalue:= (FC_value - minValue)/ 
             (maxValue-minValue)]
    
    setnames(buffer, "mmnvalue", "value")
    buffer[, FC_value:=NULL]
    buffer[, V1:=NULL]
    
    test<- rbind(dataIn,buffer)
    return(test)  
  }

cytotox.d <- normFunCS(dataIn = cytotox.d, normVar =  tox.vars[1])
cytotox.d <- normFunCS(dataIn = cytotox.d, normVar =  tox.vars[2])
cytotox.d <- normFunCS(dataIn = cytotox.d, normVar =  tox.vars[3])
# not necrosis positive
cytotox.d.sel <- cytotox.d

# add cmax values to cytotox.d
cytotox.d.sel[, replID:=NA]


cytotox.d.sel[ cytotox.d.sel$plateID == "2013_07_10_Srxn1_z1", "replID"] <-   "replicate_1"  #Srxn1 100cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_07_10_Srxn1_z1", "cmax"] <-   "100cmax"  #Srxn1 100cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_07_24_Srxn1", "replID"] <-      "replicate_2"  #Srxn1 100cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_07_24_Srxn1", "cmax"] <-      "100cmax"  #Srxn1 100cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_07_29_Srxn1", "replID"] <-      "replicate_3"  #Srxn1 100cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_07_29_Srxn1", "cmax"] <-      "100cmax"  #Srxn1 100cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_01_Srxn1", "replID"] <-      "replicate_4"  #Srxn1 100cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_01_Srxn1", "cmax"] <-      "100cmax"  #Srxn1 100cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_14_Srxn1", "replID"] <-      "replicate_5"  #Srxn1 100cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_14_Srxn1", "cmax"] <-      "100cmax"  #Srxn1 100cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_15_Srxn1", "replID"] <-      "replicate_1"  #Srxn1 50cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_15_Srxn1", "cmax"] <-      "50cmax"  #Srxn1 50cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_21_Srxn1", "replID"] <-      "replicate_2"  #Srxn1 50cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_21_Srxn1", "cmax"] <-      "50cmax"  #Srxn1 50cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_22_Srxn1", "replID"] <-      "replicate_3"  #Srxn1 50cmax     
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_22_Srxn1", "cmax"] <-      "50cmax"  #Srxn1 50cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_28_nii_Srxn1", "replID"] <-  "replicate_1"  #Srxn1 10cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_28_nii_Srxn1", "cmax"] <-  "10cmax"  #Srxn1 10cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_28_niii_Srxn1", "replID"] <- "replicate_2"  #Srxn1 10cmax  
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_08_28_niii_Srxn1", "cmax"] <- "10cmax"  #Srxn1 10cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_04_Srxn1", "replID"] <-      "replicate_3"  #Srxn1 10cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_04_Srxn1", "cmax"] <-      "10cmax"  #Srxn1 10cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_05_Srxn1", "replID"] <-      "replicate_1"  #Srxn1 5cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_05_Srxn1", "cmax"] <-      "5cmax"  #Srxn1 5cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_09_Srxn1", "replID"] <-      "replicate_2"  #Srxn1 5cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_09_Srxn1", "cmax"] <-      "5cmax"  #Srxn1 5cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_16_Srxn1", "replID"] <-      "replicate_1"  #Srxn1 1cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_16_Srxn1", "cmax"] <-      "1cmax"  #Srxn1 1cmax


cytotox.d.sel$plateID <- gsub("25-9-2013_Srxn1", "2013_09_25_Srxn1", cytotox.d.sel$plateID)

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_25_Srxn1", "replID"] <-       "replicate_2"  #Srxn1 1cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_25_Srxn1", "cmax"] <-       "1cmax"  #Srxn1 1cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_26_DDIT3", "replID"] <-      "replicate_1"  #DDIT3 100cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_09_26_DDIT3", "cmax"] <-      "100cmax"  #DDIT3 100cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_03_DDIT3", "replID"] <-      "replicate_2"  #DDIT3 100cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_03_DDIT3", "cmax"] <-      "100cmax"  #DDIT3 100cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_04_DDIT3", "replID"] <-      "replicate_1"  #DDIT3 50cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_04_DDIT3", "cmax"] <-      "50cmax"  #DDIT3 50cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_09_niii_p21", "replID"] <-   "replicate_1"  #p21 100cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_09_niii_p21", "cmax"] <-   "100cmax"  #p21 100cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_09_nii_p21", "replID"] <-    "replicate_2"  #p21 100cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_09_nii_p21", "cmax"] <-    "100cmax"  #p21 100cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_17_nii_DDIT3", "replID"] <-  "replicate_2"  #DDIT3 50cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_17_nii_DDIT3", "cmax"] <-  "50cmax"  #DDIT3 50cma

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_17_niii", "replID"] <-       "replicate_3"  #DDIT3 50 cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_17_niii", "cmax"] <-       "50cmax"  #DDIT3 50 cm

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_24_DDIT3", "replID"] <-      "replicate_1"  #DDIT3 10 cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_10_24_DDIT3", "cmax"] <-      "10cmax"  #DDIT3 10 cm

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_11_20_DDIT3", "replID"] <-      "replicate_2"  #DDIT3 10 cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_11_20_DDIT3", "cmax"] <-      "10cmax"  #DDIT3 10 cm

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_11_27_DDIT3", "replID"] <-      "replicate_1"  #DDIT3 5 cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_11_27_DDIT3", "cmax"] <-      "5cmax"  #DDIT3 5 cma

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_12_04_DDIT3_niii", "replID"] <-      "replicate_1"  #DDIT3 1 cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_12_04_DDIT3_niii", "cmax"] <-      "1cmax"  #DDIT3 

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_12_04_DDIT3_nii", "replID"] <-  "replicate_2"  #DDIT3 1 cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_12_04_DDIT3_nii", "cmax"] <-  "1cmax"  #DDIT3 1 cma

cytotox.d.sel[ cytotox.d.sel$plateID == "2013_12_11_DDIT3", "replID"] <-  "replicate_2"  #DDIT3 5 cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2013_12_11_DDIT3", "cmax"] <-  "5cmax"  #DDIT3 5 cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_05_nii_p21", "replID"] <-    "replicate_1" #p21 10cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_05_nii_p21", "cmax"] <-    "10cmax" #p21 10cma

cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_05_p21_niii", "replID"] <-   "replicate_2" #p21 10cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_05_p21_niii", "cmax"] <-   "10cmax" #p21 10cma

cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_11_p21", "replID"] <-        "replicate_1" #p21 50cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_11_p21", "cmax"] <-        "50cmax" #p21 50cma

cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_18_p21", "replID"] <-        "replicate_1" #p21 5cmax 
cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_18_p21", "cmax"] <-        "5cmax" #p21 5cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_19_p21", "replID"] <-        "replicate_2" #p21 5cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_19_p21", "cmax"] <-        "5cmax" #p21 5cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_25_p21", "replID"] <-        "replicate_1" #p21 1cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2014_02_25_p21", "cmax"] <-        "1cmax" #p21 1cmax

cytotox.d.sel[ cytotox.d.sel$plateID == "2014_03_05_p21", "replID"] <-        "replicate_2" #p21 1cmax
cytotox.d.sel[ cytotox.d.sel$plateID == "2014_03_05_p21", "cmax"] <-        "1cmax" #p21 1cmax



as.data.frame(unique(cytotox.d.sel[treatment%in%"zimelidine", list(treatment, dose_uM, cmax, plateID, cell_line)]))[
  order(unique(cytotox.d.sel[treatment%in%"zimelidine", list(treatment, dose_uM, cmax, plateID, cell_line)])$treatment),]


cytotox.d.selc <- cytotox.d.sel[, list(treatment, dose_uM, cmax, fingerprints, value)]

cytotox.d.sel.summmarized <- cytotox.d.selc[, lapply(.SD, mean, na.rm=TRUE), 
                                                by = c("treatment", "dose_uM", "cmax", "fingerprints")]


unique(cytotox.d.sel.summmarized[, fingerprints]) # final data for cytotox


cytotox.d.sel.summmarized[, log10dose:= log10(dose_uM+1)]

write.table(cytotox.d.sel.summmarized, file = "cytotox_data_originalReporters.txt", sep ="\t", col.names= TRUE)

cytotox.d.sel.summmarized <- read.delim(file = "cytotox_data_originalReporters.txt", sep ="\t", header = TRUE)


# now for the final data-set! the new cmax cytotox data.

cyto.list = list()
path.to.newTox <- "E:/DILI screen/meta analyse DILI screen/nieuwe cmax data/CDsummary"
file.tox <- dir(path.to.newTox)
for( i in seq_along(dir(path.to.newTox)) ) {
  cyto.list[[i]] <- read.delim(paste(path.to.newTox, file.tox[i], sep ="/") ) 
    
}

cytotox.d <- do.call('rbind', cyto.list)
head(cytotox.d)

mean(cytotox.d.sel.summmarized[ cytotox.d.sel.summmarized$fingerprints == "necrosis positive", "value"])# old pi CD fraction 0.02600508
mean( as.data.frame(cytotox.d)[ cytotox.d$variable == "count_obj_final_PI_AreaShape_Area.DIV.obj_nuclei_AreaShape_Area_larger_0.02" , "value"]) #0.1888792
PI.scale.factor <- 0.1888792 -  0.02600508 
# same plate full dose range: scale per plate
cytotox.d <- as.data.table(cytotox.d)
mean( cytotox.d[ variable == "count_obj_final_PI_AreaShape_Area.DIV.obj_nuclei_AreaShape_Area_larger_0.02" , value]) #0.1888792
cytotox.d[ variable == "count_obj_final_PI_AreaShape_Area.DIV.obj_nuclei_AreaShape_Area_larger_0.02" , value:= value - PI.scale.factor, ]
cytotox.d[ variable == "count_obj_final_PI_AreaShape_Area.DIV.obj_nuclei_AreaShape_Area_larger_0.02" & value < 0, value:=0]

dataIn <- cytotox.d
tox.vars <- unique(cytotox.d$variable)
setwd( "E:/DILI screen/meta analyse DILI screen/nieuwe cmax data/summaryFilesNieuwCmax" )
dir()

# normFunction has been made (however check it)
# fix variable names, add cmax. then normalize


#fix dosages: should be exactly the same digits if same concentration
cytotox.d<- as.data.table(cytotox.d)

cytotox.d$comp_dose <- paste(cytotox.d$treatment, cytotox.d$dose_uM, sep="_")

cytotox.d$dose_uM <- round(cytotox.d$dose_uM, digits = 4)

ind.more10 <- cytotox.d$dose_uM >= 10
cytotox.d$dose_uM[ind.more10] <- round(cytotox.d$dose_uM[ind.more10], digits = 1)
ind.more1 <- cytotox.d$dose_uM >= 1
cytotox.d$dose_uM[ind.more1] <- round(cytotox.d$dose_uM[ind.more1], digits = 2)

count(as.data.frame(unique(cytotox.d[, list(treatment, dose_uM) ]))[  order(as.data.frame(unique(cytotox.d[, list(treatment, dose_uM) ]))$treatment)     , ]$treatment)

# create cmax' s

compdose.table<-as.data.frame(unique(cytotox.d[, list(treatment, dose_uM) ]))
compdose.table<-compdose.table[ order(compdose.table$treatment, compdose.table$dose_uM), ]
compdose.table <- compdose.table[ compdose.table$treatment !="DMSO",]

compdose.table$cmax <- rep( c("1cmax", "5cmax", "10cmax", "25cmax", "50cmax", "100cmax"),length(unique(compdose.table$treatment)) )
setkeyv(cytotox.d, c("treatment", "dose_uM"))
compdose.table<- as.data.table(compdose.table)
setkeyv(compdose.table, c("treatment", "dose_uM"))

cytotox.d <-compdose.table[cytotox.d]
unique(cytotox.d[is.na(cmax), treatment])
# excellent =)(=
cytotox.d$comp_dose <- NULL


require(stringr)
# add cell line:
unique(cytotox.d[, plateID])
cytotox.d[, cell_line:= str_extract(plateID, "_([ICAM1CHOPP21SRXN]{3,10})$")]
cytotox.d[, cell_line := gsub("_","",  cell_line)]
unique(cytotox.d[,cell_line])


# now unify and select variable names:




cytotox.d[ ,variable:= gsub("imageCountTracked", "imageCountParentObj", variable) ]
cytotox.d[ , variable:= gsub("count_PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.02", "pi_pos", variable) ]
cytotox.d[ , variable:=gsub("count_obj_final_PI_AreaShape_Area.DIV.obj_nuclei_AreaShape_Area_larger_0.02", "pi_pos", variable) ]
cytotox.d[ , variable:=gsub("count_AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.02", "Anx_pos", variable) ]
cytotox.d[ ,variable:= gsub("count_obj_final_AnV_AreaShape_Area.DIV.obj_nuclei_AreaShape_Area_larger_0.02", "Anx_pos", variable) ]
cytotox.d[ ,variable:= gsub("obj_nuclei_AreaShape_Area", "Nuclei_AreaShape_Area", variable) ]
cytotox.d[ ,variable:= gsub("obj_nuclei_Intensity_MeanIntensity_image_hoechst", "Nuclei_Intensity_MeanIntensity_image_hoechst", variable) ]

cytotox.d <- cytotox.d[  !variable %in% c( "obj_final_PI_AreaShape_Area", "obj_final_AnV_AreaShape_Area" ), ]
cytotox.d <- cytotox.d[  !grepl("min_maxNorm", variable),  ]

# better idea...
cytotox.d <- cytotox.d[  variable %in% c( "imageCountParentObj", "Nuclei_AreaShape_Area", "Nuclei_Intensity_MeanIntensity_image_hoechst",
                                         "pi_pos", "Anx_pos" ), ]
unique(cytotox.d[, variable])
#dataIn <- cytotox.d
#normVar <- "imageCountParentObj"
cytotox.d[, varCell:=NULL]
setnames(cytotox.d, "variable", "fingerprints")
normFunCS_P <- # function for normalizing single time point data with first a fold change wrt plate specific DMSO. also includes outlier removal
  function(dataIn, normVar) {

        buffer <- dataIn[ fingerprints%in% normVar] # select 1 variable
    Checkcontrol <- buffer[ treatment %in% "DMSO"] # select the DMSO control
    controlV<-Checkcontrol[ ,mean(value, na.rm=TRUE), by = plateID]   # only single time point exist for this data
    setkey(controlV, "plateID")
    setkey(buffer, "plateID")
    buffer<- buffer[controlV]
    # remove outliers that are < 1.5 X DMSO values
    buffer <- buffer[ value < ( 1.5 * V1 ) ]
    
    buffer[, FC_value:= value/V1 ]
    buffer[, value:=NULL]
    buffer[, fingerprints:= paste("mmnFC", normVar, sep ="_")]
    
    buffer[, minValue:= min(FC_value, na.rm = TRUE), by = "plateID"]
    buffer[, maxValue:= max(FC_value, na.rm = TRUE), by = "plateID"]
    
    buffer[, mmnvalue:= (FC_value - minValue)/ 
             (maxValue-minValue)]
    
    setnames(buffer, "mmnvalue", "value")
    buffer[, FC_value:=NULL]
    buffer[, V1:=NULL]
    buffer[, maxValue:=NULL]
    buffer[, minValue:=NULL]
    
    output<- rbind(dataIn,buffer)
    return(output)  
  }


unique(cytotox.d[, fingerprints])

cytotox.d <- normFunCS_P(dataIn = cytotox.d, normVar = "imageCountParentObj")
cytotox.d <- normFunCS_P(dataIn = cytotox.d, normVar = "Nuclei_AreaShape_Area")
cytotox.d <- normFunCS_P(dataIn = cytotox.d, normVar = "Nuclei_Intensity_MeanIntensity_image_hoechst")




write.table( cytotox.d, file = "nieuweCmaxCytotox21072016.txt", sep ="\t", col.names = TRUE)
dir()
getwd() #"E:/DILI screen/meta analyse DILI screen/nieuwe cmax data/summaryFilesNieuwCmax"
#nieuw cmax data = klaar. alle data is nu schoon.
# nu in twee stukken opdelen: quality control en reporter + PI
# variabelnamen tussen afzonderlijke datasets gelijk maken
# oude data van nieuwe cmax compounds weg


# paden van data:
#ICAM1-eerst: I:\suz icam1\Summary files CD
#niewe cmax: I:\DILI screen\meta analyse DILI screen\nieuwe cmax data\CDsummary
# origninele reporters: I:/DILI screen/all cytotoxicity output/cytotox Summary Data.txt"

# toxiciteit data orig reporters: E:\DILI screen\meta analyse DILI screen\nieuwe cmax data\summaryFilesNieuwCmax (cytotox_data_originalReporters.txt)

rm(list=ls())
setwd("E:/DILI screen/meta analyse DILI screen/nieuwe cmax data/summaryFilesNieuwCmax")
oldCytotox <- as.data.table(read.delim(file = "cytotox_data_originalReporters.txt", sep ="\t", header =TRUE))
#including ICAM1:


newCmaxCytotox <- as.data.table(read.delim( file = "E:/DILI screen/meta analyse DILI screen/nieuwe cmax data/summaryFilesNieuwCmax/nieuweCmaxCytotox21072016.txt", sep ="\t", header = TRUE) )# modify treatment if TNF


# first process the reporter data.. combine then treat as before
oldReporter <-as.data.table(read.delim( file = "E:/suz icam1/norm_orig_reporters.txt", sep ="\t", header = TRUE))
oldICAM <-  as.data.table(read.delim( file = "E:/suz icam1/ICAM1_orig.all.mean.data.long.txt", sep = "\t", header = TRUE))
setnames(oldICAM, "variable", "fingerprints")
#E:/DILI screen/meta analyse DILI screen/nieuwe cmax data/norm_nieuweCMAX_reporters21072016.txt", col.names= TRUE, sep ="\t")
newComps <- unique(newReporter$treatment) 
newReporter <- as.data.table( read.delim( file = "E:/DILI screen/meta analyse DILI screen/nieuwe cmax data/norm_nieuweCMAX_reporters21072016.txt", sep = "\t", header = TRUE))
newICAM <- as.data.table(read.delim( file = "E:/suz icam1/new.cmaxICAM1.long.txt", sep = "\t", header = TRUE))


##==## loading as of 26/07/2016

setwd("C:\Users\steve_000\Documents\work\DILI paper\DILIpaper")
test <- as.data.table( read.delim( file = "data/norm_origandnew_reporters.txt", sep = "\t", header = TRUE))
nieuwPlates
oldReporter <- test[ !plateID %in% nieuwPlates, ]
newReporter <- test[ plateID %in% nieuwPlates, ]
newComps <- unique(newReporter$treatment) 
newICAM <- as.data.table(read.delim( file = "data/new.cmaxICAM1.long.txt", sep = "\t", header = TRUE))
oldICAM <-  as.data.table(read.delim( file = "data/ICAM1_orig.all.mean.data.long.txt", sep = "\t", header = TRUE))
setnames(oldICAM, "variable", "fingerprints")
oldCytotox <- as.data.table(read.delim(file = "data/cytotox_data_originalReporters.txt", sep ="\t", header =TRUE))
newCmaxCytotox <- as.data.table(read.delim( file = "data/nieuweCmaxCytotox21072016.txt", sep ="\t", header = TRUE) )# modify treatment if TNF

##==##
# combine reporter data. 
# 1) make consistent
# 2) remove the old...
# 3) rbind

dim(oldReporter)
dim(newReporter)
unique(oldReporter[ , cell_line])

newReporter[, variable:=NULL]
unique(oldReporter[, fingerprints])
unique(newReporter[, fingerprints])
#what to keep for main and what to use for QC:
oldReporter.qc <- oldReporter[ !fingerprints %in% c("GFP_pos.2m_Srxn1", "GFP_pos.2m_DDIT3", "GFP_pos.m3sd_p21", "GFP_pos.3m_p21",
                                                    "mmnDMSO.sub_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1",
                                                    "mmnDMSO.sub_Nuclei_Intensity_IntegratedIntensity_Image_GFP_DDIT3",
                                                    "mmnDMSO.sub_Nuclei_Intensity_IntegratedIntensity_Image_GFP_p21",
                                                    "GFP_pos.m3sd_Srxn1", "GFP_pos.3m_Srxn1", "GFP_pos.m3sd_DDIT3",
                                                    "GFP_pos.3m_DDIT3","GFP_pos.2m_p21","mmnDMSO.sub_Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1",
                                                    "mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3",
                                                    "mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_p21"), ]

oldReporter.resp <- oldReporter[ fingerprints %in% c("GFP_pos.2m_Srxn1", "GFP_pos.2m_DDIT3", "GFP_pos.m3sd_p21", "GFP_pos.3m_p21",
                                  "mmnDMSO.sub_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1",
                                  "GFP_pos.m3sd_Srxn1", "GFP_pos.3m_Srxn1", "GFP_pos.m3sd_DDIT3",
                                  "GFP_pos.3m_DDIT3","GFP_pos.2m_p21","mmnDMSO.sub_Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1",
                                  "mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3",
                                  "mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_p21"), ]
unique(oldReporter.resp[,fingerprints])
oldReporter.resp[, plateID:=NULL]



# change names to be same as in old. remove plate ID.
unique(newReporter[, fingerprints])
newReporter[, plateID:= NULL]
# mis integrated intensity voor nuclei gfp metingen, haal dus weg uit oude data
newReporter.qc <- newReporter[ !fingerprints %in% c("GFP_pos.2m_Srxn1", "GFP_pos.2m_DDIT3", "GFP_pos.m3sd_p21", "GFP_pos.3m_p21",
                                                   "mmnDMSO.sub_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1",
                                                   "GFP_pos.m3sd_Srxn1", "GFP_pos.3m_Srxn1", "GFP_pos.m3sd_DDIT3",
                                                   "GFP_pos.3m_DDIT3","GFP_pos.2m_p21","mmnDMSO.sub_Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1",
                                                   "mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3",
                                                   "mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_p21"), ]

newReporter.resp <- newReporter[ fingerprints %in% c("GFP_pos.2m_Srxn1", "GFP_pos.2m_DDIT3", "GFP_pos.m3sd_p21", "GFP_pos.3m_p21",
                                                    "mmnDMSO.sub_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1",
                                                    "GFP_pos.m3sd_Srxn1", "GFP_pos.3m_Srxn1", "GFP_pos.m3sd_DDIT3",
                                                    "GFP_pos.3m_DDIT3","GFP_pos.2m_p21","mmnDMSO.sub_Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1",
                                                    "mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3",
                                                    "mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_p21"), ]










# plot and check for outliers



# remove old cmax

oldReporter.resp
unique(oldReporter.resp[, treatment])

newReporter.resp <- newReporter.resp[ !cmax %in% "25cmax", ]

#remove updated compounds

newReporter.resp[, treatment:=  tolower(treatment)]
newCmaxComps <- unique(newReporter.resp[, treatment])
newCmaxComps <- gsub("dmso", "DMSO" , newCmaxComps)
newCmaxComps <- gsub("dextromethorphan hbr", "dextromethorphan HBr" , newCmaxComps)


indinSet <- newCmaxComps %in% unique(oldReporter.resp[, treatment])
newCmaxComps[ !indinSet ]


indRemoveOld <- oldReporter.resp[, treatment] %in% newCmaxComps & !oldReporter.resp[, treatment] %in% "DMSO"
sum(indRemoveOld)
unique(oldReporter.resp[ indRemoveOld, treatment] )

oldReporter.resp <- oldReporter.resp[ !indRemoveOld, ] 

# add new

dim(oldReporter.resp)

combined.resp <- rbind(oldReporter.resp, newReporter.resp)

combined.resp[is.na(cmax) & treatment %in% "dmso", cmax:= "100cmax"]
combined.resp[, treatment:= gsub("dmso", "DMSO", treatment)]

# remove old cmax from ICAM dataset then add ICAM dataset
# add cmax
# add correct replID
dim(oldICAM)
dim(newICAM)

# use lookup table to add real dose to oldICAM
unique(oldICAM$dose_uM)

oldICAM[ !treatment %in% "DMSO", cmax := paste(round(dose_uM, digit=0), "cmax", sep="")]

sort(unique(oldICAM[ treatment %in% "DMSO", cmax]))
# wetf...

oldICAM[ treatment %in% "DMSO" & dose_uM <=0.002, cmax := "1cmax"]
oldICAM[ treatment %in% "DMSO" & dose_uM > 0.002  & dose_uM <= 0.025 , cmax:= "5cmax"]
oldICAM[ treatment %in% "DMSO" & dose_uM > 0.025  & dose_uM <= 0.05 , cmax:= "10cmax"]
oldICAM[ treatment %in% "DMSO" & dose_uM > 0.05  & dose_uM <= 0.1 , cmax:= "50cmax"]
oldICAM[ treatment %in% "DMSO" & dose_uM > 0.1  , cmax:= "100cmax"]

oldICAM[ treatment %in% "DMSO" & cmax == "1cmax"  , dose_uM:= 0.002 ]
oldICAM[ treatment %in% "DMSO" & cmax == "5cmax"  , dose_uM:= 0.01 ]
oldICAM[ treatment %in% "DMSO" & cmax == "10cmax"  , dose_uM:= 0.02 ]
oldICAM[ treatment %in% "DMSO" & cmax == "50cmax"  , dose_uM:= 0.1 ]
oldICAM[ treatment %in% "DMSO" & cmax == "100cmax"  , dose_uM:= 0.2 ]


# remove old cmax
newcmaxCOmps <- unique(newICAM$treatment)
newcmaxCOmps <- newcmaxCOmps[newcmaxCOmps!="DMSO"]
newcmaxCOmps <- tolower(newcmaxCOmps)
newcmaxCOmps <- gsub("dextromethorphan hbr" , "dextromethorphan HBr"  ,newcmaxCOmps)


newcmaxCOmps %in% oldICAM[, treatment]

indRemove <- oldICAM[, treatment] %in% newcmaxCOmps



oldICAM <- oldICAM[!indRemove, ]

# add dose_uM (DMSO alrdy done)

compdosecmax <- unique(oldReporter.resp[, list(treatment, cmax, dose_uM)])
oldICAM[, dose_uM:=NULL]
setkeyv(compdosecmax, c( "treatment", "cmax"))
setkeyv(oldICAM, c( "treatment", "cmax"))

idd<-unique(oldICAM$treatment) %in% unique(compdosecmax$treatment)
unique(oldICAM$treatment)[!idd]
unique(compdosecmax$treatment) %in% unique(oldICAM$treatment)

oldICAM <- compdosecmax[oldICAM]



oldICAM[ treatment == "DMEM", dose_uM:=0]
oldICAM[ treatment == "aflatoxin B1" & cmax == "100cmax", dose_uM:=100]

# fix cmax in newICAM

lapply(newICAM, class)

  

newICAM[, cmax:=as.character(cmax)]
ind <- newICAM$cmax %in% "0.2"
newICAM[ ind, ]

newICAM[ !ind  , cmax:= paste( cmax, "cmax", sep ="" ) ]
newICAM[ ind, cmax:= "100cmax"]

# add rep ID's:
unique(oldICAM[, list(plateID, cell_line)])

oldICAM[ plateID %in% "20150830 - 100cMax_wells", replID:= "replicate_1" ]
oldICAM[ plateID %in% "20160113 - 100cMax_wells", replID:= "replicate_2" ]
oldICAM[ plateID %in% "20150604 - 10cMax_wells", replID:= "replicate_1" ]
oldICAM[ plateID %in% "20150608 - 10cMax_wells", replID:= "replicate_2" ]
oldICAM[ plateID %in% "20150829 - 1cMax_wells", replID:= "replicate_1" ]
oldICAM[ plateID %in% "20150914 - 1cMax_wells", replID:= "replicate_2" ]
oldICAM[ plateID %in% "20140701 - 50cMax_wells", replID:= "replicate_1" ]
oldICAM[ plateID %in% "20150828 - 50cMax_wells", replID:= "replicate_2" ]
oldICAM[ plateID %in% "20150601 - 5cMax_wells", replID:= "replicate_1" ]
oldICAM[ plateID %in% "20151231 - 5cMax_wells", replID:= "replicate_2" ]
oldICAM[ , plateID := NULL ]

unique(newICAM[, list(plateID, cell_line, cmax)])



newICAM<- newICAM[!cmax %in% "25cmax"]

newICAM[ plateID %in% "20150923_ICAM1_1", replID:= "replicate_1" ]
newICAM[ plateID %in% "20150923_ICAM1_0.2", replID:= "replicate_1" ]
newICAM[ plateID %in% "20150930_ICAM1_1", replID:= "replicate_2" ]
newICAM[ plateID %in% "20150930_ICAM1_0.2", replID:= "replicate_2" ]
newICAM[ plateID %in% "20150923_ICAM1_5", replID:= "replicate_1" ]
newICAM[ plateID %in% "20150930_ICAM1_5", replID:= "replicate_2" ]
newICAM[ plateID %in% "20150923_ICAM1_10", replID:= "replicate_1" ]
newICAM[ plateID %in% "20150930_ICAM1_10", replID:= "replicate_2" ]
newICAM[ plateID %in% "20150923_ICAM1_50", replID:= "replicate_1" ]
newICAM[ plateID %in% "20150930_ICAM1_50", replID:= "replicate_2" ]
newICAM[ plateID %in% "20150923_ICAM1_100", replID:= "replicate_1" ]
newICAM[ plateID %in% "20150930_ICAM1_100", replID:= "replicate_2" ]

newICAM[, plateID:=NULL]

dim(newICAM)
dim(oldICAM)
setnames(newICAM, "variable", "fingerprints")

bothICAM<- rbind( oldICAM, newICAM)
unique(combined.resp$cell_line)
combined.resp <- rbind(combined.resp, bothICAM)

unique(combined.resp$cmax)

sort(as.character(unique(combined.resp$treatment)))
combined.resp[, treatment:=tolower(treatment)]

# remove some compounds:

rmoveCOmps <- c("cccp_c", "cddo_c", "chloroquine_c", "cisplatin_c", "dmog_c","g.f. dmem", "thapsigargin_c", "dmem")
combined.resp <- combined.resp[ !treatment %in% rmoveCOmps, ]

combined.resp[ , cmax := factor(cmax, levels = c("1cmax", "5cmax", "10cmax", "50cmax", "100cmax"))]

dir.create("E:/DILI screen/meta analyse DILI screen/outlierCheck")
writepdfs <- "E:/DILI screen/meta analyse DILI screen/outlierCheck/"
source("C:\\Users\\steve_000\\Documents\\H5CellProfiler\\H5CellProfiler\\theme_sharp.R")

cell_lines <- unique( combined.resp[ , cell_line] )

unique(combined.resp[, fingerprints])



# plot each variable seperately:


for( i in seq_along(cell_lines)){
 
  plotDatatemp <- combined.resp[ cell_line %in% cell_lines[i], ]
  fingerprintVars <- unique( plotDatatemp[, fingerprints ])
 
   for( j in seq_along( fingerprintVars )) {
    plotData <- plotDatatemp[ fingerprints %in% fingerprintVars[ j]]
    p <- ggplot(data =plotData, aes(x= timeAfterExposure, y=value)) + facet_grid(treatment ~ cmax, scales = "free_x")
    p<-p + geom_point(aes(color=replID, shape = replID)) +geom_line(aes(group=replID, color=replID)) + theme_sharp() +theme(strip.text.y = element_text(angle=0))
    pdf(file = paste(writepdfs, cell_lines[i], "_", fingerprintVars[j], ".pdf", sep =""), width =20, height = 120)
      print(p)
    dev.off()
  }
}

# de nieuwe cmax counts: erug weinig te zien. Niks aan te doen want heb ik niet gedaan...

# aan de hand van bovenstaande plotjes bepalen welke fingerprints relevant zijn.
# QC fingerprints even opslaan voordat verwijderen
### plan analyze DILI 2016:
#1) max voor oude heatmap figuur met response en PI
#2) modeling dose 
#3) modeling time dynamics
#4 modeling dose of time grid
#5) in addition to old figures, include surface plots of points 3 and 4

#save.image("freeze22-04-2016.Rdata")

load("freeze22-04-2016.Rdata")
getwd() #"E:/DILI screen/meta analyse DILI screen"
#save.image("freeze24-04-2016.Rdata")

# is iets mis: srxn1 dmsosub zijn zowat allemaal super vlak - watsgebeurt???
# counts lijken wel goed.
# denk iets met de outliers? 

# wat is die ICAM1_norm_ ??
#check voor ICAM1 of DMSO inderdaad het meest naar beneden gaat. of was het tnfa?

##===gefixt====##|<

# meteen statistiek doen/ bewaar de replicates! stop aggregaten en andere afgeleiden van de data in apparte variabelen of sla op als txt.

combined.resp <- combined.resp[ !(treatment == "carmustine" & dose_uM %in% c(0.356, 1.78, 3.56, 17.8, 35.6 ) )]
combined.resp[ treatment == "cimetidine" & dose_uM ==11.89, dose_uM:=11.9]
combined.resp[ treatment == "cimetidine" & dose_uM ==59.45, dose_uM:=59.5]
combined.resp[ treatment == "isoniazid" & dose_uM ==76.56, dose_uM:=76.6]
combined.resp[ treatment == "nefazodone" & dose_uM ==19.750, dose_uM:=19.8]
combined.resp[ treatment %in% "zimelidine" & dose_uM==1.335 & cmax == "5cmax", dose_uM:= 1.34 ]
combined.resp[ treatment %in% "zimelidine" & dose_uM==13.350 & cmax == "50cmax", dose_uM:= 13.300 ]
combined.resp <- combined.resp[ !treatment %in% c("doxycycline", "tetracycline", "dmem", "doxorubicin"), ]

combined.resp[, treatment:=gsub("carbamazapine", "carbamazepine", treatment) ]
combined.resp[, treatment:=gsub("cyclohexamine", "cycloheximide", treatment) ]
combined.resp[, treatment:=gsub("etoposide ", "etoposide", treatment) ]


# SAVED REPLICATE DATA FOR STATISTICS OF TIME COURSES (21-07-2016)
getwd() #"E:/DILI screen/meta analyse DILI screen"
write.table(file = "data/combined.resp_replicates.txt", combined.resp, col.names= TRUE, sep = "\t")
unique(combined.resp[, cell_lines])
# statistics: to keep it simple: only need significance of all time courses with respect to DMSO.
combined.resp <- read.table( file = "combined.resp_replicates.txt", header = TRUE, sep ="\t")
#require(data.table)
combined.resp <- as.data.table(combined.resp)




# so after alot of considerations... the statistics method decision:
# 1) model response time courses
# 2) sample time for equal time intervals for all data
# 3) perform fda tests
# 4) would be cool to shade the significance in the 3d plots!
# 5) although.. fit the surface, not first time then dose... didnt get any satisfactory fits:
# model per dose de time courses
# maak een grid: per compound-replicate een matrix
# maak een matrix met t-test waarden tussen de replicates
# bereken gemiddelde van de replicate gemiddelde responsen
# plot met 'persp' en kleur met de p-waarden

# maar eerst de max in time heatmap figuur:
# max of time followed by mean of replicates:
# first select relevant fingerprints. and combine the bellow and aboveICAM counts, and add PI
unique(combined.resp[, fingerprints])


mean1above <- combined.resp[ fingerprints %in% "mean1above_obj_cyto_only_Intensity_IntegratedIntensity_image_GFP"]
mean1bellow <- combined.resp[ fingerprints %in% "mean1below_obj_cyto_only_Intensity_IntegratedIntensity_image_GFP"]

GFP_dif.m1.5_ICAM1 <- cbind( mean1above, mean1bellow[, value])
GFP_dif.m1.5_ICAM1[, GFP_dif.m1.5_ICAM1 := value - V2 ]
summary(GFP_dif.m1.5_ICAM1$GFP_dif.m1.5_ICAM1)
GFP_dif.m1.5_ICAM1$fingerprints <- "GFP_dif.m1.5_ICAM1"
GFP_dif.m1.5_ICAM1[, V2:=NULL]
GFP_dif.m1.5_ICAM1[, value:=NULL]
setnames(GFP_dif.m1.5_ICAM1, "GFP_dif.m1.5_ICAM1", "value")


mean2above <- combined.resp[ fingerprints %in% "mean2above_obj_cyto_only_Intensity_IntegratedIntensity_image_GFP"]
mean2bellow <- combined.resp[ fingerprints %in% "mean2below_obj_cyto_only_Intensity_IntegratedIntensity_image_GFP"]
GFP_dif.m2_ICAM1 <- cbind( mean2above, mean2bellow[, value])
GFP_dif.m2_ICAM1[, GFP_dif.m2_ICAM1 := value - V2 ]
summary(GFP_dif.m2_ICAM1$GFP_dif.m2_ICAM1)
GFP_dif.m2_ICAM1$fingerprints <- "GFP_dif.m2_ICAM1"
GFP_dif.m2_ICAM1[, V2:=NULL]
GFP_dif.m2_ICAM1[, value:=NULL]
setnames(GFP_dif.m2_ICAM1, "GFP_dif.m2_ICAM1", "value")



mean3above <- combined.resp[ fingerprints %in% "mean3above_obj_cyto_only_Intensity_IntegratedIntensity_image_GFP"]
mean3bellow <- combined.resp[ fingerprints %in% "mean3below_obj_cyto_only_Intensity_IntegratedIntensity_image_GFP"]
GFP_dif.m3_ICAM1 <- cbind( mean3above, mean3bellow[, value])
GFP_dif.m3_ICAM1[, GFP_dif.m3_ICAM1 := value - V2 ]
summary(GFP_dif.m3_ICAM1$GFP_dif.m3_ICAM1)
GFP_dif.m3_ICAM1$fingerprints <- "GFP_dif.m3_ICAM1"
GFP_dif.m3_ICAM1[, V2:=NULL]
GFP_dif.m3_ICAM1[, value:=NULL]
setnames(GFP_dif.m3_ICAM1, "GFP_dif.m3_ICAM1", "value")


dim(GFP_dif.m1.5_ICAM1)
dim(GFP_dif.m2_ICAM1)


lapply(GFP_dif.m1.5_ICAM1, class)
all.three <- rbind.data.frame(GFP_dif.m1.5_ICAM1,GFP_dif.m2_ICAM1,GFP_dif.m3_ICAM1)

combined.resp <- rbind(combined.resp, all.three)




# add min max norm norm_ data (between -1 and 1)



keepFP <- c("GFP_pos.m3sd_Srxn1", "GFP_pos.2m_Srxn1", "GFP_pos.3m_Srxn1","GFP_pos.m3sd_DDIT3","GFP_pos.2m_DDIT3",
            "GFP_pos.3m_DDIT3", "GFP_pos.m3sd_p21", "GFP_pos.2m_p21", "GFP_pos.3m_p21", "mmnDMSO.sub_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1",
            "mmnDMSO.sub_Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1","mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3",
            "mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_p21", "GFP_dif.m1.5_ICAM1", "GFP_dif.m2_ICAM1", "GFP_dif.m3_ICAM1")

GFP_data <- combined.resp[ fingerprints %in% keepFP  , list(treatment, cmax, timeAfterExposure, fingerprints, replID, value, dose_uM)]
GFP_data


dir.create("E:/DILI screen/meta analyse DILI screen/GFP_data")
writepdfs <- "E:/DILI screen/meta analyse DILI screen/GFP_data/"
source("C:\\Users\\steve_000\\Documents\\H5CellProfiler\\H5CellProfiler\\theme_sharp.R")


unique(GFP_data[, keepFP])








#getwd()"C:/Users/steve_000/Documents/work/DILI paper/DILIpaper"

save.image("freeze27072016.RData") 



# first get index of negatives: 
dmso.data <- GFP_data[treatment =="dmso"]
GFP_data <-  GFP_data[treatment !="dmso"]

# there is alot of duplicate wells of DMSO adding noise tot analysis:
dmso.data <- dmso.data[, mean(value), by = c("treatment", "cmax", "timeAfterExposure", "fingerprints", "replID", "dose_uM")]
setnames(dmso.data, "V1", "value")
GFP_data <- rbind(GFP_data, dmso.data)



# plot each variable seperately:
GFP_data[, cmax :=factor(cmax, levels = c("1cmax", "5cmax", "10cmax", "50cmax", "100cmax"))]

for( j in seq_along( keepFP )) {
  plotData <- GFP_data[ fingerprints %in% keepFP[ j]]
  p <- ggplot(data =plotData, aes(x= timeAfterExposure, y=value)) + facet_grid(treatment ~ cmax, scales = "free_x")
  p<-p + geom_point(aes(color=replID, shape = replID)) +geom_line(aes(group=replID, color=replID)) + theme_sharp() +theme(strip.text.y = element_text(angle=0))
  pdf(file = paste(writepdfs, "_", keepFP[j], ".pdf", sep =""), width =20, height = 120)
  print(p)
  dev.off()
}



# selecting max requires indexing the negative values from ICAM1 (where the max is defined as the min)

max.vals <- GFP_data[, lapply( .SD, function(x) { max(as.numeric(x), na.rm = TRUE)} # first max over time
                                ), by = c("treatment",
                                          "cmax",
                                          "fingerprints",
                                          "replID"),
                       .SDcols = "value"]

min.vals <- GFP_data[, lapply( .SD, function(x) { min(as.numeric(x), na.rm = TRUE)} # first max over time
), by = c("treatment",
          "cmax",
          "fingerprints",
          "replID"),
.SDcols = "value"]


ind.neg <- abs(min.vals[, value]) > max.vals[ , value]



GFP_data.s <- GFP_data[, lapply( .SD, function(x) { max(as.numeric(abs(x)), na.rm = TRUE)} # first max over time
), by = c("treatment",
          "cmax",
          "fingerprints",
          "replID"),
.SDcols = "value"]



GFP_data.s[ind.neg, value := -value] 
GFP_data.s[treatment=="dmso"]

max.vals[treatment=="dmso"]
min.vals[treatment=="dmso"]

GFP_data.sc <- GFP_data.s[, list(treatment, cmax, fingerprints, value)]
GFP_data.s <- GFP_data.sc[, lapply(.SD, # this is all the response data final data
                                   mean, na.rm = TRUE  # then mean over reps
                                   ), by = c("treatment",
                                             "cmax",
                                             "fingerprints")]




newCmaxCytotox <- read.delim(file = "nieuwe cmax data/summaryFilesNieuwCmax/nieuweCmaxCytotox21072016.txt",  sep ="\t")
newCmaxCytotox <- as.data.table(newCmaxCytotox)
oldCmaxCytotox <- read.delim(file = "nieuwe cmax data/summaryFilesNieuwCmax/cytotox_data_originalReporters.txt",  sep ="\t")
oldCmaxCytotox <- as.data.table(oldCmaxCytotox)
newCmaxCytotox$treatment <- tolower(newCmaxCytotox$treatment)
oldCmaxCytotox$treatment <- tolower(oldCmaxCytotox$treatment)

newCOmps <- unique(newCmaxCytotox$treatment)
newCOmps<-newCOmps[newCOmps!="dmso"]
newCOmps %in% unique(oldCmaxCytotox$treatment)
ind.newComps <- oldCmaxCytotox[, treatment] %in% newCOmps
oldCmaxCytotox <- oldCmaxCytotox[ !ind.newComps, ] # remove comps of new cmax

unique(oldCmaxCytotox[, fingerprints ])


oldAndNewCmaxCytotox <- rbind( oldCmaxCytotox[, list(treatment, dose_uM, cmax, fingerprints, value)], 
                               newCmaxCytotox[ ,list(treatment, dose_uM, cmax, fingerprints, value) ] )
oldAndNewCmaxCytotox[, fingerprints := gsub( "pi_pos", "necrosis_pos", fingerprints) ]
oldAndNewCmaxCytotox[, fingerprints := gsub( "necrosis positive", "necrosis_pos", fingerprints) ]


#write.table(oldAndNewCmaxCytotox, file = "nieuwe cmax data/summaryFilesNieuwCmax/oldAndNewCmaxCytotox21072016.txt", col.names= TRUE, sep ="\t")
oldAndNewCmaxCytotox <- as.data.table( read.delim( file = "data/oldAndNewCmaxCytotox21072016.txt", sep ="\t", header = TRUE) )
pi.data <- oldAndNewCmaxCytotox[fingerprints %in% "necrosis_pos",]



# add the PI data:

GFP_data.s <- rbind(GFP_data.s, pi.data[, list(treatment, cmax, fingerprints, value) ])
GFP_data.s[, fingerprints := factor(fingerprints)]
# bar plotje proberen:



GFP_data.s[ cmax== "0cmax",]
unique(GFP_data.s[, cmax])


GFP_data.s <- GFP_data.s[ cmax != "25cmax",]
GFP_data.s[, cmax:= factor(cmax, levels = c("1cmax","5cmax","10cmax", "50cmax", "100cmax"))]
# data molding complete, now for clustering and plotting
# reformat data shape


GFP_data.s <- GFP_data.s[!treatment %in% tolower(c("CCCP_C", "CDDO_C", "chloroquine_C", "cisplatin_C", "DMOG_C","g.f. DMEM", "thapsigargin_C"))]
#remove control compounds

which(as.data.frame(GFP_data.s[, length(value), by = treatment])$V1!=145)


contrl.data<- GFP_data.s[treatment %in% c("dmso")]
GFP_data.s<- GFP_data.s[!treatment %in% c("dmso")]

contrl.data<- contrl.data[, mean(value, na.rm=TRUE), by = list(treatment, cmax, fingerprints)]

setnames(contrl.data, "V1", "value")
GFP_data.s <- rbind(GFP_data.s, contrl.data)

GFP_data.s<-GFP_data.s[, lapply(.SD, mean(value), na.rm=TRUE), by = list(treatment, cmax, fingerprints)]

GFP_data.s[, treatment:=factor(treatment)]


GFP_data.s[, fingerprints := gsub("mmnDMSO.sub_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1", 
                                  "mmnDMSO.sub.IntCyto_Srxn1", fingerprints)]
GFP_data.s[, fingerprints := gsub("mmnDMSO.sub_Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1", 
                                  "mmnDMSO.sub.MeanCyto_Srxn1", fingerprints)]
GFP_data.s[, fingerprints := gsub("mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3", 
                                  "mmnDMSO.sub.MeanNuclei_DDIT3", fingerprints)]
GFP_data.s[, fingerprints := gsub("mmnDMSO.sub_Nuclei_Intensity_MeanIntensity_Image_GFP_p21", 
                                  "mmnDMSO.sub.MeanNuclei_p21", fingerprints)]

GFP_data.s <- GFP_data.s[ treatment != "doxycycline"] # autofl
GFP_data.s <- GFP_data.s[ treatment != "doxorubicin"] # autofl
GFP_data.s <- GFP_data.s[ treatment != "dmem"]
GFP_data.s[, treatment:=factor(treatment)]

all.data.w <- dcast(data= GFP_data.s,  treatment + fingerprints ~cmax, value.var = "value") # later I changed all.data.w: so dont use this GFP_data.s
all.data.w <- data.table(all.data.w)
head(all.data.w)
require(NMF)
all.data.w$fingerprints <- factor(all.data.w$fingerprints)
unique(all.data.w$fingerprints)

# cluster the compounds based on similarity of fingerprint dose response vectors
# then generate a heatmap 17 features 5 doses is 85 columns by 177 compounds



all.data.w[, fingerprints := factor(fingerprints)]




buffer <- all.data.w[ , (unique(fingerprints)), by ="treatment"]
as.data.frame(all.data.w[ , length(unique(fingerprints)), by ="treatment"])

setdiff(buffer[ treatment == "acarbose"], buffer[ treatment == "acetaminophen"])
buffer[ treatment == "acarbose"]
buffer[ treatment == "acetaminophen"]

unique(all.data.w$treatment)

colnames(all.data.w)[3:7] <- c("cmax1", "cmax5", "cmax10", "cmax50", "cmax100")



all.data.w[, treatment:=factor(treatment)]


# remove compounds: only fda labels compounds and gold compounds should remain:
treatmentAnot <- read.table(file = "E:/DILI screen/meta analyse DILI screen/metaanalaysis/DILI anotation.txt", sep ="\t", header = TRUE)
treatmentAnot <- na.omit(treatmentAnot)

treatmentAnot <- treatmentAnot[, c( "treatment", "updatedDILIpaper")]

# remove compounds that are not in treatmentAnot file:
# first check if all treatment in anotation file are found in data
length(treatmentAnot$treatment) == length(unique(treatmentAnot$treatment))
treatmentAnot$treatment <- tolower(treatmentAnot$treatment)
indNot <- which(!treatmentAnot$treatment %in% all.data.w[, treatment])
treatmentAnot[indNot, ]
indYes <- which(treatmentAnot$treatment %in% all.data.w[, treatment])
treatmentAnot <- treatmentAnot[indYes, ]

indKeep <- all.data.w[, treatment] %in% treatmentAnot$treatment 

all.data.w <- all.data.w[indKeep,]

all(all.data.w[, treatment] %in% treatmentAnot$treatment)
length(unique(all.data.w$treatment)) == length(treatmentAnot$treatment)


save.image("27-07-2016-woe.Rdata")
load("27-07-2016-woe.Rdata")

# voor het clusteren de compounds verwijderen die niet in annotatie bestand zijn. .

# which compounds to remove from dataset?

# all.data[, treatment:=gsub("carbamazapine", "carbamazepine", treatment) ]
# all.data[, treatment:=gsub("cyclohexamine", "cycloheximide", treatment) ]
# all.data[, treatment:=gsub("etoposide ", "etoposide", treatment) ]

all.data.w[ , treatment := factor(treatment)]
longVecsPerComp <- split(all.data.w, all.data.w$treatment)
length(unique(all.data.w$treatment))

length(treatmentAnot$treatment)

head(longVecsPerComp)
gathervecs = list()
gatherCompVecs = list()


for( j in seq_along(longVecsPerComp)){
  for(i in 1:nrow(longVecsPerComp[[j]])){

    gathervecs[[i]] <- longVecsPerComp[[j]][i,  list(cmax1, cmax5, cmax10, cmax50, cmax100)]

  }
  gatherCompVecs[[j]] <- unlist(gathervecs)
names(gatherCompVecs)[j] <- as.character(unique(longVecsPerComp[[j]]$treatment))
print(j)
}





gatherCompVecsDF<-do.call("rbind", gatherCompVecs)
head(gatherCompVecsDF)

fingerprintAnot <- longVecsPerComp[[1]]$fingerprints
fingerprintAnot<- rep(fingerprintAnot, each = 5)
# now distance measures:


dist.gatherCompVecsDF <- dist(gatherCompVecsDF, method = "manhattan")

gatherCompVecsDF.clust<-hclust(dist.gatherCompVecsDF, method = "complete")

plot(gatherCompVecsDF.clust)

names(gatherCompVecsDF.clust)
gatherCompVecsDF.clust$labels
gatherCompVecsDF.clust$order

gatherCompVecsDF.clust.order <- gatherCompVecsDF.clust$labels[gatherCompVecsDF.clust$order]


head(gatherCompVecsDF)

all.data.w
all.data <- melt(all.data.w, id.vars = c("treatment","fingerprints"), variable.name = "cmax")

#all.data <- GFP_data.s
rm("GFP_data.s")
GFP_data.s <- all.data


all.data[, treatment:=gsub("carbamazapine", "carbamazepine", treatment) ]
all.data[, treatment:=gsub("cyclohexamine", "cycloheximide", treatment) ]
all.data[, treatment:=gsub("etoposide ", "etoposide", treatment) ]


rownames(gatherCompVecsDF)<- gsub("carbamazapine", "carbamazepine", rownames(gatherCompVecsDF))
rownames(gatherCompVecsDF)<- gsub("cyclohexamine", "cycloheximide", rownames(gatherCompVecsDF))
rownames(gatherCompVecsDF)<- gsub("etoposide ", "etoposide", rownames(gatherCompVecsDF))

any(rownames(gatherCompVecsDF)== "etoposide ")
#==#


#write.table(all.comps, file = "all.comps.txt", sep = "\tab", row.names = F)

dir()
getwd()
#treatmentAnot <- read.table(file = "E:/DILI screen/meta analyse DILI screen/metaanalaysis/DILI anotation.txt", sep ="\t", header = TRUE)

all_anot <- read.delim( file = "../allDILIanotated.csv")
head(all_anot)

treatmentAnot$oldDILI <- NULL
treatmentAnot$updatedDILI <- NULL
head(treatmentAnot)
treatmentAnot <-  na.omit(treatmentAnot)



# haal eruit (deze zitten niet in dili anotatie van fda en zijn ook geen gold compounds: (verwijder ze in annotatie bestand)
# azidothymidine
# bisphenol A
# carmustine
# cromolyn
# cycloheximide
# ethionine
# ochratoxin
# oligomycin B
# phenacetin
# sulfamethoxazole
# tacrine
# valacyclovir


colnames(treatmentAnot) <- c("treatment", "DILI_Type")



unique(treatmentAnot$DILI_Type)
treatmentAnot$DILI_Type<-gsub("No-DILI-Concern\\?", "No-DILI-Concern", treatmentAnot$DILI_Type)
treatmentAnot$DILI_Type<-gsub("No-DILI-concern", "No-DILI-Concern",treatmentAnot$DILI_Type)
treatmentAnot$DILI_Type<-gsub("Most-DILI-Concern\\?", "Most-DILI-Concern",treatmentAnot$DILI_Type)
rownames(treatmentAnot) <- treatmentAnot$treatment


rownames(treatmentAnot)<- tolower(rownames(treatmentAnot))

indAnot1 <- match(rownames(treatmentAnot), rownames(gatherCompVecsDF)) # in annotation file compounds to be kept remain


treatmentAnot<-treatmentAnot[!is.na(indAnot1), ]


indAnot1 <- match(rownames(treatmentAnot), rownames(gatherCompVecsDF))

gatherCompVecsDF <- gatherCompVecsDF[indAnot1, ]

treatmentAnot$treatment<-NULL

any(is.na(indAnot1)) # FALSE

indAnot1 <- match(rownames(treatmentAnot), rownames(gatherCompVecsDF))


#treatmentAnot1 <-treatmentAnot$DILI_Type[indAnot1]

gatherCompVecsDF.clust.order <- gsub("carbamazapine", "carbamazepine", gatherCompVecsDF.clust.order)
gatherCompVecsDF.clust.order <- gsub("cyclohexamine", "cycloheximide", gatherCompVecsDF.clust.order)
gatherCompVecsDF.clust.order <- gsub("etoposide ", "etoposide", gatherCompVecsDF.clust.order)


indOrder <- match( gatherCompVecsDF.clust.order,rownames(gatherCompVecsDF))

indOrderAnot <- match( gatherCompVecsDF.clust.order, rownames(treatmentAnot))

treatmentAnot1 <- treatmentAnot[indOrderAnot,, drop= F]



require(RColorBrewer)
DILI_Type <- c(brewer.pal(9,"YlOrRd")[c(1,3,6)],"pink", "light blue")
hist(1:100 , col=DILI_Type)
names(DILI_Type) <- c("No-DILI-Concern","Less-DILI-Concern", "Most-DILI-Concern", "Mech-control", "Neg-control")

myColors <- list(DILI_Type = DILI_Type ) 

myColors$cmax <- brewer.pal(5, "BuGn" )
myColors$fingerprint = "black"
names(myColors$cmax) <- c("cmax1", "cmax5","cmax10", "cmax50", "cmax100")



colnames(gatherCompVecsDF) <- paste(fingerprintAnot, colnames(gatherCompVecsDF))

treatmentAnot1<-data.frame(DILI_Type=treatmentAnot1)
treatmentAnotCol <- list(data.frame(fingerprint = fingerprintAnot))

#write.table(colnames(gatherCompVecsDF), file = "column names.txt", sep ="\t")

reordCols <- read.table(file= "E:/DILI screen/meta analyse DILI screen/metaanalaysis/reordered column names.txt", sep="\t")

reordAnotCols <- c(
  rep("Srxn1", 5), rep(c("Chop", "p21"), each = 4), rep("ICAM1", 3), "CellDeath"
)

reordAnotCols<- rep(reordAnotCols, each = 5)

indCols <- match(reordCols$V1, colnames(gatherCompVecsDF))


annCol <- list(cmax=
  rep(
    c("cmax1","cmax5", "cmax10","cmax50", "cmax100"),  17
    ),
  fingerprint = reordAnotCols
  )

lapply(annCol, length)
lapply(myColors,length)  


head(reordAnotCols)


myColors$fingerprint <- c(rep(c("light blue", "light green", "#CC6633", "purple", "black"), each = 1))


annCol$fingerprint <- factor(annCol$fingerprint, levels=
                               c("Srxn1", "Chop", "p21", "ICAM1", "CellDeath"), ordered = FALSE)


#colBreaks <-c(seq(-0.01, 0.3,length.out=15), seq(0.301,1.01, length.out = 14))


colBreaks <-c(  -seq(1.01, 0.411, length.out = 8), -seq( 0.4, 0.01, length.out=16)  , 
                c(seq(0.01, 0.4,length.out=16), seq(0.411,1.01, length.out = 8)) )

29*1.5
43.5/4
length(colBreaks)
length(my.colors)
plot(colBreaks)

my.colors <- c(rev(brewer.pal(9, "YlGnBu")) , brewer.pal(9, "YlOrRd"))


pie(rep(1, length(my.colors)), labels = sprintf("%d (%s)", seq_along(my.colors), 
                                                my.colors), col = my.colors)

head(gatherCompVecsDF)
# doel 1: deze heatmap klaarmaken voor publicatie
# stap 1: fix normalizatie probleem met nieuwe cmax:
#a) phenobarbital dipt enorm, waardoor alles te hoog is. normalizatie zonder phenobarbital?
# stap 2: verwijder de mechanistic control compounds en check bobs mail. hou wel gold compounds er in.

getwd()
nrow(gatherCompVecsDF)
(treatmentAnot1)
pdf("figures/sel_anot_max_of_time_heatmap_manh.complete.pdf", width= 12, height = 18)
aheatmap(gatherCompVecsDF[indOrder, indCols], Rowv = NA, Colv=NA, 
         color = my.colors, breaks = colBreaks,
         fontsize = 10, width = 10, height=10, legend = TRUE,
         annRow = treatmentAnot1 , annCol = annCol   ,annColors = myColors  #,txt =gatherCompVecsDF[indOrder, indCols]
)
dev.off()
write.table(treatmentAnot, file = "finalDILIpaper-treatmentAnot.txt", sep ="\t", col.names = TRUE)

#save.image("data/freeze30-07-2016.RData")
load("data/freeze30-07-2016.RData")

#load("data/freeze27-07-2016.RData")


# next step: remove/ modify annotations according to feedback. Of mechanistic controls only keep gold compounds.
# manh + complete lijkt best tot nu toe. Deze of allen na annotatie verwerking.

#save.image("freeze25-04-2016.tmp2.RData")



# vanavond vervolg:
# modelen time courses. etc zie boven mbt statistiek enzo


# plot dose response curves

#make annotation table for adding concentrations to cmax conveniently (niet mnet cytotox verzin wat anders....)

# GFP_data.s has been processed. Use combined.resp to add dose_uM



indKeep <- combined.resp$treatment %in% all.data.w$treatment
combined.resp<-combined.resp[indKeep,]



unique(combined.resp[ treatment == "dextromethorphan hbr", list(cmax,dose_uM)])
# old dose for this compound sneaked in...?
1: 100cmax  0.7700
2:  50cmax  0.3850
3:  10cmax  0.0770
4:   5cmax  0.0385
5:   1cmax  0.0077

1:   1cmax  0.0220
7:   5cmax  0.1100
8:  10cmax  0.2200
9:  50cmax  1.1000
10: 100cmax 2.2000


combined.resp[ treatment == "dextromethorphan hbr" & dose_uM == 0.0077, dose_uM := 0.0220]
combined.resp[ treatment == "dextromethorphan hbr" & dose_uM == 0.0385, dose_uM := 0.1100] 
combined.resp[ treatment == "dextromethorphan hbr" & dose_uM == 0.0770, dose_uM := 0.2200] 
combined.resp[ treatment == "dextromethorphan hbr" & dose_uM == 0.3850, dose_uM := 1.1000] 
combined.resp[ treatment == "dextromethorphan hbr" & dose_uM == 0.7700, dose_uM := 2.2000] 



cmaxDoseTreatment <- combined.resp[, list(treatment, dose_uM, cmax)]
cmaxDoseTreatment <- unique(cmaxDoseTreatment)

length(unique(na.omit(cmaxDoseTreatment)$treatment))
cmaxDoseTreatment <- na.omit(cmaxDoseTreatment)

#fix dmso and DMEM concentrations
cmaxDoseTreatment[treatment%in%"dmso"]



cmaxDoseTreatment[ treatment%in%"dmso" & cmax == "100cmax" , dose_uM:= 0.2]
cmaxDoseTreatment <- unique(cmaxDoseTreatment)
cmaxDoseTreatment[treatment%in%"dmso"]



test <-as.data.frame(unique(cmaxDoseTreatment[, list(treatment,dose_uM, cmax)])[order(unique(cmaxDoseTreatment[, list(treatment,dose_uM, cmax)])$treatment)])
test
test<-as.data.table(test)
test[, .N,by = list(treatment)][N!=5]


newCmaxCheck <- tolower(c(
"Dextromethorphan Hbr",
"Ximelagatran",
"Phenobarbital",
"Zimelidine",
"Rotenone",
"Diltiazem",
"Carmustine",
"Acetaminophen",
"Aspirin",
"Diclofenac",
"Metformin" ,
"Bosentan",
"Nefazodone",
"Cromolyn",
"Isoproterenol",
"Verapamil",
"Clofibrate",
"Clozapine",
"Isoniazid",
"Cimetidine",
"Cyclophosphamide",
"Fialuridine"
))


test[ treatment %in% newCmaxCheck & cmax == "1cmax"]

test[treatment =="dextromethorphan hbr"]

# add concentrations to GFP_data.s
setkeyv(GFP_data.s, c("treatment","cmax"))

cmaxDoseTreatment[ , cmax := gsub( "1cmax", "cmax1", x = cmax)]
cmaxDoseTreatment[ , cmax := gsub( "5cmax", "cmax5", x = cmax)]
cmaxDoseTreatment[ , cmax := gsub( "10cmax", "cmax10", x = cmax)]
cmaxDoseTreatment[ , cmax := gsub( "50cmax", "cmax50", x = cmax)]
cmaxDoseTreatment[ , cmax := gsub( "100cmax", "cmax100", x = cmax)]

setkeyv(cmaxDoseTreatment, c("treatment","cmax"))

all.data.withDose <- cmaxDoseTreatment[GFP_data.s ]
inddd <- which(cmaxDoseTreatment[, lapply(.SD, length), by=treatment]$cmax >5)
as.data.frame(cmaxDoseTreatment[, lapply(.SD, length), by=treatment])


all.data.withDose[, dose_uM:=as.numeric(dose_uM)]
all.data.withDose[, log10Dose:= log10(dose_uM)]
all.data.withDose[is.infinite(log10Dose), log10Dose:= 0 ]

all.data.withDose[, cmax:= factor(cmax, levels = c("cmax1", "cmax5" , "cmax10","cmax50","cmax100"), ordered = TRUE )]
# 
# cur.fingers<- unique(all.data.withDose$fingerprints)
# for(i in seq_along(cur.fingers)){
#   
#   pdf(paste(cur.fingers[i],"dose response DILI.pdf"), height=30, width=30)
#   p<- ggplot(data=all.data.withDose[fingerprints%in% cur.fingers[i]] , aes(x = log10Dose, y = value)) +
#     geom_point(aes(color= cmax),size = 2) + facet_wrap( ~treatment, scales = "free_x" ) +
#     theme_sharp() +
#     theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.4, size = 12 , 
#                                       colour = "grey50") ) + theme( strip.text.x = element_text( )) +
#     ggtitle(   paste(cur.fingers[i],"Dose Response DILI")) + 
#     theme(plot.title = element_text(lineheight=.8, size = 14 )) + theme(legend.text=element_text(size=12))
#   p <-p+geom_line(aes(group=treatment))  + theme(legend.position = "bottom")
#   p <- p + coord_cartesian( ylim =c(0,1.1))
#   print(p)
#   dev.off()
# }
#

# plot dose response voor boekje: in de trend van time course format. subplotjes voor 3 cellijnen, 3 kolomen om alle compounds in te verdelen (4 voor paper)
all.data.withDose
# select relevent fingerprints

dosePlot <- all.data.withDose[ fingerprints%in% c("GFP_pos.2m_Srxn1", "GFP_pos.2m_DDIT3", "GFP_pos.2m_p21", "GFP_dif.m2_ICAM1")]
dosePlot[, fingerprints := factor(fingerprints, levels = c("GFP_pos.2m_Srxn1", "GFP_pos.2m_DDIT3", "GFP_pos.2m_p21", "GFP_dif.m2_ICAM1"))]
all.treats <- unique(dosePlot$treatment)
length(all.treats)/ 4
3*35 + 1* 36 # 141

105+36
all.treats[1:35]
all.treats[36:70]
all.treats[71:105]
all.treats[106:141]

dir()
#save.image("freeze2504_2016.RData")

getwd()
#setwd("E:/DILI screen/meta analyse DILI screen/")

dir()
getwd()

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
rcols <- gg_color_hue(4)

pdf("figures/dose plots/dose response DILI 1 of 4_line.pdf", height=55, width=30)
p<- ggplot(data=dosePlot[treatment%in% all.treats[1:35]] , aes(x = cmax, y = value)) +
  geom_point(aes(color= fingerprints),size = 2) + facet_grid( treatment~fingerprints,  scales = "free_x" ) +
  
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.4, size = 12 , 
                                    colour = "grey50") ) + theme( strip.text.x = element_text( )) +
  
  theme(plot.title = element_text(lineheight=.8, size = 14 )) + theme(legend.text=element_text(size=12))

#p+geom_smooth(aes(group=treatment, color = fingerprints, size =6), method = "loess")
p = p + geom_line(aes(group=treatment, color = fingerprints, size =6)) #+ theme(legend.position = "bottom"))
p <- p + coord_cartesian( ylim =c(-1.1,1.1))  + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + geom_hline(yintercept=0, alpha = 0.5, linetype = 2) + 
  scale_colour_manual(values =rcols[c(3,2,1,4)])
print(p)
dev.off()


warnings()
save.image("31072016.Rdata")


# plot PoD graph

all.data.withDose
# use the response data
unique(all.data.withDose[, fingerprints])
PoD <- all.data.withDose[ fingerprints%in% c("GFP_pos.2m_Srxn1", "GFP_pos.2m_DDIT3", "GFP_pos.2m_p21")]

# 50% of cells above 2X DMSO level

PoDtreatAnot<- data.table(treatmentAnot)
setkey(PoDtreatAnot,"treatment")
setkey(PoD, "treatment") 
PoD<-PoD[PoDtreatAnot]

PoD$DILI_Type <- factor( PoD$DILI_Type, levels = c("No-DILI-Concern", "Less-DILI-Concern", "Most-DILI-Concern", "control"), ordered=TRUE)
PoD$fingerprints <- factor( PoD$fingerprints, levels = c("GFP_pos.2m_Srxn1", "GFP_pos.2m_DDIT3", "GFP_pos.2m_p21"), ordered=TRUE)

PoDsel <- PoD[ value >= 0.0]

PoDMin<-PoDsel[ , min(log10Dose), by = c("fingerprints"  ,"treatment", "DILI_Type")]
whichMultiple <- PoDMin[ ,length(V1) > 1 , by = list(treatment)]
whichMultiple <- whichMultiple[ V1 == TRUE]
PoDMinM <- PoDMin[ treatment%in% whichMultiple$treatment]
setnames(PoDMinM, "V1", "PoD")
setnames(PoDMin, "V1", "PoD")


p <- ggplot(data = PoDMin, aes(x = fingerprints, y = PoD,color=DILI_Type)) +  geom_jitter() + 
  scale_color_manual(values = as.character(myColors$DILI_Type))
p + theme_sharp() 


p <- ggplot(data = PoDMinM, aes(x = fingerprints, y = PoD,color=DILI_Type)) +  geom_jitter() + 
  scale_color_manual(values = as.character(myColors$DILI_Type))
p + theme_sharp() + geom_line(aes(group=treatment), alpha=0.5) +geom_text(aes(label=treatment), size=4)





# plot Srxn1 vs Chop activation with DILI annotation


SrVSCh <- PoD[!fingerprints%in%c("GFP_pos.2m_p21")]
SrVSChMax<- SrVSCh[ , max(value), by = list(treatment, fingerprints,  DILI_Type)]
require(reshape2)
setnames(SrVSChMax, "V1", "value")

SrVSChMax <- dcast( data=SrVSChMax, treatment + DILI_Type ~ fingerprints,  value.var = "value")
SrVSCh <- dcast( data=SrVSCh, treatment+dose_uM+cmax+log10Dose+DILI_Type ~ fingerprints,  value.var = "value")


p <- ggplot(data = SrVSCh, aes(x = GFP_pos.2m_DDIT3, y = GFP_pos.2m_Srxn1,color=DILI_Type, shape=cmax) ) +geom_point(size =6) +
  scale_color_manual(values = as.character(myColors$DILI_Type)) + geom_line(aes(group=treatment), alpha=0.5)
p + theme_sharp() 



head(SrVSChMax)
p <- ggplot(data = SrVSChMax, aes(x = GFP_pos.2m_DDIT3, y = GFP_pos.2m_Srxn1,color=DILI_Type)) +geom_point(size = 6) +
  scale_color_manual(values = as.character(myColors$DILI_Type))
p + theme_sharp() 


# TGP PHH figure

TGP_d <- read.table(file = "H:/DILI screen/meta analyse DILI screen/DILI paper 3 genes.txt", sep ="\t", header =T)
#select 24 hours high
TGP_d<-as.data.table(TGP_d)
TGP_d <- TGP_d[dose%in%"High" & time=="24hr"]
TGP_d[, symbols:=factor(symbols, levels= c("SRXN1","DDIT3","CDKN1A"), ordered=TRUE) ]
TGP_d <- TGP_d[, list(treatment, symbols, logFC)]

# add imaging data: per treatment the max over dose for each gene
PoDmax <- PoD[, max(value), by = list(treatment, fingerprints, DILI_Type)]
setnames(PoDmax, "V1","GFP.2m")
PoDmax[, symbols:=fingerprints]
PoDmax[, fingerprints:=NULL]
PoDmax[, symbols:=gsub("GFP_pos.2m_DDIT3", "DDIT3", symbols)]
PoDmax[, symbols:=gsub("GFP_pos.2m_Srxn1", "SRXN1", symbols)]
PoDmax[, symbols:=gsub("GFP_pos.2m_p21", "CDKN1A", symbols)]

# manually change names to imaging convention
indConv <- unique(PoDmax$treatment) %in% unique(TGP_d$treatment)
unique(PoDmax$treatment)[!indConv]

unique(TGP_d$treatment) # veranderen

TGP_d$treatment <- gsub("WY-14643", "WY14643",TGP_d$treatment)
TGP_d$treatment <- gsub("buthionine sulfoximine", "buthionine sulfoxamine",TGP_d$treatment)
TGP_d$treatment <- gsub("cyclosporine A", "cyclosporin A",TGP_d$treatment)
TGP_d$treatment <- gsub("erythromycin ethylsuccinate", "erythromycin",TGP_d$treatment)
TGP_d$treatment <- gsub("diethyl maleate", "DEM",TGP_d$treatment)

setkeyv(PoDmax, c("treatment" , "symbols" ))
setkeyv(TGP_d, c("treatment", "symbols"))
TGPImage<- PoDmax[ TGP_d ]
TGPImage<- TGPImage[!is.na(treatment)]
TGPImage<- TGPImage[ !is.na(DILI_Type)]

p <- ggplot(data = TGPImage, aes(x=GFP.2m, y= logFC)) + geom_point(aes(color=DILI_Type), size = 5) + 
              facet_wrap(~symbols) + theme_sharp()
p + scale_color_manual(values = as.character(myColors$DILI_Type))


TGPImage.l <- melt(TGPImage, measure.vars = c("GFP.2m", "logFC"))

TGPImage.l$symbols <- factor(TGPImage.l$symbols,  levels= c("SRXN1","DDIT3","CDKN1A"), ordered=TRUE) 
p <- ggplot(data = TGPImage.l, aes(x= reorder(treatment, -value), y= value, group = variable, fill= DILI_Type, color=variable))  +
  geom_bar(stat= "identity", position = "dodge")  
p <- p + facet_wrap(~symbols, ncol=1) +
   theme_sharp()
p + scale_fill_manual(values = as.character(myColors$DILI_Type))  + theme( axis.text.x = element_text(angle = 45, hjust = 1, size = 14, colour = "grey50") ) + 
  theme( strip.text.x = element_text( size = 14)) + scale_color_manual(values = c("blue","black")) +scale_size_manual(values=c(222,222))



## model for cubic spline heatmap (only for time course data)##
################


# loop over each cell line treatment dose combination for b-spline model. store and plot for qc. use model data for plotting means and for heatmap


# compounds verwijderen die ook in max-of-time heatmap zijn verwijderd

 indKeep <- GFP_data[ , treatment] %in% all.data[, treatment]

GFP_data <- GFP_data[indKeep,]

test <-as.data.frame(unique(GFP_data[, list(treatment,dose_uM, cmax)])[order(unique(GFP_data[, list(treatment,dose_uM, cmax)])$treatment)])
test
test<-as.data.table(test)

(na.omit(test)[, .N,by = treatment])[N!=5]

test[treatment =="dmso"]

# a fix for GFP_data for dose of dmso

GFP_data[ treatment %in% "dmso" & cmax == "100cmax", dose_uM:=0.2 ]
GFP_data[ treatment %in% "dmso" & cmax == "50cmax", dose_uM:=0.1 ]
GFP_data[ treatment %in% "dmso" & cmax == "10cmax", dose_uM:=0.02 ]
GFP_data[ treatment %in% "dmso" & cmax == "5cmax", dose_uM:=0.001 ]
GFP_data[ treatment %in% "dmso" & cmax == "1cmax", dose_uM:=0.002 ]

unique(GFP_data[ treatment %in% "dmso" & cmax == "5cmax", dose_uM ])



GFP_data[ treatment == "dextromethorphan hbr" & dose_uM == 0.0077, dose_uM := 0.0220]
GFP_data[ treatment == "dextromethorphan hbr" & dose_uM == 0.0385, dose_uM := 0.1100] 
GFP_data[ treatment == "dextromethorphan hbr" & dose_uM == 0.0770, dose_uM := 0.2200] 
GFP_data[ treatment == "dextromethorphan hbr" & dose_uM == 0.3850, dose_uM := 1.1000] 
GFP_data[ treatment == "dextromethorphan hbr" & dose_uM == 0.7700, dose_uM := 2.2000] 


# redo test above must be empty


# fix problem with bundled in 3m new cmax for Srxn1 (ICAM1 on this plate is fine)

newCmaxComps


GFP_pos.3m_Srxn1
GFP_pos.2m_Srxn1
GFP_pos.m3sd_Srxn1

# find the max and minimum values to change the fingerprint to correct one. (bundled in 3m atm)
GFP_data[ treatment %in% newCmaxComps & fingerprints %in% "GFP_pos.3m_Srxn1" & replID == "replicate_2" & timeAfterExposure < 6,
          min.val := min(value) == value,
#          .N,
          by = list(treatment, cmax, timeAfterExposure) ]

all(GFP_data[ min.val ==TRUE, .N, by = list(treatment, cmax, timeAfterExposure) ]$min.val == 1) # 380 rows
# assign min rows to GFP_pos.m3sd_Srxn1
GFP_data[ min.val ==TRUE, fingerprints := "GFP_pos.m3sd_Srxn1"]


GFP_data[ treatment %in% newCmaxComps & fingerprints %in% "GFP_pos.3m_Srxn1" & replID == "replicate_2" & timeAfterExposure < 6,
          max.val := max(value) == value,
          #          .N,
          by = list(treatment, cmax, timeAfterExposure) ]

all(GFP_data[ max.val ==TRUE, .N, by = list(treatment, cmax, timeAfterExposure) ]$max.val == 1) # 380 rows
GFP_data[ max.val == TRUE] # oops 381 rows...

GFP_data[ max.val == TRUE & min.val == TRUE]
GFP_data[ treatment %in% newCmaxComps & fingerprints %in% "GFP_pos.3m_Srxn1" & replID == "replicate_2" & timeAfterExposure < 6 
          & max.val == FALSE, ] # 379
 as.data.frame(GFP_data[ max.val == TRUE, .N, by = list(treatment, cmax, timeAfterExposure) ])
# cyclophosphamide   5cmax               1.0 2

 # add insignificant amount to unbreak even...
 GFP_data[treatment %in% "cyclophosphamide" & cmax == "5cmax" & fingerprints %in% "GFP_pos.3m_Srxn1" & replID == "replicate_2" & timeAfterExposure < 6
          & c(TRUE, FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE), value := 0.114660115 ]
 # two exactly the same values:  0.001592357
 

 GFP_data[treatment %in% "cyclophosphamide" & cmax == "5cmax" & fingerprints %in% "GFP_pos.3m_Srxn1" & replID == "replicate_2" & timeAfterExposure < 6    ]$value[1] <- 0.001592367
 GFP_data[treatment %in% "cyclophosphamide" & cmax == "5cmax" & fingerprints %in% "GFP_pos.3m_Srxn1" & replID == "replicate_2" & timeAfterExposure < 6    ] 
 
 GFP_data[ treatment %in% newCmaxComps & fingerprints %in% "GFP_pos.3m_Srxn1" & replID == "replicate_2" & timeAfterExposure < 6,
           max.val := max(value) == value,
           #          .N,
           by = list(treatment, cmax, timeAfterExposure) ]
 
 all(GFP_data[ max.val ==TRUE, .N, by = list(treatment, cmax, timeAfterExposure) ]$max.val == 1) # 380 rows
 GFP_data[ max.val == TRUE] # now the correct 380 rows
 GFP_data[ max.val ==TRUE, fingerprints := "GFP_pos.2m_Srxn1"]
 
 
GFP_data[ , min.val:=NULL]
GFP_data[ , max.val:=NULL]



GFP_data.list =list()

GFP_data[, fingerprints:=factor(fingerprints)]

cmaxDoseTreatment[treatment %in% naDoseTreats]

cmaxDoseTreatment[ cmax == "cmax1", cmax := "1cmax"]
cmaxDoseTreatment[ cmax == "cmax5", cmax := "5cmax"]
cmaxDoseTreatment[ cmax == "cmax10", cmax := "10cmax"]
cmaxDoseTreatment[ cmax == "cmax50", cmax := "50cmax"]
cmaxDoseTreatment[ cmax == "cmax100", cmax := "100cmax"]

naDoseTreats <- unique(GFP_data[ is.na(dose_uM), list(treatment, cmax)])
cmaxDoseTreatment
setkey(naDoseTreats, treatment, cmax)
setkey(cmaxDoseTreatment, treatment, cmax)

toFix <- cmaxDoseTreatment[naDoseTreats]



for( i in 1 : nrow(toFix)){
  
  GFP_data[ is.na(dose_uM) & treatment == toFix[i , treatment ] &
            cmax == toFix[i, cmax],
               dose_uM := toFix[ i, dose_uM]
            ]
  
}



GFP_data$splitL <- paste(GFP_data$treatment, GFP_data$dose_uM, GFP_data$replID, GFP_data$fingerprints )



GFP_data$splitL[[1]]

GFP_data.list <- split(GFP_data, c(GFP_data$splitL))
GFP_data.list[[1]]



#dir.create("figures/modelGraphsAugust2016")
#graphDir <- "figures/modelGraphsAugust2016/"

#fit the crashed low time point data with less degrees of freedom
# fit the normal ones with 8
ind<- lapply(GFP_data.list, function(x) { length(unique(x$timeAfterExposure)) }  )
head(ind)

indd <- lapply(ind, function(x) x < 12)
inddd<- unlist(indd)
length(inddd)
length(GFP_data.list)
sum(inddd)

GFP_data.list.lowTP <- GFP_data.list[inddd]
GFP_data.list <- GFP_data.list[!inddd]
length(GFP_data.list.lowTP) # 780 (should stay the same after the dose_uM NA fix)
length(GFP_data.list) # 20764 (should become more after the dose_uM NA fix)

model.results.GFP_data_lowTP = list()

#include the fitted parameters to perform a F-test
# perform t-test for fitted/sampled values after model fit




pdf(file = paste(graphDir, "model fit graphs_lowTPwithDMSO.pdf", sep ="/"), height = 6, width = 6)


for( i in seq_along(GFP_data.list.lowTP) ){ #)
 
  # acarbose 0.1502 replicate_1 GFP_pos.m3sd_DDIT3
  fm1 <- lm(value ~ ns(timeAfterExposure, df =4), data = GFP_data.list.lowTP[[i]])
  model.results.GFP_data_lowTP[[i]] <-  predict(fm1, dtTime<-data.frame(
    timeAfterExposure = round(seq(0.6,23, length.out=24), digits=2)))
  model.results.GFP_data_lowTP[[i]] <- as.data.frame(model.results.GFP_data_lowTP[[i]] )
  colnames(model.results.GFP_data_lowTP[[i]])<- "mod"
  
  model.results.GFP_data_lowTP[[i]]$treatment <- unique(GFP_data.list.lowTP[[i]]$treatment)
  model.results.GFP_data_lowTP[[i]]$dose_uM <- unique(GFP_data.list.lowTP[[i]]$dose_uM)
  model.results.GFP_data_lowTP[[i]]$cmax <- unique(GFP_data.list.lowTP[[i]]$cmax)
  model.results.GFP_data_lowTP[[i]]$replID <- unique(GFP_data.list.lowTP[[i]]$replID)
  model.results.GFP_data_lowTP[[i]]$timeAfterExposure <- dtTime$timeAfterExposure
  model.results.GFP_data_lowTP[[i]]$fingerprints <- unique(GFP_data.list.lowTP[[i]]$fingerprints)
  
   p <- ggplot(data =  model.results.GFP_data_lowTP[[i]], aes(x=timeAfterExposure, y = mod)) + geom_point() + 
      geom_point(data= GFP_data.list.lowTP[[i]], aes(x=timeAfterExposure, y = value, color = replID)) + 
   ggtitle(unique(GFP_data.list.lowTP[[i]]$splitL)) + ylim(c(-1.2,1.2))
  
  print(p)
  
  
}
dev.off()

# new compounds counts went wrong, need to fix first. plate part which is early all went into pos3m, not 2m and 3msde
# looks like only for srxn1 repl2

# plan:
# re-run this mess of a project from start: No (learn from this and organize analyses in smaller blocks with numbered scripts)
# so far the max in time was used, since only the crashed plate start time points are bundled in the 3m, this should not effect the
# max of time for the 2m since if a response was present for the new cmax compounds, this was found in the final time points with the correct data

# solution: tape-fix: split the triple points per time point for this set. perform the modeling. continue.

# next: modeling, wait with the time-course heatmaps (its ugly and too much data). perform fda statistics for time course plots, make time course plots.
# then work on the 3D plots (worked it out alrdy)
# then the machine learning: keep it simple. use the time info in simple ways: prob not early/late peek (because dynamics is different), maybe long format.

# a problem with dose_uM being NA


save.image("E:/DILI screen/meta analyse DILI screen/03082016bkup.Rdata")
save.image("03082016bkup.Rdata")


model.results.GFP_data = list()

pdf(file = paste(graphDir, "model fit graphs with DMSO.pdf", sep ="/"), height = 6, width = 8)
#nog een keer runnen (zonder pdf) - om fit p-waarden te extraheren
for( i in seq_along(GFP_data.list)){ #
  
  # acarbose 0.1502 replicate_1 GFP_pos.m3sd_DDIT3
  fm1 <- lm(value ~ ns(timeAfterExposure, df =8), data = GFP_data.list[[i]])
  model.results.GFP_data[[i]] <-  predict(fm1, dtTime<-data.frame(
    timeAfterExposure = round(seq(0.6,23, length.out=24), digits=2)))
  model.results.GFP_data[[i]] <- as.data.frame(model.results.GFP_data[[i]] )
  colnames(model.results.GFP_data[[i]])<- "mod"
  
  
  model.results.GFP_data[[i]]$treatment <- unique(GFP_data.list[[i]]$treatment)
  model.results.GFP_data[[i]]$dose_uM <- unique(GFP_data.list[[i]]$dose_uM)
  model.results.GFP_data[[i]]$cmax <- unique(GFP_data.list[[i]]$cmax)
  model.results.GFP_data[[i]]$replID <- unique(GFP_data.list[[i]]$replID)
  model.results.GFP_data[[i]]$timeAfterExposure <- dtTime$timeAfterExposure
  model.results.GFP_data[[i]]$fingerprints <- unique(GFP_data.list[[i]]$fingerprints)
  
  p <- ggplot(data =  model.results.GFP_data[[i]], aes(x=timeAfterExposure, y = mod)) + geom_line() + 
    geom_point(data= GFP_data.list[[i]], aes(x=timeAfterExposure, y = value, color = replID)) + 
    ggtitle(unique(GFP_data.list[[i]]$splitL)) + ylim(c(-1.2,1.2))
  
  print(p)
  
  
  if((i%%1000)==0){
   write.table(file = paste(graphDir, "progres", i ,".txt", sep =""), i/28521 )
  }
  
}
dev.off()


# dont forgot to combine low TP with high TP

save.image("03082016aftermodelrun.Rdata")
head(model.results.GFP_data[[3]])
model.results.GFP_data.df <- do.call("rbind", model.results.GFP_data)

model.results.GFP_data_lowTP.df <- do.call("rbind", model.results.GFP_data_lowTP)

model.results.GFP_data.df <- rbind(model.results.GFP_data.df, model.results.GFP_data_lowTP.df)

model.results.GFP_data.df <- as.data.table(model.results.GFP_data.df)


# statistics: DILI_fda.R
save(model.results.GFP_data.df, file = "data/model.results.GFP_data.df.Rdata")



# for the 3D plots: 24 time points (~ 1 per hour). do t-tests: pure for coloring the 3D plots
model.results.GFP_data_3D <- model.results.GFP_data.df[ timeAfterExposure %in% selTP]

testje <- model.results.GFP_data_3D[ fingerprints %in% "GFP_pos.3m_Srxn1" & cmax == "100cmax" & treatment %in% c("nefazodone","fenofibrate") ]

control100  <- model.results.GFP_data_3D[ fingerprints %in% "GFP_pos.3m_Srxn1" & treatment %in% "dmso" ] # where the flying fuck is my dmso???



testje[,  , by = c(mod, treatment, dose_uM, cmax, timeAfterExposure,fingerprints)]

output<- t.test(1:5, 2:6)


model.results.GFP_data.df[, ]



require(plyr)
model.results.GFP_data.df.m <- ddply(model.results.GFP_data.df, .(treatment, dose_uM, fingerprints, timeAfterExposure),
                            summarize, meanR = mean(mod, na.rm = TRUE))

head(model.results.GFP_data.df.m)
dim(model.results.GFP_data.df.m)

GFP_data$splitL <- paste(GFP_data$fingerprints , GFP_data$treatment, GFP_data$dose_uM )
GFP_data.list <- split(GFP_data, c(GFP_data$splitL))      

model.results.GFP_data.df.m$splitL <- paste(model.results.GFP_data.df.m$fingerprints , model.results.GFP_data.df.m$treatment, model.results.GFP_data.df.m$dose_uM )
model.results.GFP_data.l.m <- split(model.results.GFP_data.df.m, model.results.GFP_data.df.m$splitL)



#pdf(file = paste(graphDir, "model fit graphs means.pdf", sep ="/"), height = 10, width = 12) 
#for(i in seq_along(model.results.GFP_data.l.m)){
#  p <- ggplot(data =  model.results.GFP_data.l.m[[i]], aes(x=timeAfterExposure, y = meanR)) + geom_point() + 
 #   geom_point(data= GFP_data.list[[i]], aes(x=timeAfterExposure, y = value, color = replID)) + 
 #   ggtitle(unique(GFP_data.list[[i]]$splitL)) + ylim(c(0,1))
 # print(p)
#}
#dev.off()


# results of model stored in:
head(model.results.GFP_data.df.m)

# create dose levels
head(model.results.GFP_data.df.m)


model.results.GFP_data.df.m <- model.results.GFP_data.df.m[  order( model.results.GFP_data.df.m[ , "fingerprints"],
                                                          model.results.GFP_data.df.m[ , "treatment"],
                                                          model.results.GFP_data.df.m[ , "dose_uM"]),   ]

counts.d <- ddply(model.results.GFP_data.df.m, .(fingerprints, treatment,timeAfterExposure ), summarize,count.d.l = length(dose_uM))
head(counts.d)
head(counts.d)
model.results.GFP_data.df.m <- as.data.table(model.results.GFP_data.df.m)

model.results.GFP_data.df.m$dose.f <- NA
model.results.GFP_data.df.m$dose.f <- as.integer(model.results.GFP_data.df.m$dose.f)
setkeyv(model.results.GFP_data.df.m, c("treatment", "timeAfterExposure", "fingerprints"))

for (i in 1 : nrow(counts.d))
{

model.results.GFP_data.df.m[ list(counts.d$treatment[i], counts.d$timeAfterExposure[i], counts.d$fingerprints[i]),
                                 dose.f:=gl(counts.d$count.d.l[i], 1)]
  
  
#     model.results.GFP_data.df.m[ treatment == counts.d$treatment[i] &
#                                 timeAfterExposure == counts.d$timeAfterExposure[i] &
#                                 fingerprints == counts.d$fingerprints[i] ,
#                               dose.f:= gl(counts.d$count.d.l[i], 1)
#                                 ]
#   model.results.GFP_data.df.m$dose.f[ model.results.GFP_data.df.m$treatment == counts.d$treatment[i] & 
#                                        model.results.GFP_data.df.m$timeAfterExposure == counts.d$timeAfterExposure[i] &
#                                        model.results.GFP_data.df.m$fingerprints == counts.d$fingerprints[i]] <- gl(counts.d$count.d.l[i], 1)
  
 
  
 
 if((i%%1000)==0){
    model.results.GFP_data.df.m[1:100,]
  }
  if((i%%100)==0){
  print(paste(i/703000, "  i=:", i))
}
}


as.data.frame(model.results.GFP_data.df.m[, length(unique(dose_uM)), by ="treatment"])

# remove control compounds
controlC <- c("CCCP_C", "CDDO_C", "chloroquine_C","cisplatin_C","DMOG_C","thapsigargin_C","g.f. DMEM")
model.results.GFP_data.df.m<- model.results.GFP_data.df.m[!treatment %in% controlC]
# check if there are 5 doses per treatment for non control compounds
any(is.na(model.results.GFP_data.df.m$dose.f))
head(model.results.GFP_data.df.m)

model.results.GFP_data.df.m[190:210,]

model.results.GFP_data.df.m$dose.f<- as.factor(model.results.GFP_data.df.m$dose.f)
# mean of reps

all.fingerprints<- unique(model.results.GFP_data.df.m[, fingerprints])
unique(model.results.GFP_data.df.m[, dose.f])



# fix the dmso concentration
unique(model.results.GFP_data.df.m[ treatment %in% "DMSO", dose_uM])

model.results.GFP_data.df.m[, dose_uM:=as.character(dose_uM)]
model.results.GFP_data.df.m[ treatment%in%"DMSO" , dose_uM:= gsub("0.004","0.002",dose_uM)]

model.results.GFP_data.df.m[ treatment%in%"DMSO" , dose_uM:= gsub("0.006","0.010",dose_uM)]
model.results.GFP_data.df.m[ treatment%in%"DMSO" , dose_uM:= gsub("0.008","0.010",dose_uM)]

model.results.GFP_data.df.m[ treatment%in%"DMSO" , dose_uM:= gsub("0.03","0.02",dose_uM)]
model.results.GFP_data.df.m[ treatment%in%"DMSO" , dose_uM:= gsub("0.04","0.02",dose_uM)]

model.results.GFP_data.df.m[ treatment%in%"DMSO" , dose_uM:= gsub("0.06","0.1",dose_uM)]
model.results.GFP_data.df.m[ treatment%in%"DMSO" , dose_uM:= gsub("0.08","0.1",dose_uM)]
#remove 0.300 0.400 0.600 0.800
model.results.GFP_data.df.m<- model.results.GFP_data.df.m[ !(treatment%in%"DMSO" & dose_uM %in% c("0.3", "0.4", "0.6", "0.8"))]

model.results.GFP_data.df.m[ treatment%in%"DMSO" & dose_uM %in% "0.002", dose.f:=as.factor(1)]
model.results.GFP_data.df.m[ treatment%in%"DMSO" & dose_uM %in% "0.01", dose.f:=as.factor(2)]
model.results.GFP_data.df.m[ treatment%in%"DMSO" & dose_uM %in% "0.02", dose.f:=as.factor(3)]
model.results.GFP_data.df.m[ treatment%in%"DMSO" & dose_uM %in% "0.1", dose.f:=as.factor(4)]
model.results.GFP_data.df.m[ treatment%in%"DMSO" & dose_uM %in% "0.2", dose.f:=as.factor(5)]
model.results.GFP_data.df.m[ , dose.f:= factor(dose.f)]

unique(model.results.GFP_data.df.m[treatment%in%"DMSO", dose.f])

# calculate mean of DMSO's
DMSOpart<- model.results.GFP_data.df.m[treatment%in%"DMSO"]
model.results.GFP_data.df.m<-model.results.GFP_data.df.m[!treatment%in%"DMSO"]
DMSOpart[, treat.dose:=NULL]
DMSOpart[, splitL:=NULL]
DMSOpart[, dose_uM:=as.numeric(dose_uM)]
DMSOpartS<-DMSOpart[, lapply(.SD, mean(meanR), na.rm=TRUE), by =list(treatment,fingerprints,timeAfterExposure,dose.f)]
model.results.GFP_data.df.m[, treat.dose:=NULL]
model.results.GFP_data.df.m[, splitL:=NULL]

model.results.GFP_data.df.m<- rbind(model.results.GFP_data.df.m, DMSOpartS)

#for( i in seq_along(all.fingerprints)) {
 # sel.data <- model.results.GFP_data.df.m[ fingerprints %in% all.fingerprints[i]]
 # pdf(paste(" time course", all.fingerprints[i], ".pdf", sep =""), height = 40, width = 40)  
#p<-ggplot( sel.data, aes(x=timeAfterExposure, y = meanR, color = dose.f)) +
#  geom_line(aes(group = dose.f )) + theme_sharp()

#p <- p + facet_wrap( ~treatment )
#print(p)
#dev.off()
#
#}

model.results.GFP_data.df.m
#save.image(file = "22_06_2015.RData")
dir()
load("22_06_2015.RData")
getwd()


#### sensitivity plot
# plot sensitivity data with time course data
all.data.withDose
unique(model.results.GFP_data.df.m[, fingerprints])
unique(model.results.GFP_data.df.m[, treatment])

treatment %in% "diclofenac"  & dose.f %in% 5 &

sel.sens <- model.results.GFP_data.df.m[ 
                                          fingerprints %in%
                                            c("GFP_pos.2m_p21","GFP_pos.3m_p21","GFP_pos.m3sd_p21",
                                              "mmnFC_Nuclei_Intensity_IntegratedIntensity_Image_GFP_p21",
                                              "mmnFC_Nuclei_Intensity_MeanIntensity_Image_GFP_p21") ]



c( "GFP_pos.2m_Srxn1",   "GFP_pos.3m_Srxn1", "GFP_pos.m3sd_Srxn1", 
   "mmnFC_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1", 
   "mmnFC_Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1")

c("GFP_pos.2m_DDIT3", "GFP_pos.3m_DDIT3","GFP_pos.m3sd_DDIT3",  
  "mmnFC_Nuclei_Intensity_IntegratedIntensity_Image_GFP_DDIT3",
  "mmnFC_Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3") 



sel.sens
# all time courses again with nice grid and thicker lines
all.feats <- unique(model.results.GFP_data.df.m[, fingerprints])

sel.feats <- c("GFP_pos.2m_Srxn1", "GFP_pos.2m_DDIT3", "GFP_pos.2m_p21")


length(all.treats)/3

120+58
61+58
model.results.GFP_data.df.m$treatment<- as.character(model.results.GFP_data.df.m$treatment)
model.results.GFP_data.df.m$treatment <- factor(model.results.GFP_data.df.m$treatment, levels =
                                                    sort(as.character(unique(model.results.GFP_data.df.m$treatment))), ordered= TRUE    )

all.treats <- sort(as.character(unique(model.results.GFP_data.df.m$treatment)))

1:60
61:119
120:length(all.treats)

sel.comps <- sort(all.treats)[120:length(all.treats)]
head(model.results.GFP_data.df.m)
myColors$DILI_Type

for(i in seq_along(all.feats)) {

  sel.sens <- model.results.GFP_data.df.m[ 
    fingerprints %in%
      sel.feats & treatment %in% sel.comps ]
  sel.sens$fingerprints <- factor(sel.sens$fingerprints, levels = sel.feats, ordered = TRUE)

pdf(file = "allCells 3 of 3 figure.pdf", height = 100, width = 12)
p<- ggplot(data=sel.sens, aes(x = timeAfterExposure, y = meanR, color = fingerprints)) + facet_grid(  treatment ~ dose.f) +
  geom_line( size = 2)  + scale_color_manual(values=myColors$fingerprint[1:3]) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
p <- p + coord_cartesian( ylim =c(0,1.1))
print(p)
dev.off()
}



theme_sharp() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.4, size = 12 , 
                                    colour = "grey50") ) + theme( strip.text.x = element_text( )) +
  theme(plot.title = element_text(lineheight=.8, size = 8 )) + theme(legend.text=element_text(size=6))




# sensitivity figure:
sel.sens.srxn1 <- model.results.GFP_data.df.m[ 
  treatment %in% c("DMSO" ,"DEM", "hydroxyurea") & dose.f %in% 4 & fingerprints %in% 
    c( "GFP_pos.2m_Srxn1",   "GFP_pos.3m_Srxn1", "GFP_pos.m3sd_Srxn1", 
       "mmnFC_Cytoplasm_Intensity_IntegratedIntensity_Image_GFP_Srxn1", 
       "mmnFC_Cytoplasm_Intensity_MeanIntensity_Image_GFP_Srxn1")
     ]

sel.sens.ddits  <- model.results.GFP_data.df.m[
  treatment %in% c( "DMSO", "thapsigargin","nitrofurantoin" ) & dose.f %in% 4  & fingerprints %in%
    c("GFP_pos.2m_DDIT3", "GFP_pos.3m_DDIT3","GFP_pos.m3sd_DDIT3",  
      "mmnFC_Nuclei_Intensity_IntegratedIntensity_Image_GFP_DDIT3",
      "mmnFC_Nuclei_Intensity_MeanIntensity_Image_GFP_DDIT3") ]

sel.sens.p21 <- model.results.GFP_data.df.m[
  treatment %in% c("DMSO", "etoposide " ,"clozapine") & dose.f %in% 5 & fingerprints %in%
    c("GFP_pos.2m_p21","GFP_pos.3m_p21","GFP_pos.m3sd_p21",
      "mmnFC_Nuclei_Intensity_IntegratedIntensity_Image_GFP_p21",
      "mmnFC_Nuclei_Intensity_MeanIntensity_Image_GFP_p21")]

unique(sel.sens.p21[, treatment])

sel.sens<- rbind(sel.sens.srxn1, sel.sens.ddits, sel.sens.p21)

indNot <- sel.sens$treatment %in% "DMSO"

sel.sens[grepl("Srxn1", fingerprints), treatment:= paste(treatment, "Srxn1")]
sel.sens[grepl("p21", fingerprints), treatment:= paste(treatment, "p21")]
sel.sens[grepl("DDIT3", fingerprints), treatment:= paste(treatment, "DDIT3")]
sel.sens[!indNot, treatment:= paste(treatment, dose_uM)]


unique(sel.sens[, fingerprints])

sel.sens[, fingerprints:= gsub("_Srxn1","", fingerprints)]
sel.sens[, fingerprints:= gsub("_DDIT3", "", fingerprints)]
sel.sens[, fingerprints:= gsub("_p21", "", fingerprints)]
sel.sens[, fingerprints:= gsub("_Cytoplasm_Intensity", "", fingerprints)]
sel.sens[, fingerprints:= gsub("_Nuclei_Intensity", "", fingerprints)]

sel.sens[, treatment:= factor(treatment, levels =
                               c("DMSO Srxn1", "hydroxyurea Srxn1 17421.6","DEM Srxn1 500",
                                "DMSO DDIT3", "nitrofurantoin DDIT3 300","thapsigargin DDIT3 5",  # 50 cmax
                                "DMSO p21","clozapine p21 98", "etoposide  p21 404.3"), ordered = TRUE )]  # 100cmax


p<- ggplot(data=sel.sens, aes(x = timeAfterExposure, y = meanR)) + facet_wrap(  ~treatment, ncol=1 ) +
  geom_line( aes(color=fingerprints),size=1)  +
  theme_sharp() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.4, size = 12 , 
                                    colour = "grey50") ) + theme( strip.text.x = element_text( )) +
  theme(plot.title = element_text(lineheight=.8, size = 8 )) + theme(legend.text=element_text(size=8))

p <- p + coord_cartesian( ylim =c(0,1.1))
print(p)


# create row-clustered heatmaps per fingerprint. Rows are treatment and columns are ordered time points
# head(model.results.df.m)
# dim(model.results.df.m)
# require(reshape2)

model.results.GFP_data.df.m$timeAfterExposure <- round(model.results.GFP_data.df.m$timeAfterExposure, digits= 2)


#write.table(all.comps, file = "all.comps.txt", sep = "\tab", row.names = F)

#dit is al gedaan (naar bovne gekopieerd)
treatmentAnot <- read.table(file = "H:/DILI screen/meta analyse DILI screen/DILI anotation.txt", sep ="\t", header = FALSE)
colnames(treatmentAnot) <- c("treatment", "DILI_Type")

unique(treatmentAnot$DILI_Type)
treatmentAnot$DILI_Type<-gsub("No-DILI-Concern\\?", "No-DILI-Concern", treatmentAnot$DILI_Type)
treatmentAnot$DILI_Type<-gsub("No-DILI-concern", "No-DILI-Concern",treatmentAnot$DILI_Type)
treatmentAnot$DILI_Type<-gsub("Most-DILI-Concern\\?", "Most-DILI-Concern",treatmentAnot$DILI_Type)
rownames(treatmentAnot) <- treatmentAnot$treatment
treatmentAnot$treatment<-NULL

DILI_Type <- c(brewer.pal(9,"YlOrRd")[c(1,3,6)],"pink")
names(DILI_Type) <- c("No-DILI-Concern","Less-DILI-Concern", "Most-DILI-Concern", "control")

myColors <- list(DILI_Type = DILI_Type ) # color for dili classification side bar



my.colors <- brewer.pal(9, "YlOrRd")
colBreaks <-c(seq(-0.01, 0.2,length.out=20), seq(0.201,1.01, length.out = 9))
length(colBreaks)



head(model.results.GFP_data.df.m)
model.results.GFP_data.df.m[ , treatment:=gsub("carbamazapine", "carbamazepine", treatment)]
model.results.GFP_data.df.m[ , treatment:=gsub("cyclohexamine", "cycloheximide", treatment)]      
model.results.GFP_data.df.m[ , treatment:=gsub("etoposide ", "etoposide", treatment)]    

model.results.GFP_data.df.m[, dose_uM:= as.character(dose_uM)]

model.results.GFP_data.df.m[, treat.dose:=paste(treatment, dose_uM)]


testCell<- model.results.GFP_data.df.m[ fingerprints %in% all.fingerprints[1]  , list(treatment, dose_uM, treat.dose, meanR, timeAfterExposure, dose.f)]
require(reshape2)

testCell.w <- dcast(data =testCell, treat.dose+treatment~timeAfterExposure, value.var = "meanR")
dim(testCell.w)

rownames(testCell.w) <- testCell.w$treat.dose
testCell.w$treat.dose<- NULL

head(testCell.w)
#names(treatmentAnotL) <- c("top10_CDKN1A" , "top10_DDIT3", "top10_HSPA5" , "top10_SRXN1" )
ind <- match(testCell.w$treatment, rownames(treatmentAnot))
testCell.w$treatment <- NULL

treatmentAnot <- treatmentAnot[ind, ,drop = FALSE]
treatmentAnot<- data.frame(DILI_Type = treatmentAnot)

require(NMF)
pdf("test.pdf", width= 12, height = 22)
aheatmap(as.matrix(testCell.w), Rowv = TRUE, Colv=NA, 
         color = my.colors, breaks = colBreaks,
         fontsize = 10, width = 10, height=10, legend = TRUE,
         distfun = "euclidean", hclustfun ="ward",annRow = treatmentAnot,annColors = myColors
)
dev.off()


?hclust
# 1) maak afstand tussen paarsgewijze fingerprints gebaseerd op alle compound time courses, doe dit per compound.dose en dan gemiddelde en visa versa

dist.treat.dose = list()

norm.data = list()


for(i in  seq_along(all.fingerprints)) {
 
  buffer.cells <- model.results.GFP_data.df.m[ fingerprints %in% all.fingerprints[i]]
  time.vectors.perCell <-  dcast(data = buffer.cells, treat.dose~timeAfterExposure, value.var = "meanR" )
  
  rownames(time.vectors.perCell)<- time.vectors.perCell$treat.dose
  time.vectors.perCell$treat.dose <- NULL
  # noise is leading to clustering problems -->
  # per-timepoint background normalize data using DMSO 75% (closest to high concentration DMSO level)
  time.vectors.perCell<- as.matrix(time.vectors.perCell)
  #backgrDMSO75 <- time.vectors.perCell[ rownames(time.vectors.perCell) %in% "DMSO 75",]
  #time.vectors.perCell <- t(apply(time.vectors.perCell, MARGIN= 1,function(x) x - backgrDMSO75))
  
  dist.treat.dose[[i]] <- dist(time.vectors.perCell, method = "euclidean", diag=TRUE)
  
  
  norm.data[[i]] <- melt(time.vectors.perCell)
  norm.data[[i]]$fingerprints <- all.fingerprints[i]
  colnames(norm.data[[i]]) <- c("treat.dose", "timeAfterExposure", "meanR", "fingerprint")
}




all.norm.data <- do.call("rbind", norm.data)
class(all.norm.data)
head(all.norm.data)




dist.treat.dose <- lapply(dist.treat.dose, function(x) x<- as.matrix(x))
length(dist.treat.dose)

my.array<- array(unlist(dist.treat.dose), c(886,886,19))

lapply(dist.treat.dose,  dim)
19*886*886-length(unlist(dist.treat.dose))  

dist.treat.dose.df <- apply(my.array, 1:2, mean)
rownames(dist.treat.dose.df) <- rownames(dist.treat.dose[[1]])
colnames(dist.treat.dose.df) <- colnames(dist.treat.dose[[1]])
dist.treat.dose.df.d <- as.dist(dist.treat.dose.df)
treat.dose.clust<-hclust(dist.treat.dose.df.d, method = "ward.D")
names(treat.dose.clust)
treat.dose.clust$labels
treat.dose.clust$order
results.treat.dose.order <- treat.dose.clust$labels[treat.dose.clust$order]
pdf("cluster comp.dose.pdf", height= 10, width = 130)
plot(treat.dose.clust)
dev.off()
?dist

# same thing for the compounds
comp.doses<- unique(model.results.GFP_data.df.m$treat.dose)


head(all.norm.data)
dist.cells = list()
for(i in  seq_along(comp.doses)) {
  
  buffer.comp.dose <- all.norm.data[ all.norm.data$treat.dose %in% comp.doses[i],]
  time.vectors.perCompound <-  dcast(data = buffer.comp.dose, fingerprint~timeAfterExposure, value.var = "meanR" )
  rownames(time.vectors.perCompound)<- time.vectors.perCompound$fingerprint
  time.vectors.perCompound$fingerprint <- NULL
  
  dist.cells[[i]] <- dist(time.vectors.perCompound, method = "manhattan",  diag=TRUE)
  
}

dist.cells <- lapply(dist.cells, function(x) x<- as.matrix(x))

length(dist.cells)
my.array<- array(unlist(dist.cells), c(19,19,886))


dist.cells.df <- apply(my.array, 1:2, mean)
rownames(dist.cells.df) <- rownames(dist.cells[[1]])
colnames(dist.cells.df) <- colnames(dist.cells[[1]])
dist.cells.df.d <- as.dist(dist.cells.df)



cells.clust<-hclust(dist.cells.df.d, method = "ward.D")
names(cells.clust)
cells.clust$labels
cells.clust$order
results.cell.order <- cells.clust$labels[cells.clust$order]
results.cell.order
plot(cells.clust)

?dist
?hclust

# 
# 
# aheatmap(as.matrix((my.data.w.all)), Rowv = NULL, Colv=NULL, 
#           annRow = treatmentAnot, annCol = reporterAnot,  color = my.colors,
#          annColors = myColors,  fontsize = 10, width = 10, height=10, legend = FALSE
#          
#   )


results.treat.dose.order
results.cell.order



min(all.norm.data$meanR)
max(all.norm.data$meanR)
mean(all.norm.data$meanR)
my.colors

display.brewer.all()
colorRampPalette(brewer.pal(9, "Set1"))
my.colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(50)
colBreaks <-c(seq(-0.32, 0.15,length.out=20), seq(0.155,1.15, length.out = 31))



# treatment annotation needed with treat.dose

treatmentAnot
treat.treat.dose <- model.results.GFP_data.df.m[, list(treatment, treat.dose)]
treat.treat.dose<- unique(treat.treat.dose)
setkey(treat.treat.dose, "treatment")
treatmentAnot<- as.data.table(treatmentAnot)
setkey(treatmentAnot, "treatment")
treatmentAnot.treat.dose <- treatmentAnot[treat.treat.dose]

write.table(treatmentAnot.treat.dose, file = "treatment.doseAnot.txt", sep ="\t", col.names=T)


Anot.treat.dose <- treatmentAnot.treat.dose[, list(DILI_Type)]
for( i in seq_along(all.fingerprints)){
  
  data.cell <- all.norm.data[ all.norm.data$fingerprint %in% all.fingerprints[i], ]
  my.data.w <-  dcast(data = data.cell, treat.dose~timeAfterExposure, value.var = "meanR" )
  
  rownames(my.data.w) <- my.data.w$treat.dose
  my.data.w$treat.dose<- NULL
  ind <- match(results.treat.dose.order, rownames(my.data.w))
  my.data.w<- my.data.w[ind, ]
  if(
    !all(rownames(my.data.w)==results.treat.dose.order)
  ){
    stop()
  }
  
  ind2 <- match( rownames(my.data.w),treatmentAnot.treat.dose$treat.dose)
  
  treatmentAnot.treat.dose <- treatmentAnot.treat.dose[ind2,]
  
  pdf(file = paste(all.fingerprints[i], ".pdf", sep =""), width = 4, height = 16)
  
  
  aheatmap(as.matrix((my.data.w)), Rowv = NA, Colv=NA, 
           annRow = Anot.treat.dose,  breaks = colBreaks, color = my.colors,
           annColors = myColors,  fontsize = 10, legend = TRUE
  )
  
  dev.off()
}







#################







++++++


  
  
  

# plot including cell tox parameters cell count and speed



time course voor n paar comps:
  unique(three_mean$treatment)
  selComps <- three_mean[treatment %in%  c( "aflatoxin B1", "doxorubicin", "iodoacetamide", "FCCP", "DMSO"),]
unique(selComps$dose_uM[ selComps$treatment == c( "aflatoxin B1", "doxorubicin", "iodoacetamide", "FCCP", "DMSO")[5]])
unique(selComps$treatment)

selComps <- selComps[ dose_uM %in% c(10,11.80, 0.100, 0.008)]


selComps <- selComps[ (cell_line %in% "Srxn1" & treatment %in% c("iodoacetamide", "FCCP", "DMSO")) |
                        (cell_line %in% "p21" & treatment %in% c("aflatoxin B1", "doxorubicin", "DMSO"))  ]


unique(selComps$plateID)
selComps <- selComps[ !plateID %in% "2013_08_28_niii_Srxn1" ]
selComps <- selComps[ !plateID %in% "2014_02_05_p21_niii" ]

selComps <- selComps[ plateID %in% c("2014_02_05_nii_p21"  ,"2013_09_04_Srxn1", "2014_03_05_p21", "2013_08_15_Srxn1")]

selComps[ treatment %in% "DMSO"]

p<- ggplot(data=selComps , aes(x = timeAfterExposure, y = value , color = cell_line)) + geom_point(size = 2) + 
  facet_wrap( ~treatment, scales = "free_x" ) +
   theme_sharp() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.4, size = 12 , 
                                    colour = "grey50") ) + theme( strip.text.x = element_text( )) +
  ggtitle(   "Time  course DDR OS") + 
  theme(plot.title = element_text(lineheight=.8, size = 14 )) + theme(legend.text=element_text(size=12))
p <-p+geom_smooth(aes( color = cell_line), method = 'loess', se = FALSE, span = 0.5)  + theme(legend.position = "bottom")
p
unique(selComps$plateID)








# select detective data (time course data & dose response data & cytotox data)

goldCompounds <- c("acetaminophen", "aflatoxin B1", "allyl alcohol", "amiodarone", "Beta-Naphthoflavone", "bosentan",
                   "CCl4" , "chlorpromazine", "DMNQ", "dirlotapide", "FCCP" , "fluoxetine" , "iodoacetamide", "methotrexate",
                   "oligomycin A", "rifampicin", "rotenone" , "TO901317", "tamoxifen", "valproic acid",
                   "carbachol", "doxorubicin")

test # for time course
all.data # for dose response and cytotox data


goldTime <- test[treatment %in% goldCompounds]
goldDose <- all.data[ treatment %in% goldCompounds]
length(unique(goldTime$treatment))
length(unique(goldDose$treatment))
length(goldCompounds)

write.table(goldTime, file = "raw gold compounds time course data.txt", sep ="\t", row.names = FALSE)
write.table(goldDose, file = "gold compounds dose response data.txt", sep ="\t", row.names = FALSE)



# comment 21 april 2016: copy and paste code that is maybe of relevance for quality control figure:


##TODO: this data (GFP_data) could be used for time based heatmaps (use code from POC screen data)
# summarizing the time dimension:
unique(GFP_data[, fingerprints])
#1) maxima of gfp based measurements (response peak can be early or later)
#2) cell number: time end/ time begin
#3) cell speed: very jumpy data when lower amoun of cells, ideal would be to capture the slope in time

# slope and mean calculation of cell speed
# first remove first time point speed zero
GFP_data<-GFP_data[ !(fingerprints %in% "mmnFC_CellSpeed" & timeID == 1)]
speedData <- GFP_data[ fingerprints %in% "mmnFC_CellSpeed" ]

speedData$splitL <- paste(speedData$plateID,speedData$cell_line , speedData$treatment, speedData$dose_uM,speedData$replID, speedData$cmax )

speedData.list <- split(speedData, c(speedData$splitL))
model.results = list()


pdf(file = paste("H:/DILI screen/meta analyse DILI screen/graphs/model speed data/speed model fit graphs.pdf", sep =""), height = 6, width = 6)

for( i in seq_along(speedData.list)){
  
  fm1 <- lm(value ~ timeAfterExposure, data = speedData.list[[i]])
  model.results[[i]] <-  predict(fm1, dtTime<-data.frame(
    timeAfterExposure = round(seq(2,23, length.out=200), digits=2))) # prediciton for qc
  
  model.results[[i]] <- as.data.frame(model.results[[i]] )
  colnames(model.results[[i]])<- "mod"
  
  model.results[[i]]$cell_line <- unique(speedData.list[[i]]$cell_line)  
  model.results[[i]]$treatment <- unique(speedData.list[[i]]$treatment)
  model.results[[i]]$dose_uM <- unique(speedData.list[[i]]$dose_uM)
  model.results[[i]]$replID <- unique(speedData.list[[i]]$replID)
  model.results[[i]]$plateID <- unique(speedData.list[[i]]$plateID)
  model.results[[i]]$cmax <- unique(speedData.list[[i]]$cmax)
  model.results[[i]]$timeAfterExposure <- dtTime$timeAfterExposure
  model.results[[i]]$slope <- fm1$coefficients[[2]]
  
  p <- ggplot(data =  model.results[[i]], aes(x=timeAfterExposure, y = mod)) + geom_point() + 
    geom_point(data= speedData.list[[i]], aes(x=timeAfterExposure, y = value, color = replID)) + 
    ggtitle(unique(speedData.list[[i]]$splitL)) + ylim(c(0,1))
  
  print(p)
}
dev.off()

model.results.t <- do.call('rbind', model.results)
model.results.t <- as.data.table(model.results.t)

# normalize model results

require(reshape2)

model.results.t.l <- melt(data = model.results.t, id.vars = c("treatment", "cmax", "replID", "plateID"), measure.vars = c( "mod","slope") )
setnames(model.results.t.l,"variable","fingerprints")
model.results.t.l<-normFun(data= model.results.t.l, normVar = "slope")
model.results.t.l<-normFun(data= model.results.t.l, normVar = "mod")
model.results.t.l<- model.results.t.l[fingerprints %in% c("mmnFC_slope", "mmnFC_mod")]
min(model.results.t.l[,value])
model.results.t.l[,fingerprints:=factor(fingerprints)]
unique(model.results.t.l[,fingerprints])
model.results.t.l[, fingerprints:= gsub("mmnFC_mod", "mmnFC_CellSpeed", fingerprints)]
model.results.t.l[, fingerprints:= gsub("mmnFC_slope", "mmnFC_CellSpeed slope", fingerprints)]
dim(model.results.t.l)
model.results.t.l <- model.results.t.l[, list(treatment, cmax, fingerprints, value)]
model.results.t.l.summarized <- model.results.t.l[  , lapply(.SD, mean, na.rm=TRUE), by = c( "treatment", "fingerprints", "cmax" )] #mean of speed slopes

max(model.results.t.l.summarized[, value])

speedData.summarized <- model.results.t.l.summarized  # completed speed data (slopes)












# now for slopes and means:


nonGFP_data <- GFP_data[ fingerprints %in% c("mmnFC_Nuclei_AreaShape_Area",
                                            "mmnFC_Nuclei_Intensity_MeanIntensity_Image_Hoechst",
                                            "mmnFC_CellNumber")]

nonGFP_data$splitL <- paste(nonGFP_data$plateID,nonGFP_data$cell_line , nonGFP_data$treatment, nonGFP_data$cmax,nonGFP_data$replID, nonGFP_data$fingerprints )

nonGFP_data.list <- split(nonGFP_data, c(nonGFP_data$splitL))

model.results = list()


pdf(file = paste("H:/DILI screen/meta analyse DILI screen/graphs/model speed data/nucleiArea hoechstInt cellNumber model fit graphs.pdf", sep =""), height = 6, width = 6)

for( i in seq_along(nonGFP_data.list)){
  
  fm1 <- lm(value ~ timeAfterExposure, data = nonGFP_data.list[[i]])
  model.results[[i]] <-  predict(fm1, dtTime<-data.frame(
    timeAfterExposure = round(seq(2,23, length.out=200), digits=2))) # prediciton for qc
  
  model.results[[i]] <- as.data.frame(model.results[[i]] )
  colnames(model.results[[i]])<- "mod"
  
  model.results[[i]]$cell_line <- unique(nonGFP_data.list[[i]]$cell_line)  
  model.results[[i]]$treatment <- unique(nonGFP_data.list[[i]]$treatment)
  model.results[[i]]$cmax <- unique(nonGFP_data.list[[i]]$cmax)
  model.results[[i]]$replID <- unique(nonGFP_data.list[[i]]$replID)
  model.results[[i]]$plateID <- unique(nonGFP_data.list[[i]]$plateID)
  model.results[[i]]$timeAfterExposure <- dtTime$timeAfterExposure
  model.results[[i]]$slope <- fm1$coefficients[[2]]
  model.results[[i]]$fingerprints <- unique(nonGFP_data.list[[i]]$fingerprints)
  
  p <- ggplot(data =  model.results[[i]], aes(x=timeAfterExposure, y = mod)) + geom_point() + 
    geom_point(data= nonGFP_data.list[[i]], aes(x=timeAfterExposure, y = value, color = replID)) + 
    ggtitle(unique(nonGFP_data.list[[i]]$splitL)) + ylim(c(0,1))
  
  print(p)
}
dev.off()

model.results.t <- do.call('rbind', model.results)
model.results.t <- as.data.table(model.results.t)

#hierhier
model.results.t <- model.results.t[, list(mod,  treatment, cmax, slope, fingerprints)]
model.results.t.s <- model.results.t[  , lapply(.SD, mean, na.rm=TRUE), by = c( "fingerprints" ,"treatment", "cmax" )] #mean of 
model.results.t.sd <- model.results.t[  , lapply(.SD, sd, na.rm=TRUE), by = c( "fingerprints" ,"treatment", "cmax" )] #sd of 



nonGFP_data.l <- melt(data = model.results.t.s, id.vars = c("treatment", "cmax", "fingerprints") )

nonGFP_data.l[, fingerprints:= paste(fingerprints, variable)]

nonGFP_data.l[, variable:=NULL]

unique(GFP_data.s[, fingerprints])

unique(nonGFP_data.l[, fingerprints])
nonGFP_data.l[, fingerprints:=gsub("mod", "mean", fingerprints)]
