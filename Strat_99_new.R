
#Load packages
library(tidyverse)
library(ape)
library(phylotools)
library(dplyr)
library(foreach)
library(phangorn)
getwd()
#Read the fish dataset in R
fish_4520 <- read.table("Fish_Realigned_multigene.phy")
# Install from CRAN
#install.packages("tidyverse")
library(tidyverse)
#Rename columns
colnames(fish_4520) <- c("species_name","seq")

#New dataframe by removing first row as it consists of number of species and the multigene seq length
fish_multigene <- fish_4520[-1,]

head(fish_multigene)

#Extract the range of COI gene to be used in EPA
fishsamplefull_co1 <- substr(fish_multigene$seq,2710,3391)
#Check the class
class(fishsamplefull_co1)
#Convert to dataframe
fishsamplefull_co1 <- data.frame(fishsamplefull_co1)
#Check the class
class(fishsamplefull_co1)
#Combine columns of COI dataframe and fishsample100 dataframe to get the sepecies name
fishsamplefull_co1_with_speciesname <- cbind(fish_multigene$species_name,fishsamplefull_co1)
#Assign column names
colnames(fishsamplefull_co1_with_speciesname) <- c("species_name","seq")
head(fishsamplefull_co1_with_speciesname)
library(phylotools)


CO1_20_1 <- read.fasta(file="CO1_20_1S.fasta")
CO1_20_2 <- read.fasta(file="CO1_20_2S.fasta")
CO1_20_3 <- read.fasta(file="CO1_20_3S.fasta")
CO1_20_4 <- read.fasta(file="CO1_20_4S.fasta")
CO1_20_5 <- read.fasta(file="CO1_20_5S.fasta")
CO1_20_6 <- read.fasta(file="CO1_20_6S.fasta")
CO1_20_7 <- read.fasta(file="CO1_20_7S.fasta")
CO1_20_8 <- read.fasta(file="CO1_20_8S.fasta")
CO1_20_9 <- read.fasta(file="CO1_20_9S.fasta")
CO1_20_10 <- read.fasta(file="CO1_20_10S.fasta")
colnames(CO1_20_1) <- c("species_name","seq")
colnames(CO1_20_2) <- c("species_name","seq")
colnames(CO1_20_3) <- c("species_name","seq")
colnames(CO1_20_4) <- c("species_name","seq")
colnames(CO1_20_5) <- c("species_name","seq")
colnames(CO1_20_6) <- c("species_name","seq")
colnames(CO1_20_7) <- c("species_name","seq")
colnames(CO1_20_8) <- c("species_name","seq")
colnames(CO1_20_9) <- c("species_name","seq")
colnames(CO1_20_10) <- c("species_name","seq")


#Remove 27 co1 species to get 1% co1 from fish dataset
set.seed(1000)
CO1_99_1<- CO1_20_1[sample(nrow(CO1_20_1), 27),]
CO1_99_2<- CO1_20_2[sample(nrow(CO1_20_2), 27),]
CO1_99_3<- CO1_20_3[sample(nrow(CO1_20_3), 27),]
CO1_99_4<- CO1_20_4[sample(nrow(CO1_20_4), 27),]
CO1_99_5<- CO1_20_5[sample(nrow(CO1_20_5), 27),]
CO1_99_6<- CO1_20_6[sample(nrow(CO1_20_6), 27),]
CO1_99_7<- CO1_20_7[sample(nrow(CO1_20_7), 27),]
CO1_99_8<- CO1_20_8[sample(nrow(CO1_20_8), 27),]
CO1_99_9<- CO1_20_9[sample(nrow(CO1_20_9), 27),]
CO1_99_10<- CO1_20_10[sample(nrow(CO1_20_10), 27),]

#Now get the 99% multi gene dataset by geting antijoin of full data set and 27 co1 data set
ST99_1S <- anti_join(fish_multigene,CO1_99_1,by="species_name")
ST99_2S <- anti_join(fish_multigene,CO1_99_2,by="species_name")
ST99_3S <- anti_join(fish_multigene,CO1_99_3,by="species_name")
ST99_4S <- anti_join(fish_multigene,CO1_99_4,by="species_name")
ST99_5S <- anti_join(fish_multigene,CO1_99_5,by="species_name")
ST99_6S <- anti_join(fish_multigene,CO1_99_6,by="species_name")
ST99_7S <- anti_join(fish_multigene,CO1_99_7,by="species_name")
ST99_8S <- anti_join(fish_multigene,CO1_99_8,by="species_name")
ST99_9S <- anti_join(fish_multigene,CO1_99_9,by="species_name")
ST99_10S <- anti_join(fish_multigene,CO1_99_10,by="species_name")
library(seqRFLP)
dataframe2fas(ST99_1S , file = "ST99_1.fasta")
dataframe2fas(ST99_2S , file = "ST99_2.fasta")
dataframe2fas(ST99_3S , file = "ST99_3.fasta")
dataframe2fas(ST99_4S , file = "ST99_4.fasta")
dataframe2fas(ST99_5S , file = "ST99_5.fasta")
dataframe2fas(ST99_6S , file = "ST99_6.fasta")
dataframe2fas(ST99_7S , file = "ST99_7.fasta")
dataframe2fas(ST99_8S , file = "ST99_8.fasta")
dataframe2fas(ST99_9S , file = "ST99_9.fasta")
dataframe2fas(ST99_10S , file = "ST99_10.fasta")

#Use the function to generate phylip files from the dataframe
dff <- function(xx){
  return(dat2phylip(xx,outfile = "s.phy"))
}
#generate phylip files
dff(ST99_10S)

#For the mafft alignment we need 1% (27 species) co1 in fasta
library(seqRFLP)
dataframe2fas(CO1_99_1 , file = "ST_CO1_1_1.fasta")
dataframe2fas(CO1_99_2 , file = "ST_CO1_1_2.fasta")
dataframe2fas(CO1_99_3 , file = "ST_CO1_1_3.fasta")
dataframe2fas(CO1_99_4 , file = "ST_CO1_1_4.fasta")
dataframe2fas(CO1_99_5 , file = "ST_CO1_1_5.fasta")
dataframe2fas(CO1_99_6 , file = "ST_CO1_1_6.fasta")
dataframe2fas(CO1_99_7 , file = "ST_CO1_1_7.fasta")
dataframe2fas(CO1_99_8 , file = "ST_CO1_1_8.fasta")
dataframe2fas(CO1_99_9 , file = "ST_CO1_1_9.fasta")
dataframe2fas(CO1_99_10 , file = "ST_CO1_1_10.fasta")

#Use the function to generate phylip files from the dataframe
#generate phylip files
dff(CO1_99_10)
