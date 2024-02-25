
#To get 99% from stratified sampling (A part of this code extracted from my "Strat_4520.R" code)
#Code for stratified samples

#Coded by M. A. Thanuja M. Fernando

#Nov 16,2023
#Last updated on JAn 14,2024
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
library(phylotools)
#Convert to phylip format
#dat2phylip(fishsamplefull_co1_with_speciesname,outfile= "fishsample4520_co1_with_speciesname.phy")

#Read the fish taxonomy spread sheet in R
fish_taxonomy <- read.csv("/home/thanu/Desktop/FishData/Stratified/PFCtaxonomy.csv")

#Extract family column and genus_species column for downstream analysis
df_fish_family_species <- fish_taxonomy[,c("family","genus.species")]
#rename columns
colnames(df_fish_family_species) <- c("family","species_name")
#Replace the space with "_" and covert the resultant list into a df
df_fish_family_species1 <- as.data.frame(gsub(' ','_',df_fish_family_species$species_name))
#Bind the colum "family" of fish_family_species dataset with new dataframe df_fish_family_species1
newdf_fish_taxon <- cbind(df_fish_family_species1,df_fish_family_species$family)
#rename columns
colnames(newdf_fish_taxon) <- c("species_name","family")
#merge fish taxon and sequences
df_fish <- merge(newdf_fish_taxon,fish_multigene,by="species_name")
colnames(df_fish)

# #remove the family column as it's not needed for phylip format
# df_fish<- df_fish[,-2]

#Checking the families with the number of species in each family
tt <- table(df_fish$family)
#Check the umber of families in fish sample
length(unique(df_fish$family))#358

#extract families where sample count is less than 5 and keep the details in a data frame
df_less5 <- df_fish[df_fish$family %in% names(tt[tt < 5]), ]
#remove the column of "family" from the df_less5 dataframe
# df_less5df <- subset( df_less5, select = -2 )
#remove the first row as it consist of unwanted data(just the number of nucleotides). remove the family column as it's not needed for phylip format
df_less5df0 <- df_less5[,-2]

#Get a data frame where the sample count is more than 5
df_fish_more5 <- anti_join(df_fish, df_less5, by='species_name')#4136

#function to generate data frames from fraction(80%,60%,40%,20%) family group. df means the data frame and fraction is the percentage
fraction <- function(df,fraction){
  df_frac <- df %>% 
    group_by(family) %>%
    sample_frac(fraction)
  return(df_frac)
}

#Get the sequentially dropped data frames for newly created function "fraction" using another function called df_frac. This is a bit tricky with the percentage. As I’ve done sequential data dropping technique, I’ve used 0.75 to get 80% from whole dataset. Likewise, 0.67 for 40% and 0.5 for 20%. (I’ve calculated this using (100*60)/80=75%)
df_frac <- function(df){
  df_99 <- fraction(df,.99)
  df_80 <- fraction(df_99,.81)
  df_60 <- fraction(df_80,.75)
  df_40 <- fraction(df_60,.67)
  df_20 <- fraction(df_40,.5)
  results <- list(df_99,df_80,df_60,df_40,df_20)
}

#call the function df_frac and replicate the sequantial dropping process 10 times
set.seed(2000)
frac_list <- replicate(10,df_frac(df_fish_more5))

#Remove the family column of all the dataframes in the list, frac_list (otherwise no way of converting the df to phylip format). Efficiently use the bracket function "[" here.
#frac_list2 <- lapply(frac_list, "[", -grep(c("family"), names(frac_list[[1]])))
frac_list2 <- lapply(frac_list, "[", -2)

#split the large list into 40 dataframes by assigning names as df_1, df_2, df_3... Here df1 to df 4 represent one cycle of sequential data dropping as 80%, 60%, 40%, 20% respectively. Then df 5 is the 80% of second cycle.
for(i in 1:length(frac_list2)) {
  assign(paste0("df_", i), frac_list2[[i]])
}

#A list of generated datarames
all_df_list <- list(df_1,df_2,df_3,df_4,df_5,df_6,df_7,df_8,df_9,df_10,df_11,df_12,df_13,df_14,df_15,df_16,df_17,df_18,df_19,df_20,df_21,df_22,df_23,df_24,df_25,df_26,df_27,df_28,df_29,df_30,df_31,df_32,df_33,df_34,df_35,df_36,df_37,df_38,df_39,df_40,df_41,df_42,df_43,df_44,df_45,df_46,df_47,df_48,df_49,df_50)

#we can use foreach loop to iterate through each data frame (There were 40 data frames)
fish_frac_list <- foreach(i=1:length(all_df_list)) %do% rbind(df_less5df0,all_df_list[[i]])

#This step is to get individual data frames form the list of data  frames (So, there will be 40 dataframes named as newdf1, newdf2 etc. note:all these are generated from randomly stratified samples). paste0 convert its arguments to character strings.
list2env(setNames(fish_frac_list, paste0("newdf", seq_along(fish_frac_list))),
         envir=.GlobalEnv)
#check the number of fish families in each dataframe. Now the new dataframes names are newdf1.... newdf40. These are dataframe with fsmilies less than 5 species plus families more than 5 species (Families more than 5 species were used to dtratify the sampling using family)

#----------------------------------------
#Attach family names to each dataframe
#Read the fish taxonomy spread sheet in R
fish_taxonomy <- read.csv("/home/thanu/Desktop/FishData/Stratified/PFCtaxonomy.csv")

#Extract family column and genus_species column for downstream analysis
df_fish_family_species <- fish_taxonomy[,c("family","genus.species")]
length(unique(df_fish_family_species$family))
#rename columns
colnames(df_fish_family_species) <- c("family","species_name")

#Replace the space with "_" and covert the resultant list into a df
df_fish_family_species1 <- as.data.frame(gsub(' ','_',df_fish_family_species$species_name))

#Bind the colum "family" of fish_family_species dataset with new dataframe df_fish_family_species1
newdf_fish_taxon <- cbind(df_fish_family_species1,df_fish_family_species$family)

# #rename columns
# colnames(newdf_fish_taxon) <- c("species_name","family")
# colnames(newdf1)<- c("species_name","seq")
# class(newdf_fish_taxon)
# #merge fish taxon and sequences
# #df_intersect <- intersect(newdf_fish_taxon,fish_multigene)
# df_fish_80 <- merge(newdf_fish_taxon,newdf1,by="species_name")
# length(unique(df_fish_80$family))
# 
# colnames(fishsample20_10) <- c("species_name","seq")
# #merge fish taxon and sequences
# #df_intersect <- intersect(newdf_fish_taxon,fish_multigene)
# df_fish_20 <- merge(newdf_fish_taxon,fishsample20_10,by="species_name")
# length(unique(df_fish_20$family))


#Use the function to generate phylip files from the dataframe
dff <- function(xx){
  return(dat2phylip(xx,outfile = "SN20_.phy"))
}
#generate phylip files
dff(newdf50)


library(seqRFLP)
#fasta files of bb datasets needed for mafft alignment
#Generate fasta files for aligned co1 dataset
dataframe2fas(newdf1,file="99_1S.fasta")
dataframe2fas(newdf6,file="99_2S.fasta")
dataframe2fas(newdf11,file="99_3S.fasta")
dataframe2fas(newdf16,file="99_4S.fasta")
dataframe2fas(newdf21,file="99_5S.fasta")
dataframe2fas(newdf26,file="99_6S.fasta")
dataframe2fas(newdf31,file="99_7S.fasta")
dataframe2fas(newdf36,file="99_8S.fasta")
dataframe2fas(newdf41,file="99_9S.fasta")
dataframe2fas(newdf46,file="99_10S.fasta")


#newdfF1 (multigene 99%) and fishsamplefull_co1_with_speciesname (Co1 100% marker) antijoin to get 1% of co1 that exist only in co1 100% not in multigene 99%
CO1_99_1S <- anti_join(fishsamplefull_co1_with_speciesname,newdf1,by="species_name")
names(newdf2)
names(fishsamplefull_co1_with_speciesname)
CO1_99_2S <- anti_join(fishsamplefull_co1_with_speciesname,newdf6,by="species_name")
CO1_99_3S <- anti_join(fishsamplefull_co1_with_speciesname,newdf11,by="species_name")
CO1_99_4S <- anti_join(fishsamplefull_co1_with_speciesname,newdf16,by="species_name")
CO1_99_5S <- anti_join(fishsamplefull_co1_with_speciesname,newdf21,by="species_name")
CO1_99_6S <- anti_join(fishsamplefull_co1_with_speciesname,newdf26,by="species_name")
CO1_99_7S <- anti_join(fishsamplefull_co1_with_speciesname,newdf31,by="species_name")
CO1_99_8S <- anti_join(fishsamplefull_co1_with_speciesname,newdf36,by="species_name")
CO1_99_9S <- anti_join(fishsamplefull_co1_with_speciesname,newdf41,by="species_name")
CO1_99_10S <- anti_join(fishsamplefull_co1_with_speciesname,newdf46,by="species_name")


#Genearte phylip files for aligned co1 dataset 
#Just change the dataset name (CO1_20_1...CO1_80_10) for the function dff (function dff leads to datarame to phylip format)
dff(CO1_99_10S)

#Generate fasta files for aligned co1 dataset
dataframe2fas(CO1_99_1S, file = "CO1_99_1S.fasta")
dataframe2fas(CO1_99_2S, file = "CO1_99_2S.fasta")
dataframe2fas(CO1_99_3S, file = "CO1_99_3S.fasta")
dataframe2fas(CO1_99_4S, file = "CO1_99_4S.fasta")
dataframe2fas(CO1_99_5S, file = "CO1_99_5S.fasta")
dataframe2fas(CO1_99_6S, file = "CO1_99_6S.fasta")
dataframe2fas(CO1_99_7S, file = "CO1_99_7S.fasta")
dataframe2fas(CO1_99_8S, file = "CO1_99_8S.fasta")
dataframe2fas(CO1_99_9S, file = "CO1_99_9S.fasta")
dataframe2fas(CO1_99_10S, file = "CO1_99_10S.fasta")


#The following need for the raxml-epa
#Convert fas to phylip format
#install.packages("seqmagick")
library(seqmagick)
library(seqinr)
library(phylotools)
#Here fish20_co1_80_1.fasta is the reference aligned sequence set which is an output of multigene 20% and co1 80% aligning by using muscle in Linux terminal
getwd()

fish20_co1_80_1S=read.fasta(file = "fish20_co1_80_1.fasta")
class(fish20_co1_80_1)

fish20_co1_80_2=read.fasta(file = "fish20_co1_80_2.fasta")
fish20_co1_80_3=read.fasta(file = "fish20_co1_80_3.fasta")
fish20_co1_80_4=read.fasta(file = "fish20_co1_80_4.fasta")
fish20_co1_80_5=read.fasta(file = "fish20_co1_80_5.fasta")
fish20_co1_80_6=read.fasta(file = "fish20_co1_80_6.fasta")
fish20_co1_80_7=read.fasta(file = "fish20_co1_80_7.fasta")
fish20_co1_80_8=read.fasta(file = "fish20_co1_80_8.fasta")
fish20_co1_80_9=read.fasta(file = "fish20_co1_80_9.fasta")
fish20_co1_80_10=read.fasta(file = "fish20_co1_80_10.fasta")

fish40_co1_60_1=read.fasta(file = "fish40_co1_60_1.fasta")
fish40_co1_60_2=read.fasta(file = "fish40_co1_60_2.fasta")
fish40_co1_60_3=read.fasta(file = "fish40_co1_60_3.fasta")
fish40_co1_60_4=read.fasta(file = "fish40_co1_60_4.fasta")
fish40_co1_60_5=read.fasta(file = "fish40_co1_60_5.fasta")
fish40_co1_60_6=read.fasta(file = "fish40_co1_60_6.fasta")
fish40_co1_60_7=read.fasta(file = "fish40_co1_60_7.fasta")
fish40_co1_60_8=read.fasta(file = "fish40_co1_60_8.fasta")
fish40_co1_60_9=read.fasta(file = "fish40_co1_60_9.fasta")
fish40_co1_60_10=read.fasta(file = "fish40_co1_60_10.fasta")

fish60_co1_40_1=read.fasta(file = "fish60_co1_40_1.fasta")
fish60_co1_40_2=read.fasta(file = "fish60_co1_40_2.fasta")
fish60_co1_40_3=read.fasta(file = "fish60_co1_40_3.fasta")
fish60_co1_40_4=read.fasta(file = "fish60_co1_40_4.fasta")
fish60_co1_40_5=read.fasta(file = "fish60_co1_40_5.fasta")
fish60_co1_40_6=read.fasta(file = "fish60_co1_40_6.fasta")
fish60_co1_40_7=read.fasta(file = "fish60_co1_40_7.fasta")
fish60_co1_40_8=read.fasta(file = "fish60_co1_40_8.fasta")
fish60_co1_40_9=read.fasta(file = "fish60_co1_40_9.fasta")
fish60_co1_40_10=read.fasta(file = "fish60_co1_40_10.fasta")

fish80_co1_20_1=read.fasta(file = "fish80_co1_20.fasta")
fish80_co1_20_2=read.fasta(file = "fish80_co1_20_2.fasta")
fish80_co1_20_3=read.fasta(file = "fish80_co1_20_3.fasta")
fish80_co1_20_4=read.fasta(file = "fish80_co1_20_4.fasta")
fish80_co1_20_5=read.fasta(file = "fish80_co1_20_5.fasta")
fish80_co1_20_6=read.fasta(file = "fish80_co1_20_6.fasta")
fish80_co1_20_7=read.fasta(file = "fish80_co1_20_7.fasta")
fish80_co1_20_8=read.fasta(file = "fish80_co1_20_8.fasta")
fish80_co1_20_9=read.fasta(file = "fish80_co1_20_9.fasta")
fish80_co1_20_10=read.fasta(file = "fish80_co1_20_10.fasta")

#Just change the dataset name (fish20_co1_80_1...fish80_co1_20_10) for the function dff (function dff leads to datarame to phylip format)
dff(fish80_co1_20_10)
