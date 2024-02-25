#Statistical test results:

#This code is to check the results of phylogenetic metrics using statistical tests.

#Coded by M. A. Thanuja M. Fernando

#Last updated on August 18, 2023

#Acknowledgement: Idea of using "foreach package" in this  code by Matthew Orton 

#Load packages
#install.packages("foreach")
library(foreach)
library(phangorn)
library(reshape) #melt
library(ez)
library(ggpubr)
library(rstatix)#outliers
library(lme4)
library(rmcorr)
#install.packages("lmerTest")
library(lmerTest)


#Read trees generated from RAxML into R. Read all of the trees at once ensuring I have the correct path to our small_dataset file. Here the path is /home/thanu/Desktop/FishData/small_dataset and the pattern is a regular expression to extract out the files I need

files <- list.files(path="/home/thanu/Desktop/FishData/New_dataset/STRATIFIED", pattern="*_[0-9]{1,2}_new", full.names=TRUE, recursive=FALSE)
#foreach package to just iterate through each file
treeList <- foreach(i=1:length(files)) %do% read.tree(files[i])

#Now I have a list with each element being its own separate tree, so treeList[[1]] is my first tree 20_1_new, treeList[[2]] is 20_2_new etc. Now I name each list element by substituting out the path from the names of each file

names(treeList) <- gsub("/home/thanu/Desktop/FishData/New_dataset/STRATIFIED/", "", files)

#This will let us keep track of which list element corresponds to which tree file

names(treeList)

#See the names of the files that correspond to each list element and can reference any list element by name, for example typing:

treeList$`20_1_new`

#Then I keep the complete full tree ("RAxML_bestTree.complete100") a separate variable since I'm using it for all the analyses below. I call it parseTreeComp

parseTreeComp <- read.tree("RAxML_bestTree.FR100_new")

#For Robinson-Foulds values, I use foreach loop to iterate through each tree (There were 40 trees)

RF_List <- foreach(i=1:length(treeList)) %do% RF.dist(parseTreeComp,treeList[[i]],normalize = TRUE)

#Same thing for path distances

P_List <- foreach(i=1:length(treeList)) %do% path.dist(parseTreeComp,treeList[[i]])

#Then reorder the lists such that the ordering goes 20_1, 40_1, 60_1, 80_1, 20_2, 40_2, 60_2, 80_2 etc.

order <- c(2, 12, 22, 32,# 20_1, 40_1, 60_1, 80_1
           3, 13, 23, 33, # 20_2, 40_2, 60_2, 80_2
           4, 14, 24, 34,
           5, 15, 25, 35,
           6, 16, 26, 36,
           7, 17, 27, 37,
           8, 18, 28, 38,
           9, 19, 29, 39,
           10, 20, 30, 40,
           1, 11, 21, 31) # 20_10, 40_10, 60_10, 80_10

#When plotting the percentage of backbone tree represent by the x axis. Assigned percentages to a variable called GRF

GRF <- c(20, 40, 60, 80)

#Here I have used a function to reorder the metric lists, group the values and assign values to variable names.From metric list (L_list) I'm going to group metric values (Robinson-Foulds values, CADM values, path distances etc.) (eg; 20_1, 40_1, 60_1, 80_1 as L_1). Then unlist it to get a vector and assign variable names for down stream analysis

metric <- function(L_List){
  #reordering
  L_List <- L_List[order]
  #group the values and assign values to variable names
  L_1 <- unlist(L_List[1:4])
  L_2 <- unlist(L_List[5:8])
  L_3 <- unlist(L_List[9:12])
  L_4 <- unlist(L_List[13:16])
  L_5 <- unlist(L_List[17:20])
  L_6 <- unlist(L_List[21:24])
  L_7 <- unlist(L_List[25:28])
  L_8 <- unlist(L_List[29:32])
  L_9 <- unlist(L_List[33:36])
  L_10 <- unlist(L_List[37:40])
  #generate a dataframe using backbone % and metric values
  results <- data.frame(GRF,L_1,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10)
}

#RF metric
#Now we can use the above function to generate dataframes in different metrics
#First, Robinson-Foulds metric,

RF_metric<- metric(RF_List)

#transpose the dataframe
trans_RF_metric <- data.frame(t(RF_metric))

#delete the 1st row of the dataframe
trans_RF_metric = trans_RF_metric[-1,]
dim(trans_RF_metric)

#adding unique id
trans_RF_metric$id <- cumsum(!duplicated(trans_RF_metric[1:2]))
colnames(trans_RF_metric) <- c("0.2","0.4","0.6","0.8","id")

#Melt the dataframe into a form suitable for easy casting
longdata <- melt(trans_RF_metric,
                 id="id",
                 measure=c("0.2","0.4","0.6","0.8"))
#assign column names
colnames(longdata)=c("id","sample_completeness","RF_distance")

longdata

#ANOVA test. dv means dependant variable, The one-way repeated measures ANOVA can be used to determine whether the means of longdata RF distances are significantly different among four backbone levels(sample completeness 20%,40%,60%,80%)
ezANOVA(data = longdata, dv = RF_distance, wid = id, within = sample_completeness)

#assumptions
#Outliers can be easily identified using box plot methods, implemented in the R function identify_outliers() [rstatix package].

longdata %>%
  group_by(sample_completeness) %>%
  identify_outliers(RF_distance)

#visualization
#Create a violin plot
#violin plot with dot plot (geom_dotplot())
violin_plot <- ggplot(longdata, aes(x=sample_completeness, y=RF_distance,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of RF distance for sampling completeness")
violin_plot

#Normality assumption
#The normality assumption can be checked by computing Shapiro-Wilk test for each sample completeness. If the data is normally distributed, the p-value should be greater than 0.05.

longdata %>%
  group_by(sample_completeness) %>%
  shapiro_test(RF_distance)

#QQ plot draws the correlation between a given data and the normal distribution. Create QQ plots for each level backbone tree(sample completeness)

ggqqplot(longdata, "RF_distance", facet.by = "sample_completeness",color ="sample_completeness")

#From the resultant plots, as all the points fall approximately along the reference line, we can assume normality

#The assumption of sphericity will be automatically checked during the computation of the ANOVA test using the R function anova_test() [rstatix package]. The Mauchly’s test is internally used to assess the sphericity assumption

#By using the function get_anova_table() [rstatix] to extract the ANOVA table, the Greenhouse-Geisser sphericity correction is automatically applied to factors violating the sphericity assumption

res.aov <- anova_test(data = longdata, dv = RF_distance, wid = id, within = sample_completeness)

get_anova_table(res.aov)

#Post-hoc tests

#perform multiple pairwise paired t-tests between the levels of the within-subjects factor (Sample_completeness). P-values are adjusted using the Bonferroni multiple testing correction method.

# pairwise comparisons

pwc <- longdata %>%
  pairwise_t_test(
    RF_distance ~ sample_completeness, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc


#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
longdat <- lapply(longdata,as.numeric)

#Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
rmcorr(id, RF_distance, sample_completeness, longdat)

#PathDistance
PathDistance_metric<- metric(P_List)

#transpose the dataframe
trans_PathDistance_metric <- data.frame(t(PathDistance_metric))
#delete the 1st row of the dataframe
trans_PathDistance_metric = trans_PathDistance_metric[-1,]

#adding unique id
trans_PathDistance_metric$id <- cumsum(!duplicated(trans_PathDistance_metric[1:2]))
colnames(trans_PathDistance_metric) <- c("0.2","0.4","0.6","0.8","id")

#Melt the dataframe into a form suitable for easy casting
longdata <- melt(trans_PathDistance_metric,
                 id="id",
                 measure=c("0.2","0.4","0.6","0.8"))
colnames(longdata)=c("id","sample_completeness","Path_Distance")
longdata

#The one-way repeated measures ANOVA can be used to determine whether the means of longdata path distances are significantly different among four backbone levels(sample completeness 20%,40%,60%,80%)
ezANOVA(data = longdata, dv = Path_Distance, wid = id, within = sample_completeness)

#check assumptions
#Outliers can be easily identified using box plot methods, implemented in the R function identify_outliers() [rstatix package]
longdata %>%
  group_by(sample_completeness) %>%
  identify_outliers(Path_Distance)
#There were no extreme outliers

#visualization
#Create a violin plot
violin_plot <- ggplot(longdata, aes(x=sample_completeness, y=Path_Distance,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of Path Distance for sampling completeness")
violin_plot

#Normality assumption
#The normality assumption can be checked by computing Shapiro-Wilk test for each time point. If the data is normally distributed, the p-value should be greater than 0.05.
longdata %>%
  group_by(sample_completeness) %>%  
  shapiro_test(Path_Distance)
#The longdata Path_distance was normally distributed at each sample completeness, as assessed by Shapiro-Wilk’s test (p > 0.05).

#QQ plot draws the correlation between a given data and the normal distribution. Create QQ plots for each level backbone tree(sample completeness)
ggqqplot(longdata, "Path_Distance", facet.by = "sample_completeness",color = "sample_completeness")

#From the resultant plots, as all the points fall approximately along the reference line, we can assume normality

#The assumption of sphericity will be automatically checked during the computation of the ANOVA test using the R function anova_test() [rstatix package]. The Mauchly’s test is internally used to assess the sphericity assumption
#By using the function get_anova_table() [rstatix] to extract the ANOVA table, the Greenhouse-Geisser sphericity correction is automatically applied to factors violating the sphericity assumption
# res.aov <- anova_test(data = longdata, dv = RF_distance, wid = id,within = sample_completeness)
# get_anova_table(res.aov)

#Post-hoc tests
#perform multiple pairwise paired t-tests between the levels of the within-subjects factor (Sample_completeness). P-values are adjusted using the Bonferroni multiple testing correction method.

# pairwise comparisons
pwc <- longdata %>%
  pairwise_t_test(
    Path_Distance ~ sample_completeness, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#All the pairwise differences are statistically significant

#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
longdat <- lapply(longdata,as.numeric)

#Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
rmcorr(id, Path_Distance, sample_completeness, longdat)

###Simplified function abort R, so let's try one by one
library(phangorn)
#Read trees into R generated from RAxML
aF <- read.tree("RAxML_bestTree.FR100_new")
bF <- read.tree("20_1_new")   
cF <- read.tree("20_2_new") 
dF <- read.tree("20_3_new") 
eF <- read.tree("20_4_new") 
fF <- read.tree("20_5_new") 
gF <- read.tree("20_6_new")   
hF <- read.tree("20_7_new")     
iF <- read.tree("20_8_new") 
jF <- read.tree("20_9_new") 
kF <- read.tree("20_10_new") 
lF <- read.tree("40_1_new") 
mF <- read.tree("40_2_new") 
nF <- read.tree("40_3_new") 
oF <- read.tree("40_4_new") 
pF <- read.tree("40_5_new") 
qF <- read.tree("40_6_new") 
rF <- read.tree("40_7_new") 
sF <- read.tree("40_8_new") 
tF <- read.tree("40_9_new") 
uF <- read.tree("40_10_new") 
vF <- read.tree("60_1_new") 
wF <- read.tree("60_2_new") 
xF <- read.tree("60_3_new") 
yF <- read.tree("60_4_new") 
zF <- read.tree("60_5_new") 
aaF <- read.tree("60_6_new") 
bbF <- read.tree("60_7_new") 
ccF <- read.tree("60_8_new")
ddF <- read.tree("60_9_new")
eeF <- read.tree("60_10_new")
ffF <- read.tree("80_1_new")
ggF <- read.tree("80_2_new")
hhF <- read.tree("80_3_new")
iiF <- read.tree("80_4_new")
jjF <- read.tree("80_5_new")
kkF <- read.tree("80_6_new")
llF <- read.tree("80_7_new")
mmF <- read.tree("80_8_new")
nnF <- read.tree("80_9_new")
ooF <- read.tree("80_10_new")

library(ape)
#Newick to distancematrix
DM_100F <- cophenetic(aF)
DM_20_1F <- cophenetic(bF)
DM_20_2F <- cophenetic(cF)
DM_20_3F <- cophenetic(dF)
DM_20_4F <- cophenetic(eF)
DM_20_5F <- cophenetic(fF)
DM_20_6F <- cophenetic(gF)
DM_20_7F <- cophenetic(hF)
DM_20_8F <- cophenetic(iF)
DM_20_9F <- cophenetic(jF)
DM_20_10F <- cophenetic(kF)
DM_40_1F <- cophenetic(lF)
DM_40_2F <- cophenetic(mF)
DM_40_3F <- cophenetic(nF)
DM_40_4F <- cophenetic(oF)
DM_40_5F <- cophenetic(pF)
DM_40_6F <- cophenetic(qF)
DM_40_7F <- cophenetic(rF)
DM_40_8F <- cophenetic(sF)
DM_40_9F <- cophenetic(tF)
DM_40_10F <- cophenetic(uF)
DM_60_1F <- cophenetic(vF)
DM_60_2F <- cophenetic(wF)
DM_60_3F <- cophenetic(xF)
DM_60_4F <- cophenetic(yF)
DM_60_5F <- cophenetic(zF)
DM_60_6F <- cophenetic(aaF)
DM_60_7F <- cophenetic(bbF)
DM_60_8F <- cophenetic(ccF)
DM_60_9F <- cophenetic(ddF)
DM_60_10F <- cophenetic(eeF)
DM_80_1F <- cophenetic(ffF)
DM_80_2F <- cophenetic(ggF)
DM_80_3F <- cophenetic(hhF)
DM_80_4F <- cophenetic(iiF)
DM_80_5F <- cophenetic(jjF)
DM_80_6F <- cophenetic(kkF)
DM_80_7F <- cophenetic(llF)
DM_80_8F <- cophenetic(mmF)
DM_80_9F <- cophenetic(nnF)
DM_80_10F <- cophenetic(ooF)
#function to generate CADM results
CADM <- function(DM){
  #add two matrics together
  AB <- rbind(DM_100F, DM)
  cadm_global <- CADM.global(AB, 2, 6796)$congruence_analysis[1,]
  return(cadm_global)
}
cadm_global <- CADM.global(AB, 2, 6796)
#cadm_AB <- CADM(DM_B)
#cadm_AB$congruence_analysis
#Use CADM function to compare backbone trees with the complete tree
CADM(DM_20_1F)
CADM(DM_20_2)
CADM(DM_20_3)
CADM(DM_20_4)
CADM(DM_20_5)
CADM(DM_20_6)
CADM(DM_20_7)
CADM(DM_20_8)
CADM(DM_20_9)
CADM(DM_20_10)
CADM(DM_40_1)
CADM(DM_40_2)
CADM(DM_40_3)
CADM(DM_40_4)
CADM(DM_40_5)
CADM(DM_40_6)
CADM(DM_40_7)
CADM(DM_40_8)
CADM(DM_40_9)
CADM(DM_40_10)
CADM(DM_60_1)
CADM(DM_60_2)
CADM(DM_60_3)
CADM(DM_60_4)
CADM(DM_60_5)
CADM(DM_60_6)
CADM(DM_60_7)
CADM(DM_60_8)
CADM(DM_60_9)
CADM(DM_60_10)
CADM(DM_80_1)
CADM(DM_80_2)
CADM(DM_80_3)
CADM(DM_80_4)
CADM(DM_80_5)
CADM(DM_80_6)
CADM(DM_80_7)
CADM(DM_80_8)
CADM(DM_80_9)
CADM(DM_80_10)













#CADM
#Newick to distance matrix. here I have used foreach loop to simplify the code 
DM_List <- foreach(i=1:length(treeList)) %do% (cophenetic(treeList[[i]])/max(cophenetic(treeList[[i]]))) 

#Distance matrix of 100% tree
DM_100 <- cophenetic(parseTreeComp)/max(cophenetic(parseTreeComp))  

#function to generate CADM results
CADM <- function(DM){
  #add two matrices together
  AB <- rbind(DM_100F, DM)
  cadm_global <- CADM.global(AB, 2, 4520)$congruence_analysis[1,]
  return(cadm_global)
}

#Use CADM function to compare backbone trees with the complete tree (foreach loop to iterate through CADM list)
CADM_List <- foreach(i=1:length(treeList)) %do% CADM(DM_List[[i]])

#Using the metric(), generate CADM dataframe
CADM_metric <- metric(CADM_List)

#transpose the dataframe
trans_CADM_metric <- data.frame(t(CADM_metric))

#delete the 1st row of the dataframe
trans_CADM_metric = trans_CADM_metric[-1,]

#adding unique id
trans_CADM_metric$id <- cumsum(!duplicated(trans_CADM_metric[1:2]))
colnames(trans_CADM_metric) <- c("0.2","0.4","0.6","0.8","id")

longdata <- melt(trans_CADM_metric,
                 id="id",
                 measure=c("0.2","0.4","0.6","0.8"))
colnames(longdata)=c("id","sample_completeness","CADM_metric")
longdata

ezANOVA(data = longdata, dv = CADM_metric, wid = id,within = sample_completeness)
#The one-way repeated measures ANOVA can be used to determine whether the means of longdata RF distances are significantly different among four backbone levels(sample completeness 20%,40%,60%,80%)

#check assumptions
#Outliers can be easily identified using box plot methods, implemented in the R function identify_outliers() [rstatix package].

longdata %>%
  group_by(sample_completeness) %>%
  identify_outliers(CADM_metric)


#visualization
#Create a violin plot
violin_plot <- ggplot(longdata, aes(x=sample_completeness, y=CADM_metric,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of CADM for sampling completeness")
violin_plot

#Normality assumption
#The normality assumption can be checked by computing Shapiro-Wilk test for each time point. If the data is normally distributed, the p-value should be greater than 0.05.
longdata %>%
  group_by(sample_completeness) %>%  
  shapiro_test(CADM_metric)


#The longdata RF_distance was normally distributed at each sample completeness, as assessed by Shapiro-Wilk’s test (p > 0.05).

#QQ plot draws the correlation between a given data and the normal distribution. Create QQ plots for each level backbone tree(sample completeness)
ggqqplot(longdata, "CADM_metric", facet.by = "sample_completeness", color = "sample_completeness")

#From the resultant plots, as all the points fall approximately along the reference line, we can assume normality
#The assumption of sphericity will be automatically checked during the computation of the ANOVA test using the R function anova_test() [rstatix package]. The Mauchly’s test is internally used to assess the sphericity assumption

#By using the function get_anova_table() [rstatix] to extract the ANOVA table, the Greenhouse-Geisser sphericity correction is automatically applied to factors violating the sphericity assumption
# res.aov <- anova_test(data = longdata, dv = RF_distance, wid = id,within = sample_completeness)
# get_anova_table(res.aov)

#Post-hoc tests
#perform multiple pairwise paired t-tests between the levels of the within-subjects factor (Sample_completeness). P-values are adjusted using the Bonferroni multiple testing correction method.

# pairwise comparisons
pwc <- longdata %>%
  pairwise_t_test(CADM_metric ~ sample_completeness, paired = TRUE, p.adjust.method = "bonferroni")
pwc

#All the pairwise differences are not statistically significant

#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
longdat <- lapply(longdata,as.numeric)

#Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
rmcorr(id, CADM_metric, sample_completeness, longdat)                                  





