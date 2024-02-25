#Stat_Code for Stratified swTools( RAxMLepa and EPAng)
#Jan 15, 2024
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
getwd()

library(ape)
library(foreach)
library(phangorn)
library(ggplot2)

#First read Stratified epang trees_for RF distance
files <- list.files(path="/home/thanu/Desktop/FishData/New_dataset/STRATIFIED", pattern="^epa[0-9]{1,2}_[0-9]{1,2}S_new.nwk", full.names=TRUE, recursive=FALSE)
treeList <- foreach(i=1:length(files)) %do% read.tree(files[i])

# Now I have a list with each element being its own separate tree, so treeList[[1]] is my first tree 20_1_new, treeList[[2]] is 20_2_new etc.

# Also we can name each list element by substituting out our path from the names of each file
names(treeList) <- gsub("/home/thanu/Desktop/FishData/New_dataset/STRATIFIED/", "", files)

# This will let us keep track of which list element corresponds to which tree file
names(treeList) 

# we can see the names of the files that correspond to each list element and we can reference any list element by name, for example typing:
#treeList$`20_1_new`

# Then I keep the complete full tree ("RAxML_parsimonyTree.complete100") a separate variable since I'm using it for all the analyses below
# I call it parseTreeComp
parseTreeComp <- read.tree("RAxML_bestTree.FR100_new")

# For Robinson-Foulds values, we can use foreach loop to iterate through each tree (There were 40 trees)
RF_List <- foreach(i=1:length(treeList)) %do% RF.dist(parseTreeComp,treeList[[i]],normalize = TRUE)

# Same thing for path distances 
#P_List <- foreach(i=1:length(treeList)) %do% path.dist(parseTreeComp,treeList[[i]])

# Then reorder the lists such that the ordering goes 20_1, 40_1, 60_1, 80_1, 20_2, 40_2, 60_2, 80_2 etc. 
order <- c(2, 12, 22, 32,42,# 20_10, 40_10, 60_10, 80_10
           3, 13, 23, 33,43, # 20_2, 40_2, 60_2, 80_2
           4, 14, 24, 34,44, 
           5, 15, 25, 35,45,
           6, 16, 26, 36,46,
           7, 17, 27, 37,47,
           8, 18, 28, 38,48,
           9, 19, 29, 39,49,
           10, 20, 30, 40,50,
           1, 11, 21, 31,41) # 20_1, 40_1, 60_1, 80_1

# When plotting we need the percentage of backbone tree as the x axis. For that we can use,
GRF <- c(20, 40, 60, 80,99)

# Here I have used a function to reorder the metric lists, group the values and assign values to variable names.
# From metric list (L_list) I'm going to group metric values (Robinson-Foulds values, CADM values, path distances etc.) (eg; 20_1, 40_1, 60_1, 80_1 as L_1). Then unlist it to get a vector and assign variable names for down stream analysis

metric <- function(L_List){
  #reordering
  L_List <- L_List[order]
  #group the values and assign values to variable names
  L_1 <- unlist(L_List[1:5])
  L_2 <- unlist(L_List[6:10])
  L_3 <- unlist(L_List[11:15])
  L_4 <- unlist(L_List[16:20])
  L_5 <- unlist(L_List[21:25])
  L_6 <- unlist(L_List[26:30])
  L_7 <- unlist(L_List[31:35])
  L_8 <- unlist(L_List[36:40])
  L_9 <- unlist(L_List[41:45])
  L_10 <- unlist(L_List[46:50])
  #generate a dataframe using backbone % and metric values
  results <- data.frame(GRF,L_1,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10)
}

# Now we can use the above function to generate dataframes in different metrics
# First, Robinson-Foulds metric,
RF_metric<- metric(RF_List)#if need line plot then use this line and run the plot directly. do not convert to transpose, stacks..
RF_metric <- t(RF_metric)
RF_metric <- data.frame(RF_metric)
dim(RF_metric)
RF_metric <- RF_metric[-1,]
colnames(RF_metric) <- c("20","40","60","80","99")
RF_metric <- stack(RF_metric)
head(RF_metric)
colnames(RF_metric) <- c("RF_distance","sample_completeness")
library(ggpubr)
bp <- ggboxplot(RF_metric,x="sample_completeness",y="RF_distance", color = "sample_completeness",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#003366")
                ,main="Boxplot of RF distance for stratified samples based on sample completeness using EPAng",add= "point", width=0.8,outlier.shape = NA)
ggpar(bp,orientation ="reverse")
#-----------------------------------
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
colnames(trans_RF_metric) <- c("0.2","0.4","0.6","0.8","0.99","id")

#Melt the dataframe into a form suitable for easy casting
longdata <- melt(trans_RF_metric,
                 id="id",
                 measure=c("0.2","0.4","0.6","0.8","0.99"))
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

#----------------

#Stratified RF distance by RAxMLEPA

RF<- read.csv("RF_dist_sample_types.csv")
#In this csv file Random,stratified and biased results are included. As I need Stratified results first I subset another dataframe from Correct_position dataset
Strati_RF <- RF[51:100,]
colnames(Strati_RF) <- c("Sample_type","RF_dist")
#As we need an id to perform ezANOVA I'm adding an id column to this dataframe. In RF random sample dataset (longdata) we have id so I combine longdata to Random_Correct_position

#id <- c(rep(1:10,5))
Strati_RF <- cbind(longdata, Strati_RF)
Strati_RF <- Strati_RF[,c("id","sample_completeness","RF_dist")]
class(Strati_RF)
colnames(Strati_RF) <- c("id","sample_completeness","RF_distance")
ezANOVA(data = Strati_RF, dv = RF_distance, wid = id, within = sample_completeness)

#check assumptions
#Outliers can be easily identified using box plot methods, implemented in the R function identify_outliers() [rstatix package]
Strati_RF %>%
  group_by(sample_completeness) %>%
  identify_outliers(RF_distance)
#There were no extreme outliers

#visualization
#Create a violin plot
violin_plot <- ggplot(Strati_RF, aes(x=sample_completeness, y=RF_distance,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of Correct position percentage for sampling completeness")
violin_plot

#Normality assumption
#The normality assumption can be checked by computing Shapiro-Wilk test for each time point. If the data is normally distributed, the p-value should be greater than 0.05.
Strati_RF %>%
  group_by(sample_completeness) %>%  
  shapiro_test(RF_distance)
#The longdata Path_distance was normally distributed at each sample completeness, as assessed by Shapiro-Wilk’s test (p > 0.05).

#QQ plot draws the correlation between a given data and the normal distribution. Create QQ plots for each level backbone tree(sample completeness)
ggqqplot(Strati_RF, "RF_distance", facet.by = "sample_completeness",color = "sample_completeness")

#From the resultant plots, as all the points fall approximately along the reference line, we can assume normality

#The assumption of sphericity will be automatically checked during the computation of the ANOVA test using the R function anova_test() [rstatix package]. The Mauchly’s test is internally used to assess the sphericity assumption
#By using the function get_anova_table() [rstatix] to extract the ANOVA table, the Greenhouse-Geisser sphericity correction is automatically applied to factors violating the sphericity assumption
# res.aov <- anova_test(data = longdata, dv = RF_distance, wid = id,within = sample_completeness)
# get_anova_table(res.aov)

#Post-hoc tests
#perform multiple pairwise paired t-tests between the levels of the within-subjects factor (Sample_completeness). P-values are adjusted using the Bonferroni multiple testing correction method.

# pairwise comparisons
pwc <- Strati_RF %>%
  pairwise_t_test(
    RF_distance ~ sample_completeness, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#All the pairwise differences are statistically significant

#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
# Strati_RF_new <- lapply(Strati_RF,as.numeric)#here SAmple type has 20S not only numeric values. So have to get only numeric values. otherwise rmcorr gives an error. As longdata has that. I combine that column from long data

Strati_RF_new <- lapply(Strati_RF,as.numeric)
#Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
rmcorr(id, RF_distance, sample_completeness, Strati_RF_new )

#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
# Stratified_Correct_position_new <- lapply(Stratified_Correct_position,as.numeric)
# 
# #Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
# rmcorr(id, correct_placement_percentage , sample_completeness, Stratified_Correct_position_new )

#non-parametric for repeated messures anova is friedman test

res.fried <- Strati_RF %>% friedman_test(RF_distance ~ sample_completeness |id)
res.fried
#The correct species percentage was statistically significantly different at the different time points during the diet, X2(2) = 30, p = 0.00000138

#effect size
#The Kendall’s W can be used as the measure of the Friedman test effect size. It is calculated as follow : W = X2/N(K-1); where W is the Kendall’s W value; X2 is the Friedman test statistic value; N is the sample size. k is the number of measurements per subject (M. T. Tomczak and Tomczak 2014).

#The Kendall’s W coefficient assumes the value from 0 (indicating no relationship) to 1 (indicating a perfect relationship).

#Kendall’s W uses the Cohen’s interpretation guidelines of 0.1 - < 0.3 (small effect), 0.3 - < 0.5 (moderate effect) and >= 0.5 (large effect). Confidence intervals are calculated by bootstap.
Strati_RF %>% friedman_effsize(RF_distance ~ sample_completeness |id)
#A large size is detected W=1

#multiple pair-wise comparisons
# From the output of the Friedman test, we know that there is a significant difference between groups, but we don’t know which pairs of groups are different.
# 
# A significant Friedman test can be followed up by pairwise Wilcoxon signed-rank tests for identifying which groups are different.
# 
# Note that, the data must be correctly ordered by the blocking variable (id) so that the first observation for time t1 will be paired with the first observation for time t2, and so on.
# 
# Pairwise comparisons using paired Wilcoxon signed-rank test. P-values are adjusted using the Bonferroni multiple testing correction method.
pwc <- Strati_RF%>%
  wilcox_test(RF_distance ~ sample_completeness, paired = TRUE, p.adjust.method = "bonferroni")
pwc


#-----------------------------------

Correct_position<- read.csv("Strati_SWtools_sister.csv")
dim(Correct_position)
#In this csv file RAxMLEPA and EPAng results are included. As I need EPAngresults first I subset another dataframe from Correct_position dataset
EPAng_Correct_position <- Correct_position[51:100,]
head(EPAng_Correct_position,11)
names(EPAng_Correct_position)
#As we need an id to perform ezANOVA I'm adding an id column to this dataframe. In RF random sample dataset (longdata) we have id so I combine longdata to Random_Correct_position

#id <- c(rep(1:10,4))
EPAng_Correct_position<- cbind(longdata, EPAng_Correct_position)
class(EPAng_Correct_position)
names(EPAng_Correct_position)
head(EPAng_Correct_position)
ezANOVA(data =EPAng_Correct_position, dv = Correct_percentage, wid = id, within = sample_completeness)

#check assumptions
#Outliers can be easily identified using box plot methods, implemented in the R function identify_outliers() [rstatix package]
EPAng_Correct_position %>%
  group_by(sample_completeness) %>%
  identify_outliers(Correct_percentage)
#There were no extreme outliers

#visualization
#Create a violin plot
violin_plot <- ggplot(EPAng_Correct_position, aes(x=sample_completeness, y=Correct_percentage,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of Correct position percentage for sampling completeness")
violin_plot

#Normality assumption
#The normality assumption can be checked by computing Shapiro-Wilk test for each time point. If the data is normally distributed, the p-value should be greater than 0.05.
EPAng_Correct_position %>%
  group_by(sample_completeness) %>%  
  shapiro_test(Correct_percentage)
#The longdata Path_distance was normally distributed at each sample completeness, as assessed by Shapiro-Wilk’s test (p > 0.05).

#QQ plot draws the correlation between a given data and the normal distribution. Create QQ plots for each level backbone tree(sample completeness)
ggqqplot(EPAng_Correct_position, "Correct_percentage", facet.by = "sample_completeness",color = "sample_completeness")

#From the resultant plots, as all the points fall approximately along the reference line, we can assume normality

#The assumption of sphericity will be automatically checked during the computation of the ANOVA test using the R function anova_test() [rstatix package]. The Mauchly’s test is internally used to assess the sphericity assumption
#By using the function get_anova_table() [rstatix] to extract the ANOVA table, the Greenhouse-Geisser sphericity correction is automatically applied to factors violating the sphericity assumption
# res.aov <- anova_test(data = longdata, dv = RF_distance, wid = id,within = sample_completeness)
# get_anova_table(res.aov)

#Post-hoc tests
#perform multiple pairwise paired t-tests between the levels of the within-subjects factor (Sample_completeness). P-values are adjusted using the Bonferroni multiple testing correction method.

# pairwise comparisons
pwc <- EPAng_Correct_position %>%
  pairwise_t_test(
    Correct_percentage ~ sample_completeness, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#All the pairwise differences are statistically significant

#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
EPAng_Correct_position_new <- lapply(EPAng_Correct_position,as.numeric)

#Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
rmcorr(id, Correct_percentage , sample_completeness, EPAng_Correct_position_new )

#non-parametric for repeated messures anova is friedman test

res.fried <- EPAng_Correct_position %>% friedman_test(Correct_percentage ~ sample_completeness |id)
res.fried
#The correct species percentage was statistically significantly different at the different time points during the diet, X2(2) = 30, p = 0.00000138

#effect size
#The Kendall’s W can be used as the measure of the Friedman test effect size. It is calculated as follow : W = X2/N(K-1); where W is the Kendall’s W value; X2 is the Friedman test statistic value; N is the sample size. k is the number of measurements per subject (M. T. Tomczak and Tomczak 2014).

#The Kendall’s W coefficient assumes the value from 0 (indicating no relationship) to 1 (indicating a perfect relationship).

#Kendall’s W uses the Cohen’s interpretation guidelines of 0.1 - < 0.3 (small effect), 0.3 - < 0.5 (moderate effect) and >= 0.5 (large effect). Confidence intervals are calculated by bootstap.
EPAng_Correct_position %>% friedman_effsize(Correct_percentage ~ sample_completeness |id)
#A large size is detected W=1

#multiple pair-wise comparisons
# From the output of the Friedman test, we know that there is a significant difference between groups, but we don’t know which pairs of groups are different.
# 
# A significant Friedman test can be followed up by pairwise Wilcoxon signed-rank tests for identifying which groups are different.
# 
# Note that, the data must be correctly ordered by the blocking variable (id) so that the first observation for time t1 will be paired with the first observation for time t2, and so on.
# 
# Pairwise comparisons using paired Wilcoxon signed-rank test. P-values are adjusted using the Bonferroni multiple testing correction method.
pwc <- EPAng_Correct_position%>%
  wilcox_test(Correct_percentage ~ sample_completeness, paired = TRUE, p.adjust.method = "bonferroni")
pwc
#In this csv file RAxMLEPA and EPAng results are included. As I need RAxMLEPA results first I subset another dataframe from Correct_position dataset
Stratified_Correct_position <- Correct_position[1:50,]
names(Stratified_Correct_position)
#As we need an id to perform ezANOVA I'm adding an id column to this dataframe. In RF random sample dataset (longdata) we have id so I combine longdata to Random_Correct_position

#id <- c(rep(1:10,4))
Stratified_Correct_position <- cbind(longdata, Stratified_Correct_position)

ezANOVA(data = Stratified_Correct_position, dv = Correct_percentage, wid = id, within = sample_completeness)

#check assumptions
#Outliers can be easily identified using box plot methods, implemented in the R function identify_outliers() [rstatix package]
Stratified_Correct_position %>%
  group_by(sample_completeness) %>%
  identify_outliers(Correct_percentage)
#There were no extreme outliers

#visualization
#Create a violin plot
violin_plot <- ggplot(Stratified_Correct_position, aes(x=sample_completeness, y=Correct_percentage,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of Correct position percentage for sampling completeness of stratified sampling")
violin_plot

#Normality assumption
#The normality assumption can be checked by computing Shapiro-Wilk test for each time point. If the data is normally distributed, the p-value should be greater than 0.05.
Stratified_Correct_position %>%
  group_by(sample_completeness) %>%  
  shapiro_test(Correct_percentage)
#The stratified correct position percentage was normally distributed at each sample completeness, as assessed by Shapiro-Wilk’s test (p > 0.05).

#QQ plot draws the correlation between a given data and the normal distribution. Create QQ plots for each level backbone tree(sample completeness)
ggqqplot(Stratified_Correct_position, "Correct_percentage", facet.by = "sample_completeness",color = "sample_completeness")

#From the resultant plots, as all the points fall approximately along the reference line, we can assume normality

#The assumption of sphericity will be automatically checked during the computation of the ANOVA test using the R function anova_test() [rstatix package]. The Mauchly’s test is internally used to assess the sphericity assumption
#By using the function get_anova_table() [rstatix] to extract the ANOVA table, the Greenhouse-Geisser sphericity correction is automatically applied to factors violating the sphericity assumption
# res.aov <- anova_test(data = longdata, dv = RF_distance, wid = id,within = sample_completeness)
# get_anova_table(res.aov)

#Post-hoc tests
#perform multiple pairwise paired t-tests between the levels of the within-subjects factor (Sample_completeness). P-values are adjusted using the Bonferroni multiple testing correction method.

# pairwise comparisons
pwc <- Stratified_Correct_position %>%
  pairwise_t_test(
    Correct_percentage ~ sample_completeness, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#All the pairwise differences are statistically significant

#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
Stratified_Correct_position_new <- lapply(Stratified_Correct_position,as.numeric)

#Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
rmcorr(id, Correct_percentage, sample_completeness, Stratified_Correct_position_new )



#---------------------------------------
