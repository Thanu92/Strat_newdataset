#Stat_Code for Stratified swTools( RAxMLepa and EPAng)
#Jan 21, 2024
#The two-way repeated measures ANOVA can be performed in order to determine whether there is a significant interaction between Sample type and software tool on the correct position percentage.
#In stratified folder
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
library(ggpubr)
#Read the stratified sw tools for 5 different sampling types csv file
Correct_position<- read.csv("Strati_SWtools_new_extened.csv")
head(Correct_position)

#Summary statistics
Correct_position %>%
  group_by(SW_Tool, Sample_Type) %>%
  get_summary_stats(Correct_percentage, type = "mean_sd")

#Create box plots of the Correct_percentage colored by SW_Tool groups:
bxp <- ggboxplot(
  Correct_position, x = "Sample_Type", y = "Correct_percentage",
  color = "SW_Tool", palette = "jco"
)
bxp

#Check assumptions
#Outliers

Correct_position %>%
  group_by(SW_Tool, Sample_Type) %>%
  identify_outliers(Correct_percentage)

#Normality assumprion
#Compute Shapiro-Wilk test for each combinations of factor levels:
  
Correct_position %>%
  group_by(SW_Tool, Sample_Type) %>%
  shapiro_test(Correct_percentage)

#The Correct position percentage was normally distributed at each sample type point (p > 0.05), except for EPAng at 40% and 60% and RAxMLEPA at 60%, as assessed by Shapiro-Wilk’s test.

#Create QQ plot for each cell of design:
  
ggqqplot(Correct_position, "Correct_percentage", ggtheme = theme_bw()) +
  facet_grid(Sample_Type ~ SW_Tool, labeller = "label_both")
#From the plot above, as all the points fall approximately along the reference line except in 60% EPAng and RAxMLEPA, can we assume normality???
#Computation
Correct_position$SW_Tool<- as.factor(Correct_position$SW_Tool)
Correct_position$id <- as.factor(Correct_position$id)
Correct_position$Sample_Type <- as.factor(Correct_position$Sample_Type)
class(Correct_position$Correct_percentage)
res.aov <- anova_test(
  data = Correct_position, dv = Correct_percentage , wid = id,
  between = c(SW_Tool, Sample_Type)
)
get_anova_table(res.aov)

#There is a statistically significant two-way interactions between SW_Tool and Sample_Type, F(4,90) = 4.731, p=2.0e-03,p < 0.05.

# Effect of SW tools at each sample completeness 
one.way <- Correct_position %>%
  group_by(Sample_Type) %>%
  anova_test(dv = Correct_percentage, wid = id, between = SW_Tool) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
df <- data.frame(one.way)
# Pairwise comparisons between SW_tools 
pwc <- Correct_position %>%
  group_by(Sample_Type) %>%
  pairwise_t_test(
    Correct_percentage ~ SW_Tool, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#Considering the Bonferroni adjusted p-value, it can be seen that the simple main effect of SW_tools was not significant at the sample completeness of 60% (p = 0.183) and 99% (p= 0.281). 
#It becomes significant at 20% (p = 0.016),40% (1.35e-11) and 80% (p = 0.029).

#Pairwise comparisons show that the mean Correct position percentage was significantly different between RAxMLEPA and EPAng sw tools at 20% 40%,60%,80% (p = 6.29e-9, 1.04e-6, 1.44e-10, 2.68e-9 respectively) but not at 99% (p = 0.151).


# Effect of sample type for each sw tool
one.way2 <- Correct_position %>%
  group_by(SW_Tool) %>%
  anova_test(dv = Correct_percentage, wid = id, between = Sample_Type) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2
# Pairwise comparisons between sample completeness
pwc2 <- Correct_position %>%
  group_by(SW_Tool) %>%
  pairwise_t_test(
    Correct_percentage ~ Sample_Type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2


#According to Bonferroni, the sample type was significant for both the sw tools for RF distance
#Pairwise comparisons show that all comparisons among sample completeness were statistically significant for both the sw tools.

#We could report the result as follow:

#A two-way repeated measures ANOVA was performed to evaluate the effect of different sw tools over sample completeness on correct position percentage.

#There was a statistically significant interaction between sw tool and sample completeness on correct position percentage, F(4, 90) = 4.731,p= 0.002 p < 0.05. 
#Therefore, the effect of sw tool variable was analyzed at each sample completeness. 
#P-values were adjusted using the Bonferroni multiple testing correction method. 
#The effect of sw tool was significant at  20% (p = 0.016),40% (1.35e-11) and 80% (p = 0.029) but not at the sample completeness of 60% (p = 0.183) and 99% (p= 0.281).

#Pairwise comparisons, using paired t-test, show that the mean correct position percentage was significantly different between RAxMLEPA and EPAng at 20% 40%,60%,80% (p = 0.028, 0.0000000604, 0.038, 0.025 respectively) and t3 (p = 0.00017) but not at 99% (p = 0.169)..

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Sample_Type")
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

#------------------------------

#The two-way repeated measures ANOVA can be performed in order to determine whether there is a significant interaction between Sample type and software tool on the RF distance.
#In stratified folder

#Read the stratified sw tools for 5 different sampling types csv file
RF<- read.csv("RF_Strati_SWtools_extended.csv")
head(RF)

#Summary statistics
RF %>%
  group_by(SW_Tool, Sample_Type) %>%
  get_summary_stats(RF_distance, type = "mean_sd")

#Create box plots of the Correct_percentage colored by SW_Tool groups:
bxp <- ggboxplot(
  RF, x = "Sample_Type", y = "RF_distance",
  color = "SW_Tool", palette = "jco"
)
bxp
bxp_rec <- ggpar(bxp,orientation ="reverse")
bxp_rec
#Check assumptions
#Outliers

RF %>%
  group_by(SW_Tool, Sample_Type) %>%
  identify_outliers(RF_distance)

#Normality assumprion
#Compute Shapiro-Wilk test for each combinations of factor levels:

RF %>%
  group_by(SW_Tool, Sample_Type) %>%
  shapiro_test(RF_distance)

#The RF distance was normally distributed at each sample type point (p > 0.05), except for  RAxMLEPA at 80%, as assessed by Shapiro-Wilk’s test.

#Create QQ plot for each cell of design:

ggqqplot(RF, "RF_distance", ggtheme = theme_bw()) +
  facet_grid(Sample_Type ~ SW_Tool, labeller = "label_both")
#From the plot above, as all the points fall approximately along the reference line except in 60% EPAng and RAxMLEPA, can we assume normality???
#Computation
RF$SW_Tool<- as.factor(RF$SW_Tool)
RF$id <- as.factor(RF$id)
RF$Sample_Type <- as.factor(RF$Sample_Type)

res.aov <- anova_test(
  data = RF, dv = RF_distance , wid = id,
  between = c(SW_Tool, Sample_Type)
)
get_anova_table(res.aov)

#There is a statistically significant two-way interactions between SW_Tool and Sample_Type, F(4,90) = 369.161, p=6.38e-55 ,p < 0.05.

# Effect of SW tools at each sample completeness 
one.way <- RF %>%
  group_by(Sample_Type) %>%
  anova_test(dv = RF_distance, wid = id, between = SW_Tool) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# Pairwise comparisons between SW_tools 
pwc <- RF %>%
  group_by(Sample_Type) %>%
  pairwise_t_test(
    RF_distance ~ SW_Tool, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#Considering the Bonferroni adjusted p-value, it can be seen that the simple main effect of SW_tools was not significant at the sample completeness of 99% (p = 0.197). 
#It becomes significant at 20%,40%,60%,80% (p = 2.51e-14,8.13e-8,6.25e-18,8.26e-16 respectively).

#Pairwise comparisons show that the mean RF distance was significantly different between RAxMLEPA and EPAng sw tools at 20% 40%,60%,80% (p = 6.29e-9, 1.04e-6, 1.44e-10, 2.68e-9 respectively) but not at 99% (p = 0.151).


# Effect of sample type for each sw tool
one.way2 <- RF %>%
  group_by(SW_Tool) %>%
  anova_test(dv = RF_distance, wid = id, between = Sample_Type) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2
# Pairwise comparisons between sample completeness
pwc2 <- RF %>%
  group_by(SW_Tool) %>%
  pairwise_t_test(
    RF_distance ~ Sample_Type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2
#According to Bonferroni, the sample type was significant for both the sw tools based on RF distance
#Pairwise comparisons show that all comparisons among RF distance were statistically significant for both the sw tools.

#We could report the result as follow:

#A two-way repeated measures ANOVA was performed to evaluate the effect of different sw tools over sample completeness on RF_distance.

#There was a statistically significant interaction between sw tool and sample completeness on RF distance, F(4, 90) = 369.161,p= 6.38e-55 p < 0.05. 
#Therefore, the effect of sw tool variable was analyzed at each sample completeness based on RF distance. 
#P-values were adjusted using the Bonferroni multiple testing correction method. 
#The effect of sw tool was significant at 20%,40%,60%,80% (p = 2.51e-14,8.13e-8,6.25e-18,8.26e-16 respectively)but not at the sample completeness of 99% (p = 0.197).

#Pairwise comparisons, using paired t-test, show that the mean correct position percentage was significantly different between RAxMLEPA and EPAng at 20% 40%,60%,80% (p = 6.29e-9, 1.04e-6, 1.44e-10, 2.68e-9 respectively) but not at 99% (p = 0.151).


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Sample_Type")
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )
#------------------------------

#The two-way repeated measures ANOVA can be performed in order to determine whether there is a significant interaction between Sample type and sample completeness on the correct position percentage based on RAxMLEPA sw tool.
#In stratified folder
library(ggpubr)
#Read the SAmple types (RAndom, stratified and biased) for 5 different sampling completeness based on RAxMLEPA csv file
Correct_position<- read.csv("RAxML_R_S_correct_extended.csv")
#Correct_position<- read.csv("RAxML_SAmple_type_extended.csv")
head(Correct_position)

#Summary statistics
Correct_position %>%
  group_by(Sample, sample_type) %>%
  get_summary_stats(correct_placement_percentage, type = "mean_sd")

#Create box plots of the Correct_percentage colored by SW_Tool groups:
bxp <- ggboxplot(
  Correct_position, x = "sample_type", y = "correct_placement_percentage",
  color = "Sample", palette = "jco"
)
bxp

#Check assumptions
#Outliers

Correct_position %>%
  group_by(Sample, sample_type) %>%
  identify_outliers(correct_placement_percentage)

#Normality assumprion
#Compute Shapiro-Wilk test for each combinations of factor levels:

Correct_position %>%
  group_by(Sample, sample_type) %>%
  shapiro_test(correct_placement_percentage)

#The Correct position percentage was normally distributed at each sample type point (p > 0.05), except for random sample at 40%,60% and 99% stratified samples at 60%, as assessed by Shapiro-Wilk’s test.

#Create QQ plot for each cell of design:

ggqqplot(Correct_position, "correct_placement_percentage", ggtheme = theme_bw()) +
  facet_grid(sample_type ~ Sample, labeller = "label_both")
#From the plot above, as all the points fall approximately along the reference line except in 60% EPAng and RAxMLEPA, can we assume normality???
#Computation
Correct_position$Sample<- as.factor(Correct_position$Sample)
Correct_position$id <- as.factor(Correct_position$id)
Correct_position$sample_type <- as.factor(Correct_position$sample_type)

res.aov <- anova_test(
  data = Correct_position, dv = correct_placement_percentage , wid = id,
  between = c(Sample,sample_type)
)
get_anova_table(res.aov)

#There is a statistically significant two-way interactions between SW_Tool and Sample_Type, F(4,90) = 8.311, p=9.43e-06,p < 0.05.

# Effect of SW tools at each sample completeness 
one.way <- Correct_position %>%
  group_by(sample_type) %>%
  anova_test(dv = correct_placement_percentage, wid = id, between = Sample) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# Pairwise comparisons between SW_tools 
pwc <- Correct_position %>%
  group_by(sample_type) %>%
  pairwise_t_test(
    correct_placement_percentage ~ Sample, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#Considering the Bonferroni adjusted p-value, it can be seen that the simple main effect of Sample type (Random, stratified) significant at all the sampling completeness, since, the p<0.05 (p = 1.82e-14, 4.22e- 8 , 1.08e- 8 , 6.84e-14, 1.38e- 4  respectively)

#Pairwise comparisons show that the mean Correct position percentage was significantly different between Random and stratified samples at all the sampling levels 20% 40%,60%,80%,99% (p = 0.0000000117, 0.0000176 , 0.0000126 , 0.00000000147, 0.0000135  respectively)


# Effect of sample completeness for each samples (RAndom and stratified)
one.way2 <- Correct_position %>%
  group_by(Sample) %>%
  anova_test(dv = correct_placement_percentage, wid = id, between = sample_type) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2
# Pairwise comparisons between sample completeness
pwc2 <- Correct_position %>%
  group_by(Sample) %>%
  pairwise_t_test(
    correct_placement_percentage ~ sample_type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2
#According to Bonferroni, the sample completeness was significant for both the sample types (RAndom, stratified) for correct position percentage
#Pairwise comparisons show that all comparisons among sample completeness were statistically significant for both the random and stratified samples.

#We could report the result as follow:

#A two-way repeated measures ANOVA was performed to evaluate the effect of different type of sample (random and stratified) over sample completeness on correct position percentage.

#There was a statistically significant interaction between sampling types (random and stratified) and sample completeness on correct position percentage,  F(4,90) = 8.311, p=9.43e-06,p < 0.05. 
#Therefore, the effect of sample type (random, stratified) variable was analyzed at each sample completeness. 
#P-values were adjusted using the Bonferroni multiple testing correction method. 
#The effect of sample types (random,stratified) was significant at  significant at all the sampling completeness, since, the p<0.05 (p = 1.82e-14, 4.22e- 8 , 1.08e- 8 , 6.84e-14, 1.38e- 4  respectively)

#Pairwise comparisons, using paired t-test, show that the mean correct position percentage was significantly different between random and stratified at all the sampling levels 20% 40%,60%,80%,99% (p = 0.0000000117, 0.0000176 , 0.0000126 , 0.00000000147, 0.0000135  respectively)

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "sample_type")
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

#-----------------------------
#The two-way repeated measures ANOVA can be performed in order to determine whether there is a significant interaction between Sample type and sample completeness on the RF DISTANCE based on RAxMLEPA sw tool.
#In stratified folder
library(ggpubr)
#Read the SAmple types (RAndom, stratified and biased) for 5 different sampling completeness based on RAxMLEPA csv file
RF<- read.csv("RAxML_R_S_RF_extended.csv")
head(RF)

#Summary statistics
RF %>%
  group_by(Sample, sample_type) %>%
  get_summary_stats(RF_dist, type = "mean_sd")

#Create box plots of the Correct_percentage colored by SW_Tool groups:
bxp <- ggboxplot(
  RF, x = "sample_type", y = "RF_dist",
  color = "Sample", palette = "jco"
)
bxp
ggpar(bxp,orientation ="reverse")
#Check assumptions
#Outliers

RF%>%
  group_by(Sample, sample_type) %>%
  identify_outliers(RF_dist)

#Normality assumprion
#Compute Shapiro-Wilk test for each combinations of factor levels:

RF %>%
  group_by(Sample, sample_type) %>%
  shapiro_test(RF_dist)

#The Correct position percentage was normally distributed at each sample type point (p > 0.05), except for Stratified sample at 80%, as assessed by Shapiro-Wilk’s test.

#Create QQ plot for each cell of design:

ggqqplot(RF, "RF_dist", ggtheme = theme_bw()) +
  facet_grid(sample_type ~ Sample, labeller = "label_both")
#From the plot above, as all the points fall approximately along the reference line except in 60% EPAng and RAxMLEPA, can we assume normality???
#Computation
RF$Sample<- as.factor(RF$Sample)
RF$id <- as.factor(RF$id)
RF$sample_type <- as.factor(RF$sample_type)

res.aov <- anova_test(
  data = RF, dv = RF_dist, wid = id,
  between = c(Sample,sample_type)
)
get_anova_table(res.aov)

#There is a statistically significant two-way interactions between SW_Tool and Sample_Type, F(4,90) = 98.651, p=4.72e-32 ,p < 0.05.

# Effect of Sample types (random,stratified) at each sample completeness 
one.way <- RF %>%
  group_by(sample_type) %>%
  anova_test(dv = RF_dist, wid = id, between = Sample) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# Pairwise comparisons between SW_tools 
pwc <- RF %>%
  group_by(sample_type) %>%
  pairwise_t_test(
    RF_dist ~ Sample, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#Considering the Bonferroni adjusted p-value, it can be seen that the simple main effect of Sample type (Random, stratified) significant at all the sampling completeness, since, the p<0.05 (p =3.26e-14, 1.57e-11, 5.34e-16 , 2.39e-14, 3.9 e-2  respectively)

#Pairwise comparisons show that the mean Correct position percentage was significantly different between Random and stratified samples at all the sampling levels 20% 40%,60%,80%,99% (p = 0.00000000946, 0.000000261, 0.00000000044 , 0.00000000738, 0.038 respectively)


# Effect of sample completeness for each samples (RAndom and stratified)
one.way2 <- RF %>%
  group_by(Sample) %>%
  anova_test(dv = RF_dist, wid = id, between = sample_type) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2
# Pairwise comparisons between sample completeness
pwc2 <- RF %>%
  group_by(Sample) %>%
  pairwise_t_test(
    RF_dist ~ sample_type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2
#According to Bonferroni, the sample completeness was significant for both the sample types (RAndom, stratified) for RF distance
#Pairwise comparisons show that all comparisons among sample completeness were statistically significant for both the random and stratified samples.

#We could report the result as follow:

#A two-way repeated measures ANOVA was performed to evaluate the effect of different type of sample (random and stratified) over sample completeness on RF distance.

#There was a statistically significant interaction between sampling types (random and stratified) and sample completeness on RF distance, F(4,90) = 98.651, p=4.72e-32 ,p < 0.05.
#Therefore, the effect of sample type (random, stratified) variable was analyzed at each sample completeness. 
#P-values were adjusted using the Bonferroni multiple testing correction method. 
#The effect of sample types (random,stratified) was significant at  all the sampling completeness, since, the p<0.05 (p =3.26e-14, 1.57e-11, 5.34e-16 , 2.39e-14, 3.9 e-2  respectively)

#Pairwise comparisons, using paired t-test, show that the mean correct position percentage was significantly different between random and stratified at all the sampling levels 20% 40%,60%,80%,99% (p = 0.00000000946, 0.000000261, 0.00000000044 , 0.00000000738, 0.038 respectively)

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "sample_type")
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )



#------------------------------
#In Random folder
#Random, correct position, sw_tools data

#Read the random sw tools for 5 different sampling types csv file
Correct_position<- read.csv("Random_correct_position_SW_tools_extended.csv")
head(Correct_position)

#Summary statistics
Correct_position %>%
  group_by(SW_Tool, Sample_Type) %>%
  get_summary_stats(Correct_percentage, type = "mean_sd")

#Create box plots of the Correct_percentage colored by SW_Tool groups:
bxp <- ggboxplot(
  Correct_position, x = "Sample_Type", y = "Correct_percentage",
  color = "SW_Tool", palette = "jco"
)
bxp

#Check assumptions
#Outliers

Correct_position %>%
  group_by(SW_Tool, Sample_Type) %>%
  identify_outliers(Correct_percentage)

#Normality assumprion
#Compute Shapiro-Wilk test for each combinations of factor levels:

Correct_position %>%
  group_by(SW_Tool, Sample_Type) %>%
  shapiro_test(Correct_percentage)

#The Correct position percentage was normally distributed at each sample type point (p > 0.05), except for EPAng at 99% and RAxMLEPA at 40%, 60%,99% as assessed by Shapiro-Wilk’s test.

#Create QQ plot for each cell of design:

ggqqplot(Correct_position, "Correct_percentage", ggtheme = theme_bw()) +
  facet_grid(Sample_Type ~ SW_Tool, labeller = "label_both")
#From the plot above, as all the points fall approximately along the reference line except at 40% and 99% in RAxMLEPA, can we assume normality???
#Computation
Correct_position$SW_Tool<- as.factor(Correct_position$SW_Tool)
Correct_position$id <- as.factor(Correct_position$id)
Correct_position$Sample_Type <- as.factor(Correct_position$Sample_Type)
class(Correct_position$Correct_percentage)
res.aov <- anova_test(
  data = Correct_position, dv = Correct_percentage , wid = id,
  between = c(SW_Tool, Sample_Type)
)
get_anova_table(res.aov)

#There is a statistically significant two-way interactions between SW_Tool and Sample completeness, F(4,90) = 7.500, p=2.92e-05,p < 0.05.

# Effect of SW tools at each sample completeness 
one.way <- Correct_position %>%
  group_by(Sample_Type) %>%
  anova_test(dv = Correct_percentage, wid = id, between = SW_Tool) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# Pairwise comparisons between SW_tools 
pwc <- Correct_position %>%
  group_by(Sample_Type) %>%
  pairwise_t_test(
    Correct_percentage ~ SW_Tool, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#Considering the Bonferroni adjusted p-value, it can be seen that the simple main effect of SW_tools for random samples was not significant at the sample completeness of 60% (p = 0.104). 
#It becomes significant at 20% (p = 2.29e-13),40% (0.042), 80% (p = 0.009) and 99% (p=0.053).

#Pairwise comparisons show that the mean Correct position percentage was significantly different between RAxMLEPA and EPAng sw tools at 20%,80%,99% (p = 0.00000000483, 0.005, 0.000403 respectively) but not at 40% (p = 0.07) and 60%(p=0.051).


# Effect of sample type for each sw tool
one.way2 <- Correct_position %>%
  group_by(SW_Tool) %>%
  anova_test(dv = Correct_percentage, wid = id, between = Sample_Type) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2
# Pairwise comparisons between sample completeness
pwc2 <- Correct_position %>%
  group_by(SW_Tool) %>%
  pairwise_t_test(
    Correct_percentage ~ Sample_Type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2
#According to Bonferroni, the sample type was significant for both the sw tools for correct position percentage
#Pairwise comparisons show that all comparisons among sample completeness were statistically significant for both the sw tools.

#We could report the result as follow:

#A two-way repeated measures ANOVA was performed to evaluate the effect of different sw tools over sample completeness on correct position percentage for Random sample.

#There was a statistically significant interaction between sw tool and sample completeness on correct position percentage, F(4,90) = 7.500, p=2.92e-05,p < 0.05. 
#Therefore, the effect of sw tool variable was analyzed at each sample completeness. 
#P-values were adjusted using the Bonferroni multiple testing correction method. 
#The effect of sw tool for random samples was not significant at the sample completeness of 60% (p = 0.104). It becomes significant at 20% (p = 2.29e-13),40% (0.042), 80% (p = 0.009) and 99% (p=0.053).

#Pairwise comparisons, using paired t-test, show that the mean correct position percentage was significantly different between RAxMLEPA and EPAng significantly different between RAxMLEPA and EPAng sw tools at 20%,80%,99% (p = 0.00000000483, 0.005, 0.000403 respectively) but not at 40% (p = 0.07) and 60%(p=0.051).



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Sample_Type")
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

#------------------------------

#The two-way repeated measures ANOVA can be performed in order to determine whether there is a significant interaction between Sample type and software tool on the RF distance for random sample.
#In RANDOM folder

#Read the stratified sw tools for 5 different sampling types csv file
RF<- read.csv("Random_RF_SW_Tools_extended.csv")
head(RF)

#Summary statistics
RF %>%
  group_by(SW_Tool, Sample_Type) %>%
  get_summary_stats(RF_dist, type = "mean_sd")

#Create box plots of the Correct_percentage colored by SW_Tool groups:
bxp <- ggboxplot(
  RF, x = "Sample_Type", y = "RF_dist",
  color = "SW_Tool", palette = "jco"
)
bxp
bxp_rev <- ggpar(bxp,orientation ="reverse")
bxp_rev
#Check assumptions
#Outliers

RF %>%
  group_by(SW_Tool, Sample_Type) %>%
  identify_outliers(RF_dist)

#Normality assumprion
#Compute Shapiro-Wilk test for each combinations of factor levels:

RF %>%
  group_by(SW_Tool, Sample_Type) %>%
  shapiro_test(RF_dist)

#The RF distance was normally distributed at each sample type point (p > 0.05), as assessed by Shapiro-Wilk’s test.

#Create QQ plot for each cell of design:

ggqqplot(RF, "RF_dist", ggtheme = theme_bw()) +
  facet_grid(Sample_Type ~ SW_Tool, labeller = "label_both")
#From the plot above, as all the points fall approximately along the reference line except in 60% EPAng and RAxMLEPA, can we assume normality???
#Computation
RF$SW_Tool<- as.factor(RF$SW_Tool)
RF$id <- as.factor(RF$id)
RF$Sample_Type <- as.factor(RF$Sample_Type)

res.aov <- anova_test(
  data = RF, dv = RF_dist , wid = id,
  between = c(SW_Tool, Sample_Type)
)
get_anova_table(res.aov)

#There is a statistically significant two-way interactions between SW_Tool and Sample_Type, F(4,90) = 144.442, p= 2.72e-38  ,p < 0.05.

# Effect of SW tools at each sample completeness 
one.way <- RF %>%
  group_by(Sample_Type) %>%
  anova_test(dv = RF_dist, wid = id, between = SW_Tool) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# Pairwise comparisons between SW_tools 
pwc <- RF %>%
  group_by(Sample_Type) %>%
  pairwise_t_test(
    RF_dist ~ SW_Tool, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#Considering the Bonferroni adjusted p-value, it can be seen that the simple main effect of SW_tools was not significant at the sample completeness of 99% (p = 0.646). 
#It becomes significant at 20%,40%,60%,80% (p = 2   e-16,5.71e-12,2.26e-11,1.61e- 5  respectively).

#Pairwise comparisons show that the mean RF distance was significantly different between RAxMLEPA and EPAng sw tools at 20% 40%,60%,80% (p = 8.33e-11, 2.20e- 9, 6.93e- 7, 6.15e- 6 respectively) but not at 99% (p = 0.064).


# Effect of sample type for each sw tool
one.way2 <- RF %>%
  group_by(SW_Tool) %>%
  anova_test(dv = RF_dist, wid = id, between = Sample_Type) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2
# Pairwise comparisons between sample completeness
pwc2 <- RF %>%
  group_by(SW_Tool) %>%
  pairwise_t_test(
    RF_dist ~ Sample_Type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2
#According to Bonferroni, the sample type was significant for both the sw tools based on RF distance in random sampling
#Pairwise comparisons show that all comparisons among RF distance were statistically significant for both the sw tools.

#We could report the result as follow:

#A two-way repeated measures ANOVA was performed to evaluate the effect of different sw tools over sample completeness on RF_distance for random sampling.

#There was a statistically significant interaction between sw tool and sample completeness on RF distance, F(4,90) = 144.442, p= 2.72e-38  ,p < 0.05.. 
#Therefore, the effect of sw tool variable was analyzed at each sample completeness based on RF distance. 
#P-values were adjusted using the Bonferroni multiple testing correction method. 
#The effect of sw tool was significant at 20%,40%,60%,80% (p = 2   e-16,5.71e-12,2.26e-11,1.61e- 5  respectively).but not at the sample completeness of  99% (p = 0.646).

#Pairwise comparisons, using paired t-test, show that the mean correct position percentage was significantly different between RAxMLEPA and EPAng at 20% 40%,60%,80% (p = 8.33e-11, 2.20e- 9, 6.93e- 7, 6.15e- 6 respectively) but not at 99% (p = 0.064).


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Sample_Type")
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )
bxp_rev + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

#------------------------------
#Updated on Feb5, 2024
#For Random_Genus_swTools comparison
#The two-way repeated measures ANOVA can be performed in order to determine whether there is a significant interaction between Sample type and software tool on the Genus_position for random sample.
#In RANDOM folder

#Read the stratified sw tools for 5 different sampling types csv file
Genus<- read.csv("R_Genus_SWTools_extended.csv")
head(Genus)

#Summary statistics
Genus %>%
  group_by(SW_Tool, Sample_type) %>%
  get_summary_stats(Genus_position, type = "mean_sd")

#Create box plots of the Correct_percentage colored by SW_Tool groups:
bxp <- ggboxplot(
  Genus, x = "Sample_type", y = "Genus_position",
  color = "SW_Tool", palette = "jco"
)
bxp

#Check assumptions
#Outliers

Genus %>%
  group_by(SW_Tool, Sample_type) %>%
  identify_outliers(Genus_position)

#Normality assumprion
#Compute Shapiro-Wilk test for each combinations of factor levels:

Genus %>%
  group_by(SW_Tool, Sample_type) %>%
  shapiro_test(Genus_position)

#The Genus_position was normally distributed at each sample type point (p > 0.05) , as assessed by Shapiro-Wilk’s test.

#Create QQ plot for each cell of design:

ggqqplot(Genus, "Genus_position", ggtheme = theme_bw()) +
  facet_grid(Sample_type ~ SW_Tool, labeller = "label_both")
#From the plot above, as all the points fall approximately along the reference line, we can assume normality
#Computation
Genus$SW_Tool<- as.factor(Genus$SW_Tool)
Genus$id <- as.factor(Genus$id)
Genus$Sample_type <- as.factor(Genus$Sample_type)

res.aov <- anova_test(
  data = Genus, dv = Genus_position , wid = id,
  between = c(SW_Tool, Sample_type)
)
get_anova_table(res.aov)

#There is a statistically significant two-way interactions between SW_Tool and Sample_Type, F(4,90) = 19.153, p=2.00e-11  ,p < 0.05.

# Effect of SW tools at each sample completeness 
one.way <- Genus %>%
  group_by(Sample_type) %>%
  anova_test(dv = Genus_position, wid = id, between = SW_Tool) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# Pairwise comparisons between SW_tools 
pwc <- Genus %>%
  group_by(Sample_type) %>%
  pairwise_t_test(
    Genus_position ~ SW_Tool, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#Considering the Bonferroni adjusted p-value, it can be seen that the simple main effect of SW_tools was not significant at the sample completeness of 99% (p = 0.14). 
#It becomes significant at 20%,40%,60%,80% (p =  1.05e-15, 4.54e-10, 2.2e-2, 3e-3 respectively).

#Pairwise comparisons show that the mean Genus_position was significantly different between RAxMLEPA and EPAng sw tools at 20% 40%,60%,80%,99% (p = 2.6 e-11, 1.83e-7, 1.64e-4, 3.61e-5, 4.3 e-2 respectively).


# Effect of sample type for each sw tool
one.way2 <- Genus %>%
  group_by(SW_Tool) %>%
  anova_test(dv = Genus_position, wid = id, between = Sample_type) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2
# Pairwise comparisons between sample completeness
pwc2 <- Genus %>%
  group_by(SW_Tool) %>%
  pairwise_t_test(
    Genus_position ~ Sample_type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2
#According to Bonferroni, the sample type was significant for both the sw tools based on Genus_position in random sampling
#Pairwise comparisons show that all comparisons among Genus_position were statistically significant for both the sw tools.

#We could report the result as follow:

#A two-way repeated measures ANOVA was performed to evaluate the effect of different sw tools over sample completeness on Genus_position for random sampling.

#There was a statistically significant interaction between sw tool and sample completeness on Genus_position, F(4,90) = 19.153, p=2.00e-11 ,p < 0.05.. 
#Therefore, the effect of sw tool variable was analyzed at each sample completeness based on FGenus_position. 
#P-values were adjusted using the Bonferroni multiple testing correction method. 
#The effect of sw tool was significant at 20%,40%,60%,80% (p =  1.05e-15, 4.54e-10, 2.2e-2, 3e-3 respectively). not significant at the sample completeness of 99% (p = 0.14)

#Pairwise comparisons show that the mean Genus_position was significantly different between RAxMLEPA and EPAng sw tools at 20% 40%,60%,80% (p = 2.6 e-11, 1.83e-7, 1.64e-4, 3.61e-5, 4.3 e-2 respectively).

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Sample_type")
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

#------------------------------------------------

#------------------------------
#For Stratified_Genus_swTools comparison
#The two-way repeated measures ANOVA can be performed in order to determine whether there is a significant interaction between Sample type and software tool on the Genus_position for Stratified sample.
#In STRATIFIED folder

#Read the stratified sw tools for 5 different sampling types csv file
Genus<- read.csv("S_Genus_SWTools_extended.csv")
head(Genus)

#Summary statistics
Genus %>%
  group_by(SW_Tool, Sample_Type) %>%
  get_summary_stats(Genus_position, type = "mean_sd")

#Create box plots of the Correct_percentage colored by SW_Tool groups:
bxp <- ggboxplot(
  Genus, x = "Sample_Type", y = "Genus_position",
  color = "SW_Tool", palette = "jco"
)
bxp

#Check assumptions
#Outliers

Genus %>%
  group_by(SW_Tool, Sample_Type) %>%
  identify_outliers(Genus_position)

#Normality assumprion
#Compute Shapiro-Wilk test for each combinations of factor levels:

Genus %>%
  group_by(SW_Tool, Sample_Type) %>%
  shapiro_test(Genus_position)

#The Genus_position was normally distributed at each sample type point (p > 0.05) , as assessed by Shapiro-Wilk’s test.

#Create QQ plot for each cell of design:

ggqqplot(Genus, "Genus_position", ggtheme = theme_bw()) +
  facet_grid(Sample_Type ~ SW_Tool, labeller = "label_both")
#From the plot above, as all the points fall approximately along the reference line, we can assume normality
#Computation
Genus$SW_Tool<- as.factor(Genus$SW_Tool)
Genus$id <- as.factor(Genus$id)
Genus$Sample_Type <- as.factor(Genus$Sample_Type)

res.aov <- anova_test(
  data = Genus, dv = Genus_position , wid = id,
  between = c(SW_Tool, Sample_Type)
)
get_anova_table(res.aov)

#There is a statistically significant two-way interactions between SW_Tool and Sample_Type, F(4,90) = 83.954, 1.54e-29  ,p < 0.05.

# Effect of SW tools at each sample completeness 
one.way <- Genus %>%
  group_by(Sample_Type) %>%
  anova_test(dv = Genus_position, wid = id, between = SW_Tool) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way
# Pairwise comparisons between SW_tools 
pwc <- Genus %>%
  group_by(Sample_Type) %>%
  pairwise_t_test(
    Genus_position ~ SW_Tool, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

#Considering the Bonferroni adjusted p-value, it can be seen that the simple main effect of SW_tools was significant at 20%,40%,60%,80%,99% (p =  4.64e-8, 4.43e-7, 2.65e-10, 1.28e-16,3.71e-9 respectively).

#Pairwise comparisons show that the mean Genus_position was significantly different between RAxMLEPA and EPAng sw tools at 20% 40%,60%,80%,99% (p = 6.99e- 5, 3.93e- 5, 7.5 e- 7, 1.54e-11, 1.23e- 6 respectively).


# Effect of sample type for each sw tool
one.way2 <- Genus %>%
  group_by(SW_Tool) %>%
  anova_test(dv = Genus_position, wid = id, between = Sample_Type) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2
# Pairwise comparisons between sample completeness
pwc2 <- Genus %>%
  group_by(SW_Tool) %>%
  pairwise_t_test(
    Genus_position ~ Sample_Type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2
#According to Bonferroni, the sample type was significant for both the sw tools based on Genus_position in random sampling
#Pairwise comparisons show that all comparisons among Genus_position were statistically significant for both the sw tools except EPAng 0.8S - 0.99S (p=1.25e- 1)

#We could report the result as follow:

#A two-way repeated measures ANOVA was performed to evaluate the effect of different sw tools over sample completeness on Genus_position for random sampling.

#There was a statistically significant interaction between sw tool and sample completeness on Genus_position, F(4,90) = 83.954, 1.54e-29  ,p < 0.05. 
#Therefore, the effect of sw tool variable was analyzed at each sample completeness based on FGenus_position. 
#P-values were adjusted using the Bonferroni multiple testing correction method. 
#The effect of sw tool was significant at 20%,40%,60%,80%,99% (p =  4.64e-8, 4.43e-7, 2.65e-10, 1.28e-16,3.71e-9 respectively).

#Pairwise comparisons show that the mean Genus_position was significantly different between RAxMLEPA and EPAng sw tools at 20% 40%,60%,80%,99% (p = 6.99e- 5, 3.93e- 5, 7.5 e- 7, 1.54e-11, 1.23e- 6 respectively).

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Sample_Type")
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

#------------------------------After this line old code---------

#Load packages
#install.packages("foreach")


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

Correct_position<- read.csv("Strati_SWtools_new.csv")

dim(Correct_position)
# #In this csv file RAxMLEPA and EPAng results are included. As I need EPAngresults first I subset another dataframe from Correct_position dataset
# EPAng_Correct_position <- Correct_position[51:100,]
head(Correct_position,11)
names(Correct_position)
#As we need an id to perform ezANOVA I'm adding an id column to this dataframe. In RF random sample dataset (longdata) we have id so I combine longdata to Random_Correct_position

id <- c(rep(1:10,10))
SW_Correct_position<- cbind(id, Correct_position)
# class(EPAng_Correct_position)
# names(EPAng_Correct_position)
# head(EPAng_Correct_position)
ezANOVA(data =SW_Correct_position, dv = Correct_percentage, wid = id, within = Sample_Type)

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
