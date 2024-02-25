#Coded by M. A. Thanuja M. Fernando

#Last updated on19,Jan 2024
#With 99% sample

library(ape)
library(foreach)
library(phangorn)
library(ggplot2)
files <- list.files(path="/home/thanu/Desktop/FishData/New_dataset/STRATIFIED", pattern="^EPA[0-9]{1,2}_[0-9]{1,2}S_new.nwk", full.names=TRUE, recursive=FALSE)
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
           1, 11, 21, 31, 41) # 20_1, 40_1, 60_1, 80_1

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

# metric <- function(L_List){
#   #reordering
#   L_List <- L_List[order]
#   #group the values and assign values to variable names
#   L_1 <- unlist(L_List[1:4])
#   L_2 <- unlist(L_List[5:8])
#   L_3 <- unlist(L_List[9:12])
#   L_4 <- unlist(L_List[13:16])
#   L_5 <- unlist(L_List[17:20])
#   L_6 <- unlist(L_List[21:24])
#   L_7 <- unlist(L_List[25:28])
#   L_8 <- unlist(L_List[29:32])
#   L_9 <- unlist(L_List[33:36])
#   L_10 <- unlist(L_List[37:40])
#   #generate a dataframe using backbone % and metric values
#   results <- data.frame(GRF,L_1,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10)
# }

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
                ,main="Boxplot of RF distance for stratified samples based on sample completeness",add= "point", width=0.8,outlier.shape = NA)
ggpar(bp,orientation ="reverse")
# Plot Robinson-Fould distances 
RF_Stratified_Plot <- ggplot(RF_metric , aes(GRF)) +
  scale_y_reverse() + xlab("The level of backbone (% species completeness)")+ylab("RF Distance")+labs(title="The plot of RF distance vs. the level of backbone")+
  geom_line(aes(y=L_1),
            colour="red") +
  geom_line(aes(y=L_2),
            colour="green") +
  geom_line(aes(y=L_3),
            colour="blue") +
  geom_line(aes(y=L_4),
            colour="yellow") +
  geom_line(aes(y=L_5),
            colour="orange") +
  geom_line(aes(y=L_6),
            colour="black")+
  geom_line(aes(y=L_7),
            colour="purple") +
  geom_line(aes(y=L_8),
            colour="pink") +
  geom_line(aes(y=L_9),
            colour="brown") +
  geom_line(aes(y=L_10),
            colour="violet")

#theme_bw() will get rid of the background #remove grid using theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#remove border lines using heme(panel.border = element_blank())
#remove background (remove backgroud colour and border lines, but does not remove grid lines) using theme(panel.background = element_blank())
#add axis line using heme(axis.line = element_line(colour = "black"))
RF_Stratified_Plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#oath distance
P_metric<- metric(P_List)

# Plot Robinson-Fould distances 
Path_Stratified_Plot <- ggplot(P_metric , aes(GRF)) +
  scale_y_reverse() + xlab("The level of backbone (% species completeness)")+ylab("Path Distance")+labs(title="The plot of Path distance vs. the level of backbone")+
  geom_line(aes(y=L_1),
            colour="red") +
  geom_line(aes(y=L_2),
            colour="green") +
  geom_line(aes(y=L_3),
            colour="blue") +
  geom_line(aes(y=L_4),
            colour="yellow") +
  geom_line(aes(y=L_5),
            colour="orange") +
  geom_line(aes(y=L_6),
            colour="black")+
  geom_line(aes(y=L_7),
            colour="purple") +
  geom_line(aes(y=L_8),
            colour="pink") +
  geom_line(aes(y=L_9),
            colour="brown") +
  geom_line(aes(y=L_10),
            colour="violet")


#theme_bw() will get rid of the background #remove grid using theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#remove border lines using heme(panel.border = element_blank())
#remove background (remove backgroud colour and border lines, but does not remove grid lines) using theme(panel.background = element_blank())
#add axis line using heme(axis.line = element_line(colour = "black"))
Path_Stratified_Plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
getwd()
