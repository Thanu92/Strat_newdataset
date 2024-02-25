#updated 27, Jan 2024
Anova <- read.csv("2way_Repeated_Measures_ANOVA.csv")
class(Anova)
df1 <- data.frame(Anova)
names(df1)
df1_new <- df1[,c("group1","group2","p")]
df1_melt <- melt(df1_new)
plot <- ggplot(df1_melt, aes(group1, group2))+ geom_tile(aes(fill=value))
plot+geom_text(aes(label = value), fontface = "bold",size = 5,color="white")+theme(panel.background = element_rect(fill='white', colour='black'))                                                

# plot <- ggplot(df1_melt, aes(group1, group2,fill=value) )+ geom_tile(color="black")+scale_fill_gradient(low = "yellow", high = "red")
# plot+geom_text(aes(label = value), fontface = "bold",size = 5)+theme(panel.background = element_rect(fill='white', colour='black'))                                                
# Anova <- read.csv("2way_Repeated_Measures_ANOVA.csv")
# class(Anova)
# df1 <- data.frame(Anova)
# names(df1)
# df1_new <- df1[c(1,3,5),c("group1","group2","p")]
# df1_melt <- melt(df1_new)
# plot <- ggplot(df1_melt, aes(group1, group2)) + geom_tile(aes(fill=value))
# plot+theme(panel.background = element_rect(fill='white', colour='black'))  

PWC_CSp <- read.csv("pairwise_comparison_CorrectSp.csv")
class(PWC_CSp)
df1 <- data.frame(PWC_CSp)
names(df1)
#df1_new <- df1[,c("group1","group2","p")]
df1_new <- df1
df1_melt <- melt(df1_new)
plot <- ggplot(df1_melt, aes(group1, group2))+ geom_tile(aes(fill=value))
plot+geom_text(aes(label = value), fontface = "bold",size = 2,color="white")+theme(panel.background = element_rect(fill='white', colour='black'))                                                
