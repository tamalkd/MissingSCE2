library(ggplot2)
library(ggExtra)
library(gridExtra)
library(readxl)

Raw <- read_xlsx("Results.xlsx", sheet = "plots")

Full <- Raw[Raw$Method == "Complete",]

Full10 <- Full
Full30 <- Full
Full50 <- Full
Full10$Missprop <- 10
Full30$Missprop <- 30
Full50$Missprop <- 50

Full2 <- rbind(Full10, Full30, Full50)

All <- Raw[Raw$Method != "Complete",]
All <- rbind(Full2, All)

All$Legend <- paste(All$Method, All$Model, sep = " - ")

All$Missprop <- factor(All$Missprop)
# All$Method <- factor(All$Method, levels = c("Complete", "RM", "MI"))
# All$Design <- factor(All$Design, levels = c("ABAB", "RBD"))

# Mode filters
mode = "mnar"

if(mode == "uni")
{
  AllF <- All[All$Model %in% c("normal", "uniform", "AR1"),]
  AllF$Legend <- factor(AllF$Legend, levels = c(
    "Complete - normal", "RM - normal", "MI - normal", 
    "Complete - uniform", "RM - uniform", "MI - uniform", 
    "Complete - AR1", "RM - AR1", "MI - AR1"
  ))
} else
if(mode == "mcar")
{
  AllF <- All[All$Model %in% c("mvn.3", "mvn.6"),]
  AllF <- AllF[AllF$Misstype %in% c("mcar", "none"),]
  AllF$Legend <- factor(AllF$Legend, levels = c(
    "Complete - mvn.3", "RM - mvn.3", "MI - mvn.3", 
    "Complete - mvn.6", "RM - mvn.6", "MI - mvn.6"
  ))
} else
if(mode == "mar")
{
  AllF <- All[All$Model %in% c("mvn.3", "mvn.6"),]
  AllF <- AllF[AllF$Misstype %in% c("mar+", "mar-", "none"),]
  AllF$Legend <- factor(AllF$Legend, levels = c(
    "Complete - mvn.3", "RM - mvn.3", "MI - mvn.3", 
    "Complete - mvn.6", "RM - mvn.6", "MI - mvn.6"
  ))
} else
if(mode == "mnar")
{
  AllF <- All[All$Model %in% c("mvn.3", "mvn.6"),]
  AllF <- AllF[AllF$Misstype %in% c("mnar+", "mnar-", "none"),]
  AllF$Legend <- factor(AllF$Legend, levels = c(
    "Complete - mvn.3", "RM - mvn.3", "MI - mvn.3", 
    "Complete - mvn.6", "RM - mvn.6", "MI - mvn.6"
  ))
}

All0 <- AllF[AllF$ES == 0,]
All1 <- AllF[AllF$ES == 1,]
All2 <- AllF[AllF$ES == 2,]

Agg0 <- aggregate(All0$Power, by=list(Missprop=All0$Missprop, Legend=All0$Legend), FUN=mean)
Agg1 <- aggregate(All1$Power, by=list(Missprop=All1$Missprop, Legend=All1$Legend), FUN=mean)
Agg2 <- aggregate(All2$Power, by=list(Missprop=All2$Missprop, Legend=All2$Legend), FUN=mean)

theme_set(theme_gray(base_size = 20))

plot1 <- ggplot(data=Agg0, aes(x=Missprop, y=x, group=Legend, colour = Legend)) + 
  geom_point(aes(shape = Legend), size = 3) +
  scale_shape_manual(values=c(15:18,21:25)) +
  geom_line(aes(linetype = Legend), size = 1.1) +
  xlab("Missing data percentage") +
  ylab("Type I error rate (%)") +
  ylim(0, 100) +
  ggtitle("Effect size = 0") +
  theme(plot.title = element_text(hjust=0.5), legend.position="bottom", legend.direction = "vertical") +
  removeGrid()

plot2 <- ggplot(data=Agg1, aes(x=Missprop, y=x, group=Legend, colour = Legend)) + 
  geom_point(aes(shape = Legend), size = 3) +
  scale_shape_manual(values=c(15:18,21:25)) +
  geom_line(aes(linetype = Legend), size = 1.1) +
  xlab("Missing data percentage") +
  ylab("Estimated power (%)") +
  ylim(0, 100) +
  ggtitle("Effect size = 1") +
  theme(plot.title = element_text(hjust=0.5), legend.position="bottom", legend.direction = "vertical") +
  removeGrid()

plot3 <- ggplot(data=Agg2, aes(x=Missprop, y=x, group=Legend, colour = Legend)) + 
  geom_point(aes(shape = Legend), size = 3) +
  scale_shape_manual(values=c(15:18,21:25)) +
  geom_line(aes(linetype = Legend), size = 1.1) +
  xlab("Missing data percentage") +
  ylab("Estimated power (%)") +
  ylim(0, 100) +
  ggtitle("Effect size = 2") +
  theme(plot.title = element_text(hjust=0.5), legend.position="bottom", legend.direction = "vertical") +
  removeGrid()

grid.arrange(plot1, plot2, plot3, ncol=3)
