"Perhaps of I changed it"
# Functions

stat_box_data <- function(y, upper_limit = max(un.gr.gen$Range) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'mean =', round(mean(y), 2), '\n',
                    'median =', round(median(y), 2), '\n')
    )
  )
}



# Specialists #
for (BBB in 2:24) {
  take.a.look <- read.csv(file = paste("./Mapping.single.sp/Mantel.0.05/Maps/",BBB,"otu.range.habitat.csv"), header = TRUE)
  take.a.look <- take.a.look[complete.cases(take.a.look), ] # remove the NA's
  take.a.look <- take.a.look[take.a.look$Range>1,] # remove garbage variograms (you can also take a look at the maps in the corresponding directory)
  #take.a.look <- take.a.look[take.a.look$Range<30,] # some ranges are really big
  print(c(BBB,wilcox.test(Range ~ Habitat, data = take.a.look, exact = FALSE)$p.value))
}

BBB=2 #min occur
take.a.look <- read.csv(file = paste("./Mapping.single.sp/Mantel.0.05/Maps/",BBB,"otu.range.habitat.csv"), header = TRUE)
take.a.look <- take.a.look[complete.cases(take.a.look), ] # remove the NA's
take.a.look <- take.a.look[take.a.look$Range>1,] # remove garbage variograms (you can also take a look at the maps in the corresponding directory)
take.a.look <- take.a.look[take.a.look$Range<10,] # there are three OTU with large ranges - This is what they look like -> take.a.look[take.a.look$Range>30,]
spec <- take.a.look
# this is the table you will want to use for publication 
#take.a.look # OTU / Range / Habitat / Max.abun / Species
take.a.look %>%
  group_by(Habitat) %>%
  summarise(
    median = median(Range),
    mean =mean(Range),
    sd = sd(Range),
    n = n(),
    min.range = min(Range),
    max.range = max(Range)
  )


# Generalists #
for (BBB in 2:30) { # which occupancy has significant p-values 
  take.a.look <- read.csv(file = paste("./Mapping.single.sp/Mantel.0.05/Maps/",BBB,"otu.range.habitat.generalist.csv"), header = TRUE)
  take.a.look <- take.a.look[take.a.look$Range>0.5,] # remove garbage variograms (you can also take a look at the maps in the corresponding directory)
  take.a.look <- take.a.look[take.a.look$Range<30,]
  print(c(BBB,wilcox.test(Range ~ Habitat, data = take.a.look, exact = FALSE)$p.value))
}


BBB=2 #min occur
take.a.look <- read.csv(file = paste("./Mapping.single.sp/Mantel.0.05/Maps/",BBB,"otu.range.habitat.generalist.csv"), header = TRUE)
take.a.look <- take.a.look[complete.cases(take.a.look), ] # remove the NA's
take.a.look <- take.a.look[take.a.look$Range>1,] # remove garbage variograms (you can also take a look at the maps in the corresponding directory)
take.a.look <- take.a.look[take.a.look$Range<10,]# there are three OTU with large ranges - This is what they look like -> take.a.look[take.a.look$Range>30,]
genr <- take.a.look
#take.a.look # OTU / Range / Habitat / Max.abun / Species
take.a.look %>%
  group_by(Habitat) %>%
  summarise(
    median = median(Range),
    mean =mean(Range),
    sd = sd(Range),
    n = n(),
    min.range = min(Range),
    max.range = max(Range)
  )
wilcox.test(Range ~ Habitat, data = take.a.look, exact = FALSE)$p.value




# now you have two datasets that are either specialists or generalists
spec$SG <- "Specialist"
genr$SG <- "Generalist"
# add a id tag and combine 
un.gr.gen <- rbind(spec,genr)


un.gr.gen %>%
  group_by(SG) %>%
  summarise(
    median = median(Range),
    mean =mean(Range),
    sd = sd(Range),
    n = n(),
    min.range = min(Range),
    max.range = max(Range)
  )



# check for normality 
shapiro.test(spec$Range)
shapiro.test(genr$Range)
# For the both datasets p < 0.05 suggesting strong evidence of non-normality and a nonparametric test should be used


wilcox.test(Range ~ Habitat, data = spec, exact = FALSE)$p.value
wilcox.test(Range ~ Habitat, data = genr, exact = FALSE)$p.value
wilcox.test(Range ~ Habitat, data = un.gr.gen, exact = FALSE)$p.value

un.gr.gen$group <- paste(un.gr.gen$Habitat,un.gr.gen$SG)
pairwise.wilcox.test(un.gr.gen$Range,un.gr.gen$group,p.adjust.method = "BH")



# set the order you'd like
un.gr.gen$group <- factor(un.gr.gen$group , levels=c("Understory Specialist", "Ground Specialist", "Understory Generalist", "Ground Generalist"))
# annotation
anno <- data.frame(x1 = 1, x2 = 2, y1 = 8, y2 = 8.5, xstar = 1.5, ystar = 8.7, lab = "***", SG = "Specialist")
# change labels to more appropriate names 
levels(un.gr.gen$Habitat)[match("Understory",levels(un.gr.gen$Habitat))] <- "Phyllosphere"
levels(un.gr.gen$Habitat)[match("Ground",levels(un.gr.gen$Habitat))] <- "Soil"



ggplot(un.gr.gen, aes(x = Habitat, y = Range, fill=Habitat)) +
  geom_violin() + 
  stat_summary(fun.data = stat_box_data, 
               geom = "text", 
               hjust = 0.5,
               vjust = 0.9) +
  stat_summary(fun.y = mean, 
               fun.ymin = mean, 
               fun.ymax = mean, 
               geom = "crossbar", 
               width = 0.3)+
  geom_jitter(shape=16, size=.3, position=position_jitter(0.3))+
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "OTU patch radius size (m)") +
  scale_fill_grey(start = 0.4, end = 0.9) +
  facet_grid(. ~ SG) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab,fill=NULL),family = "Times") +
  geom_segment(data = anno, aes(x = x1, xend = x1, 
                                y = y1, yend = y2,fill=NULL),
               colour = "black") +
  geom_segment(data = anno, aes(x = x2, xend = x2, 
                                y = y1, yend = y2,fill=NULL),
               colour = "black") +
  geom_segment(data = anno, aes(x = x1, xend = x2, 
                                y = y2, yend = y2,fill=NULL),
               colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.text.x = element_text(size=12, face="bold"))

ggsave("OTU range size.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 14, height = 8, units = "in", 
       dpi = 150, limitsize = TRUE)



