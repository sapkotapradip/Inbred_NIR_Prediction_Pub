data = read.csv("unique.hybrids.csv")
library(ggplot2)


a <- ggplot(data, aes(female, male)) +
  geom_point(size=3) +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = -90, vjust = 0.25, hjust=0.25)) +
  #scale_x_discrete(labels = setNames(female_labs$label, female_labs$female)) +
  #scale_y_discrete(labels = setNames(male_labs$label, male_labs$male))+
  labs(y = "Seed Parents", x= "Pollinator Parents")

a

jpeg("crosses.jpeg",width = 12,height =7,units = "in", res=600)
a
dev.off()

