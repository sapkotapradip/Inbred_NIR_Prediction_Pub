# Plot results
getwd()
setwd("output/results_8_29/")
###lets plot
library(plyr)
library(readr)
library(dplyr)
library(stringr)
getwd()
setwd("output/results_8_29/")


list_csv_files <- list.files("../results_8_29/")
df2 <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
df2
write.csv(df2, "summary.csv")
df2 = read.csv("summary.csv")

df1 <- as.data.frame(df2 %>%  dplyr::group_by(traits,models,CV,grain) %>% 
                       dplyr::summarise(M = mean(cor, na.rm=TRUE),
                                        SD = sd(cor, na.rm=TRUE)))

library(tidyr)
library(ggplot2)

df1_CV1 = df1[df1$CV == "CV1",]
df1_CV2 = df1[df1$CV == "CV2",]

df1_CV1 =
  df1_CV1  %>%
  mutate(
    traits = case_when(
      traits == "ash" ~ "Ash",
      traits == "da" ~ "Days to Anthesis",
      traits == "fat" ~ "Fat",
      traits == "fiber" ~ "Fiber",
      traits == "ph" ~ "Plant Height",
      traits == "protein" ~ "Protein",
      traits == "starch" ~ "Starch",
      traits == "yield" ~ "Grain Yield",
      TRUE~traits)
  )

df1_CV1$traits <- factor(df1_CV1$traits, levels =  c("Grain Yield", 
                                                     "Plant Height", 
                                                     "Days to Anthesis", 
                                                     "Starch", 
                                                     "Protein", 
                                                     "Fat", 
                                                     "Fiber", 
                                                     "Ash"))


p = ggplot(df1_CV1, aes(models, y=M)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(y = M + SD/2 + 0.05,   # push above error bar
                label = round(M, 2)), 
            vjust = 0,                 # align bottom of text to point
            color = "black",
            position = position_dodge(0.9),
            angle = 0, size = 3) +
  geom_errorbar(aes(ymin=M, ymax=M+(SD/2)), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  facet_grid(~traits~grain)+
  scale_y_continuous("Prediction Ability", expand = expansion(mult = c(0, 0.1)))+
  xlab("Prediction Models") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=5, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=20, face="bold", colour = "black"),    
    axis.title.y = element_text(size=20, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=10, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 15, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 9, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position = "none"
    #axis.text.x.bottom = element_blank()
  )

p

jpeg("Pred.ability_CV1_monochrom.jpeg",width = 12,height =10,units = "in", res=600)
p
dev.off()


# CV2

df1_CV2 =
  df1_CV2  %>%
  mutate(
    traits = case_when(
      traits == "ash" ~ "Ash",
      traits == "da" ~ "Days to Anthesis",
      traits == "fat" ~ "Fat",
      traits == "fiber" ~ "Fiber",
      traits == "ph" ~ "Plant Height",
      traits == "protein" ~ "Protein",
      traits == "starch" ~ "Starch",
      traits == "yield" ~ "Grain Yield",
      TRUE~traits)
  )

df1_CV2$traits <- factor(df1_CV1$traits, levels =  c("Grain Yield", 
                                                     "Plant Height", 
                                                     "Days to Anthesis", 
                                                     "Starch", 
                                                     "Protein", 
                                                     "Fat", 
                                                     "Fiber", 
                                                     "Ash"))


p = ggplot(df1_CV2, aes(models, y=M)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(y = M + SD/2 + 0.05,   # push above error bar
                label = round(M, 2)), 
            vjust = 0,                 # align bottom of text to point
            color = "black",
            position = position_dodge(0.9),
            angle = 0, size = 3) +
  geom_errorbar(aes(ymin=M, ymax=M+(SD/2)), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  facet_grid(~traits~grain)+
  scale_y_continuous("Prediction Ability", expand = expansion(mult = c(0, 0.1)))+
  xlab("Prediction Models") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=5, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=20, face="bold", colour = "black"),    
    axis.title.y = element_text(size=20, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=10, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 15, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 9, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position = "none"
    #axis.text.x.bottom = element_blank()
  )

p

jpeg("Pred.ability_CV2_monochrom.jpeg",width = 12,height =10,units = "in", res=600)
p
dev.off()


df2 = read.csv("summary.csv")

df2 = subset(df2, CV == "CV1")
df2 = subset(df2, grain = "20CS")
# significance letters
# Create a placeholder to store all letters

library(emmeans)

sig_letters <- data.frame()

for (grain_s in unique(df2$grain)) {
  for (trait_level in unique(df2$traits))  {
    
    subset_data <- df2 %>% filter(grain == grain_s, traits == trait_level)
    
    # Fit linear model
    mod <- lm(cor ~ models, data = subset_data)
    
    # Tukey test
    tukey <- emmeans(mod, pairwise ~ models, adjust = "tukey")
    
    # Get compact letter display
    cld_df <- multcomp::cld(tukey$emmeans, Letters = letters, adjust = "tukey")
    
    cld_df$CV <- cv_level
    cld_df$trait <- trait_level
    
    sig_letters <- bind_rows(sig_letters, cld_df)
  }
}

df = sig_letters %>% arrange(models)


df = as.data.frame(sig_letters)

df$CV <- factor(df$CV, levels =  c("CV1", "CV2"))
df$models = factor(df$models, levels = c("M1", "M2", "M3", "M4", "M5", "M6"))
df$trait = factor(df$trait, levels = c("yield", "da", "ph", "starch", "protein", "fat", "fiber", "ash"))
df$.group = trimws(df$.group)
colnames(df)[1] = "Models"

df1_CV1 = df[df$CV == "CV1",]
df1_CV2 = df[df$CV == "CV2",]

p = ggplot(df, aes(Models, y=emmean, fill=Models)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=round(emmean,2), y = emmean + 0.3), color="black",
            position = position_dodge(0.9), angle = 90,size=3, vjust = 0.4)+
  geom_text(aes(label = .group, y = emmean + SE + 0.03,
                group = Models),
            position = position_dodge(width = 0.9), 
            size = 3, 
            fontface = "bold") +
  geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width=.2,
                position=position_dodge(.9))+
  
  theme_bw()+
  facet_grid(~trait~grain)+
  scale_y_continuous("Prediction ability")+
  xlab("Prediction Models") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=5, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=17, face="bold", colour = "black"),    
    axis.title.y = element_text(size=17, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 13, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 8, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position = "none"
    #axis.text.x.bottom = element_blank()
  )

p
      