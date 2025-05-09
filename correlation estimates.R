# calculation of heterosis

pheno_yield_CS = read.csv("pheno_yield_20CSMP.csv")
pheno_da_CS = read.csv("pheno_da_20CSMP.csv")
pheno_ph_CS= read.csv("pheno_ph_20CSMP.csv")
pheno_yield_MA = read.csv("pheno_yield_20MAMP.csv")
pheno_da_MA = read.csv("pheno_da_20MAMP.csv")
pheno_ph_MA = read.csv("pheno_ph_20MAMP.csv")

cor.test(pheno_yield_CS[pheno_yield_CS$env == "20CS",]$gy_mid_parent, 
    pheno_yield_CS[pheno_yield_CS$env == "20CS",]$blue) 

cor.test(pheno_da_CS[pheno_da_CS$env == "20CS",]$dy_mid_parent, 
    pheno_da_CS[pheno_da_CS$env == "20CS",]$blue) 

cor.test(pheno_ph_CS[pheno_ph_CS$env == "20CS",]$ph_mid_parent, 
         pheno_ph_CS[pheno_ph_CS$env == "20CS",]$blue) 

cor.test(pheno_yield_MA[pheno_yield_MA$env == "20LY",]$gy_mid_parent, 
         pheno_yield_MA[pheno_yield_MA$env == "20LY",]$blue) 

cor.test(pheno_da_MA[pheno_da_MA$env == "20LY",]$dy_mid_parent, 
         pheno_da_MA[pheno_da_MA$env == "20LY",]$blue) 

cor.test(pheno_ph_MA[pheno_ph_MA$env == "20LY",]$ph_mid_parent, 
         pheno_ph_MA[pheno_ph_MA$env == "20LY",]$blue) 


LY_ph = as.data.frame(pheno_ph_MA %>% group_by(env) %>% dplyr::summarize(cor = cor(blue, ph_mid_parent , use = "complete.obs"),
                                                                         p_value = cor.test(blue, ph_mid_parent, use = "complete.obs")$p.value))

CS_ph = as.data.frame(pheno_ph_CS %>% group_by(env) %>% dplyr::summarize(cor = cor(blue, ph_mid_parent , use = "complete.obs"),
                                                                         p_value = cor.test(blue, ph_mid_parent, use = "complete.obs")$p.value))

LY_da = as.data.frame(pheno_da_MA %>% group_by(env) %>% dplyr::summarize(cor = cor(blue, dy_mid_parent , use = "complete.obs"),
                                                                         p_value = cor.test(blue, dy_mid_parent, use = "complete.obs")$p.value))
CS_da = as.data.frame(pheno_da_CS %>% group_by(env) %>% dplyr::summarize(cor = cor(blue, dy_mid_parent , use = "complete.obs"),
                                                                         p_value = cor.test(blue, dy_mid_parent, use = "complete.obs")$p.value))

LY_gy = as.data.frame(pheno_yield_MA %>% group_by(env) %>% dplyr::summarize(cor = cor(blue, gy_mid_parent , use = "complete.obs"),
                                                                            p_value = cor.test(blue, gy_mid_parent, use = "complete.obs")$p.value))

CS_gy = as.data.frame(pheno_yield_CS %>% group_by(env) %>% dplyr::summarize(cor = cor(blue, gy_mid_parent , use = "complete.obs"),
                                                                            p_value = cor.test(blue, gy_mid_parent, use = "complete.obs")$p.value))


######### estimate correlation between hybrid and inbred NIR bands ############

#### reading pheno data for hybrids
pheno_19CS = read.csv("blues_19CS_hybrids_spectra.csv", check.names = FALSE)
pheno_19TA = read.csv("blues_19TA_hybrids_spectra.csv", check.names = FALSE)
pheno_20CS = read.csv("blues_20CS_hybrids_spectra.csv", check.names = FALSE)
pheno_20MA = read.csv("blues_20MA_hybrids_spectra.csv", check.names = FALSE)

pheno_combined = rbind(pheno_19CS,pheno_19TA, pheno_20CS, pheno_20MA)
ne <- as.vector(table(pheno_combined$env)) ## counting the number of observations
ne

#### reading pheno data for inbreds
parents_20CS = read.csv("blues_20CS_inbreds_spectra.csv", check.names = FALSE)
parents_20MA = read.csv("blues_20MA_inbreds_spectra.csv", check.names = FALSE)

midparent_20CS = read.csv("mid_parent_NIR_20CS_all.csv", check.names = FALSE)
midparent_20MA = read.csv("mid_parent_NIR_20MA_all.csv", check.names = FALSE)

rownames(midparent_20CS) = midparent_20CS[,1]
rownames(midparent_20MA) = midparent_20MA[,1]

pedigrees_list = pheno_combined$pedigree
filtered_df_20CS <- midparent_20CS[rownames(midparent_20CS) %in% pedigrees_list, ]
final_pedigree_20CS <- filtered_df_20CS[match(pedigrees_list, rownames(filtered_df_20CS)), ]

filtered_df_20MA <- midparent_20MA[rownames(midparent_20MA) %in% pedigrees_list, ]
final_pedigree_20MA <- filtered_df_20MA[match(pedigrees_list, rownames(filtered_df_20MA)), ]

mid_parent_19CS_20CS = final_pedigree_20CS[1:364,]
mid_parent_19TA_20CS = final_pedigree_20CS[365:543,]
mid_parent_20CS_20CS = final_pedigree_20CS[544:907,]
mid_parent_20MA_20CS = final_pedigree_20CS[-c(1:907),]

mid_parent_19CS_20MA = final_pedigree_20MA[1:364,]
mid_parent_19TA_20MA = final_pedigree_20MA[365:543,]
mid_parent_20CS_20MA = final_pedigree_20MA[544:907,]
mid_parent_20MA_20MA = final_pedigree_20MA[-c(1:907),]


rownames(pheno_20MA) = pheno_20MA$pedigree
pheno_20MA = pheno_20MA[,-c(1:7)]
dim(pheno_20MA); dim(mid_parent_20MA_20MA)

rownames(mid_parent_20MA_20MA) = mid_parent_20MA_20MA[,1]
mid_parent_20MA_20MA = mid_parent_20MA_20MA[,-1]
dim(pheno_20MA); dim(mid_parent_20MA_20MA)

cor = list()
p_value = list()

for (i in 1:4200) {
  cor[i] <- cor(pheno_20MA[,i], mid_parent_20MA_20MA[,i], use = "complete.obs")
  p_value[i] <- cor.test(pheno_20MA[,i], mid_parent_20MA_20MA[,i], use = "complete.obs")$p.value
}



cor = as.data.frame(cor, check.names = FALSE)
a = t(cor)
p_value = as.data.frame(p_value, check.names = FALSE)
p_value = t(p_value)

res = cbind(a, p_value)
rownames(res) = NULL

range(res[,2])

#########################
rownames(pheno_20CS) = pheno_20CS$pedigree
pheno_20CS = pheno_20CS[,-c(1:7)]
dim(pheno_20CS); dim(mid_parent_20CS_20CS)

rownames(mid_parent_20CS_20CS) = mid_parent_20CS_20CS[,1]
mid_parent_20CS_20CS = mid_parent_20CS_20CS[,-1]
dim(pheno_20CS); dim(mid_parent_20CS_20CS)

cor = list()
p_value = list()

for (i in 1:4200) {
  cor[i] <- cor(pheno_20CS[,i], mid_parent_20CS_20CS[,i], use = "complete.obs")
  p_value[i] <- cor.test(pheno_20CS[,i], mid_parent_20CS_20CS[,i], use = "complete.obs")$p.value
}



cor = as.data.frame(cor, check.names = FALSE)
a = t(cor)
p_value = as.data.frame(p_value, check.names = FALSE)
p_value = t(p_value)

res = cbind(a, p_value)
rownames(res) = NULL

range(res[,2])
#################################
# do correlation with NIR bands between 20CS and 20LY with other environments as well

# 19CS and 20MA
mid_parent_19CS_20MA[1:10,1:10]
mid_parent_19CS_20MA = mid_parent_19CS_20MA[,-1]
mid_parent_20MA_20MA





write.csv(a, "mp_cor.NIR.20LY.csv")

x = as.data.frame(colnames(pheno_20CS))
a$sepc = x$`colnames(pheno_20CS)`
a = as.data.frame(a, check.names = FALSE)

#### plotting mid-parent correlation
data = read.csv("mp_cor.NIR_CS_LY.csv")


a = ggplot(data, aes(x=wave, y=cor))+
  geom_line()+
  scale_y_continuous("Correlation of Mid-parent Reflectance with Hybrid Reflectance")+
  scale_x_continuous(bquote("NIRS Bands"~(nm)),guide = guide_axis(n.dodge=1) )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  facet_grid(~env)+
  guides(color = TRUE)+
  #ggtitle("A) Genotypic values of parents")+
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=10, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=15, face="bold", colour = "black"),    
    axis.title.y = element_text(size=8, face="bold", colour = "black"),    
    axis.text.x = element_text(size=10, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=8, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 15, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.4),
    axis.line.y = element_line(color="black", size = 0.4),
    panel.border = element_rect(colour = "black", fill=NA, size=0.4),
    #axis.text.x.bottom = element_blank()
  )

a

jpeg("corbands_mid_parents.jpeg",width = 9,height =4,units = "in", res=600)
a
dev.off()


dim(data)
data_CS = subset(data, env == "20LY")
which.min(data_CS$cor)
data_CS[1426,]
