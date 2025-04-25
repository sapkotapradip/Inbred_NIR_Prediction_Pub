#############DPAC
library(adegenet)
dapc1 = dapc(pheno_combined[,-c(1:7)], pheno_combined$env)
p = scatter(dapc1)

jpeg("dapc1.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()

inbreds_combined = rbind(parents_20CS, parents_20MA)
dapc2 = dapc(inbreds_combined[,-c(1:5)], inbreds_combined$env)
scatter(dapc2)

inbreds_20CS = parents_20CS
parents_20CS$group = c(rep("A",45), rep("R",44))

dapc2 = dapc(parents_20CS[,-c(1:5,4206)], parents_20CS$group)
scatter(dapc2)

inbreds_20MA = parents_20MA
parents_20MA$group = c(rep("A",45), rep("R",44))

dapc2 = dapc(parents_20MA[,-c(1:5,4206)], parents_20MA$group)
scatter(dapc2)



jpeg("pred_acc.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()
