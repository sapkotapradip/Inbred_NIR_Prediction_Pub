# Inbred_NIR_Prediction
NIR of parents and hybrids to predict performance for publication

### pipeline of data analysis for inbred NIR data....
Step 1: calculated blues for each hybrids and inbreds at each environments
Step 2: Use blues to calculate mid-parent heterosis, high parent heterosis for each wavelengths (reflectance)
Step 3: fit into six models that are presented in proposal
Step 4: use lasso, deep learning, or AI models at first; try machine learning models such as svm, random forest, elastic net and lasso
Step 5: Fit multi- trait phenomic model using DA and PH from parents to compare predictability


Step 2: done perfectly; may revise when got time later
Step 3: Fit the model
write a code for MP heterosis, high-parent heterosis, 

# Daniel and Noah inputs 

GCA from CS and other locations 
prediction of new cv schemes leave certain % inbred !
GCA from
leave one parent, both parents left out !
line x tester: one evaluated one left out !

# This was loaded once to create hybrid relationship but no more needed so commented here
# K3=kronecker(K1,K2)
# dim(K3)
# 
# ##### extract only available hybrid combinations from marker matrix
# 
# ################################
# K3.data = as.data.frame(K3)
# 
# complete.list = row.names(G)
# numF = dim(K1)[1]
# Fnames = complete.list[1:numF]
# numM= dim(K2)[1]
# Mnames = complete.list[numF+1:numM]
# 
# total.hybrids = c()
# #loop
# for (f in Fnames) {
#   for (m in Mnames) {
#     total.hybrids = c(total.hybrids, paste(f,"/",m, sep = ""))
#   }
# }
# 
# final.k3 = as.data.frame(matrix(1, ncol = 366, nrow = 366))
# colnames(final.k3) <- available.hybrid
# rownames(final.k3) <- available.hybrid
# 
# 
# colnames(K3.data) <- total.hybrids
# rownames(K3.data) <- total.hybrids
# 
# 
# for(av.hyb1 in available.hybrid){
#   for(av.hyb2 in available.hybrid){
#     # print( K3.data[av.hyb1, av.hyb2])
#     final.k3[av.hyb1, av.hyb2] = K3.data[av.hyb1, av.hyb2]
#   }
# }
# 
# final.k3
# K3 = as.matrix(final.k3)
# # load hybrid_matrix
# save(K3, file = "K3.Rdata")


