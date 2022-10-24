#################### Phylogenetic signal OGTL ####################

library(caper)
library(ggplot2)
library(vegan)
library(phangorn)

# create an ideal TEA tree with neighbours
tea_tree <- read.tree(text = "((Korean,(Japanese,((((Ikema,Tarama),(Hateruma,Yonaguni)),Ogami),((Okinoerabu,Yuwan),Shuri)))),(((((EvenD,EvenB),(Negidal,(Evenki,Solon)))),(((Nanai,Ulch),Orok),(Oroch,Udihe))),((((Buriat,Khalkha),Kalmyk),(Baoan,Mangghuer)),(Chuvash,(Khalaj,((Yakut,(Tuvan,(Khakas, Shor))),(((Bashkir,Tatar),CrimTatar),(((Turkish,Gagauz),Azerbaijani),Turkmen))))))));")

par(mar=c(rep(1,4)))
ideal_topology <- plot(tea_tree,font = 1, label.offset = 0.4, cex = 0.8)
ggsave("./OGTL/plots/Figure_1.pdf", ideal_topology, width=7,height=7)

# simulate branch lengths for the ideal TEA tree with neighbours
tea_tree_br.len <- compute.brlen(tea_tree, method = "Grafen", power = 1)

# binarize the tree
tea_tree_br.len_di <- di2multi(tea_tree_br.len)

# read in raw data

tea_data_raw <- read.csv("./OGTL/data/data_OGTL.csv", sep=";") 

# prepare the data: drop the features with no variation
drops <- names(tea_data_raw[, sapply(tea_data_raw, function(v) var(v, na.rm=TRUE)==0)])
tea_with_variation <- tea_data_raw[ , !(names(tea_data_raw) %in% drops)]

# calculate D for an ideal topology
# the code strikes at 59, thinking there is no variation in that feature

for (i in 1:59) {
  small_dataframe <- tea_with_variation[,c(1,i)]
  colnames(small_dataframe) <- c("ID","V1")
  print(colnames(tea_with_variation)[i])
  #print(small_dataframe)
  result <- phylo.d(data = small_dataframe, 
                    names.col = ID, phy = tea_tree_br.len_di, 
                    binvar = V1, permut = 10000)
  sink("fritz_d.txt",append=TRUE)
  print(colnames(tea_with_variation)[i])
  print(result)
  sink()
}

for (i in 60:152) {
  small_dataframe <- tea_with_variation[,c(1,i)]
  colnames(small_dataframe) <- c("ID","V1")
  print(colnames(tea_with_variation)[i])
  #print(small_dataframe)
  result <- phylo.d(data = small_dataframe, 
                    names.col = ID, phy = tea_tree_br.len_di, 
                    binvar = V1, permut = 10000)
  sink("fritz_d.txt",append=TRUE)
  print(colnames(tea_with_variation)[i])
  print(result)
  sink()
}

rownames(tea_with_variation) <- tea_with_variation[,1]
tea_with_variation<- tea_with_variation[,-1]
tea_with_variation_random <- tea_with_variation

# randomize the data
for (i in 1:38) {
  tea_with_variation_random[i,] <- sample(tea_with_variation)
}

tea_data_raw <- read.csv("data_OGTL.csv", sep=";") 
tea_with_variation_random[,152] <- tea_data_raw[,1]
tea_with_variation_random[,152]
colnames(tea_with_variation_random)[152] <- "ID"
colnames(tea_with_variation_random)[152]
rownames(tea_with_variation_random) <- tea_with_variation_random[,152]
tea_with_variation_random <- tea_with_variation_random[,-1]

write.csv(tea_with_variation_random[,-151],"./OGLT/data/tea_with_variation_random.csv")

# calculate D for randomized data

for (i in 1:151) {
  small_dataframe <- tea_with_variation_random[,c(151,i)]
  colnames(small_dataframe) <- c("ID","V1")
  print(colnames(tea_with_variation_random)[i])
  #print(small_dataframe)
  result <- phylo.d(data = small_dataframe, 
                    names.col = ID, phy = tea_tree_br.len_di, 
                    binvar = V1, permut = 10000)
  sink("fritz_d_random.txt",append=TRUE)
  print(colnames(tea_with_variation_random)[i])
  print(result)
  sink()
}

# Plot Brownian phylogenetic structure for both ideal and cov rel topologies

fritz_d <- read.csv("./OGTL/data/fritz_d_table.txt", sep = "\t")
fritz_d_random <- read.csv("fritz_d_random_table.txt", sep = "\t")

ggplot() + 
  geom_histogram(data=fritz_d_random, aes(fritz_d_random$EstimatedD), 
                 fill="#D53E4F", alpha = 0.3)

sd(fritz_d_random$EstimatedD)
mean(fritz_d_random$EstimatedD)
cut_off_point <- mean(fritz_d_random$EstimatedD) - 2*sd(fritz_d_random$EstimatedD)

mean(fritz_d$EstimatedD)
sd(fritz_d$EstimatedD)
summary(fritz_d$EstimatedD)
length(fritz_d$EstimatedD[fritz_d$EstimatedD<cut_off_point])/length(fritz_d$EstimatedD)*100

plot_random_estimated <- ggplot() + 
  geom_histogram(data=fritz_d, aes(fritz_d$EstimatedD), 
                 fill="#9EBCDA", alpha = 0.8) +
  geom_histogram(data=fritz_d_random, aes(fritz_d_random$EstimatedD),
                 fill="#D53E4F", alpha = 0.3) +
  labs(x="Estimated D", y="Frequency")+
  theme(text = element_text(size=12),legend.text=element_text(size=12)) +
  geom_vline(xintercept = cut_off_point, col="red")

ggsave("./OGTL/plots/Figure_2.pdf", plot_random_estimated, width=7,height=5)

# extract all the rows with the stable features
unstable_d <- subset(fritz_d, fritz_d$EstimatedD>cut_off_point)
stable_d <- subset(fritz_d, fritz_d$EstimatedD<cut_off_point)
length(stable_d$ID)

# save the data frame with stable features
tea_with_variation_stable <- tea_with_variation[ , -which(names(tea_with_variation) %in% unstable_d$ID)]
tea_with_variation_stable <- tea_with_variation_stable[,-1]
tea_with_variation_stable <- tea_with_variation_stable[,-40]

write.csv(tea_with_variation_stable,"./OGTL/data/stable_dataframe.csv")

fritz_d <- read.csv("./OGTL/data/fritz_d_table.txt", sep = "\t")
features <- read.csv("./OGTL/data/feature_description.txt", sep="\t")
total <- merge(fritz_d,features,by="ID")

write.csv(total,"total.csv")

stable_dataframe_d <- total[total$EstimatedD<cut_off_point,]

# read in results on D
fritz_d_feature <- read.csv("./OGTL/results/fritz_d_feature.txt", sep = "\t")

ggplot(data = fritz_d_feature, aes(x=fritz_d_feature$sort, y=fritz_d_feature$EstimatedD)) + 
  geom_boxplot(aes(fill=sort)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot D per part of speech

pos <- read.csv("./OGTL/data/pos.txt", sep = "\t")

plot_pos_d <- ggplot(data = pos, aes(x=PartOfSpeech, y=EstimatedD)) + 
  geom_boxplot(aes(fill=PartOfSpeech)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="Estimated D", y="Frequency") +
  geom_hline(yintercept = min(fritz_d_random$EstimatedD), col="red")

ggsave("./OGTL/plots/Figure_3.pdf", plot_pos_d, width=7,height=5)

table(pos$PartOfSpeech)

# Compare topologies of trees with the whole feature set and with stable
# features only

tea_mcct <- read.nexus("./OGTL/models/OGTL_cov_rel_final/OGTL_cov_rel_final_mcct.trees")
tea_mcct_stable <- read.nexus("./OGTL/models/OGTL_cov_rel_stable_final/OGTL_cov_rel_stable_final_mcct.trees")

par(mar=c(rep(1,4)))
par(mfrow=c(1,2))
plot(tea_mcct,font = 1, label.offset = 0.4, cex = 0.8)
plot(tea_mcct_stable,font = 1, label.offset = 0.4, cex = 0.8)
