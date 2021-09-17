################################################ 
# PCA over chr15 region for pops...
# Purpose of this script is to visualise complex multi-parallelism at a candidate region...
################################################ 

lib <- c("adegenet","vcfR","GenotypePlot","pbmcapply","tidyverse","ggplot2","parallel","regioneR","bedr","cowplot")
lapply(lib,library,character.only=T)

# What region and VCF
region <- "chr15:5028361-5066375"
chr="chr15"
start=5028361
end=5066375
vcf_path <- "~/Exeter/VCFs/five_aside_STAR_3033083_allchr_shapeit_beagle.vcf.gz"
popmap <- data.frame(ind=system(paste0("bcftools query -l ",vcf_path),intern = T))
popmap$pop <- gsub('[0-9]+', '', popmap$ind)
popmap$pop <- gsub('_M', '', popmap$pop)
popmap$pop <- gsub('_F', '', popmap$pop)

# Do genotype plot and keep the PCA
geno_plot <- genotype_plot(vcf_path,
                           chr= 'chr15',
                           start=start,
                           end=end,
                           popmap = popmap,
                           cluster=T)
geno_pca <- geno_plot$cluster_pca

# Fetch centroids
scores <- geno_pca$li[,1:2]
scores$ind <- rownames(scores)
for(i in 1:nrow(scores)){
  scores$pop[i] <- popmap[popmap$ind == scores$ind[i],"pop"]
}
score_centroids <- data.frame(scores %>% group_by(pop) %>% summarise(PC1=mean(Axis1),
                                                                     PC2=mean(Axis2)))
                    
# Add river and pred
pops <- c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")
river <- rep(c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas"),each=2)    
pred <- rep(c("HP","LP"),5)
for(i in 1:nrow(score_centroids)){
  score_centroids$river[i] <- river[which(pops == score_centroids$pop[i])]
  score_centroids$pred[i] <- pred[which(pops == score_centroids$pop[i])]
}

# Now plot...
ggplot(score_centroids,aes(PC1,PC2,shape=pred,colour=river))+
  geom_point(size=5)+
  geom_path(aes(group=river),arrow = arrow(),colour="black")









