# Analysis of Drosophila experimental parallelism from Kelly and Hughes 2019
lib <- c("regioneR","tidyverse","ggplot2","parallel","data.table","pbmcapply","ggVennDiagram","cowplot")
lapply(lib,library,character.only=T)
devtools::load_all("~/Exeter/afvaper/")
n_cores = 6

# Window size
window_snps = 25

# Prepare original results... ---------------------------------------------
original_results <- "data/kelly_hughes_supp_data.xlsx"
original_res <- data.table(na.omit(readxl::read_xlsx(original_results,skip = 1)))
colnames(original_res) <- c("chr","pos","ref","alt",
                            "A0_reads","A0_pR",
                            "B0_reads","B0_pR",
                            "C0_reads","C0_pR",
                            "A7_reads","A7_pR",
                            "B7_reads","B7_pR",
                            "C7_reads","C7_pR",
                            "LRT-parallel")
original_res$chr <- gsub("Scf_","",original_res$chr)

significant_snps <- data.table(na.omit(readxl::read_xlsx("data/kelly_hughes_supp_data_S8.xlsx")))
significant_snps$plot_y = 4.2
significant_snps$Chrom <- gsub("Scf_","",significant_snps$Chrom)
significant_snps$chr <- significant_snps$Chrom
significant_snps_LRT <- significant_snps[significant_snps$`LRT sig` != "ns",]
significant_snps_CMH <- significant_snps[significant_snps$CMHsig != "ns",]

# Also assemble the '30' unique loci
significant_snps_LRT_unique <- significant_snps_LRT
significant_snps_LRT_unique$snp_flank3 <- significant_snps_LRT_unique$snp - 500000
significant_snps_LRT_unique$snp_flank5 <- significant_snps_LRT_unique$snp + 500000
significant_snps_LRT_unique$snp_flank3[significant_snps_LRT_unique$snp_flank3 < 0] <- 0

snp_ranges <- toGRanges(significant_snps_LRT_unique[,c("Chrom","snp_flank3","snp_flank5"),with=F])
significant_snps_LRT_unique <- data.frame(mergeRegions(snp_ranges,snp_ranges))
significant_snps_LRT_unique$window_id <- paste0(significant_snps_LRT_unique$seqnames,":",significant_snps_LRT_unique$start,"-",significant_snps_LRT_unique$end)
max_locus_LRT <- rbindlist(lapply(1:nrow(significant_snps_LRT_unique),function(x){
  
  # Subset for snps in the window
  tmp <- significant_snps_LRT[significant_snps_LRT$Chrom == significant_snps_LRT_unique$seqnames[x] &
                                significant_snps_LRT$snp >= significant_snps_LRT_unique$start[x] & 
                                significant_snps_LRT$snp <= significant_snps_LRT_unique$end[x],]
  # Return
  return(data.frame(max_snp = paste0(tmp[tmp$LRT == max(tmp$LRT),"Chrom"],":",tmp[tmp$LRT == max(tmp$LRT),"snp"]),
                    max_LRT = max(tmp$LRT)))
}))
significant_snps_LRT_unique <- cbind(significant_snps_LRT_unique,max_locus_LRT)

# Calculate per window linkage for our AF-vapeR windows based on the AF frequencies...
original_res_winds <- rbindlist(lapply(unique(original_res$chr),function(chr_tmp){
  tmp <- original_res[chr == chr_tmp,]
  window_N <- floor(nrow(tmp)/window_snps)
  tmp <- tmp[1:(window_snps*window_N),]
  tmp$window <- paste0(chr_tmp,"-",rep(1:window_N,each = window_snps))
  tmp[,c("chr","pos","window",grep("pR",colnames(tmp),value = T)),with=F]
}))

# Now do correlations per window...
original_res_winds_linkage <- original_res_winds[,
                                                 .(LD=mean(cor(t(.SD[,grep("pR",colnames(.SD),value = T),with=F]))^2),
                                                   window_id=paste0(unique(.SD$chr),":",min(.SD$pos),"-",max(.SD$pos))),
                                                 by=window]

# Plot genome-wide linkage
original_res_winds_linkage <- separate(original_res_winds_linkage,col="window_id",sep = ":",into = c("chr","pos")) %>%
  separate(col="pos",sep = "-",into = c("start","end"))
ggplot(original_res_winds_linkage,aes(as.integer(start),LD))+
  geom_point()+
  facet_grid(.~chr, scales = "free", space='free',switch="x")+
  theme_minimal()+
  theme(axis.text=element_text(size=16),
        axis.title = element_text(size=18),
        title = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        strip.text = element_text(size=14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(y=expression(Allele~Frequency~R^2),x="Chromosome")

# Process AFvapeR results -------------------------------------------------
drosophila_res <- grep(paste0(window_snps,".rds"),list.files("outputs/",pattern = "drosophila_test_"),value=T)

# Read them all back in...
all_eigen_res <- merge_eigen_res(lapply(drosophila_res,function(res){
  tmp <- readRDS(paste0("outputs/",res))
  tmp[[3]]
}))

all_AFV <- merge_eigen_res(lapply(drosophila_res,function(res){
  tmp <- readRDS(paste0("outputs/",res))
  tmp[[1]]
}))

all_null_AFV <- merge_eigen_res(lapply(drosophila_res,function(res){
  tmp <- readRDS(paste0("outputs/",res))
  tmp[[2]]
}))

null_cutoffs <- afvaper:::find_null_cutoff(all_null_AFV,c(0.95,0.99,0.999))

# What's the average window size???
window_info <- afvaper:::window2pos_df(names(all_eigen_res))
window_info$wind_size <- as.integer(window_info$end) - as.integer(window_info$start)
median(window_info$wind_size)
window_info$window_id <- names(all_eigen_res)

# Plot
eig1_fig <- afvaper::eigenval_plot(all_eigen_res,null_vectors = all_null_AFV,plot.pvalues = T)[[1]]

# Histogram of pvalues
hist(-log10(drosophila_pvals$Eigenvalue_1))

# Fetch significant windows...
afvaper_signif <- data.frame(window_id = signif_eigen_windows(all_eigen_res,cutoffs = null_cutoffs[,1])[[1]])
afvaper_signif <- cbind(afvaper_signif,window2pos_df(afvaper_signif$window_id))
afvaper_signif$start <- as.integer(afvaper_signif$start)
afvaper_signif$end <- as.integer(afvaper_signif$end)

# Summarise these
afvaper_signif_summary <- afvaper:::summarise_window_parallelism(afvaper_signif$window_id,all_eigen_res,eigenvector = 1)

# How many of these overlap with significant LRT1 SNPs... -----------------
afvaper_signif_list <- lapply(3:1,function(x) {
  tmp = data.frame(window_id = signif_eigen_windows(all_eigen_res,cutoffs = null_cutoffs[,x])[[1]])
  tmp_info = afvaper:::summarise_window_parallelism(tmp$window_id,all_eigen_res)
  return(data.frame(window_id = tmp_info[tmp_info$parallel_pops == "A,B,C","window_id"]))
})
afvaper_signif_dd <- cbind(rbindlist(afvaper_signif_list),
                           signif = c(rep("<0.001",nrow(afvaper_signif_list[[1]])),
                                      rep("<0.01",nrow(afvaper_signif_list[[2]])),
                                      rep("<0.05",nrow(afvaper_signif_list[[3]]))))
afvaper_signif_dd <- afvaper_signif_dd[!duplicated(afvaper_signif_dd$window_id),]

# Intersect these
signif_snp_overlap <- rbindlist(pbmclapply(1:nrow(significant_snps_LRT),function(snp){
  
  tmp <- window_info[window_info$chr == significant_snps_LRT$chr[snp] &
                       as.integer(window_info$start) <= significant_snps_LRT$snp[snp] &
                       as.integer(window_info$end) >= significant_snps_LRT$snp[snp],]
  data.frame(chr=significant_snps_LRT$chr[snp],
             pos=significant_snps_LRT$snp[snp],
             snp = paste0(significant_snps_LRT$chr[snp],":",significant_snps_LRT$snp[snp]),
             window_id = paste0(tmp$chr,":",tmp$start,"-",tmp$end),
             is_overlap = ifelse(paste0(tmp$chr,":",tmp$start,"-",tmp$end) %in% afvaper_signif_dd$window_id,"Overlapping","Not Overlapping"),
             LRT = significant_snps_LRT$LRT[snp])
},mc.cores=n_cores))

# Add the level of significance...
signif_snp_overlap <- merge(signif_snp_overlap,afvaper_signif_dd,by="window_id",all=T)
signif_snp_overlap <- signif_snp_overlap[!is.na(signif_snp_overlap$snp),]
signif_snp_overlap$signif[is.na(signif_snp_overlap$signif)] <- "None"

# Add the plotting y value...
signif_snp_overlap$plot_y <- 4
plot_y_vals <- c(4.2,4.35,4.5)
signif_vals <- c("<0.05","<0.01","<0.001")
for(i in 1:length(plot_y_vals)){
  signif_snp_overlap[signif_snp_overlap$signif == signif_vals[i],"plot_y"] <- plot_y_vals[i]
}

# How many windows have hit LRT SNPs?
table(signif_snp_overlap$is_overlap)/sum(table(signif_snp_overlap$is_overlap))

# Fetch the pvals and plot again overlapped with significant LRT snps
drosophila_pvals <- data.frame(afvaper:::eigen_pvals(all_eigen_res,all_null_AFV))
drosophila_pvals <- cbind(drosophila_pvals,afvaper:::window2pos_df(rownames(drosophila_pvals)))
LRT_eig1_overlap <- ggplot(drosophila_pvals,aes(y=-log10(Eigenvalue_1),x=as.integer(start)))+
  geom_step()+
  facet_grid(.~chr, scales = "free", space='free',switch="x")+
  geom_point(data=signif_snp_overlap,aes(y=plot_y,x=pos,colour=signif))+
  theme_minimal()+
  theme(axis.text=element_text(size=16),
        axis.title = element_text(size=18),
        title = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        strip.text = element_text(size=14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size=16),
        legend.text= element_text(size=14))+
  scale_colour_manual(values=c("None"="gray50",
                               "<0.05" = "red2",
                               "<0.01" = "orange2",
                               "<0.001" = "gold2"))+
  labs(y=expression(-log[10](p)),x="Chromosome",colour="")

# # How many of these overlap with UNIQUE LRT1 SNPs... -----------------
# signif_snp_overlap <- rbindlist(pbmclapply(1:nrow(significant_snps_LRT_unique),function(snp){
#   
#   tmp <- window_info[window_info$chr == significant_snps_LRT_unique$seqnames[snp] &
#                        as.integer(window_info$start) <= significant_snps_LRT$snp[snp] &
#                        as.integer(window_info$end) >= significant_snps_LRT$snp[snp],]
#   data.frame(chr=significant_snps_LRT$chr[snp],
#              pos=significant_snps_LRT$snp[snp],
#              snp = paste0(significant_snps_LRT$chr[snp],":",significant_snps_LRT$snp[snp]),
#              afvaper_window = paste0(tmp$chr,":",tmp$start,"-",tmp$end),
#              is_overlap = ifelse(paste0(tmp$chr,":",tmp$start,"-",tmp$end) %in% afvaper_signif$window_id,"Overlapping","Not Overlapping"))
# },mc.cores=n_cores))
# signif_snp_overlap$plot_y <- ifelse(signif_snp_overlap$is_overlap == "Overlapping",4.2,4)
# 
# # How many windows have hit LRT SNPs?
# table(signif_snp_overlap$is_overlap)/sum(table(signif_snp_overlap$is_overlap))
# 
# # Fetch the pvals and plot again overlapped with significant LRT snps
# drosophila_pvals <- data.frame(afvaper:::eigen_pvals(all_eigen_res,all_null_AFV))
# drosophila_pvals <- cbind(drosophila_pvals,afvaper:::window2pos_df(rownames(drosophila_pvals)))
# ggplot(drosophila_pvals,aes(y=-log10(Eigenvalue_1),x=as.integer(start)))+
#   geom_step()+
#   facet_grid(.~chr, scales = "free", space='free',switch="x")+
#   geom_point(data=signif_snp_overlap,aes(y=plot_y,x=pos,colour=is_overlap))+
#   theme_minimal()+
#   theme(axis.text=element_text(size=16),
#         axis.title = element_text(size=18),
#         title = element_text(size=20),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.spacing.x=unit(0.1, "lines"),
#         strip.text = element_text(size=14),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())+
#   scale_colour_manual(values=c("Not Overlapping"="gray50",
#                                "Overlapping" = "red2"))+
#   labs(y=expression(Eig~1~-log[10](p)),x="Chromosome",colour="")



# How many windows overlap with significant CMH SNPs ----------------------
signif_snp_overlap_CMH <- rbindlist(pbmclapply(1:nrow(significant_snps_CMH),function(snp){
  
  tmp <- window_info[window_info$chr == significant_snps_CMH$chr[snp] &
                       as.integer(window_info$start) <= significant_snps_CMH$snp[snp] &
                       as.integer(window_info$end) >= significant_snps_CMH$snp[snp],]
  data.frame(chr=significant_snps_CMH$chr[snp],
             pos=significant_snps_CMH$snp[snp],
             snp = paste0(significant_snps_CMH$chr[snp],":",significant_snps_CMH$snp[snp]),
             window_id = paste0(tmp$chr,":",tmp$start,"-",tmp$end),
             is_overlap = ifelse(paste0(tmp$chr,":",tmp$start,"-",tmp$end) %in% afvaper_signif_dd$window_id,"Overlapping","Not Overlapping"),
             CMH = significant_snps_CMH$CMH[snp])
},mc.cores=n_cores))

# Add the level of significance...
signif_snp_overlap_CMH <- merge(signif_snp_overlap_CMH,afvaper_signif_dd,by="window_id",all=T)
signif_snp_overlap_CMH <- signif_snp_overlap_CMH[!is.na(signif_snp_overlap_CMH$snp),]
signif_snp_overlap_CMH$signif[is.na(signif_snp_overlap_CMH$signif)] <- "None"

# How many windows have hit LRT SNPs?
table(signif_snp_overlap_CMH$is_overlap)/sum(table(signif_snp_overlap_CMH$is_overlap))

# Fetch the pvals and plot again overlapped with significant LRT snps
drosophila_pvals <- data.frame(afvaper:::eigen_pvals(all_eigen_res,all_null_AFV))
drosophila_pvals <- cbind(drosophila_pvals,afvaper:::window2pos_df(rownames(drosophila_pvals)))
ggplot(drosophila_pvals,aes(y=-log10(Eigenvalue_1),x=as.integer(start)))+
  geom_step()+
  facet_grid(.~chr, scales = "free", space='free',switch="x")+
  geom_point(data=signif_snp_overlap_CMH,aes(y=plot_y,x=pos,colour=is_overlap))+
  theme_minimal()+
  theme(axis.text=element_text(size=16),
        axis.title = element_text(size=18),
        title = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        strip.text = element_text(size=14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  scale_colour_manual(values=c("Not Overlapping"="gray50",
                               "Overlapping" = "red2"))+
  labs(y=expression(Eig~1~-log[10](p)),x="Chromosome",colour="")


# Plot intersection of all LRT, CMH, and AF-vapeR outliers ---------------------------
# Repeat the plotting over three cutoffs...
LRT_winds <- unique(signif_snp_overlap$window_id)
CMH_winds <- unique(signif_snp_overlap_CMH$window_id)

# Venn diagram for each set of significant afvaper windows...
vaper_labs <- c("<0.001","<0.01","<0.05")
signif_venns <- lapply(length(afvaper_signif_list):1,function(vaper_i){
  ggVennDiagram(list(LRT_winds,afvaper_signif_list[[vaper_i]][,1],CMH_winds),
                label_alpha=0,colour="black",
                category.names = c("LRT",paste0("AF-vapeR ",vaper_labs[vaper_i]),"CMH"))+
    scale_fill_gradient(low="white",high = "red2")+
    scale_color_manual(values=c("gray50","gray50","gray50"))+
    theme(legend.position = "none")
})
pdf("figs/FigureX_drosophila_LRT_CMH_overlap.pdf",width=14,height=7)
plot_grid(plot_grid(plotlist = signif_venns,ncol=3,labels=c("A","B","C"),label_size = 24),
          LRT_eig1_overlap,
          ncol=1,
          rel_heights=c(1.2,1),
          labels=c("","D"),label_size = 24,
          align = "v",axis = "tblr")
dev.off()

# # Replace windows by SNPs that intersect the two tests...
# afvaper_signif_window_snps <- afvaper_signif$window_id
# for(wind in afvaper_signif_window_snps){
#   if(wind %in% signif_snp_overlap$afvaper_window){
#     afvaper_signif_window_snps <- afvaper_signif_window_snps[afvaper_signif_window_snps != wind]
#     afvaper_signif_window_snps <- c(afvaper_signif_window_snps,signif_snp_overlap[signif_snp_overlap$afvaper_window == wind,snp])
#   } 
#   if (wind %in% signif_snp_overlap_CMH$afvaper_window){
#     afvaper_signif_window_snps <- afvaper_signif_window_snps[afvaper_signif_window_snps != wind]
#     afvaper_signif_window_snps <- c(afvaper_signif_window_snps,signif_snp_overlap_CMH[signif_snp_overlap_CMH$afvaper_window == wind,snp])
#   }
# }
afvaper_signif_window_snps <- unique(afvaper_signif_window_snps)

# Plot Venn...
venn_list <- list(LRT_snps,CMH_snps,afvaper_signif_window_snps)
names(venn_list) <- c("LRT","CMH","AF-vapeR")
ggvenn(venn_list)

# How does the overlapping of windows vs parallel SNPs change by l --------
# Compare to window size
signif_snp_overlap$window_id <- signif_snp_overlap$afvaper_window
signif_snp_overlap <- merge(signif_snp_overlap,window_info[,c("window_id","wind_size")],by="window_id")
ggplot(signif_snp_overlap,aes(x=is_overlap,y=wind_size))+
  geom_violin(draw_quantiles = c(0.5))

# Compare to LD - Clear effect here...
original_res_winds_linkage$window_id <- paste0(original_res_winds_linkage$chr,":",original_res_winds_linkage$start,"-",original_res_winds_linkage$end)
signif_snp_overlap <- merge(signif_snp_overlap,original_res_winds_linkage[,c("LD","window_id"),with=F],by='window_id')
ggplot(signif_snp_overlap,aes(x=is_overlap,y=LD))+
  geom_violin(draw_quantiles = c(0.5))

# Compare to the LRT - Similar but marginal loss of signal...
signif_snp_overlap[,median(LRT),by=is_overlap]
ggplot(signif_snp_overlap,aes(x=is_overlap,y=log10(LRT)))+
  geom_violin(draw_quantiles = c(0.5))

# Compare eig1 pvals to LD ------------------------------------------------
drosophila_pvals$LD <- original_res_winds_linkage$LD
cor.test(-log10(drosophila_pvals$Eigenvalue_1),drosophila_pvals$LD,method = "kendall")
ggplot(drosophila_pvals,aes(LD,-log10(Eigenvalue_1)))+
  geom_point()+
  stat_density2d_filled(alpha=0.5)

# Estimate linkage vs eig1 for null AFV... --------------------------------
null_windows <- afvaper::window2pos_df(names(all_null_AFV))
null_windows$start <- as.integer(null_windows$start)
null_windows$end <- as.integer(null_windows$end)
null_windows$window_id <- names(all_null_AFV)

# Loop over and calculate linkage...
null_windows_linkage <- rbindlist(pbmclapply(1:nrow(null_windows),function(x){
  
  tmp <- original_res[original_res$chr == null_windows$chr[x] &
                        original_res$pos >= as.integer(null_windows$start[x]) &
                        original_res$pos <= as.integer(null_windows$end[x]),grep("pR",colnames(original_res),value=T),with=F]
  cor_res <- mean(cor(t(tmp))^2)
  return(data.frame(window_id=null_windows$window_id[x],
                    LD=cor_res))
},mc.cores=n_cores))

# Merge with eigenvalues and plot...
null_eigen_res <- lapply(all_null_AFV,eigen_analyse_vectors)
null_eigenvals <- lapply(null_eigen_res,'[[',1)
null_empP <- afvaper:::eigen_pvals(null_eigen_res,all_null_AFV)

# Add them on
null_windows_linkage$eig1 <- sapply(null_eigenvals,function(x) x[1])
null_windows_linkage$eig1_p <- null_empP[,1]

# Compare correlation
cor.test(null_windows_linkage$LD,-log10(null_windows_linkage$eig1_p),method = "kendall")

# Plot
ggplot(null_windows_linkage,aes(LD,-log10(eig1_p)))+
  geom_point()+
  stat_density2d_filled(alpha=0.5)


# Split null windows into deciles on the basis of LD and re-estimate pvals --------
LD_deciles <- c(0,quantile(null_windows_linkage$LD,probs = seq(0.1,1,0.1)))
LD_based_null_cutoffs <- lapply(2:length(LD_deciles),function(x){
  
  # Get the windows for this decile
  windows_to_keep <- null_windows_linkage[null_windows_linkage$LD <= LD_deciles[x] &
                                            null_windows_linkage$LD > LD_deciles[x-1],window_id]
  
  # Get null cutoffs for this window
  cutoff_tmp <- find_null_cutoff(all_null_AFV[windows_to_keep],cutoffs = c(0.95,0.99,0.999))
  return(list(cutoffs = cutoff_tmp,
              null_AFV = all_null_AFV[windows_to_keep]))
})
names(LD_based_null_cutoffs) <- seq(0.1,1,0.1)

# Inspect the eig1 cutoffs based on different deciles of linkage...
decile_eig1_cutoffs <- lapply(LD_based_null_cutoffs,'[[',1)
to_plot <- data.frame(quantile = names(LD_based_null_cutoffs),
                      eig1_cutoff = sapply(decile_eig1_cutoffs,function(x) x[1,2]))
ggplot(to_plot,aes(quantile,eig1_cutoff))+
  geom_point()

# Assign each of the windows to the quantile it is represents in the null...
original_res_winds_linkage$null_LD_quantile <- sapply(original_res_winds_linkage$LD,function(LD) names(LD_based_null_cutoffs)[max(which(LD_deciles < LD))])
original_res_winds_linkage$window_id <- paste0(original_res_winds_linkage$chr,":",original_res_winds_linkage$start,"-",original_res_winds_linkage$end)

# Calculate new empirical pvalues on the basis of LD-decile null distirbutions. Note that these will have a reduced maximum to a power of 10 because each one only has 1,000 null windows...
LD_decile_empPvals <- rbindlist(lapply(names(LD_based_null_cutoffs),function(quant){
  
  # Find the windows...
  windows_to_test <- original_res_winds_linkage[original_res_winds_linkage$null_LD_quantile == quant,window_id]
  
  # Test the windows...
  test_window_pvals <- data.table(afvaper:::eigen_pvals(all_eigen_res[windows_to_test],
                                                        LD_based_null_cutoffs[[as.character(quant)]]$null_AFV))
  
  # Add back in linkage info and window ids...
  test_window_pvals <- cbind(test_window_pvals,
                             original_res_winds_linkage[original_res_winds_linkage$null_LD_quantile == quant,c("LD","window_id"),with=F])
  
  return(test_window_pvals)
}))

# Plot these new results...
LD_decile_empPvals <- cbind(LD_decile_empPvals,afvaper:::window2pos_df(LD_decile_empPvals$window_id))
ggplot(LD_decile_empPvals,aes(as.integer(start),-log10(Eigenvalue_1)))+
  geom_point()+
  facet_grid(.~chr, scales = "free", space='free',switch="x")+
  theme_minimal()+
  theme(axis.text=element_text(size=16),
        axis.title = element_text(size=18),
        title = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        strip.text = element_text(size=14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(y="Eig1 (-log10(p))",x="Chromosome")

# Take significant windows from these and compare against the LRT-1 SNPs...
LD_based_signif_windows <- LD_decile_empPvals[-log10(LD_decile_empPvals$Eigenvalue_1) > 2,window_id]
signif_snp_overlap_LD <- rbindlist(pbmclapply(1:nrow(significant_snps_LRT),function(snp){
  
  tmp <- window_info[window_info$chr == significant_snps_LRT$chr[snp] &
                       as.integer(window_info$start) <= significant_snps_LRT$snp[snp] &
                       as.integer(window_info$end) >= significant_snps_LRT$snp[snp],]
  data.frame(chr=significant_snps_LRT$chr[snp],
             pos=significant_snps_LRT$snp[snp],
             snp = paste0(significant_snps_LRT$chr[snp],":",significant_snps_LRT$snp[snp]),
             afvaper_window = paste0(tmp$chr,":",tmp$start,"-",tmp$end),
             is_overlap = ifelse(paste0(tmp$chr,":",tmp$start,"-",tmp$end) %in% LD_based_signif_windows,"Overlapping","Not Overlapping"),
             LRT = significant_snps_LRT$LRT[snp])
},mc.cores=n_cores))
signif_snp_overlap_LD$plot_y <- ifelse(signif_snp_overlap_LD$is_overlap == "Overlapping",4.2,4)

# How many windows have hit LRT SNPs?
table(signif_snp_overlap_LD$is_overlap)/sum(table(signif_snp_overlap_LD$is_overlap))

# Plot again...
ggplot(LD_decile_empPvals,aes(as.integer(start),-log10(Eigenvalue_1)))+
  geom_point()+
  facet_grid(.~chr, scales = "free", space='free',switch="x")+
  theme_minimal()+
  theme(axis.text=element_text(size=16),
        axis.title = element_text(size=18),
        title = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        strip.text = element_text(size=14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(y="Eig1 (-log10(p))",x="Chromosome")+
  geom_point(data=signif_snp_overlap_LD,aes(y=plot_y,x=pos,colour=is_overlap))+
  scale_colour_manual(values=c("Not Overlapping"="gray50","Overlapping" = "red2"))                                                                                                         

# What are we missing...
signif_snp_overlap_LD$window_id <- signif_snp_overlap_LD$afvaper_window
signif_snp_overlap_LD <- merge(signif_snp_overlap_LD,original_res_winds_linkage[,c("LD","window_id"),with=F],by='window_id')
ggplot(signif_snp_overlap_LD,aes(x=is_overlap,y=LD))+
  geom_violin(draw_quantiles = c(0.5))
ggplot(signif_snp_overlap_LD,aes(x=is_overlap,y=LRT))+
  geom_violin(draw_quantiles = c(0.5))
