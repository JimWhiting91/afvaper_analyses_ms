################################################ 
# Genome-wide Test of AF-vapeR analysis for HP-LP guppies...
################################################ 

lib <- c("tidyverse","ggplot2","parallel")
lapply(lib,library,character.only=T)
devtools::load_all("~/Exeter/afvaper/")
source("~/Exeter/five_aside/five_aside_functions.R")

# Test inputs
five_aside_test <- "~/Exeter/VCFs/five_aside_STAR_3033083_final.vcf.gz"
five_aside_test_chr20 <- "~/Exeter/VCFs/five_aside_STAR_chr20_scaf94_merged_shapeit_beagle.vcf.gz"
CORES <- 6
permutations <- 10000
wind_sizes <- c(25,50,100,200)

# Make the vectors
test_vectors=list(c("LT","UT"),c("GH","GL"),c("APHP","APLP"),c("LO","UQ"),c("LMD","UMD"))
names(test_vectors) <- c("Tac","Guanapo","Aripo","Oropouche","Madamas")

# Get per chrom perms for a total of N permutations
chr_sizes <- read.table("~/Exeter/Genomes/STAR.chromosomes.release.fasta.fai")
chr_sizes <- chr_sizes[chr_sizes$V2 > 1000000,]
chr_props <- chr_sizes$V2/sum(chr_sizes$V2)
chr_perms <- ceiling(chr_props * permutations)

# Run over the following chr
chrs <- chr_sizes$V1
names(chr_perms) <- chrs 

# Merge chr20 and scaf94
chr_perms["chr20"] <- chr_perms["chr20"] + chr_perms["000094F_0"]
chr_perms <- chr_perms[which(names(chr_perms) != "000094F_0")]
chrs <- names(chr_perms)

#chr_res <- lapply(chrs,function(chr){
for(chr in chrs){
  
  message(paste0("STARTING CHR ",which(chrs == chr)," of ",length(chrs)))
  
  # Make the input - we merge chr20 and scaf94 as described in Whiting et al. 2021 Plos Genetics
  if(chr != "chr20"){
    system(paste0("bcftools view -r ",chr," ",five_aside_test," > outputs/tmp.vcf"))
    vcf_in <- read.vcfR("outputs/tmp.vcf")
    system("rm -f outputs/tmp.vcf")
  } else {
    vcf_in <- read.vcfR(five_aside_test_chr20)
  }
  
  # Make the popmaps
  test_popmap <- data.frame(inds=colnames(vcf_in@gt)[2:ncol(vcf_in@gt)])
  test_popmap$pop <- gsub('[0-9]+', '', test_popmap$inds)
  test_popmap$pop <- gsub('_F', '', test_popmap$pop)
  test_popmap$pop <- gsub('_M', '', test_popmap$pop)
  
  # Loop over window sizes
  for(window_snps in wind_sizes){
    
    # Calculate allele frequency matrices
    AF_input <- calc_AF_vectors(vcf = vcf_in,
                                window_size = window_snps,
                                popmap = test_popmap,
                                vectors = test_vectors,
                                n_cores = CORES)
    
    # Calculate null frequencies
    null_input <- calc_AF_vectors(vcf = vcf_in,
                                  window_size = window_snps,
                                  popmap = test_popmap,
                                  vectors = test_vectors,
                                  n_cores = CORES,
                                  null_perms = chr_perms[chr])
    
    # Calculate Eigen statistics
    eigen_res <- lapply(AF_input,eigen_analyse_vectors)
    
    # Get null matrix
    null_cutoffs <- find_null_cutoff(null_res = null_input,
                                     cutoffs = c(0.95,0.99,0.999))

    # Save both
    saveRDS(list(AF_input,null_input,eigen_res,null_cutoffs),
            paste0("outputs/five_aside_test_",chr,"_AF_eigen_res_windsize_",window_snps,".rds"))
  }
}