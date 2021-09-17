################################################ 
# Genome-wide Test of AF-vapeR analysis for Finch data
################################################ 

source("R/angle_grinder_functions.R")
source("~/Exeter/five_aside/five_aside_functions.R")
library(tidyverse)
library(ggplot2)
library(vcfR)

# Input VCF path
vcf_path <- "data/111_sample_maf0.05.vcf.gz"

# Get species groups...
C.pallidus_Z <-  c("PAL1.bam_PAL1.bam","PAL2.bam_PAL2.bam","PAL3.bam_PAL3.bam","PAL4.bam_PAL4.bam","PAL5.bam_PAL5.bam")					
C.parvulus_Z	<-  c("Sample_STF1","Sample_STF13","Sample_STF17","Sample_STF18","Sample_STF2","Sample_STF4","Sample_STF5","Sample_STF6","Sample_STF7","Sample_STF8","PARV1.bam_PARV1.bam","PARV2.bam_PARV2.bam")
P.crassirostris_Z	<-  c("PL15.bam","PL16.bam","PL4.bam","PL7.bam","PL9.bam")	
P.crassirostris_Z <- paste0(P.crassirostris_Z,"_",P.crassirostris_Z)
G.propinqua_G	<-  c("CG1.bam","CG10.bam","CG2.bam","CG3.bam","CG4.bam","CG5.bam","CG6.bam","CG7.bam","CG8.bam","CG9.bam")
G.propinqua_G	<-  paste0(G.propinqua_G,"_",G.propinqua_G)
G.acutirostris_G <- c("DG101.bam","DG103.bam","DG107.bam","DG108.bam")	
G.acutirostris_G <- paste0(G.acutirostris_G,"_",G.acutirostris_G)
G.magnirostris_G <-	c("MG2.bam","MG3.bam","MG4.bam","MG7.bam","mg11.bam")					
G.magnirostris_G <- paste0(G.magnirostris_G,"_",G.magnirostris_G)
C.psittacula_P <- c("Sample_LTF1","Sample_LTF22","Sample_LTF27","Sample_LTF28","Sample_LTF3","Sample_LTF380","Sample_LTF4","Sample_LTF5","Sample_LTF7","Sample_LTF8")
G.difficilis_P <-  c("DP10.bam","DP11.bam","DP13.bam","DP2.bam","DP3.bam","DP4.bam","DP5.bam","DP6.bam","DP7.bam","DP9.bam")
G.difficilis_P <-  paste0(G.difficilis_P,"_",G.difficilis_P)
C.pauper_F <-	 c("Sample_MTF1","Sample_MTF11","Sample_MTF13","Sample_MTF15","Sample_MTF16","Sample_MTF17","Sample_MTF19","Sample_MTF2","Sample_MTF4","Sample_MTF5","PAU1.bam_PAU1.bam","PAU2.bam_PAU2.bam")
G.fuliginosa_S <-  c("Sample_ful1","Sample_ful12","Sample_ful13","Sample_ful16","Sample_ful19","Sample_ful2","Sample_ful22","Sample_ful23","Sample_ful3","Sample_ful8","ful1.bam_ful1.bam","ful2.bam_ful2.bam")
G.magnirostris_M <-  c("Sample_9287","Sample_9309","Sample_9379","Sample_9387","Sample_9398","Sample_9703","Sample_9709","Sample_9737","Sample_9773","Sample_9779")
P.inornata_C <-  c("PIN1.bam","PIN2.bam","PIN3.bam","PIN4.bam","PIN5.bam","PIN6.bam","PIN7.bam","PIN8.bam")
P.inornata_C <- paste0(P.inornata_C,"_",P.inornata_C)
L.noctis <- c("L7.bam","L1.bam","L2.bam","L3.bam","L6.bam")			
L.noctis <- paste0(L.noctis,"_",L.noctis)
#T.bicolor	<- c("T1.bam","	T2.bam","	T3.bam")

# List out
species <- list(C.pallidus_Z,C.parvulus_Z,P.crassirostris_Z,G.propinqua_G,G.acutirostris_G,G.magnirostris_G,C.psittacula_P,G.difficilis_P,C.pauper_F,G.fuliginosa_S,G.magnirostris_M,P.inornata_C,L.noctis)
names(species) <- c("C.pallidus_Z","C.parvulus_Z","P.crassirostris_Z","G.propinqua_G","G.acutirostris_G","G.magnirostris_G","C.psittacula_P","G.difficilis_P","C.pauper_F","G.fuliginosa_S","G.magnirostris_M","P.inornata_C","L.noctis")
# Get the inds
inds <- read.table("data/finch_inds.txt")

# Write some popmaps
for(i in 1:length(species)){
  write.table(data.frame(ind=species[[i]]),
              paste0("data/finch_species_",names(species)[i],".popmap"),
              row.names = F,col.names = F,quote = F)
}

# Check for missing inds...
inds[!(inds[,1] %in% unlist(species)),1] # All species are present

# Build the vectors
finch_vectors <- list(# Allopatry (Cocos)
  c("P.inornata_C","G.propinqua_G"),
  c("P.inornata_C","G.magnirostris_M"),
  c("P.inornata_C","G.difficilis_P"),
  # Allopatry
  c("C.pauper_F","C.psittacula_P"),
  c("G.magnirostris_G","G.magnirostris_M"),
  c("G.fuliginosa_S","C.parvulus_Z"),
  c("G.magnirostris_M","C.parvulus_Z"),
  c("G.magnirostris_M","C.pallidus_Z"),
  # Sympatry
  c("G.magnirostris_G","G.propinqua_G"),
  c("G.magnirostris_G","G.acutirostris_G"),
  c("P.crassirostris_Z","C.parvulus_Z"),
  c("P.crassirostris_Z","C.pallidus_Z"))

# Separate vector pairs just for HMGA2 haplogroup - specifcally small > large
hmga2_vectors <-  list(# Allopatry (Cocos)
  c("P.inornata_C","G.propinqua_G"),
  c("P.inornata_C","G.magnirostris_M"),
  # Allopatry
  c("C.pauper_F","C.psittacula_P"),
  c("C.parvulus_Z","G.magnirostris_M"),
  # Sympatry
  c("G.acutirostris_G","G.magnirostris_G"),
  c("C.parvulus_Z","P.crassirostris_Z"))

# ALX1 haplotypes, pointed > blunt
alx1_vectors <- list(# Allopatry (Cocos)
  c("P.inornata_C","G.magnirostris_M"),
  # Allopatry
  c("C.parvulus_Z","G.magnirostris_M"),
  c("C.pallidus_Z","G.magnirostris_M"),
  # Sympatry
  c("G.propinqua_G","G.magnirostris_G"),
  c("G.acutirostris_G","G.magnirostris_G"))

names(finch_vectors) <- paste0("Comparison",1:length(finch_vectors))
names(hmga2_vectors) <- paste0("Comparison",1:length(hmga2_vectors))
names(alx1_vectors) <- paste0("Comparison",1:length(alx1_vectors))

# Get popmap
species_to_keep2 <- unique(unlist(finch_vectors))
popmap_finch <- data.frame(rbindlist(lapply(species_to_keep2,function(x){
  out <- data.frame(ind=species[x],
                    pop=x)
  colnames(out) <- c("ind","species")
  return(out)
})))

# Get the chr
chrs <- unique(read.table("data/finch_chrs.txt")[,1])
which(chrs == "JH739900")
# HMGA-2 is on:
chr_of_interest <- "JH739900"
which(chrs == chr_of_interest)
# ALX=1 is on:
chr_of_interest2 <- "JH739921"
which(chrs == chr_of_interest2)

# Let's only keep reasonable contigs
finch_fai <- read.table("data/finch_genome.fai")
large_chr <- finch_fai[finch_fai$V2 > 5000000,1]

# Get null windows proportional to chr size
perms = 10000
prop_chr_sizes <- finch_fai[finch_fai[,1] %in% large_chr,2]/sum(finch_fai[finch_fai[,1] %in% large_chr,2])
chr_perms <- as.integer(perms * prop_chr_sizes)
names(chr_perms) <- large_chr

######################################################################## Run ####
# What dataset are we looking at?
run_names <- c("alx1_finch_final","hmga2_finch_final")
vector_run_list <- list("alx1_finch_final"=alx1_vectors,
                        "hmga2_finch_final"=hmga2_vectors)

for(run in run_names){
  
  run_name <- run
  vector_run <- vector_run_list[[run]]
  
  # Run per chromosome
  for(chr in large_chr){
    
    # Subset popmap if needed
    popmap_run <- popmap_finch[popmap_finch[,2] %in% unlist(vector_run),]
    
    # Read in the VCF
    system(paste0("bcftools view -r ",chr," ",vcf_path," > outputs/tmp.vcf"))
    vcf_in <- read.vcfR("outputs/tmp.vcf")
    system("rm -f outputs/tmp.vcf")
    
    message(paste0("Starting chr ",which(large_chr == chr)," of ",length(large_chr)))
    
    # Run over these chromosome windows
    for(wind_size in c(50,100,200,500)){

      # Run over this chromosome...
      AF_vectors <- calc_AF_vectors(vcf = vcf_in,popmap = popmap_run,vectors = vector_run,n_cores = 6,window_size = wind_size,end_cutoff = wind_size/2)
      
      # Get our nulls
      #null_vectors <- calc_null_AF_vectors(vcf = vcf_in,popmap = popmap_finch,vectors = finch_vectors,n_cores = 6,window_size = wind_size,perm = chr_perms[chr])
      null_vectors <- calc_AF_vectors(vcf = vcf_in,popmap = popmap_run,vectors = vector_run,n_cores = 6,window_size = wind_size,end_cutoff = wind_size/2,null_perms = chr_perms[chr])
      cutoffs <- find_null_cutoff(null_vectors,0.99)
      
      # Get Eigen
      eigen_res <- lapply(AF_vectors,eigen_analyse_vectors)
      
      plots <- eigenval_plot(eigen_res,cutoffs = cutoffs[,1])
      plots[[1]]+geom_vline(xintercept=c(411351,428776),linetype="dotted",colour="blue2")
      
      saveRDS(list(AF_vectors,null_vectors,eigen_res),
              paste0("outputs/",run_name,"_chr",chr,"_windowsize_",wind_size,".rds"))
    }
  }
}
