##################################################
# Script visualises AFvapeR results over simulations from different recombination rates
##################################################

lib <- c("ggplot2","data.table","pbmcapply","tidyverse","cowplot")
lapply(lib,library,character.only=T)

# Name of sim run
sim_name <- "21,07,26_recombination_intervals_REC_"
res_dir <- list.files("outputs/slim/",pattern="21,07,26_recombination_intervals_REC_")
rec_rate_N <- length(res_dir)

# Fetch the rec rates
rec_rates <- read.table("data/sim_recombination_rates.txt")

# What window sizes?
wind_sizes <- c(50,200,500)

# Read in each and save to file all eigen_res plus avg
all_rec_res <- data.frame(rbindlist(pbmclapply(1:rec_rate_N,function(x){
  
  # Fetch the eigen_res
  vaper_res_files <- list.files(paste0("outputs/slim/",sim_name,x,"/final_results/"),pattern = "EIGEN_RES")
  
  # Read in all files from a given window
  all_window_res <- data.frame(rbindlist(lapply(wind_sizes,function(wind){
    
    # Subset for window size files
    tmp_wind_files <- grep(paste0("wind",wind,".EIGEN"),vaper_res_files,value=T)
    
    # Loop over these and read in
    tmp_window_res <- data.frame(rbindlist(lapply(tmp_wind_files,function(input){
      
      # Get res
      tmp_res <- readRDS(paste0("outputs/slim/",sim_name,x,"/final_results/",input))
      
      # Reformat...
      eigen_res <- tmp_res[[2]]
      pos <- sapply(strsplit(names(eigen_res),":"),'[[',2)
      start <- as.integer(sapply(strsplit(pos,"-"),'[[',1))
      end <- as.integer(sapply(strsplit(pos,"-"),'[[',2))
      all_eig1 <- sapply(eigen_res,function(i){
        return(i$eigenvals[1])
      })
      out <- data.frame(start=start,
                        end=end,
                        eig1=all_eig1,
                        sim_group=input)
    })))
    
    # Add in window size and rec_rate
    tmp_window_res$window_size <- wind
    tmp_window_res$rec_rate <- rec_rates$V1[x]
    
    # We also want to have separately, the averaged value...
    tmp_window_res$mid <- rowMeans(tmp_window_res[,c("start","end")])
    tmp_window_res$window <- cut_interval(tmp_window_res$mid,length = 50000)
    return(tmp_window_res)
  })))
  
  # Return all
  return(all_window_res)
},mc.cores = 4)))

# We also want to have separately, the averaged per rec_rate, wind_size and window
avg_eig1 <- data.frame(all_rec_res %>% group_by(window_size,rec_rate,window) %>%
                         summarise(eig1=mean(eig1)))

# Get the start of windows
avg_eig1$window_starts <- sapply(strsplit(as.character(avg_eig1$window),","),'[[',1)
avg_eig1$window_starts <- gsub("[()]","",avg_eig1$window_starts)
avg_eig1$window_starts <- as.numeric(gsub("\\[|\\]","",avg_eig1$window_starts))

# Fix rec rate
avg_eig1$rec_rate_trans <- round(log10(avg_eig1$rec_rate),3)

# And plot
rec_peaks <- ggplot(avg_eig1,aes(x=window_starts+25000,y=eig1))+
  geom_line()+
  facet_grid(window_size~rec_rate_trans)+
  geom_vline(xintercept = 1e6,linetype="dotted")+
  geom_hline(yintercept = 4,linetype="solid",colour="red2")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=15),
        strip.text = element_text(size=14))+
  labs(y="Eigenvalue 1",x="Position")

# Visualise the output from a whole chromosome ----------------------------
# Use parameter combo 15
# s = 0.05
# mig = 0.0025
# DP = 10000

# Get recombination groups from log files
get_rec_groups <- function(log_dir){
  logs <- list.files(log_dir)
  
  # Fetch IDs
  recs <- c("low","med","high")
  rec_sim_ids <- lapply(recs,function(rec){
    log_tmp <- grep(rec,logs,value=T)
    
    sim_outs <- sapply(log_tmp,function(log){
      tmp_in <- read.table(paste0(log_dir,"/",log),fill=T)
      return(tmp_in[2,1])
    })
    return(sim_outs)
  })
  names(rec_sim_ids) <- recs
  return(rec_sim_ids)
}

input_dir = "outputs/slim/21,07,27_Variable_recombination_along_chromosome_fixed_generations_15/final_results/"
window_size = 200
eigenvec = 1
n_cores = 4
variable_recomb = TRUE
iterations = 50

# Get the recombination groups if needs be
if(variable_recomb){
  rec_groups <- get_rec_groups(gsub("/final_results","",input_dir))
}

# Get the inputs
tmp_inputs <- list.files(input_dir)
tmp_inputs <- grep("EIGEN",tmp_inputs,value=T)
tmp_inputs <- grep(paste0("wind",window_size,"."),tmp_inputs,value=T)

# Fetch results and merge
all_res <- data.frame(rbindlist(mclapply(1:length(tmp_inputs),function(iter){
  
  # Get raw res
  tmp <- readRDS(paste0(input_dir,"/",tmp_inputs[iter]))
  eigen_res <- tmp[[2]]
  
  # Reformat
  plot_format <- data.frame(winds=names(eigen_res))
  plot_format$eigen <- sapply(eigen_res,function(x){return(x[[1]][eigenvec])})
  
  # Get mids
  positions <- data.frame(rbindlist(lapply(plot_format$winds,function(wind){
    pos <- strsplit(wind,":")[[1]][2]
    start <- as.integer(strsplit(pos,"-")[[1]][1])
    end <- as.integer(strsplit(pos,"-")[[1]][2])
    mids <- mean(c(start,end))
    data.frame(start,end,mids)
  })))
  
  plot_format <- cbind(plot_format,positions)
  
  # Add iteration ID
  plot_format$iteration <- iter
  
  # Add recombination identifier
  tmp_id <- strsplit(tmp_inputs[iter],"_")[[1]][2]
  plot_format$rec <- c("low","med","high")[grep(tmp_id,rec_groups)]
  
  return(plot_format)
},mc.cores=n_cores)))

# Get rec averages...
rec_rates <- c("low","med","high")
window_avg <- data.frame(rbindlist(lapply(rec_rates,function(rate){
  print(rate)
  wind_size=50000
  tmp_plot <- all_res[all_res$rec == rate,]
  
  # Windowise
  start <- seq(0,max(tmp_plot$mids),wind_size)
  end <-   start+wind_size
  
  # Fetch the avg
  window_calcs <- sapply(start,function(pos){
    mean(tmp_plot[tmp_plot$end >= pos &
                    tmp_plot$start <= pos+wind_size,"eigen"])
  })
  
  # output
  out <- data.frame(mids = start+(wind_size/2),
                    eigen = window_calcs,
                    rec = rate)
  return(out)
})))

# Visualise for each recombination
window_avg$rec_F <- factor(window_avg$rec,levels=rec_rates)
all_res$rec_F <- factor(all_res$rec,levels=rec_rates)
base_plot <- ggplot(window_avg,aes(x=mids,y=eigen))+
  facet_wrap(~rec_F,ncol=1,strip.position = "right")+
  theme_minimal()+
  labs(y="Eigenvalue 1", x="Chromosome Position")+
  scale_x_continuous(breaks=seq(0,max(window_avg$mids),2e6),
                     labels=seq(0,max(window_avg$mids),2e6)/1e6)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=14))

# Now add all the replicates
for(i in unique(all_res$iteration)){
  base_plot <- base_plot +
    geom_line(data=all_res[all_res$iteration == i,],
              aes(x=mids,y=eigen),colour="black",alpha=0.1)
}

# Add back in avg
base_plot <- base_plot +
  geom_line(colour="green3")

# Add in recombination regions
recombination_along_chr <- base_plot +
  scale_x_continuous(breaks=seq(1,max(window_avg$mids),2e6)+1e6,
                     labels=c(seq(1,10,1.8)*1e-9,seq(2.8,10,1.8)*1e-8))+
  labs(x="Recombination Rate/Chromosome Position")


# Final recombination figure ----------------------------------------------
pdf("figs/FigureSX_effect_of_recombination.pdf",width=14,height=10)
plot_grid(rec_peaks,
          recombination_along_chr,
          ncol=1,nrow=2,labels="AUTO",label_size=28,
          align="v",axis = "tblr")
dev.off()




