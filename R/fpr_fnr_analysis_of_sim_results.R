# Analysis of SLiM simulations for AF-vapeR simulations
rm(list=ls())
lib <- c("pbmcapply","cowplot","data.table","ggplot2","dplyr","viridis")
sapply(lib,library,character.only=T)
devtools::load_all("~/Exeter/afvaper/")

# Get our results
sim_results <- "22,04,20_Rerunning_full_parallel_with_new_v2_null"
metadata <- read.table("data/sim_parameter_metadata_v3.txt")
colnames(metadata) <- c("selection","pop2","migration")

# Remember to update this!
sim_perms <- 1:27

# Eigen SNPs - update accordingly
eigen_snps <- c(50,200,500)

###### FUNCTION LIBRARY #######
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

# Get FST fpr and fnr
fst_fpr_fnr <- function(input_dir,variable_recomb=FALSE,iterations=1:100){
  
  # these are daughters
  daughter_pops <- paste0("p",2:5)
  
  # Get the recombination groups if needs be
  if(variable_recomb){
    rec_groups <- get_rec_groups(gsub("/final_results","",input_dir))
    iterations <- 1:length(unlist(rec_groups))
  }
  
  # Get the inputs
  tmp_inputs <- list.files(input_dir)
  tmp_inputs <- grep(".fst",tmp_inputs,value=T)
  sim_IDs <- unique(sapply(strsplit(tmp_inputs,"_finished.trees"),function(i){return(i[[1]][1])}))
  
  # Fetch results and merge
  all_res <- lapply(iterations,function(iter){
    
    sim_ID_tmp <- sim_IDs[iter]
    
    # Get recombination type if needs be
    if(variable_recomb){
      rec_tmp <- c("low","med","high")[grep(gsub("slim_","",sim_ID_tmp),rec_groups)]
    }
    
    # Read in all 4 fsts 
    daughter_res <- data.frame(rbindlist(lapply(2:5,function(daughter){
      
      # Fetch specific results
      to_read <- grep(sim_ID_tmp,tmp_inputs,value=T)
      to_read <- grep(paste0("pop",daughter),to_read,value=T)
      tmp_in <- data.frame(fread(paste0(input_dir,"/",to_read)))
      
      # Add run information
      tmp_in$iteration <- iter
      tmp_in$daughter <- paste0("p",daughter)
      
      # Mark outliers for each alpha
      tmp_in$outlier <- "No"
      for(alpha in alphas){
        tmp_in[tmp_in$WEIGHTED_FST > quantile(tmp_in$WEIGHTED_FST,alpha),"outlier"] <- alpha
      }
      return(tmp_in)
    })))
    
    # Based on all daughters, did we detect the selected site?
    covergence_tiers <- 4:2
    positive_mat <- matrix(nrow=3,ncol=length(alphas))
    negative_mat <- matrix(nrow=3,ncol=length(alphas))
    for(alpha in 1:length(alphas)){
      
      signif <- daughter_res[daughter_res$outlier %in% alphas[alpha:length(alphas)],]
      
      # Count up all the windows
      window_counts <- table(signif$BIN_START)
      
      # Regard false positives as those outside 1 Mb from selected site
      if(!(variable_recomb)){
        false_positives <- window_counts[as.integer(names(window_counts)) < 9500000 |
                                           as.integer(names(window_counts)) > 10500000]
      } else {
        # Or if recombination is variable
        if(rec_tmp == "low"){
          false_positives <- window_counts[as.integer(names(window_counts)) > 2000000]
          
        } else if (rec_tmp == "med"){
          false_positives <- window_counts[as.integer(names(window_counts)) < 10900000 |
                                             as.integer(names(window_counts)) > 11100000]
          
        } else if (rec_tmp == "high"){
          false_positives <- window_counts[as.integer(names(window_counts)) < 20090000 |
                                             as.integer(names(window_counts)) > 21010000]
          
        }
      }
      # False positive count
      positive_mat[1,alpha] <- length(false_positives[false_positives==4])/length(unique(daughter_res$BIN_START))
      positive_mat[2,alpha] <- length(false_positives[false_positives %in% c(3,4)])/length(unique(daughter_res$BIN_START))
      positive_mat[3,alpha] <- length(false_positives[false_positives %in% c(2,3,4)])/length(unique(daughter_res$BIN_START))
      
      # Are there false negatives?
      if(!(variable_recomb)){
        true_positives <- window_counts[as.integer(names(window_counts)) > 9950000 &
                                          as.integer(names(window_counts)) < 10050000]
      } else {
        # Or if recombination is variable
        if(rec_tmp == "low"){
          true_positives <- window_counts[as.integer(names(window_counts)) >= 900000 & 
                                            as.integer(names(window_counts)) <= 1100000]
          
        } else if (rec_tmp == "med"){
          true_positives <- window_counts[as.integer(names(window_counts)) >= 10090000 |
                                            as.integer(names(window_counts)) <= 11010000]
          
        } else if (rec_tmp == "high"){
          true_positives <- window_counts[as.integer(names(window_counts)) >= 20090000 |
                                            as.integer(names(window_counts)) <= 21010000]
          
        }
      }
      
      # False positive count
      negative_mat[1,alpha] <- ifelse(any(true_positives==4),1,0)
      negative_mat[2,alpha] <- ifelse(any(true_positives %in% c(4,3)),1,0)
      negative_mat[3,alpha] <- ifelse(any(true_positives %in% c(2,3,4)),1,0)
    }
    colnames(positive_mat) <- alphas
    rownames(positive_mat) <- covergence_tiers
    colnames(negative_mat) <- alphas
    rownames(negative_mat) <- covergence_tiers
    
    # Return these
    return(list(positive_mat,negative_mat))
  })
  
  # Label with IDs
  names(all_res) <- sim_IDs
  
  # To get the final false positive matrix, average through and retain standard error...
  all_fpr_mats <- lapply(all_res,function(iter){
    return(iter[[1]])
  })
  
  # Handle separately if rec rate varies
  if(!(variable_recomb)){
    fpr_avg <- apply(simplify2array(all_fpr_mats), 1:2, mean)
    fpr_sd <- apply(simplify2array(all_fpr_mats), 1:2, sd)
    fpr_res <- list(avg=fpr_avg,
                    sd=fpr_sd)
  } else {
    names(all_fpr_mats) <- sim_IDs
    fpr_res <- lapply(c("low","med","high"),function(rate){
      
      # Get IDs
      rec_IDs <- paste0("slim_",rec_groups[[rate]])
      rec_fpr_mats <- all_fpr_mats[rec_IDs]
      fpr_avg <- apply(simplify2array(rec_fpr_mats), 1:2, mean)
      fpr_sd <- apply(simplify2array(rec_fpr_mats), 1:2, sd)
      fpr_res <- list(avg=fpr_avg,
                      sd=fpr_sd)
    })
    names(fpr_res) <- c("low","med","high")
  }
  
  # For final false negative, we just get a single value for all runs
  all_fnr_mats <- lapply(all_res,function(iter){
    return(iter[[2]])
  })
  
  if(!(variable_recomb)){
    fnr_sum <- apply(simplify2array(all_fnr_mats), 1:2, sum)    
    fnr_rates <- 1-(fnr_sum/length(iterations))
  } else {
    
    fnr_rates <- lapply(c("low","med","high"),function(rate){
      
      # Get IDs
      rec_IDs <- paste0("slim_",rec_groups[[rate]])
      rec_fnr_mats <- all_fnr_mats[rec_IDs]
      
      fnr_sum <- apply(simplify2array(rec_fnr_mats), 1:2, sum)    
      out_rates <- 1-(fnr_sum/length(iterations))
      return(out_rates)
    })
  }
  names(fnr_rates) <- c("low","med","high")
  
  # Return both
  return(list(FPR=fpr_res,
              FNR=fnr_rates))
}

# Get eigen_fpr_fnr
eigen_fpr_fnr <- function(input_dir,window_size,eigenvec=1,variable_recomb=FALSE){
  
  # Get the recombination groups if needs be
  if(variable_recomb){
    rec_groups <- get_rec_groups(gsub("/final_results","",input_dir))
    iterations <- 1:length(unlist(rec_groups))
  }
  
  # Get the inputs
  tmp_inputs <- list.files(input_dir)
  tmp_inputs <- grep("EIGEN",tmp_inputs,value=T)
  tmp_inputs <- grep(paste0("wind",window_size,".EIGEN"),tmp_inputs,value=T)
  
  if(!variable_recomb){
    iterations = 1:length(tmp_inputs)
  }
  
  # Set alphas
  alphas <- c(0.95,0.99,0.999)
  
  # Fetch results and merge
  all_res <- lapply(iterations,function(iter){
    
    # Get raw res
    tmp <- readRDS(paste0(input_dir,"/",tmp_inputs[iter]))
    af <- tmp[[1]]
    eigen_res <- tmp[[2]]
    null <- tmp[[3]]
    
    # Get cutoffs
    null_cuts <- find_null_cutoff(null,cutoffs = alphas)
    
    # Get all the signif windows
    signif_list <- lapply(1:length(alphas),function(alpha){
      signif_eigen_windows(eigen_res,cutoffs = null_cuts[,alpha])
    })
    names(signif_list) <- alphas
    
    # Return the signif windows
    out_list <- lapply(1:length(alphas),function(z){
      out = data.frame(window=names(eigen_res),
                       eigenvector=eigenvec,
                       iteration=iter,
                       alpha=names(signif_list)[z])
      out$signif <- "No"
      out[out$window %in% signif_list[[z]][[eigenvec]],"signif"] <- "Yes"
      return(out)
    })
    return(out_list)
    
  })
  
  # Summarise each
  alpha_res <- lapply(1:length(alphas),function(alpha){
    tmp2 <- data.frame(rbindlist(lapply(all_res,function(tmp){
      return(tmp[[alpha]])
    })))
    
    # Transform to analyse
    pos <- sapply(strsplit(tmp2$window,":"),function(x){return(x[[2]])})
    start <- as.integer(sapply(strsplit(pos,"-"),function(x){return(x[[1]])}))
    end <- as.integer(sapply(strsplit(pos,"-"),function(x){return(x[[2]])}))
    tmp2$start <- start
    tmp2$end <- end
    tmp2$mid <- rowMeans(tmp2[,c("start","end")])
    
    return(tmp2)
  })
  
  # Get all the FPR and FNR rates
  if(!(variable_recomb)){
    out_mat <- matrix(ncol=length(alphas),nrow=2)
    rownames(out_mat) <- c("FPR","FNR")
    colnames(out_mat) <- alphas
    
    # Loop to fill - Here we define false positives that are > 0.5 Mb either side of the selected site
    for(i in 1:length(alphas)){
      res_tmp <- alpha_res[[i]]
      
      # Get all FPR (per iterations)
      iter_fpr <- sapply(iterations,function(iter){
        res_tmp2 <- res_tmp[res_tmp$iteration==iter,]
        non_selected <- res_tmp2[res_tmp2$mid < 9500000 | res_tmp2$mid > 10500000,]
        non_selected_counts <- table(non_selected$signif)
        fpr_tmp <- ifelse(non_selected_counts["No"] == nrow(non_selected),0,round(non_selected_counts["Yes"]/sum(non_selected_counts),4))
        return(fpr_tmp)
      })
      
      # Set FPR
      out_mat[1,i] <- mean(iter_fpr)
      
      # Get FNR
      windows_with_selected <- res_tmp[res_tmp$start <= 10050000 &
                                         res_tmp$end >= 9950000,]
      
      # Do we detect a significant window around the selected site in every iteration?
      successful_iters <- unique(windows_with_selected[windows_with_selected$signif=="Yes","iteration"])
      out_mat[2,i] <- round(1-(length(successful_iters)/length(iterations)),4)
    }
  } else { 
    out_mat <- lapply(c("low","med","high"),function(rate){
      out_mat2 <- matrix(ncol=length(alphas),nrow=2)
      rownames(out_mat2) <- c("FPR","FNR")
      colnames(out_mat2) <- alphas
      
      # Loop to fill - Here we define false positives that are > 0.5 Mb either side of the selected site
      for(i in 1:length(alphas)){
        res_tmp <- alpha_res[[i]]
        
        # Get relevant iterations
        iterations_sub <- sapply(rec_groups[[rate]],function(x){
          grep(x,tmp_inputs)
        })
        
        # Get all FPR (per iterations)
        iter_fpr <- sapply(iterations_sub,function(iter){
          res_tmp2 <- res_tmp[res_tmp$iteration==iter,]
          if(rate == "low"){
            non_selected <- res_tmp2[res_tmp2$mid > 2000000,]
          } else if (rate == "med"){
            non_selected <- res_tmp2[res_tmp2$mid < 10900000 | res_tmp2$mid > 11100000,]
          } else if (rate == "high"){
            non_selected <- res_tmp2[res_tmp2$mid < 20090000 | res_tmp2$mid > 21010000,]
          }
          non_selected_counts <- table(non_selected$signif)
          fpr_tmp <- ifelse(non_selected_counts["No"] == nrow(non_selected),0,round(non_selected_counts["Yes"]/sum(non_selected_counts),4))
          return(fpr_tmp)
        })
        
        # Set FPR
        out_mat2[1,i] <- mean(iter_fpr)
        
        # Get FNR
        res_tmp2 <- res_tmp[res_tmp$iteration %in% iterations_sub,]
        if(rate == "low"){
          windows_with_selected <- res_tmp2[res_tmp2$start <= 1100000 &
                                              res_tmp2$end >= 900000,]
        } else if (rate == "med"){
          windows_with_selected <- res_tmp2[res_tmp2$start <= 11010000 &
                                              res_tmp2$end >= 10090000,] 
        } else if (rate == "high"){
          windows_with_selected <- res_tmp2[res_tmp2$start <= 21010000 &
                                              res_tmp2$end >= 20090000,] 
        }
        # Do we detect a significant window around the selected site in every iteration?
        successful_iters <- unique(windows_with_selected[windows_with_selected$signif=="Yes","iteration"])
        out_mat2[2,i] <- round(1-(length(successful_iters)/length(iterations_sub)),4)
      }
      out_mat2  
    })
    
  }
  
  # Return
  return(out_mat)
}

# Visualise eigen results with variable recombination
eigen_visualise <- function(input_dir,window_size,eigenvec=1,n_cores=1){
  
  # Get the recombination groups if needs be
  if(variable_recomb){
    rec_groups <- get_rec_groups(gsub("/final_results","",input_dir))
  }
  
  # Get the inputs
  tmp_inputs <- list.files(input_dir)
  tmp_inputs <- grep("EIGEN",tmp_inputs,value=T)
  tmp_inputs <- grep(paste0("wind",window_size),tmp_inputs,value=T)
  
  # Fetch results and merge
  all_res <- data.frame(rbindlist(mclapply(iterations,function(iter){
    
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
    wind_size=10000
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
    labs(y="Eigenvalue", x="Chromosome Position")+
    scale_x_continuous(breaks=seq(0,max(window_avg$mids),2e6),
                       labels=seq(0,max(window_avg$mids),2e6)/1e6)+
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          strip.text = element_text(size=14))
  
  # Now add all the replicates
  for(i in iterations){
    base_plot <- base_plot +
      geom_line(data=all_res[all_res$iteration == i,],
                aes(x=mids,y=eigen),colour="black",alpha=0.05)
  }
  
  # Add back in avg
  base_plot <- base_plot +
    geom_line(colour="green3")
  
  return(base_plot)
}

# Visualise null distribution
null_distribution_visualise <- function(input_dir,rec_rate="low",window_size=200,n_cores=1){
  
  # Get the recombination groups if needs be
  if(rec_rate != "fixed"){
    rec_groups <- get_rec_groups(gsub("/final_results","",input_dir))
  }
  
  # Get the inputs
  tmp_inputs <- list.files(input_dir)
  tmp_inputs <- grep("EIGEN",tmp_inputs,value=T)
  tmp_inputs <- grep(paste0("wind",window_size),tmp_inputs,value=T)
  
  # Filter for rec_rate
  tmp_inputs2 <- sapply(rec_groups[[rec_rate]],function(id){
    return(grep(id,tmp_inputs,value=T))
  })
  
  # Fetch results and merge
  all_res <- data.frame(rbindlist(mclapply(1:length(tmp_inputs2),function(iter){
    
    # Get raw res
    tmp <- readRDS(paste0(input_dir,"/",tmp_inputs2[iter]))
    eigen_res <- tmp[[2]]
    null_res <- tmp[[3]]
    
    # Get null eigens
    null_eigens <- lapply(null_res,eigen_analyse_vectors)
    
    # Get eigenvalue distributions and reformat
    plot_format <- data.frame(rbindlist(lapply(null_eigens,function(eigen){
      
      out <- data.frame(eigenvalue=eigen$eigenvals,
                        eigenvector=paste0("Eigenvector ",1:length(eigen$eigenvals)))
      
      # Sum them as well...
      for(i in 2:nrow(out)){
        out[i,"eigenvalue"] <- out[i,"eigenvalue"] + out[i-1,"eigenvalue"]
      }
      
      # Output
      return(out)
    })))
    
  },mc.cores=n_cores)))
  
  # Get quantiles
  quantiles <- data.frame(all_res %>% group_by(eigenvector) %>%
                            summarise(mean=mean(eigenvalue),
                                      lower=quantile(eigenvalue,probs = 0.05),
                                      upper=quantile(eigenvalue,probs = 0.95)))
  
  # Plot these
  quantiles$eigenvector <- gsub("Eigenvector ","",quantiles$eigenvector)
  ggplot(quantiles,aes(eigenvector,y=mean))+
    geom_pointrange(aes(ymin=lower,ymax=upper),shape=15,size=1)+
    labs(y=expression(sum(lambda[i], i=1, m)),x="Eigenvector (i)")+
    theme_minimal()+
    theme(axis.title=element_text(size=20),
          axis.text = element_text(size=18))
}

# Get AFs from logs
fetch_log_AF <- function(input_dir){
  
  # List logs
  log_files <- grep(".log",list.files(input_dir,pattern = "run_"),value=T)
  
  # For each log file we want to find the final AFs and save as a matrix
  AF_mat <- matrix(ncol=5,nrow=length(log_files))
  for(i in 1:nrow(AF_mat)){
    # Find the AF
    tmp_in <- read.table(paste0(input_dir,"/",log_files[i]),fill=T)
    # Fetch the final freqs from these ugly files
    to_keep <- grep("p1",tmp_in[,1])
    if(length(to_keep) > 0){
      freqs <- tmp_in[to_keep:(to_keep+10),]
      AF_mat[i,] <- as.numeric(freqs[seq(2,10,2),1])
    } else {
      AF_mat[i,] <- NA
    }
  }
  
  # Average across columns
  class(AF_mat) <- "numeric"
  colnames(AF_mat) <- paste0("p",1:5)
  return(AF_mat)
}

###### ANALYSIS #######
# We want to summarise false positive rate and false negative rate for eigen analyses first
for(i in 1:3){
  eigen_snp_window <- eigen_snps[i]
  print(paste0("STARTING SNP WINDOW SIZE:",eigen_snp_window))
  if(!(file.exists(paste0("outputs/slim/",sim_results,"_all_eigen_window",eigen_snp_window,".rds")))){
    all_fpr_fnr <- pbmclapply(sim_perms,function(sim){
      # tmp_dir <- paste0("outputs/slim/",sim_results,"_",sim,"/final_results")
      tmp_dir <- paste0("/Volumes/jimwhiting_external/afvaper_sims/",sim_results,"_",sim,"/final_results")
      eigen_fpr_fnr(input_dir = tmp_dir,window_size = eigen_snp_window,eigenvec = 1)
    },mc.cores=6)
    saveRDS(all_fpr_fnr,
            paste0("outputs/slim/",sim_results,"_all_eigen_window",eigen_snp_window,".rds"))
  }
}

# Fetch all the eigen results for different window sizes and merge
all_fpr_fnr <- data.frame(rbindlist(lapply(1:length(eigen_snps),function(x){
  
  print(x)
  # Read
  tmp <- readRDS(paste0("outputs/slim/",sim_results,"_all_eigen_window",eigen_snps[x],".rds"))
  
  # rbind and add metadata
  tmp2 <- data.frame(rbindlist(lapply(1:length(tmp),function(y){
    tmp_df <- data.frame(tmp[[y]])
    tmp_df <- reshape2::melt(tmp_df)
    colnames(tmp_df) <- c("alpha","rate")
    tmp_df$type <- rep(c("FPR","FNR"),nrow(tmp_df)/2)
    tmp_df$selection <- metadata[y,1]
    tmp_df$pop2 <- metadata[y,2]
    tmp_df$migration <- metadata[y,3]
    return(tmp_df)
  })))
  
  tmp2$window_size <- eigen_snps[x]
  return(tmp2)
})))

# What explains what...
# Selection
ggplot(all_fpr_fnr,aes(x=rate))+
  geom_histogram()+
  facet_grid(selection~type,scales = "free_y")

# Pop2
ggplot(all_fpr_fnr,aes(x=rate))+
  geom_histogram()+
  facet_grid(pop2~type,scales = "free_y")

# Migration
ggplot(all_fpr_fnr,aes(x=rate))+
  geom_histogram()+
  facet_grid(migration~type,scales = "free_y")

################################################################
# We want to summarise false positive rate and false negative rate for Fst next
if(!(file.exists(paste0("outputs/slim/",sim_results,"_all_fst_results.rds")))){
  all_fpr_fnr_fst <- pbmclapply(sim_perms,function(sim){
    # tmp_dir <- paste0("outputs/slim/",sim_results,"_",sim,"/final_results")
    tmp_dir <- paste0("/Volumes/jimwhiting_external/afvaper_sims/",sim_results,"_",sim,"/final_results")
    fst_fpr_fnr(tmp_dir)
  },mc.cores=6)
  saveRDS(all_fpr_fnr_fst,
          paste0("outputs/slim/",sim_results,"_all_fst_results.rds"))
}
# all_paramater_fst_res <- readRDS(paste0("outputs/slim/",sim_results,"_all_fst_results.rds"))
all_paramater_fst_res <- readRDS(paste0("outputs/slim/",sim_results,"_all_fst_results.rds"))

# Merge to a comparable data.frame
all_fst_error <- data.frame(rbindlist(lapply(1:length(all_paramater_fst_res),function(x){
  
  # Get the FPR results
  fpr_res <- reshape2::melt(all_paramater_fst_res[[x]]$FPR$avg)
  colnames(fpr_res) <- c("overlap","alpha","rate")
  fpr_res$type <- "FPR"
  
  # Get the FNR results
  fnr_res <- reshape2::melt(all_paramater_fst_res[[x]]$FNR)
  colnames(fnr_res) <- c("overlap","alpha","rate")
  fnr_res$type <- "FNR"
  
  # Merge
  fpr_fnr <- rbind(fpr_res,fnr_res)
  
  # Add parameters
  fpr_fnr$selection <- metadata[x,1]
  fpr_fnr$pop2 <- metadata[x,2]
  fpr_fnr$migration <- metadata[x,3]
  
  return(fpr_fnr)
})))

ggplot(all_fst_error,aes(x=rate))+
  geom_histogram()+
  facet_grid(selection~type,scales = "free")

# Pop2
ggplot(all_fst_error,aes(x=rate))+
  geom_histogram()+
  facet_grid(pop2~type,scales = "free_y")

# Migration
ggplot(all_fst_error,aes(x=rate))+
  geom_histogram()+
  facet_grid(migration~type,scales = "free_y")

# Compare the two data structures
head(all_fpr_fnr)
head(all_fst_error)

# Tidy up
all_fpr_fnr$alpha <- as.numeric(gsub("X","",all_fpr_fnr$alpha))

# Now merge these together
colnames(all_fpr_fnr)[colnames(all_fpr_fnr)=="window_size"] <- "Window Size/Overlap"
colnames(all_fst_error)[colnames(all_fst_error)=="overlap"] <- "Window Size/Overlap"

# And in the darkness, bind them
plot_rates <- rbind(all_fpr_fnr,all_fst_error)
plot_rates$`Window Size/Overlap` <- factor(plot_rates$`Window Size/Overlap`,levels=c(50,200,500,4,3,2))
plot_rates$selection <- factor(plot_rates$selection,levels=c(0.01,0.05,0.1))
plot_rates$pop2 <- factor(plot_rates$pop2,levels=c(400,2000,10000))
plot_rates$migration <- factor(plot_rates$migration,levels=c(0,0.0025,0.01))

# For each of the parameters, we'll make a 4x2 grid showing results and mean/se
parameters <- c("selection","pop2","migration")
parameter_labs <- c("Selection Coefficient","Daughter Pop Size","Migration Rate")
names(parameter_labs) <- parameters

param_plots <- lapply(parameters,function(param){
  
  # Make new dummy data
  tmp <- plot_rates[,c("alpha","rate","type","Window Size/Overlap",param)]
  colnames(tmp)[ncol(tmp)] <- "Var"
  
  # Get summaries
  sum_plot <- data.frame(tmp %>% 
                           group_by(`Window Size/Overlap`,alpha,type,Var) %>%
                           summarise(mean_rate=mean(rate),
                                     se=sd(rate)/sqrt(length(rate))))
  
  # Add AF-VapeR vs FST flags
  sum_plot$analysis <- "FST"
  sum_plot[sum_plot$Window.Size.Overlap %in% c(50,200,500),"analysis"] <- "AF-vapeR"
  
  # Make param factor again...
  if(class(sum_plot$Var) != "factor"){
    sum_plot$Var <- factor(sum_plot$Var,levels = sort(as.numeric(unique(sum_plot$Var))))
  }
  
  # And plot...
  ggplot(sum_plot[sum_plot$alpha != 0.9,],aes(x=Var,y=mean_rate,colour=Window.Size.Overlap,group=Window.Size.Overlap,shape=analysis))+
    geom_point(position=position_dodge(width = 0.5))+
    geom_linerange(aes(ymin = mean_rate-se, ymax = mean_rate+se),position=position_dodge(width = 0.5),show.legend = F)+
    facet_grid(type~alpha,scales="free_y")+
    labs(y="Error Rate",x=parameter_labs[param],colour="Window Size/Overlap",shape="Analysis")+
    ylim(0,1)
  
})

# Combine
legend <- get_legend(
  # create some space to the left of the legend
  param_plots[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)
pdf("figs/FigureSX_simulation_effects_fprfnr.pdf",width=8,height=8)
plot_grid(
  plot_grid(param_plots[[1]] + theme(legend.position = "none"),
            param_plots[[2]] + theme(legend.position = "none"),
            param_plots[[3]] + theme(legend.position = "none"),labels="AUTO",label_size = 32,ncol=1,align='v',axis="tblr",hjust = 0.1),
  legend,ncol=2,rel_widths = c(3,1))
dev.off()

# Save these to a suppp table...
write.table(plot_rates,"tables/TableSX_afvaper_simulation_error_rates.tsv",
            sep="\t",quote=F,row.names = F)

#### Plot FNR heatmap for each analysis subset ####
analyses <- unique(plot_rates$`Window Size/Overlap`)

# Now plot over different alphas...
# FNR + FPR Matrices for alpha of 0.99
alphas <- c(0.95,0.99,0.999)
for(error_rate in c("FPR","FNR")){
  alpha_heats <- lapply(alphas,function(alpha){
    # Get our tmp data
    tmp <- plot_rates[plot_rates$alpha == alpha &
                        plot_rates$type == error_rate,]
    
    # Make a selection label
    tmp$selection_lab <- paste0("s:",tmp$selection)
    
    # Make new facet labels
    labs <- c(50,200,500,4,3,2)
    new_labs <- c("AFvapeR:50","AFvapeR:200","AFvapeR:500",
                  "F[ST]:4","F[ST]:3","F[ST]:2")
    tmp$`Window Size/Overlap` <- as.character(tmp$`Window Size/Overlap`)
    for(i in 1:length(labs)){
      tmp[tmp$`Window Size/Overlap`== labs[i],"Window Size/Overlap"] <- new_labs[i]
    }
    tmp$`Window Size/Overlap` <- factor( tmp$`Window Size/Overlap`,levels=new_labs)
    
    # Make heatmap...
    p1 <- ggplot(tmp,aes(x=factor(pop2),y=factor(migration),fill=rate))+
      facet_grid(`Window Size/Overlap`~selection_lab,labeller = label_parsed)+
      theme_minimal()+
      theme(panel.grid = element_blank(),
            strip.text = element_text(size = 11),
            axis.text = element_text(size=12),
            axis.title = element_text(size=14),
            strip.background = element_rect(colour="black"))+
      geom_tile()+
      labs(y="Migration",x="Daughter Population Size",fill=error_rate)
    
    if(error_rate == "FPR"){
      p1 +  scale_fill_viridis(option = "A",limits = c(0, 0.2))
    } else {
      p1 +  scale_fill_viridis(option = "A",limits = c(0, 1))
    }
  })
  
  # Also get legend
  # Combine
  legend <- get_legend(
    # create some space to the left of the legend
    alpha_heats[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12),
                             legend.title = element_text(size=18),
                             legend.text = element_text(size=16))
  )
  
  pdf(paste0("figs/FigureX_",error_rate,"_heats_all_alpha.pdf"),width=14,height=8)
  print(cowplot::plot_grid(alpha_heats[[1]]+theme(legend.position = "none",
                                                  axis.text.x = element_text(angle=45,hjust=1))+
                             ggtitle(expression(alpha:0.05)),
                           alpha_heats[[2]]+theme(legend.position = "none",
                                                  axis.text.x = element_text(angle=45,hjust=1))+
                             ggtitle(expression(alpha:0.01)),
                           alpha_heats[[3]]+theme(legend.position = "none",
                                                  axis.text.x = element_text(angle=45,hjust=1))+
                             ggtitle(expression(alpha:0.001)),
                           legend,
                           ncol=4,labels=c("A","B","C",""),label_size = 32,rel_widths = c(2,2,2,0.7)))
  dev.off()
}

# Fetch max FPR for all afvaper and fst sets
for(analysis in unique(plot_rates$`Window Size/Overlap`)){
  print(analysis)
  print(max(plot_rates[plot_rates$type=="FPR" & 
                         plot_rates$`Window Size/Overlap` == analysis &
                         plot_rates$alpha == 0.99,"rate"]))
}

# Fetch max FNR for eigen=50 under weak selection
plot_rates[plot_rates$type == "FNR" &
             plot_rates$selection == 0.01 &
             plot_rates$`Window Size/Overlap` == "50" &
             plot_rates$alpha == 0.950,"rate"]

plot_rates[plot_rates$type == "FNR" &
             plot_rates$selection == 0.01 &
             plot_rates$`Window Size/Overlap` == "4" &
             plot_rates$alpha == 0.950,"rate"]

plot_rates[plot_rates$type == "FNR" &
             plot_rates$selection == 0.01 &
             plot_rates$`Window Size/Overlap` == "2" &
             plot_rates$alpha == 0.950,"rate"]

################################################################
#### Read in logs for AFs ####
AF_list <- pbmclapply(sim_perms,function(iter){
  print(iter)
  # AF_mat <- data.frame(fetch_log_AF(paste0("outputs/slim/",sim_results,"_",iter)))
  AF_mat <- data.frame(fetch_log_AF(paste0("/Volumes/jimwhiting_external/afvaper_sims/",sim_results,"_",iter)))
  
  # Merge this with information about the metadata
  AF_mat$selection <- metadata[iter,1]
  AF_mat$pop2 <- metadata[iter,2]
  AF_mat$migration <- metadata[iter,3]
  return(AF_mat)
},mc.cores=6)

# Combine
AF_dd <- data.frame(rbindlist(AF_list))

# For each row, we want variance of AFD and min AFD...
AF_dd$p2_AFD <- AF_dd$p2 - AF_dd$p1
AF_dd$p3_AFD <- AF_dd$p3 - AF_dd$p1
AF_dd$p4_AFD <- AF_dd$p4 - AF_dd$p1
AF_dd$p5_AFD <- AF_dd$p5 - AF_dd$p1
AF_dd$AFD_mean <- apply(AF_dd[,paste0("p",2:5,"_AFD")],1,mean)
AF_dd$AFD_var <- apply(AF_dd[,paste0("p",2:5,"_AFD")],1,var)
AF_dd$AFD_min <- apply(AF_dd[,paste0("p",2:5,"_AFD")],1,min)

# Summarise within treatments...
AF_dd_treat <- data.frame(AF_dd %>% group_by(selection,migration,pop2) %>% summarise(mean_AFD=mean(AFD_mean,na.rm=T),
                                                                                     var_AFD=mean(AFD_var,na.rm=T),
                                                                                     min_AFD=mean(AFD_min,na.rm=T)))

# Heatmap of AF
AF_dd_treat$selection_lab <- paste0("s = ",AF_dd_treat$selection)
AF_heat <- ggplot(AF_dd_treat,aes(x=factor(pop2),y=factor(migration),fill=mean_AFD))+
  facet_wrap(~selection_lab,nrow=1)+
  theme_minimal()+
  geom_tile()+
  scale_fill_viridis(option = "A",limits = c(0, 1))+
  labs(y="Migration",x="Daughter Population Size",fill="Mean AFD")+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=14),
        axis.text.x = element_text(size=14,angle=45,hjust=1),
        strip.text = element_text(size=14),
        strip.background = element_rect(colour='black'),
        legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        axis.title = element_text(size=14))
AF_heat

# Combine with the FNR rates for each of the six test metrics to demonstrate relationship with AF
FNR_rates <- plot_rates[plot_rates$type=="FNR" &
                          plot_rates$alpha==0.99,]
plot_AF <- NULL
for(i in 1:27){
  tmp_FNR <- FNR_rates[FNR_rates$selection == AF_dd_treat$selection[i] &
                         FNR_rates$pop2 == AF_dd_treat$pop2[i] &
                         FNR_rates$migration == AF_dd_treat$migration[i],]
  tmp_FNR$mean_AFD <- AF_dd_treat$mean_AFD[i]
  tmp_FNR$var_AFD <- AF_dd_treat$var_AFD[i]
  tmp_FNR$min_AFD <- AF_dd_treat$min_AFD[i]
  
  plot_AF <- rbind(plot_AF,tmp_FNR)
}

# And plot...
AF_FNR <- ggplot(plot_AF,aes(x=mean_AFD,y=rate,colour=`Window Size/Overlap`,shape=selection))+
  geom_point(size=2)+
  geom_line(aes(group=`Window Size/Overlap`))+
  facet_grid(migration~pop2,scales = "free")+
  labs(x="Mean AFD between founding\npop and all daughters",y="FNR",
       shape="Selection Coefficient")+
  theme_bw()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        axis.text.x = element_text(size=12,angle=45,hjust=1),
        strip.text = element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))

# Plot these together for a supp fig
pdf("figs/FigureSX_FNR_vs_AFD.pdf",width=8,height=8)
plot_grid(AF_heat,
          AF_FNR,
          rel_heights = c(1,2),ncol=1,nrow=2,align="v",axis='tblr',
          labels = "AUTO",label_size=32)
dev.off()

##### Visualise representative null for figure #####
null_fig <- null_distribution_visualise(input_dir = "outputs/slim/20,12,07_Variable_recombination_along_chromosome_15/final_results/",rec_rate = "med",window_size = 200,n_cores = 6)
pdf("figs/eigen_null_distributions.pdf",width=4,height=2.5)
null_fig
dev.off()

#### Two Parallel Results #####
sim_results2 <- "22,04,20_Rerunning_two_parallel_with_new_v2_null"
# We want to summarise false positive rate and false negative rate for eigen analyses first
for(i in 1:3){
  eigen_snp_window <- eigen_snps[i]
  print(paste0("STARTING TWO PARALLEL FOR SNP SIZE:",eigen_snp_window))
  if(!(file.exists(paste0("outputs/slim/",sim_results2,"_all_eigen_window",eigen_snp_window,".rds")))){
    all_fpr_fnr <- pbmclapply(sim_perms,function(sim){
      # tmp_dir <- paste0("outputs/slim/",sim_results2,"_",sim,"/final_results")
      tmp_dir <- paste0("/Volumes/jimwhiting_external/afvaper_sims/",sim_results2,"_",sim,"/final_results")
      eigen_fpr_fnr(tmp_dir,eigen_snp_window,2)
    },mc.cores=6)
    saveRDS(all_fpr_fnr,
            paste0("outputs/slim/",sim_results2,"_all_eigen_window",eigen_snp_window,".rds"))
  }
}

# Fetch all the eigen results for different window sizes and merge
all_fpr_fnr <- data.frame(rbindlist(lapply(1:length(eigen_snps),function(x){
  
  print(x)
  # Read
  tmp <- readRDS(paste0("outputs/slim/",sim_results2,"_all_eigen_window",eigen_snps[x],".rds"))
  
  # rbind and add metadata
  tmp2 <- data.frame(rbindlist(lapply(1:length(tmp),function(y){
    tmp_df <- data.frame(tmp[[y]])
    tmp_df <- reshape2::melt(tmp_df)
    colnames(tmp_df) <- c("alpha","rate")
    tmp_df$type <- rep(c("FPR","FNR"),nrow(tmp_df)/2)
    tmp_df$selection <- metadata[y,1]
    tmp_df$pop2 <- metadata[y,2]
    tmp_df$migration <- metadata[y,3]
    return(tmp_df)
  })))
  
  tmp2$window_size <- eigen_snps[x]
  return(tmp2)
})))

# What explains what...
# Selection
ggplot(all_fpr_fnr,aes(x=rate))+
  geom_histogram()+
  facet_grid(selection~type,scales = "free_y")

# Pop2
ggplot(all_fpr_fnr,aes(x=rate))+
  geom_histogram()+
  facet_grid(pop2~type,scales = "free_y")

# Migration
ggplot(all_fpr_fnr,aes(x=rate))+
  geom_histogram()+
  facet_grid(migration~type,scales = "free_y")

################################################################
# We want to summarise false positive rate and false negative rate for Fst next
if(!(file.exists(paste0("outputs/slim/",sim_results2,"_all_fst_results.rds")))){
  all_fpr_fnr_fst <- pbmclapply(sim_perms,function(sim){
    # tmp_dir <- paste0("outputs/slim/",sim_results2,"_",sim,"/final_results")
    tmp_dir <- paste0("/Volumes/jimwhiting_external/afvaper_sims/",sim_results2,"_",sim,"/final_results")
    fst_fpr_fnr(tmp_dir)
  },mc.cores=6)
  saveRDS(all_fpr_fnr_fst,
          paste0("outputs/slim/",sim_results2,"_all_fst_results.rds"))
}
all_paramater_fst_res <- readRDS(paste0("outputs/slim/",sim_results2,"_all_fst_results.rds"))

# Merge to a comparable data.frame
all_fst_error <- data.frame(rbindlist(lapply(1:length(all_paramater_fst_res),function(x){
  
  # Get the FPR results
  fpr_res <- reshape2::melt(all_paramater_fst_res[[x]]$FPR$avg)
  colnames(fpr_res) <- c("overlap","alpha","rate")
  fpr_res$type <- "FPR"
  
  # Get the FNR results
  fnr_res <- reshape2::melt(all_paramater_fst_res[[x]]$FNR)
  colnames(fnr_res) <- c("overlap","alpha","rate")
  fnr_res$type <- "FNR"
  
  # Merge
  fpr_fnr <- rbind(fpr_res,fnr_res)
  
  # Add parameters
  fpr_fnr$selection <- metadata[x,1]
  fpr_fnr$pop2 <- metadata[x,2]
  fpr_fnr$migration <- metadata[x,3]
  
  return(fpr_fnr)
})))

ggplot(all_fst_error,aes(x=rate))+
  geom_histogram()+
  facet_grid(selection~type,scales = "free")

# Pop2
ggplot(all_fst_error,aes(x=rate))+
  geom_histogram()+
  facet_grid(pop2~type,scales = "free_y")

# Migration
ggplot(all_fst_error,aes(x=rate))+
  geom_histogram()+
  facet_grid(migration~type,scales = "free_y")

# Compare the two data structures
head(all_fpr_fnr)
head(all_fst_error)

# Tidy up
all_fpr_fnr$alpha <- as.numeric(gsub("X","",all_fpr_fnr$alpha))

# Now merge these together
colnames(all_fpr_fnr)[colnames(all_fpr_fnr)=="window_size"] <- "Window Size/Overlap"
colnames(all_fst_error)[colnames(all_fst_error)=="overlap"] <- "Window Size/Overlap"

# And in the darkness, bind them
plot_rates <- rbind(all_fpr_fnr,all_fst_error)
plot_rates$`Window Size/Overlap` <- factor(plot_rates$`Window Size/Overlap`,levels=c(50,200,500,4,3,2))
plot_rates$selection <- factor(plot_rates$selection,levels=c(0.01,0.05,0.1))
plot_rates$pop2 <- factor(plot_rates$pop2,levels=c(400,2000,10000))
plot_rates$migration <- factor(plot_rates$migration,levels=c(0,0.0025,0.01))

# For each of the parameters, we'll make a 4x2 grid showing results and mean/se
parameters <- c("selection","pop2","migration")
parameter_labs <- c("Selection Coefficient","Daughter Pop Size","Migration Rate")
names(parameter_labs) <- parameters

param_plots <- lapply(parameters,function(param){
  
  # Make new dummy data
  tmp <- plot_rates[,c("alpha","rate","type","Window Size/Overlap",param)]
  colnames(tmp)[ncol(tmp)] <- "Var"
  
  # Get summaries
  sum_plot <- data.frame(tmp %>% 
                           group_by(`Window Size/Overlap`,alpha,type,Var) %>%
                           summarise(mean_rate=mean(rate),
                                     se=sd(rate)/sqrt(length(rate))))
  
  # Add AF-VapeR vs FST flags
  sum_plot$analysis <- "FST"
  sum_plot[sum_plot$Window.Size.Overlap %in% c(50,200,500),"analysis"] <- "AF-vapeR"
  
  # Make param factor again...
  if(class(sum_plot$Var) != "factor"){
    sum_plot$Var <- factor(sum_plot$Var,levels = sort(as.numeric(unique(sum_plot$Var))))
  }
  
  # And plot...
  ggplot(sum_plot[sum_plot$alpha != 0.9,],aes(x=Var,y=mean_rate,colour=Window.Size.Overlap,group=Window.Size.Overlap,shape=analysis))+
    geom_point(position=position_dodge(width = 0.5))+
    geom_linerange(aes(ymin = mean_rate-se, ymax = mean_rate+se),position=position_dodge(width = 0.5),show.legend = F)+
    facet_grid(type~alpha,scales="free_y")+
    labs(y="Error Rate",x=parameter_labs[param],colour="Window Size/Overlap",shape="Analysis")+
    ylim(0,1)
  
})

# Combine
legend <- get_legend(
  # create some space to the left of the legend
  param_plots[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)
pdf("figs/FigureSX_simulation_effects_fprfnr_TWO_PARALLEL.pdf",width=8,height=8)
plot_grid(
  plot_grid(param_plots[[1]] + theme(legend.position = "none"),
            param_plots[[2]] + theme(legend.position = "none"),
            param_plots[[3]] + theme(legend.position = "none"),labels="AUTO",label_size = 32,ncol=1,align='v',axis="tblr",hjust = 0.1),
  legend,ncol=2,rel_widths = c(3,1))
dev.off()

#### Plot FNR heatmap for each analysis subset ####
analyses <- unique(plot_rates$`Window Size/Overlap`)

# FNR + FPR Matrices for alpha of 0.99
error_heats <- lapply(c("FNR","FPR"),function(false_type){
  # Get our tmp data
  tmp <- plot_rates[plot_rates$alpha == 0.99 &
                      plot_rates$type == false_type,]
  
  # Make a selection label
  tmp$selection_lab <- paste0("s = ",tmp$selection)
  
  # Make heatmap...
  ggplot(tmp,aes(x=factor(pop2),y=factor(migration),fill=rate))+
    facet_grid(`Window Size/Overlap`~selection_lab)+
    theme_minimal()+
    geom_tile()+
    scale_fill_viridis(option = "A",limits = c(0, 1))+
    labs(y="Migration",x="Daughter Population Size",fill=false_type)
})

# Now plot over different alphas...
# FNR + FPR Matrices for alpha of 0.99
alphas <- c(0.95,0.99,0.999)
for(error_rate in c("FPR","FNR")){
  alpha_heats <- lapply(alphas,function(alpha){
    # Get our tmp data
    tmp <- plot_rates[plot_rates$alpha == alpha &
                        plot_rates$type == error_rate,]
    
    # Make a selection label
    tmp$selection_lab <- paste0("s:",tmp$selection)
    
    # Make new facet labels
    labs <- c(50,200,500,4,3,2)
    new_labs <- c("AFvapeR:50","AFvapeR:200","AFvapeR:500",
                  "F[ST]:4","F[ST]:3","F[ST]:2")
    tmp$`Window Size/Overlap` <- as.character(tmp$`Window Size/Overlap`)
    for(i in 1:length(labs)){
      tmp[tmp$`Window Size/Overlap`== labs[i],"Window Size/Overlap"] <- new_labs[i]
    }
    tmp$`Window Size/Overlap` <- factor( tmp$`Window Size/Overlap`,levels=new_labs)
    
    # Make heatmap...
    p1 <- ggplot(tmp,aes(x=factor(pop2),y=factor(migration),fill=rate))+
      facet_grid(`Window Size/Overlap`~selection_lab,labeller = label_parsed)+
      theme_minimal()+
      theme(panel.grid = element_blank(),
            strip.text = element_text(size = 11),
            axis.text = element_text(size=12),
            axis.title = element_text(size=14),
            strip.background = element_rect(colour="black"))+
      geom_tile()+
      scale_fill_viridis(option = "A")+
      labs(y="Migration",x="Daughter Population Size",fill=error_rate)
    
    if(error_rate == "FPR"){
      p1 +  scale_fill_viridis(option = "A",limits = c(0, 0.2))
    } else {
      p1 +  scale_fill_viridis(option = "A",limits = c(0, 1))
    }
    
  })
  
  # Also get legend
  # Combine
  legend <- get_legend(
    # create some space to the left of the legend
    alpha_heats[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12),
                             legend.title = element_text(size=18),
                             legend.text = element_text(size=16))
  )
  
  pdf(paste0("figs/FigureX_",error_rate,"_heats_all_alpha_TWO_PARALLEL.pdf"),width=14,height=8)
  print(cowplot::plot_grid(alpha_heats[[1]]+theme(legend.position = "none",
                                                  axis.text.x = element_text(angle=45,hjust=1))+
                             ggtitle(expression(alpha:0.05)),
                           alpha_heats[[2]]+theme(legend.position = "none",
                                                  axis.text.x = element_text(angle=45,hjust=1))+
                             ggtitle(expression(alpha:0.01)),
                           alpha_heats[[3]]+theme(legend.position = "none",
                                                  axis.text.x = element_text(angle=45,hjust=1))+
                             ggtitle(expression(alpha:0.001)),
                           legend,
                           ncol=4,labels=c("A","B","C",""),label_size = 32,rel_widths = c(2,2,2,0.7)))
  dev.off()
}


