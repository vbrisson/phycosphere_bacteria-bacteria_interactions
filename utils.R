## Modules used frequently throughout
get_df <- function(fileloc = 'data/flow2_bact.csv'){
  df <- read.csv(file=fileloc)
  # df <- df[df$Time!=0,]
  df$Treatment[df$Treatment==1] <- 'Devosia'
  df$Treatment[df$Treatment==2] <- 'none'
  df$Treatment[df$Treatment==3] <- 'Alcanivorax'
  df$Treatment[df$Treatment==4] <- 'Marinobacter'
  
  return(df)
}


get_df_alg <- function(fileloc = 'data/flow2_pt.csv'){
  df_alg <- read.csv(file='data/flow2_pt.csv') # algae
  df_alg$Treatment[df_alg$Treatment==1] <- 'Devosia'
  df_alg$Treatment[df_alg$Treatment==2] <- 'none'
  df_alg$Treatment[df_alg$Treatment==3] <- 'Alcanivorax'
  df_alg$Treatment[df_alg$Treatment==4] <- 'Marinobacter'
  
  return(df_alg)
}


get_df_fold <- function(df, df_alg){  # FOLD INCREASE day 5 -- 14
  df_fold <- data.frame(matrix(nrow=0, ncol=7))
  colnames(df_fold) <- c('Strain', 'Treatment', 'Microplate', 'Ring', 'Direction',
                         'fold', 'fold_adj')
  conds <- c('Treatment', 'Microplate', 'Ring', 'Direction')
  df_cond <- df[,conds]
  
  for(i in 1:nrow(unique(df_cond))){
    ind <- (df$Treatment==df_cond[i,1]) & (df$Microplate==df_cond[i,2]) & 
      (df$Ring==df_cond[i,3]) & (df$Direction==df_cond[i,4])
    df_ext <- df[ind,]
    ind2 <- (df_alg$Treatment==df_cond[i,1]) & (df_alg$Microplate==df_cond[i,2])
    df_alg_ext <- df_alg[ind2,]
    if (length(df_ext[,1])>1){
      val <- df_ext[2,8] / df_ext[1,8]
      val2 <- df_alg_ext[2,5] / df_alg_ext[1,5]
      df_app <- data.frame(Strain=df_ext[1,2], 
                           Treatment=df_ext[1,4], 
                           Microplate=df_ext[1,5], 
                           Ring=df_ext[1,6], 
                           Direction=df_ext[1,7], 
                           fold=val,
                           fold_adj=val/val2)
      df_fold[nrow(df_fold)+1, ] = df_app
    }
  }
  
  df_fold$Strain <- factor(df_fold$Strain, 
                           levels=c('Devosia EAB7WZ', 'Alcanivorax EA2', 
                                    'None', 'Marinobacter 3-2'))
  df_fold$fold_log <- log(df_fold$fold)  # natural log
  df_fold$fold_adj_log <- log(df_fold$fold_adj)
  df_fold <- df_fold[(df_fold$Ring!=1 | df_fold$Treatment!='none'),]
  df_fold <- df_fold[-nrow(df_fold),]
  df_fold$rate <- df_fold$fold_log / 9  # bacteria growth rate
  
  return(df_fold)
}


get_df14 <- function(df){
  ## SUBSET BY TIMEPOINT
  df14 <- df[df$Time==14,]
  df_alg14 <- df_alg[df_alg$Time==14,]
  df14$Abd_per_alga <- matrix(0, nrow=nrow(df14))
  conds <- c('Treatment', 'Microplate')
  df14_cond <- df14[,conds]
  df14_cond_unique <- unique(df14_cond)
  
  for(i in 1:nrow(df14_cond_unique)){
    ind <- (df14$Treatment==df14_cond_unique[i,1]) & 
      (df14$Microplate==df14_cond_unique[i,2])
    
    ind2 <- (df_alg14$Treatment==df14_cond_unique[i,1]) & 
      (df_alg14$Microplate==df14_cond_unique[i,2])
    
    df14[ind, 9] <- df14[ind, 8] / df_alg14[ind2, 5]
  }
  df14$Strain <- factor(df14$Strain, levels=c('Alcanivorax EA2', 'Devosia EAB7WZ',
                                              'Marinobacter 3-2', 'None'))
  df14$Treatment <- factor(df14$Treatment, levels=c(3,1,4,2))
  df14_se <- summarySE(df14, measurevar = 'Abd_per_alga', 
                       groupvars = c('Treatment', 'Ring'))
  
  # export day14
  write.csv(df14, "data/flow2_day14_bact.csv", row.names=FALSE)
  write.csv(df_alg14, "data/flow2_day14_pt.csv", row.names=FALSE)
  
  return(c(df14, df_alg14))
}


get_df5 <- function(df){
  ## subset day 5
  df5 <- df[df$Time==5,]
  df5$Strain <- factor(df5$Strain,
                       levels=c('Alcanivorax EA2', 'Devosia EAB7WZ',
                                'Marinobacter 3-2', 'None'))
  df5$Treatment <- factor(df5$Treatment, levels=c(3,1,4,2))
  
  return(df5)
}


append_xnet <- function(xnet, xnet_info){
  xnet_app <- xnet
  xnet_app$microplate <- NaN
  xnet_app$ring <- NaN
  xnet_app$strain <- NaN
  xnet_app$treatment <- NaN
  list_sample_name <- unique(xnet$sample_name)
  for (s in list_sample_name){
    mi <- xnet_info$microplate[xnet_info$sample_name==s]
    ri <- xnet_info$ring[xnet_info$sample_name==s]
    st <- xnet_info$strain[xnet_info$sample_name==s]
    tr <- xnet_info$treatment[xnet_info$sample_name==s]
    
    xnet_app$microplate[xnet_app$sample_name==s] <- mi
    xnet_app$ring[xnet_app$sample_name==s] <- ri
    xnet_app$strain[xnet_app$sample_name==s] <- st
    xnet_app$treatment[xnet_app$sample_name==s] <- tr
  }
  xnet_app <- xnet_app[(xnet_app$treatment!='none' | xnet_app$ring!='inner'),]
  return(xnet_app)
}


## merge Cnet and bacterial count for plotting fig 4e
summarize_cnet <- function(cnet, cnet_info){
  # summarize cnet
  cnet$microplate <- NA
  cnet$ring <- NA
  cnet$treatment <- NA
  cnet$strain <- NA
  for (s in unique(cnet_info$sample_name)){
    m <- cnet_info$microplate[cnet_info$sample_name==s]
    r <- cnet_info$ring[cnet_info$sample_name==s]
    t <- cnet_info$treatment[cnet_info$sample_name==s]
    st <- cnet_info$strain[cnet_info$sample_name==s]
    
    cnet$microplate[cnet$sample_name==s] <- m
    cnet$ring[cnet$sample_name==s] <- r
    cnet$treatment[cnet$sample_name==s] <- t
    cnet$strain[cnet$sample_name==s] <- st
  }
  
  cnet_summ <- cnet %>%
    group_by(treatment, microplate, ring) %>%
    summarize(cnet_q25 = quantile(Cnet, probs = 0.25),
              cnet_q50 = quantile(Cnet, probs = 0.5),
              cnet_q75 = quantile(Cnet, probs = 0.75),
              cnet_n = n()
              # cnet_mean = mean(Cnet),
              # cnet_sd = sd(Cnet),
              # replaced because we look for total incorp rate
    )
  return(list('df' = cnet, 'df_summ' = cnet_summ))
}


merge_count_cnet <- function(count_bact_stat, cnet_stat){
  count_bact_stat$Ring[count_bact_stat$Ring==1] <- 'inner'
  count_bact_stat$Ring[count_bact_stat$Ring==2] <- 'outer'
  df <- cnet_stat
  colnames(df) <- c("treatment", "microplate", "ring", "cnet_q25", "cnet_q50", 
                    "cnet_q75", "n_cnet")
  df$count_mean <- NA
  df$count_sd <- NA
  
  for (row in 1:nrow(df)){
    t <- df$treatment[row]
    r <- df$ring[row]
    m <- df$microplate[row]
    df[row, 'count_mean'] <- 
      count_bact_stat$mean[count_bact_stat$Treatment==t 
                           & count_bact_stat$Ring==r
                           & count_bact_stat$Microplate==m]
    df[row, 'count_sd'] <- 
      count_bact_stat$sd[count_bact_stat$Treatment==t 
                         & count_bact_stat$Ring==r
                         & count_bact_stat$Microplate==m]
  }
  
  df <- df[(df$treatment!='none' | df$ring!='inner'),]
  
  return(df)
}