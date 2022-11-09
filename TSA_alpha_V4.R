rm(list = ls())

packages <- c("dplyr", "netmeta", "parallel", "foreach", "doParallel", "doRNG","DEoptimR")
lapply(packages, library, character.only = TRUE)
#setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/data/V2')
BRD <- read.csv("./pairwise_dataset.csv", 
                stringsAsFactors = F)
wide_long_pair <- function(MTCdata){
  N <- max(MTCdata$Number.of.arms)
  study <- rep(MTCdata$Study.number,N)
  t <- NULL
  n <- NULL
  r <- NULL
  for(i in c(2,1)){
    r <- c(r, eval(parse(text = paste0("MTCdata$Number.of.Event.in.arm.",i, sep = ""))))
    n <- c(n, eval(parse(text = paste0("MTCdata$Total.number.in.arm.",i, sep = ""))))
    t <- c(t, eval(parse(text = paste0("MTCdata$Arm.",i, sep = ""))))
  }
  res <- data.frame(id = study, t = t, r = r, n = n)
  res <- res %>% dplyr::filter(!is.na(n)) %>% arrange(id)
  res
}

bio_equal <- function(p_baseline, p_trt2, sigma, re_sigma = 0, pig_alloc, data_prev){
  #res_TE <- numeric(1)
  #res_seTE <- numeric(1)
  #power <- numeric(length(p_trt2_vector))
  new_s <- new_study(p_baseline, p_trt2, sigma = sigma, pig_alloc = pig_alloc)
  data_final <- rbind(data_prev, new_s)
  BRD_new <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = data_final, allstudies = T, sm = "OR")
  nma_res <- netmeta::netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new,sm="OR",
                              comb.fixed =F,comb.random = T, tol.multiarm.se = 0.01)
  #z <- true_lor/se
  #power <- pnorm(z-qnorm(0.975))+pnorm(-z-qnorm(0.975))
  lor_hat <- nma_res$TE.fixed[1,2]
  se <- nma_res$seTE.fixed[1,2]
  p_value <- nma_res$pval.fixed[1,2]
  
  return(c(p_value,abs(lor_hat)/se))
}

bioeq_single_s <- function(p_baseline, p_trt2, pig_alloc){
  #power <- numeric(1)
  
  # random part in the simulation, rbinom
  BRD_s <- new_study(p_baseline, p_trt2, sigma = 0, pig_alloc = pig_alloc)
  BRD_new_s <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_s, allstudies = T, sm = "OR")
  # TE: Estimate of treatment effect (log odds ratio, mean difference)
  # seTE: S.E. of TE
  # sm: summary measure,
  # comb.fixed: whether a fixed effects (common effects) network meta-analysis should be conducted.
  # comb.random: whether a random effects network meta-analysis should be conducted.
  nma_s <- netmeta::netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new_s,sm="OR",
                            comb.fixed =F,comb.random = T, tol.multiarm.se = 0.01)
  #z <- true_lor/se
  #power <- pnorm(z-qnorm(0.975))+pnorm(-z-qnorm(0.975))
  # lor_TILD_2_CEFTP
  p_value <- nma_s$pval.fixed[1,2]
  p_value
}

# the type of the new_s is the same as what wide2long() generate above.
new_study <- function(p_baseline, p_trt2, sigma, pig_alloc = c(50,50)){
  r_vec <- rbinom(2, size = pig_alloc, prob = c(p_baseline, p_trt2))
  new_s <- data.frame(id = rep(100,2), t = c("Florfenicol", "Gamithromycin"),
                      r = r_vec, n = pig_alloc)
  return(new_s)
}

c1 <- 2.16
c2 <- 2.04

nrep <- 100000
p <- 0.29
sample_size_new <- c(51,52)

set.seed(20211031, kind = "L'Ecuyer-CMRG")
registerDoParallel(cores = 16)
table_with_prev <- foreach (rep=1:nrep, .combine = rbind)%dopar%{
  
  ### simulation the old study###
  BRD$Number.of.Event.in.arm.1 <- rbinom(nrow(BRD), 
                                         size = BRD$Total.number.in.arm.1, 
                                         prob = p)
  BRD$Number.of.Event.in.arm.2 <- rbinom(nrow(BRD), 
                                         size = BRD$Total.number.in.arm.2, 
                                         prob = p)
  
  BRD_long <- wide_long_pair(BRD)
  BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
  nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_pair,sm="OR",comb.fixed = T,comb.random = F)
  #lor_hat <- nma_old$TE.fixed[1,2]
  #se <- nma_old$seTE.fixed[1,2]
  p_value_1 <- nma_old$pval.fixed[1,2]
    
  #z1 <- abs(lor_hat)/se
  
  #res_2 <- bio_equal(p, p, sigma = 0, re_sigma = 0, sample_size_new, data_prev = BRD_long)
  p_value_2 <- bioeq_single_s(p,p, sample_size_new)
  #z2 <- res_2[2]
  
  sig_single <- ifelse(p_value_2 < 0.05,1,0)
  #sig <- ifelse(z1>c1 | z2>c2,1,0)
  return(c(p_value_1,sig_single))
}

# only conduct the new study when result is promising 
BE_1 <- table_with_prev[table_with_prev[,1]<0.1,] # < 0.1
BE_2 <- BE_1[BE_1[,1]>0.05,] # >0.05 and <0.1

BE_1 <- BE_1[,-1]
BE_2 <- BE_2[,-1]
BE_3 <- table_with_prev[,-1]

error_rate_1 <- mean(BE_1) # single, < 0.1
error_rate_2 <- mean(BE_2) # single, < 0.1, >0.05
error_rate_3 <- mean(BE_3) # single, indep

sink(paste0("TSA_alpha",'.txt'), append = TRUE)
cat("\n", "simulation result \n")
cat( "=====================================================","\n")
cat("simulation time: ", nrep, "\n")
cat("c1 = ", c1, "\n")
cat("c2 = ", c2, "\n")
cat( "When p-value < 0.1 is the criterion for promising result","\n")
cat("Total number of new study:", length(BE_1), "\n")
cat("Type I error when analyzing in isolation:", error_rate_1, "\n")
cat( "When 0.05 < p-value < 0.1 is the criterion for promising result","\n")
cat("Total number of new study:", length(BE_2), "\n")
cat("Type I error when analyzing in isolation:", error_rate_2, "\n")
cat( "=====================================================","\n")
cat( "When the new trial will be conducted regardless of the previous result","\n")
cat("Total number of new study:", length(BE_3), "\n")
cat("Type I error when analyzing in isolation:", error_rate_3, "\n")
cat( "=====================================================","\n")
sink()


