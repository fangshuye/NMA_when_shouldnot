rm(list = ls())
library(ldbounds)
#setwd("/Users/fangshu/OneDrive - Iowa State University/Research/Project_II_When_not/Simulation/")
source(file = "./functions_when_not.R")
option_list <- list(
  optparse::make_option("--r", default = 1,
                        help = "different row [default %default].")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

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
  lor_hat <- nma_res$TE.fixed[12,1]
  se <- nma_res$seTE.fixed[12,1]
  
  abs(lor_hat)/se
}


# setwd("/Users/fangshu/OneDrive - Iowa State University/Research/Project_II_When_not/")

BRD <- read.csv("./data/updated_dataset.csv", stringsAsFactors = F)
BRD_r <- BRD 
BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] <- BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] - 0.5
BRD_long <- wide2long(BRD_r)
BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_pair,sm="OR",comb.fixed = T,comb.random = F)

p_nac <- BRD_long %>% filter(t == "No active control") %>% summarise(p = sum(r)/sum(n)) %>% pull
#lor_nac_2_TILD <- nma_old$TE.fixed[8, 10]
#p <- lor2prob(p_nac,lor_nac_2_TILD)

lor_nac_2_CEFTH <- nma_old$TE.fixed[8, 1]
p <- lor2prob(p_nac,lor_nac_2_CEFTH)

# get the risk table
trt1 <- 'Trimethoprim'
trt2 <- 'Ceftiofur hydrochloride'
lor_nac_2_all <- nma_old$TE.fixed[8, ]
risk_all <- lor2prob(p_nac,lor_nac_2_all)
risk_all <- as.data.frame(risk_all)
risk_all$trt <- rownames(risk_all)
rownames(risk_all) <- rep(1:nrow(risk_all))
colnames(risk_all)[1] <- "p"

risk_all[risk_all$trt==trt1,"p"] <- risk_all[risk_all$trt==trt2,"p"]

# need to change every time; from 1 to 4
r=opt$r
n <- c(50,100,150,200) # total sample size
n_r <- n[r]
sample_size_new <- rep(n_r/2,2)
# alpha-spending function
alpha_function <- function(IF, alpha = 0.05){
  quantile <- qnorm(1-alpha/2)/sqrt(IF)
  2-2*pnorm(quantile)
}

BRD_subset <- BRD_r %>%
  filter(Arm.1 == trt1 | 
           Arm.1 == trt2 |
           Arm.2 == trt1 | 
           Arm.2 == trt2 )

old_n_trt1 <- 0
old_n_trt2 <- 0
for (i in 1:nrow(BRD_subset)) {
  if(BRD_subset$Arm.1[i] == trt1){
    old_n_trt1 <- old_n_trt1 + BRD_subset$Total.number.in.arm.1[i]
  }else if(BRD_subset$Arm.1[i] == trt2){
    old_n_trt2 <- old_n_trt2 + BRD_subset$Total.number.in.arm.1[i]
  }
  
  if(BRD_subset$Arm.2[i] == trt1){
    old_n_trt1 <- old_n_trt1 + BRD_subset$Total.number.in.arm.2[i]
  }else if(BRD_subset$Arm.2[i] == trt2){
    old_n_trt2 <- old_n_trt2 + BRD_subset$Total.number.in.arm.2[i]
  }
}

old_n <- old_n_trt1 + old_n_trt2
IF <-  old_n + n_r
t <- c(old_n/IF,1)
obf.bd <- bounds(t,iuse=c(1,1),alpha=c(0.025,0.025))
thresholds <- summary(obf.bd)$bounds[,3]
c1 <- thresholds[1]
c2 <- thresholds[2]
alpha1 <- alpha_function(old_n/IF)
alpha2 <- alpha_function(1) - alpha1


BRD_extend <- BRD
risk_1 <- risk_all
colnames(risk_1) <- c("p1","Arm.1")
risk_2 <- risk_all
colnames(risk_2) <- c("p2","Arm.2")
risk_3<- risk_all
colnames(risk_3) <- c("p3","Arm.3")
BRD_extend <- merge(BRD_extend,risk_1,by="Arm.1",all.x  = T)
BRD_extend <- merge(BRD_extend,risk_2,by="Arm.2",all.x = T)
BRD_extend <- merge(BRD_extend,risk_3,by="Arm.3",all.x = T)

nrep <- 100000
set.seed(20211031, kind = "L'Ecuyer-CMRG")
registerDoParallel(cores = 16)
table_with_prev <- foreach (rep=1:nrep, .combine = rbind)%dopar%{
  
  ### simulation the whole network###
  BRD_extend$Number.of.Event.in.arm.1=rbinom(nrow(BRD_extend), 
                                             size =BRD_extend$Total.number.in.arm.1, 
                                             prob = BRD_extend$p1)
  BRD_extend$Number.of.Event.in.arm.2=rbinom(nrow(BRD_extend), 
                                             size =BRD_extend$Total.number.in.arm.2, 
                                             prob = BRD_extend$p2)
  BRD_extend$Number.of.Event.in.arm.3=rbinom(nrow(BRD_extend), 
                                             size =BRD_extend$Total.number.in.arm.3, 
                                             prob = BRD_extend$p3)
  
  BRD_new_long <- wide2long(BRD_extend) 
  BRD_new_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_new_long, allstudies = T, sm = "OR")
  nma_sim_prev <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new_pair,sm="OR",comb.fixed = T,comb.random = F)
  
  lor_hat <- nma_sim_prev$TE.fixed[12,1]
  se <- nma_sim_prev$seTE.fixed[12,1]
  p_value_1 <- nma_sim_prev$pval.fixed[1,12]
    
  z1 <- abs(lor_hat)/se
  z2 <- bio_equal(p, p, sigma = 0, re_sigma = 0, sample_size_new, data_prev = BRD_new_long)
  
  sig_tsa <- ifelse(z1>c1 | z2>c2,1,0)
  return(c(p_value_1,sig_tsa))
}

# only conduct the new study when result is promising 
BE_1 <- table_with_prev[table_with_prev[,1]<0.1,] # < 0.1
BE_2 <- BE_1[BE_1[,1]>0.05,] # >0.05 and <0.1

BE_1 <- BE_1[,-1]
BE_2 <- BE_2[,-1]
BE_3 <- table_with_prev[,-1]

error_rate_1 <- mean(BE_1) # tsa, < 0.1
error_rate_2 <- mean(BE_2) # tsa, < 0.1, >0.05
error_rate_3 <- mean(BE_3) # tsa, indep

# output 3 values for each total sample size n and summarize the alpha-spending
sink(paste0("TSA_alpha_",r,'.txt'), append = TRUE)
cat("\n", "simulation result \n")
cat( "=====================================================","\n")
cat("simulation time: ", nrep, "\n")
cat("trt1 = ", trt1, "\n")
cat("trt2 = ", trt2, "\n")
cat("IF =  ", t, "\n")
cat("alpha1 = ", alpha1, "\n")
cat("alpha2 = ", alpha2, "\n")
cat("c1 = ", c1, "\n")
cat("c2 = ", c2, "\n")
cat( "=====================================================","\n")
cat( "When p-value < 0.1 is the criterion for promising result","\n")
cat("Total number of new study:", length(BE_1), "\n")
cat("Type I error when analyzing with TSA:", error_rate_1, "\n")
cat( "=====================================================","\n")
cat( "When 0.05 < p-value < 0.1 is the criterion for promising result","\n")
cat("Total number of new study:", length(BE_2), "\n")
cat("Type I error when analyzing with TSA:", error_rate_2, "\n")
cat( "=====================================================","\n")
cat( "When the new trial will be conducted regardless of the previous result","\n")
cat("Total number of new study:", length(BE_3), "\n")
cat("Type I error when analyzing with TSA:", error_rate_3, "\n")
cat( "=====================================================","\n")
sink()


