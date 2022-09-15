source(file = "./functions_when_not.R")

option_list <- list(
  optparse::make_option("--r", default = 1,
                        help = "different row [default %default].")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

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
lor_nac_2_all <- nma_old$TE.fixed[8, ]
risk_all <- lor2prob(p_nac,lor_nac_2_all)
risk_all <- as.data.frame(risk_all)
risk_all$trt <- rownames(risk_all)
rownames(risk_all) <- rep(1:nrow(risk_all))
colnames(risk_all)[1] <- "p"

risk_all[risk_all$trt=="Trimethoprim","p"] <- risk_all[risk_all$trt=="Ceftiofur hydrochloride","p"]

nrep <- 100000

# need to change every time; from 1 to 4
r=opt$r
n <- c(50,100,150,200)
n_r <- n[r]
dat_para <- as.data.frame(n_r)
############ with prev info##########################
pig_alloc_c <- rep(n_r/2,2)

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
#BRD$Number.of.Event.in.arm.2[BRD$Study.number == 2] <- BRD$Number.of.Event.in.arm.2[BRD$Study.number == 2] - 0.5

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
  
  # lor_hat <- nma_sim_prev$TE.fixed[10,2]
  # se <- nma_sim_prev$seTE.fixed[10,2]
  # critical_vale <- abs(lor_hat)/se
  # pvalue <- 2*pnorm(-critical_vale)
  # if the simulate previous network shows the p-value is less than 0.1, then we conduct the new study, else not;
  return(c(nma_sim_prev$pval.fixed[12,1],
           bio_equal(p, p, sigma = 0, re_sigma = 0, pig_alloc_c, data_prev = BRD_new_long)))
}

BE <- table_with_prev[table_with_prev[,1]<0.1,]
BE <- BE[BE[,1]>0.05,] # >0.05 and <0.1
BE <- BE[,-1]

dat_para[1,2:3] <- pig_alloc_c
dat_para[1,4:6] <- apply(BE, 2, mean)
lor_hat_all <- BE[,2]
dat_para[1,6] <- sd(lor_hat_all)
dat_para[1,7] <- nrow(BE)

write.csv(dat_para, file = paste0("./res_error_prev_dep_row_",r,".csv"), row.names = F)

