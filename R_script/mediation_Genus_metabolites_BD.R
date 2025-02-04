#20240508 中介分析：加上用药史的协变量-菌+代谢物~BD
setwd("G:\\BD_Genues_metabolites_MRE\\BD_code\\mediation_1210\\20140114mediation")
library(mediation)
library(MASS)
#COV_fea <- c('Gender',	'Age',	'BMI', 'Smoking', 'Drinking')
COV_fea <- c('Gender',	'Age',	'BMI', 'Smoking', 'Drinking','Biologics','Antibiotics','Corticosteroids','Probiotics','Immunomodulator')

group_df <- read.csv('treat.csv',row.names = 1)
group_df <- group_df[group_df$Group_BD %in% c('BD2','BD1'),]

#LI_score <- read.csv('Limann.csv',row.names = 1)


otu      <- as.data.frame(t(read.csv('12_combined_lefse_wilcox_otu.csv',row.names = 1)))

Q300     <- as.data.frame((read.csv('select_Q300_18_0.05_P_1_VIP.csv')))
Q300      <- as.data.frame(t(Q300[seq(4,ncol(Q300))]))
colnames(Q300) <- Q300[1, ]# 将第一行作为列名
Q300 <- Q300[-1, ]

Q600     <- as.data.frame((read.csv('select_Q600_67_0.05_P_1_VIP.csv')))
Q600      <- as.data.frame(t(Q600[seq(4,ncol(Q600))]))
colnames(Q600) <- Q600[1, ]# 将第一行作为列名
Q600 <- Q600[-1, ]

mre      <- as.data.frame((read.csv('MRE_7_sig_features.csv',row.names = 1)))
mre      <- mre[seq(5,ncol(mre))]

#cov <- as.data.frame((read.csv('6_clinical_covariables.csv',row.names = 1)))[COV_fea]
cov <- as.data.frame((read.csv('10_COV_treat.csv',row.names = 1)))[COV_fea]
samples <- intersect(intersect(rownames(otu),rownames(group_df)), intersect(rownames(Q600),rownames(Q300)))

group <- group_df[,'Group_BD',drop=FALSE]

#LI_score <- LI_score[,'LI',drop=FALSE]

colnames(otu) <- gsub("/", "_", colnames(otu))
colnames(Q300) <- gsub("/", "_", colnames(Q300))
colnames(Q600) <- gsub("/", "_", colnames(Q600))
colnames(mre) <- gsub("/", "_", colnames(mre))

colnames(otu) <- gsub(" ", "_", colnames(otu))
colnames(Q300) <- gsub(" ", "_", colnames(Q300))
colnames(Q600) <- gsub(" ", "_", colnames(Q600))
colnames(mre) <- gsub(" ", "_", colnames(mre))

colnames(otu) <- gsub("-", "_", colnames(otu))
colnames(Q300) <- gsub("-", "_", colnames(Q300))
colnames(Q600) <- gsub("-", "_", colnames(Q600))
colnames(mre) <- gsub("-", "_", colnames(mre))

colnames(otu) <- gsub(":", "_", colnames(otu))
colnames(Q300) <- gsub(":", "_", colnames(Q300))
colnames(Q600) <- gsub(":", "_", colnames(Q600))
colnames(mre) <- gsub(":", "_", colnames(mre))

colnames(otu) <- gsub('\\(', "_", colnames(otu))
colnames(Q300) <- gsub("\\(", "_", colnames(Q300))
colnames(Q600) <- gsub("\\(", "_", colnames(Q600))
colnames(mre) <- gsub("\\(", "_", colnames(mre))

colnames(otu) <- gsub(")", "_", colnames(otu))
colnames(Q300) <- gsub(")", "_", colnames(Q300))
colnames(Q600) <- gsub(")", "_", colnames(Q600))
colnames(mre) <- gsub(")", "_", colnames(mre))

colnames(otu) <- gsub("\\[", "", colnames(otu))
colnames(Q300) <- gsub("\\[", "", colnames(Q300))
colnames(Q600) <- gsub("\\[", "", colnames(Q600))
colnames(mre) <- gsub("\\[", "", colnames(mre))

colnames(otu) <- gsub("]", "_", colnames(otu))
colnames(Q300) <- gsub("]", "_", colnames(Q300))
colnames(Q600) <- gsub("]", "_", colnames(Q600))
colnames(mre) <- gsub("]", "_", colnames(mre))


colnames_otu  <- lapply(colnames(otu),function(x){strsplit(x,"g__")[[1]][2]})
colnames_Q300 <- lapply(colnames(Q300),function(x){paste("Fecal_", x , sep = '')})
colnames_Q600 <- lapply(colnames(Q600),function(x){paste("Serum_", x , sep = '')})
colnames_mre  <- colnames(mre)

colnames(otu)  <- colnames_otu
colnames(Q300) <- colnames_Q300
colnames(Q600) <- colnames_Q600
colnames(mre)  <- colnames_mre


# 确保所有数据框的行按'samples'排序
otu <- otu[match(samples, rownames(otu)), ]
Q600 <- Q600[match(samples, rownames(Q600)), ]
Q300 <- Q300[match(samples, rownames(Q300)), ]
mre <- mre[match(samples, rownames(mre)), ]
cov <- cov[match(samples, rownames(cov)), ]
group <- group[match(samples, rownames(group)),,drop=FALSE]
#LI_score <- LI_score[match(samples, rownames(LI_score)),,drop=FALSE]


rownames(otu) <- samples
rownames(Q300) <- samples
rownames(Q600) <- samples
rownames(mre) <- samples
rownames(cov) <- samples
rownames(group) <- samples
#rownames(LI_score) <- samples


otu <- data.frame(lapply(otu, function(x) as.numeric(as.character(x))))
otu<- scale(otu)
rownames(otu) <- samples
Q300 <- data.frame(lapply(Q300, function(x) as.numeric(as.character(x))))
Q300<- scale(Q300)
Q600 <- data.frame(lapply(Q600, function(x) as.numeric(as.character(x))))
Q600 <- scale(Q600)
rownames(Q300) <- samples
rownames(Q600) <- samples
mre <- data.frame(lapply(mre, function(x) as.numeric(as.character(x))))
rownames(mre) <- samples
cov <- data.frame(lapply(cov, function(x) as.numeric(as.character(x))))
rownames(cov) <- samples
#LI_score <- data.frame(lapply(LI_score, function(x) as.numeric(as.character(x))))
#rownames(LI_score) <- samples

cov <- cov
for (i in 1:ncol(cov)) {
  cov[is.na(cov[, i]), i] <- median(cov[, i], na.rm = TRUE)
}

BD_merge_matrix_Q300 <- cbind(otu, Q300, group, cov)
BD_merge_matrix_Q600 <- cbind(otu, Q600, group, cov)

BD_merge_matrix_Q300$Group_BD <- as.factor(BD_merge_matrix_Q300$Group_BD)
BD_merge_matrix_Q600$Group_BD <- as.factor(BD_merge_matrix_Q600$Group_BD)



#####group-Q600
BD_combined_result <- data.frame()
for (otu_n in colnames(otu)) {
  #print(otu_n)
  for (Q600_n in colnames(Q600)){
    #print(Q600_n)
    for (group_n in colnames(group)){
      #print(paste(otu_n,Q600_n,group_n))
      
      #定义表达式
      my_formula_direct <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_direct <- as.formula(paste(group_n,"~",my_formula_direct))
      
      my_formula_mediate <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_mediate <- as.formula(paste(Q600_n,"~",my_formula_mediate))
      
      my_formula_combine <- paste(append(append(COV_fea, otu_n), Q600_n), collapse = '+')
      my_formula_combine <- as.formula(paste(group_n,"~",my_formula_combine))
      
      
      fit.direct <- glm(my_formula_direct, data=BD_merge_matrix_Q600, family=binomial())
      fit.combine <- glm(my_formula_combine, data=BD_merge_matrix_Q600, family=binomial())
    }
    fit.mediator <- lm(my_formula_mediate, data=BD_merge_matrix_Q600)
    
    
    #计算coe和p值
    p_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 4], n=1) 
    coe_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 1], n=1) 
    p_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 4], n=1) 
    coe_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 1], n=1) 
    p_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 4], n=1) 
    coe_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 1], n=1) 
    if (p_value_direct < 0.05 & p_value_mediate < 0.05 & p_value_combine < 0.05)
      print(paste(otu_n,Q600_n,group_n))
    
    set.seed(123456)
    if (p_value_combine < 0.05){
      #计算中介分析结果
      results <- mediate(fit.mediator, fit.combine, sims=1000, treat=otu_n, mediator=Q600_n, boot=TRUE)
      results_sum <- summary(results)
      
      if (length(unique(BD_merge_matrix_Q600[[group_n]])) == 2) {
        ACME_value <- results_sum$d.avg
        ACME_lower_value <- results_sum$d.avg.ci[1]
        ACME_upper_value <- results_sum$d.avg.ci[2]
        ACME_p_value <- results_sum$d.avg.p
        total_value <- results_sum$tau.coef
        total_lower_value <- results_sum$tau.ci[1]
        total_upper_value <- results_sum$tau.ci[2]
        total_p_value <- results_sum$tau.p
      }
      ADE_value <- results_sum$z0
      ADE_lower_value <- results_sum$z0.ci[1]
      ADE_upper_value <- results_sum$z0.ci[2]
      ADE_p_value <- results_sum$z0.p
      Prop_value <- results_sum$n0
      Prop_lower_value <- results_sum$n0.ci[1]
      Prop_upper_value <- results_sum$n0.ci[2]
      Prop_p_value <- results_sum$n0.p
      
      if (ACME_p_value < 0.05 & total_p_value < 0.05 )
        print(results_sum )
      result_df <- data.frame(otu_name=otu_n,
                              Q600_name=Q600_n, group_name=group_n,
                              Pvalue_direct=p_value_direct,
                              COEvalue_direct=coe_value_direct,
                              Pvalue_mediate=p_value_mediate, 
                              COEvalue_mediate=coe_value_mediate,
                              Pvalue_combine=p_value_combine,
                              COEvalue_combine=coe_value_combine,
                              ACME=ACME_value, ACME_lower=ACME_lower_value,
                              ACME_upper=ACME_upper_value, ACME_p=ACME_p_value,
                              Total_value=total_value, Total_lower=total_lower_value,
                              Total_upper=total_upper_value, Total_p=total_p_value,
                              ADE=ADE_value, ADE_lower=ADE_lower_value,
                              ADE_upper=ADE_upper_value, ADE_p=ADE_p_value,
                              Prop_value=Prop_value, Prop_lower=Prop_lower_value,
                              Prop_upper=Prop_upper_value, Prop_p=Prop_p_value)
      BD_combined_result <- rbind(BD_combined_result, result_df)
    }
  }
}
#fdr矫正
BD_combined_result$FDR_ACME_p=p.adjust(BD_combined_result$ACME_p,method = "fdr")
BD_combined_result$FDR_Total_p=p.adjust(BD_combined_result$Total_p,method = "fdr")


write.csv(BD_combined_result, file = '20240508_mediation_12_genus_67_Q600_BD_NBD.csv', row.names=TRUE)



####group-Q300
BD_combined_result <- data.frame()
for (otu_n in colnames(otu)) {
  #print(otu_n)
  for (Q300_n in colnames(Q300)){
    #print(Q600_n)
    for (group_n in colnames(group)){
      #print(paste(otu_n,Q600_n,group_n))
      
      #定义表达式
      my_formula_direct <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_direct <- as.formula(paste(group_n,"~",my_formula_direct))
      
      my_formula_mediate <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_mediate <- as.formula(paste(Q300_n,"~",my_formula_mediate))
      
      my_formula_combine <- paste(append(append(COV_fea, otu_n), Q300_n), collapse = '+')
      my_formula_combine <- as.formula(paste(group_n,"~",my_formula_combine))
      
      
      fit.direct <- glm(my_formula_direct, data=BD_merge_matrix_Q300, family=binomial())
      fit.combine <- glm(my_formula_combine, data=BD_merge_matrix_Q300, family=binomial())
    }
    fit.mediator <- lm(my_formula_mediate, data=BD_merge_matrix_Q300)
    
    
    #计算coe和p值
    p_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 4], n=1) 
    coe_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 1], n=1) 
    p_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 4], n=1) 
    coe_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 1], n=1) 
    p_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 4], n=1) 
    coe_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 1], n=1) 
    if (p_value_direct < 0.05 & p_value_mediate < 0.05 & p_value_combine < 0.05)
      print(paste(otu_n,Q600_n,group_n))
    
    set.seed(123456)
    if (p_value_combine < 0.05){
      #计算中介分析结果
      results <- mediate(fit.mediator, fit.combine, sims=1000, treat=otu_n, mediator=Q300_n, boot=TRUE)
      results_sum <- summary(results)
      
      if (length(unique(BD_merge_matrix_Q300[[group_n]])) == 2) {
        ACME_value <- results_sum$d.avg
        ACME_lower_value <- results_sum$d.avg.ci[1]
        ACME_upper_value <- results_sum$d.avg.ci[2]
        ACME_p_value <- results_sum$d.avg.p
        total_value <- results_sum$tau.coef
        total_lower_value <- results_sum$tau.ci[1]
        total_upper_value <- results_sum$tau.ci[2]
        total_p_value <- results_sum$tau.p
      }
      ADE_value <- results_sum$z0
      ADE_lower_value <- results_sum$z0.ci[1]
      ADE_upper_value <- results_sum$z0.ci[2]
      ADE_p_value <- results_sum$z0.p
      Prop_value <- results_sum$n0
      Prop_lower_value <- results_sum$n0.ci[1]
      Prop_upper_value <- results_sum$n0.ci[2]
      Prop_p_value <- results_sum$n0.p
      
      if (ACME_p_value < 0.05 & total_p_value < 0.05 )
        print(results_sum )
      result_df <- data.frame(otu_name=otu_n,
                              Q300_name=Q300_n, group_name=group_n,
                              Pvalue_direct=p_value_direct,
                              COEvalue_direct=coe_value_direct,
                              Pvalue_mediate=p_value_mediate, 
                              COEvalue_mediate=coe_value_mediate,
                              Pvalue_combine=p_value_combine,
                              COEvalue_combine=coe_value_combine,
                              ACME=ACME_value, ACME_lower=ACME_lower_value,
                              ACME_upper=ACME_upper_value, ACME_p=ACME_p_value,
                              Total_value=total_value, Total_lower=total_lower_value,
                              Total_upper=total_upper_value, Total_p=total_p_value,
                              ADE=ADE_value, ADE_lower=ADE_lower_value,
                              ADE_upper=ADE_upper_value, ADE_p=ADE_p_value,
                              Prop_value=Prop_value, Prop_lower=Prop_lower_value,
                              Prop_upper=Prop_upper_value, Prop_p=Prop_p_value)
      BD_combined_result <- rbind(BD_combined_result, result_df)
    }
  }
}
#fdr矫正
BD_combined_result$FDR_ACME_p=p.adjust(BD_combined_result$ACME_p,method = "fdr")
BD_combined_result$FDR_Total_p=p.adjust(BD_combined_result$Total_p,method = "fdr")


write.csv(BD_combined_result, file = '20240508_p_combine_mediation_12_genus_18_Q300_BD_NBD.csv', row.names=TRUE)





##202402 中介分析-菌群+代谢物+BD
setwd("G:\\BD_Genues_metabolites_MRE\\BD_code\\mediation_1210\\20140114mediation")
library(mediation)
library(MASS)
COV_fea <- c('Gender',	'Age',	'BMI', 'Smoking', 'Drinking')

group_df <- read.csv('treat.csv',row.names = 1)
group_df <- group_df[group_df$Group_BD %in% c('BD2','BD1'),]

LI_score <- read.csv('Limann.csv',row.names = 1)


otu      <- as.data.frame(t(read.csv('12_combined_lefse_wilcox_otu.csv',row.names = 1)))

Q300     <- as.data.frame((read.csv('select_Q300_18_0.05_P_1_VIP.csv')))
Q300      <- as.data.frame(t(Q300[seq(4,ncol(Q300))]))
colnames(Q300) <- Q300[1, ]# 将第一行作为列名
Q300 <- Q300[-1, ]

Q600     <- as.data.frame((read.csv('select_Q600_67_0.05_P_1_VIP.csv')))
Q600      <- as.data.frame(t(Q600[seq(4,ncol(Q600))]))
colnames(Q600) <- Q600[1, ]# 将第一行作为列名
Q600 <- Q600[-1, ]

mre      <- as.data.frame((read.csv('MRE_7_sig_features.csv',row.names = 1)))
mre      <- mre[seq(5,ncol(mre))]

cov <- as.data.frame((read.csv('6_clinical_covariables.csv',row.names = 1)))[COV_fea]

samples <- intersect(intersect(rownames(otu),rownames(group_df)), intersect(rownames(Q600),rownames(Q300)))

group <- group_df[,'Group_BD',drop=FALSE]

LI_score <- LI_score[,'LI',drop=FALSE]

colnames(otu) <- gsub("/", "_", colnames(otu))
colnames(Q300) <- gsub("/", "_", colnames(Q300))
colnames(Q600) <- gsub("/", "_", colnames(Q600))
colnames(mre) <- gsub("/", "_", colnames(mre))

colnames(otu) <- gsub(" ", "_", colnames(otu))
colnames(Q300) <- gsub(" ", "_", colnames(Q300))
colnames(Q600) <- gsub(" ", "_", colnames(Q600))
colnames(mre) <- gsub(" ", "_", colnames(mre))

colnames(otu) <- gsub("-", "_", colnames(otu))
colnames(Q300) <- gsub("-", "_", colnames(Q300))
colnames(Q600) <- gsub("-", "_", colnames(Q600))
colnames(mre) <- gsub("-", "_", colnames(mre))

colnames(otu) <- gsub(":", "_", colnames(otu))
colnames(Q300) <- gsub(":", "_", colnames(Q300))
colnames(Q600) <- gsub(":", "_", colnames(Q600))
colnames(mre) <- gsub(":", "_", colnames(mre))

colnames(otu) <- gsub('\\(', "_", colnames(otu))
colnames(Q300) <- gsub("\\(", "_", colnames(Q300))
colnames(Q600) <- gsub("\\(", "_", colnames(Q600))
colnames(mre) <- gsub("\\(", "_", colnames(mre))

colnames(otu) <- gsub(")", "_", colnames(otu))
colnames(Q300) <- gsub(")", "_", colnames(Q300))
colnames(Q600) <- gsub(")", "_", colnames(Q600))
colnames(mre) <- gsub(")", "_", colnames(mre))

colnames(otu) <- gsub("\\[", "", colnames(otu))
colnames(Q300) <- gsub("\\[", "", colnames(Q300))
colnames(Q600) <- gsub("\\[", "", colnames(Q600))
colnames(mre) <- gsub("\\[", "", colnames(mre))

colnames(otu) <- gsub("]", "_", colnames(otu))
colnames(Q300) <- gsub("]", "_", colnames(Q300))
colnames(Q600) <- gsub("]", "_", colnames(Q600))
colnames(mre) <- gsub("]", "_", colnames(mre))


colnames_otu  <- lapply(colnames(otu),function(x){strsplit(x,"g__")[[1]][2]})
colnames_Q300 <- lapply(colnames(Q300),function(x){paste("Fecal_", x , sep = '')})
colnames_Q600 <- lapply(colnames(Q600),function(x){paste("Serum_", x , sep = '')})
colnames_mre  <- colnames(mre)

colnames(otu)  <- colnames_otu
colnames(Q300) <- colnames_Q300
colnames(Q600) <- colnames_Q600
colnames(mre)  <- colnames_mre


# 确保所有数据框的行按'samples'排序
otu <- otu[match(samples, rownames(otu)), ]
Q600 <- Q600[match(samples, rownames(Q600)), ]
Q300 <- Q300[match(samples, rownames(Q300)), ]
mre <- mre[match(samples, rownames(mre)), ]
cov <- cov[match(samples, rownames(cov)), ]
group <- group[match(samples, rownames(group)),,drop=FALSE]
LI_score <- LI_score[match(samples, rownames(LI_score)),,drop=FALSE]

  
rownames(otu) <- samples
rownames(Q300) <- samples
rownames(Q600) <- samples
rownames(mre) <- samples
rownames(cov) <- samples
rownames(group) <- samples
rownames(LI_score) <- samples


otu <- data.frame(lapply(otu, function(x) as.numeric(as.character(x))))
otu<- scale(otu)
rownames(otu) <- samples
Q300 <- data.frame(lapply(Q300, function(x) as.numeric(as.character(x))))
Q300<- scale(Q300)
Q600 <- data.frame(lapply(Q600, function(x) as.numeric(as.character(x))))
Q600 <- scale(Q600)
rownames(Q300) <- samples
rownames(Q600) <- samples
mre <- data.frame(lapply(mre, function(x) as.numeric(as.character(x))))
rownames(mre) <- samples
cov <- data.frame(lapply(cov, function(x) as.numeric(as.character(x))))
rownames(cov) <- samples
LI_score <- data.frame(lapply(LI_score, function(x) as.numeric(as.character(x))))
rownames(LI_score) <- samples

cov <- cov
for (i in 1:ncol(cov)) {
  cov[is.na(cov[, i]), i] <- median(cov[, i], na.rm = TRUE)
}

BD_merge_matrix_Q300 <- cbind(otu, Q300, group, cov)
BD_merge_matrix_Q600 <- cbind(otu, Q600, group, cov)

score_merge_matrix_Q300 <- cbind(otu, Q300, LI_score, cov)
score_merge_matrix_Q600 <- cbind(otu, Q600, LI_score, cov)

BD_merge_matrix_Q300$Group_BD <- as.factor(BD_merge_matrix_Q300$Group_BD)
BD_merge_matrix_Q600$Group_BD <- as.factor(BD_merge_matrix_Q600$Group_BD)


merge_matrix_Q600 <- cbind(otu, Q600, mre, cov)
merge_matrix_Q300 <- cbind(otu, Q300, mre, cov)
merge_matrix_Q600$Gender <- as.factor(merge_matrix_Q600$Gender)
merge_matrix_Q600$Length <- as.factor(merge_matrix_Q600$Length)
merge_matrix_Q600$Perianal.diseases <- as.factor(merge_matrix_Q600$Perianal.diseases)
merge_matrix_Q600$Comb.sign <- as.factor(merge_matrix_Q600$Comb.sign)
merge_matrix_Q600$Bowel.wall.T2WI.sinal <- as.factor(merge_matrix_Q600$Bowel.wall.T2WI.sinal)
merge_matrix_Q600$Penetration <- factor(merge_matrix_Q600$Penetration, levels = c(0, 1, 2), ordered = TRUE)
merge_matrix_Q600$Perienteric.effusion <- factor(merge_matrix_Q600$Perienteric.effusion, levels = c(0, 1, 2), ordered = TRUE)
str(merge_matrix_Q600)

merge_matrix_Q300$Gender <- as.factor(merge_matrix_Q300$Gender)
merge_matrix_Q300$Length <- as.factor(merge_matrix_Q300$Length)
merge_matrix_Q300$Perianal.diseases <- as.factor(merge_matrix_Q300$Perianal.diseases)
merge_matrix_Q300$Comb.sign <- as.factor(merge_matrix_Q300$Comb.sign)
merge_matrix_Q300$Bowel.wall.T2WI.sinal <- as.factor(merge_matrix_Q300$Bowel.wall.T2WI.sinal)
merge_matrix_Q300$Penetration <- factor(merge_matrix_Q300$Penetration, levels = c(0, 1, 2), ordered = TRUE)
merge_matrix_Q300$Perienteric.effusion <- factor(merge_matrix_Q300$Perienteric.effusion, levels = c(0, 1, 2), ordered = TRUE)
str(merge_matrix_Q300)

#####score
score_combined_result <- data.frame()
for (otu_n in colnames(otu)) {
  #print(otu_n)
  for (Q300_n in colnames(Q300)){
    #print(Q600_n)
      
      #定义表达式
      my_formula_direct <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_direct <- as.formula(paste(LI_score,"~",my_formula_direct))
      
      my_formula_mediate <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_mediate <- as.formula(paste(Q300_n,"~",my_formula_mediate))
      
      my_formula_combine <- paste(append(append(COV_fea, otu_n), Q300_n), collapse = '+')
      my_formula_combine <- as.formula(paste(LI_score,"~",my_formula_combine))
      
      
      fit.direct <- lm(my_formula_direct, data=score_merge_matrix_Q300)
      fit.combine <- lm(my_formula_combine, data=score_merge_matrix_Q300)
    }
      fit.mediator <- lm(my_formula_mediate, data=score_merge_matrix_Q300)
    
    
    #计算coe和p值
    p_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 4], n=1) 
    coe_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 1], n=1) 
    p_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 4], n=1) 
    coe_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 1], n=1) 
    p_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 4], n=1) 
    coe_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 1], n=1) 
    if (p_value_combine < 0.05)
      print(paste(otu_n,Q300_n,LI_score))
    
    set.seed(123456)
    if (p_value_combine < 0.05){
      #计算中介分析结果
      results <- mediate(fit.mediator, fit.combine, sims=1000, treat=otu_n, mediator=Q300_n, boot=TRUE)
      results_sum <- summary(results)
      
      if (length(unique(score_merge_matrix_Q300[[group_n]])) == 2) {
        ACME_value <- results_sum$d.avg
        ACME_lower_value <- results_sum$d.avg.ci[1]
        ACME_upper_value <- results_sum$d.avg.ci[2]
        ACME_p_value <- results_sum$d.avg.p
        total_value <- results_sum$tau.coef
        total_lower_value <- results_sum$tau.ci[1]
        total_upper_value <- results_sum$tau.ci[2]
        total_p_value <- results_sum$tau.p
      }
      ADE_value <- results_sum$z0
      ADE_lower_value <- results_sum$z0.ci[1]
      ADE_upper_value <- results_sum$z0.ci[2]
      ADE_p_value <- results_sum$z0.p
      Prop_value <- results_sum$n0
      Prop_lower_value <- results_sum$n0.ci[1]
      Prop_upper_value <- results_sum$n0.ci[2]
      Prop_p_value <- results_sum$n0.p
      
      if (ACME_p_value < 0.05 & total_p_value < 0.05 )
        print(results_sum )
      result_df <- data.frame(otu_name=otu_n,
                              Q300_name=Q300_n, group_name=group_n,
                              Pvalue_direct=p_value_direct,
                              COEvalue_direct=coe_value_direct,
                              Pvalue_mediate=p_value_mediate, 
                              COEvalue_mediate=coe_value_mediate,
                              Pvalue_combine=p_value_combine,
                              COEvalue_combine=coe_value_combine,
                              ACME=ACME_value, ACME_lower=ACME_lower_value,
                              ACME_upper=ACME_upper_value, ACME_p=ACME_p_value,
                              Total_value=total_value, Total_lower=total_lower_value,
                              Total_upper=total_upper_value, Total_p=total_p_value,
                              ADE=ADE_value, ADE_lower=ADE_lower_value,
                              ADE_upper=ADE_upper_value, ADE_p=ADE_p_value,
                              Prop_value=Prop_value, Prop_lower=Prop_lower_value,
                              Prop_upper=Prop_upper_value, Prop_p=Prop_p_value)
      score_combined_result <- rbind(score_combined_result, result_df)
    }
  }

#fdr矫正
score_combined_result$FDR_ACME_p=p.adjust(score_combined_result$ACME_p,method = "fdr")
score_combined_result$FDR_Total_p=p.adjust(score_combined_result$Total_p,method = "fdr")


write.csv(score_combined_result, file = '20240207_mediation_12_genus_18_Q300_Score.csv', row.names=TRUE)









#####group-Q600
BD_combined_result <- data.frame()
for (otu_n in colnames(otu)) {
  #print(otu_n)
  for (Q600_n in colnames(Q600)){
    #print(Q600_n)
    for (group_n in colnames(group)){
      #print(paste(otu_n,Q600_n,group_n))
      
      #定义表达式
      my_formula_direct <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_direct <- as.formula(paste(group_n,"~",my_formula_direct))
      
      my_formula_mediate <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_mediate <- as.formula(paste(Q600_n,"~",my_formula_mediate))
      
      my_formula_combine <- paste(append(append(COV_fea, otu_n), Q600_n), collapse = '+')
      my_formula_combine <- as.formula(paste(group_n,"~",my_formula_combine))
      
    
      fit.direct <- glm(my_formula_direct, data=BD_merge_matrix_Q600, family=binomial())
      fit.combine <- glm(my_formula_combine, data=BD_merge_matrix_Q600, family=binomial())
      }
      fit.mediator <- lm(my_formula_mediate, data=BD_merge_matrix_Q600)
      
      
      #计算coe和p值
      p_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 4], n=1) 
      coe_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 1], n=1) 
      p_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 4], n=1) 
      coe_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 1], n=1) 
      p_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 4], n=1) 
      coe_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 1], n=1) 
      if (p_value_direct < 0.05 & p_value_mediate < 0.05 & p_value_combine < 0.05)
        print(paste(otu_n,Q600_n,group_n))
      
      set.seed(123456)
      if (p_value_combine < 0.05){
        #计算中介分析结果
        results <- mediate(fit.mediator, fit.combine, sims=1000, treat=otu_n, mediator=Q600_n, boot=TRUE)
        results_sum <- summary(results)
        
        if (length(unique(BD_merge_matrix_Q600[[group_n]])) == 2) {
          ACME_value <- results_sum$d.avg
          ACME_lower_value <- results_sum$d.avg.ci[1]
          ACME_upper_value <- results_sum$d.avg.ci[2]
          ACME_p_value <- results_sum$d.avg.p
          total_value <- results_sum$tau.coef
          total_lower_value <- results_sum$tau.ci[1]
          total_upper_value <- results_sum$tau.ci[2]
          total_p_value <- results_sum$tau.p
        }
        ADE_value <- results_sum$z0
        ADE_lower_value <- results_sum$z0.ci[1]
        ADE_upper_value <- results_sum$z0.ci[2]
        ADE_p_value <- results_sum$z0.p
        Prop_value <- results_sum$n0
        Prop_lower_value <- results_sum$n0.ci[1]
        Prop_upper_value <- results_sum$n0.ci[2]
        Prop_p_value <- results_sum$n0.p
        
        if (ACME_p_value < 0.05 & total_p_value < 0.05 )
          print(results_sum )
        result_df <- data.frame(otu_name=otu_n,
                                Q600_name=Q600_n, group_name=group_n,
                                Pvalue_direct=p_value_direct,
                                COEvalue_direct=coe_value_direct,
                                Pvalue_mediate=p_value_mediate, 
                                COEvalue_mediate=coe_value_mediate,
                                Pvalue_combine=p_value_combine,
                                COEvalue_combine=coe_value_combine,
                                ACME=ACME_value, ACME_lower=ACME_lower_value,
                                ACME_upper=ACME_upper_value, ACME_p=ACME_p_value,
                                Total_value=total_value, Total_lower=total_lower_value,
                                Total_upper=total_upper_value, Total_p=total_p_value,
                                ADE=ADE_value, ADE_lower=ADE_lower_value,
                                ADE_upper=ADE_upper_value, ADE_p=ADE_p_value,
                                Prop_value=Prop_value, Prop_lower=Prop_lower_value,
                                Prop_upper=Prop_upper_value, Prop_p=Prop_p_value)
        BD_combined_result <- rbind(BD_combined_result, result_df)
      }
  }
}
#fdr矫正
BD_combined_result$FDR_ACME_p=p.adjust(BD_combined_result$ACME_p,method = "fdr")
BD_combined_result$FDR_Total_p=p.adjust(BD_combined_result$Total_p,method = "fdr")


write.csv(BD_combined_result, file = '20240207_mediation_12_genus_67_Q600_BD_NBD.csv', row.names=TRUE)



####group-Q300
BD_combined_result <- data.frame()
for (otu_n in colnames(otu)) {
  #print(otu_n)
  for (Q300_n in colnames(Q300)){
    #print(Q600_n)
    for (group_n in colnames(group)){
      #print(paste(otu_n,Q600_n,group_n))
      
      #定义表达式
      my_formula_direct <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_direct <- as.formula(paste(group_n,"~",my_formula_direct))
      
      my_formula_mediate <- paste(append(COV_fea, otu_n), collapse = '+')
      my_formula_mediate <- as.formula(paste(Q300_n,"~",my_formula_mediate))
      
      my_formula_combine <- paste(append(append(COV_fea, otu_n), Q300_n), collapse = '+')
      my_formula_combine <- as.formula(paste(group_n,"~",my_formula_combine))
      
      
      fit.direct <- glm(my_formula_direct, data=BD_merge_matrix_Q300, family=binomial())
      fit.combine <- glm(my_formula_combine, data=BD_merge_matrix_Q300, family=binomial())
    }
    fit.mediator <- lm(my_formula_mediate, data=BD_merge_matrix_Q300)
    
    
    #计算coe和p值
    p_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 4], n=1) 
    coe_value_direct <- tail(as.data.frame(summary(fit.direct)$coefficients)[, 1], n=1) 
    p_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 4], n=1) 
    coe_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 1], n=1) 
    p_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 4], n=1) 
    coe_value_combine <- tail(as.data.frame(summary(fit.combine)$coefficients)[, 1], n=1) 
    if (p_value_direct < 0.05 & p_value_mediate < 0.05 & p_value_combine < 0.05)
      print(paste(otu_n,Q600_n,group_n))
    
    set.seed(123456)
    if (p_value_combine < 0.05){
      #计算中介分析结果
      results <- mediate(fit.mediator, fit.combine, sims=1000, treat=otu_n, mediator=Q300_n, boot=TRUE)
      results_sum <- summary(results)
      
      if (length(unique(BD_merge_matrix_Q300[[group_n]])) == 2) {
        ACME_value <- results_sum$d.avg
        ACME_lower_value <- results_sum$d.avg.ci[1]
        ACME_upper_value <- results_sum$d.avg.ci[2]
        ACME_p_value <- results_sum$d.avg.p
        total_value <- results_sum$tau.coef
        total_lower_value <- results_sum$tau.ci[1]
        total_upper_value <- results_sum$tau.ci[2]
        total_p_value <- results_sum$tau.p
      }
      ADE_value <- results_sum$z0
      ADE_lower_value <- results_sum$z0.ci[1]
      ADE_upper_value <- results_sum$z0.ci[2]
      ADE_p_value <- results_sum$z0.p
      Prop_value <- results_sum$n0
      Prop_lower_value <- results_sum$n0.ci[1]
      Prop_upper_value <- results_sum$n0.ci[2]
      Prop_p_value <- results_sum$n0.p
      
      if (ACME_p_value < 0.05 & total_p_value < 0.05 )
        print(results_sum )
      result_df <- data.frame(otu_name=otu_n,
                              Q300_name=Q300_n, group_name=group_n,
                              Pvalue_direct=p_value_direct,
                              COEvalue_direct=coe_value_direct,
                              Pvalue_mediate=p_value_mediate, 
                              COEvalue_mediate=coe_value_mediate,
                              Pvalue_combine=p_value_combine,
                              COEvalue_combine=coe_value_combine,
                              ACME=ACME_value, ACME_lower=ACME_lower_value,
                              ACME_upper=ACME_upper_value, ACME_p=ACME_p_value,
                              Total_value=total_value, Total_lower=total_lower_value,
                              Total_upper=total_upper_value, Total_p=total_p_value,
                              ADE=ADE_value, ADE_lower=ADE_lower_value,
                              ADE_upper=ADE_upper_value, ADE_p=ADE_p_value,
                              Prop_value=Prop_value, Prop_lower=Prop_lower_value,
                              Prop_upper=Prop_upper_value, Prop_p=Prop_p_value)
      BD_combined_result <- rbind(BD_combined_result, result_df)
    }
  }
}
#fdr矫正
BD_combined_result$FDR_ACME_p=p.adjust(BD_combined_result$ACME_p,method = "fdr")
BD_combined_result$FDR_Total_p=p.adjust(BD_combined_result$Total_p,method = "fdr")


write.csv(BD_combined_result, file = '20240208_p_combine_mediation_12_genus_18_Q300_BD_NBD.csv', row.names=TRUE)

