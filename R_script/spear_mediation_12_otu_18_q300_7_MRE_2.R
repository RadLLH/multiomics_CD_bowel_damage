setwd("G:\\BD_Genues_metabolites_MRE\\BD_code\\mediation_1210\\20140114mediation")
library(mediation)
COV_fea <- c('Gender',	'Age')

group_df <- read.csv('treat.csv',row.names = 1)
group_df <- group_df[group_df$Group_BD=='BD2',]


otu      <- as.data.frame(t(read.csv('filtered_otus.csv',row.names = 1)))

Q300     <- as.data.frame(read.csv('filtered_Q300.csv',row.names = 1))

Q600     <- as.data.frame(read.csv('filtered_Q600.csv',row.names = 1))

mre      <- as.data.frame((read.csv('MRE_7_sig_features.csv',row.names = 1)))
mre      <- mre[seq(5,ncol(mre))]

cov <- as.data.frame((read.csv('6_clinical_covariables.csv',row.names = 1)))[COV_fea]

samples <- intersect(intersect(rownames(otu),rownames(group_df)), intersect(rownames(Q600),rownames(Q300)))



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

##如果用28个菌属csv，不需要修改菌属名，注释掉下行代码
colnames_otu  <- lapply(colnames(otu),function(x){strsplit(x,"g__")[[1]][2]})

colnames_Q300 <- lapply(colnames(Q300),function(x){paste("Fecal_", x , sep = '')})
colnames_Q600 <- lapply(colnames(Q600),function(x){paste("Serum_", x , sep = '')})
colnames_mre  <- colnames(mre)


otu  <- as.data.frame(lapply(otu[samples,],as.numeric))
#对代谢物数据进行标准化，需去掉以下两行的注释
Q300 <- as.data.frame(scale(as.data.frame(lapply(Q300[samples,],as.numeric))))
Q600 <- as.data.frame(scale(as.data.frame(lapply(Q600[samples,],as.numeric))))

#对代谢物数据进行标准化需注释掉以下两行
#Q300 <- as.data.frame(lapply(Q300[samples,],as.numeric))
#Q600 <- as.data.frame(lapply(Q600[samples,],as.numeric))


mre  <- as.data.frame(lapply(mre[samples,],as.numeric))
cov <- as.data.frame(lapply(cov[samples,],as.numeric))
cov_imputed <- cov
for (i in 1:ncol(cov)) {
  cov_imputed[is.na(cov[, i]), i] <- median(cov[, i], na.rm = TRUE)
}

##如果用28个菌属csv，不需要修改菌属名,注释掉下行代码
colnames(otu)  <- colnames_otu

colnames(Q300) <- colnames_Q300
colnames(Q600) <- colnames_Q600
colnames(mre)  <- colnames_mre


merge_matrix <- cbind(otu,Q300,mre,cov) 

combined_result <- data.frame()
for (otu_n in colnames(otu)) {
  #print(otu_n)
  for (Q300_n in colnames(Q300)){
    #print(Q300_n)
    for (mre_n in colnames(mre)){
      print(paste(otu_n,Q300_n,mre_n))
      
      #定义表达式
      my_formula_direct <- paste(append(cov_imputed, otu_n), collapse = '+')
      my_formula_direct <- as.formula(paste(mre_n,"~",my_formula_direct))
      
      my_formula_mediate <- paste(append(cov_imputed, otu_n), collapse = '+')
      my_formula_mediate <- as.formula(paste(Q300_n,"~",my_formula_mediate))
      
      my_formula_combine <- paste(append(append(cov_imputed, otu_n), Q300_n), collapse = '+')
      my_formula_combine <- as.formula(paste(mre_n,"~",my_formula_combine))
      
      #定义lm模型
      if (length(unique(merge_matrix[[mre_n]])) > 2) {
        merge_matrix[, mre_n] <- as.numeric(merge_matrix[[mre_n]])
        tmp.model_mre_otu <- lm(my_formula_direct, data=merge_matrix)
        fit.dv <- lm(my_formula_combine, data=merge_matrix)
      }
      if (length(unique(merge_matrix[[mre_n]])) == 2) {
        merge_matrix[, mre_n] <- as.factor(merge_matrix[[mre_n]])
        tmp.model_mre_otu <- glm(my_formula_direct, data=merge_matrix, family=binomial())
        fit.dv <- glm(my_formula_combine, data=merge_matrix, family=binomial())
      }
      fit.mediator <- lm(my_formula_mediate, data=merge_matrix)
      
      #计算coe和p值
      p_value_direct <- tail(as.data.frame(summary(tmp.model_mre_otu)$coefficients)[, 4], n=1) 
      coe_value_direct <- tail(as.data.frame(summary(tmp.model_mre_otu)$coefficients)[, 1], n=1) 
      p_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 4], n=1) 
      coe_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 1], n=1) 
      p_value_combine <- tail(as.data.frame(summary(fit.dv)$coefficients)[, 4], n=1) 
      coe_value_combine <- tail(as.data.frame(summary(fit.dv)$coefficients)[, 1], n=1) 
      print(p_value_combine)
      
      set.seed(666)
      if (p_value_combine < 0.1){
      #计算中介分析结果
      results <- mediate(fit.mediator, fit.dv, sims=1000, treat=otu_n, mediator=Q300_n, boot=TRUE)
      results_sum <- summary(results)
      
      if (length(unique(merge_matrix[[mre_n]])) > 2) {
        ACME_value <- results_sum$d0
        ACME_lower_value <- results_sum$d0.ci[1]
        ACME_upper_value <- results_sum$d0.ci[2]
        ACME_p_value <- results_sum$d0.p
        total_value <- results_sum$tau.coef
        total_lower_value <- results_sum$tau.ci[1]
        total_upper_value <- results_sum$tau.ci[2]
        total_p_value <- results_sum$tau.p
      }
      if (length(unique(merge_matrix[[mre_n]])) == 2) {
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
      
      print(results_sum)
      result_df <- data.frame(otu_name=otu_n,
                              q300_name=Q300_n, mre_name=mre_n,
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
      combined_result <- rbind(combined_result, result_df)
      }
    }
  }
}
#fdr矫正
combined_result$FDR_ACME_p=p.adjust(combined_result$ACME_p,method = "fdr")
combined_result$FDR_Total_p=p.adjust(combined_result$Total_p,method = "fdr")


write.csv(combined_result, file = '20240126_spear_mediation_12_otu_18_q300_7_MRE_2.csv', row.names=TRUE)