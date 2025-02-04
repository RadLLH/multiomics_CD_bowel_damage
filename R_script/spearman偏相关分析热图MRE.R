###spearman偏相关分析--画热图MRE
# 加载 stats 包
library(stats)
# 加载包
library(igraph)
library(psych)
library(ppcor)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(grid)

setwd('G:\\BD_Genues_metabolites_MRE\\BD_code\\heatmap\\热图聚类和棒棒糖图\\bangbangtangtu')

#COV_fea <- c('Gender',	'Age',	'BMI', 'Smoking', 'Drinking')
COV_fea <- c('Gender',	'Age',	'BMI', 'Smoking', 'Drinking','Biologics','Antibiotics','Corticosteroids','Probiotics','Immunomodulator')

group_df <- read.csv('treat.csv',row.names = 1)
group_df <- group_df[group_df$Group_BD == 'BD2',]

#用49个菌属
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


colnames_otu  <- lapply(colnames(otu),function(x){strsplit(x,"g__")[[1]][2]})
colnames_Q300 <- lapply(colnames(Q300),function(x){paste("Fecal_", x , sep = '')})
colnames_Q600 <- lapply(colnames(Q600),function(x){paste("Serum_", x , sep = '')})
colnames_mre  <- colnames(mre)

group_df <- group_df[samples,]
otu  <- as.data.frame(lapply(otu[samples,],as.numeric))

Q300 <- as.data.frame(scale(as.data.frame(lapply(Q300[samples,],as.numeric))))
Q600 <- as.data.frame(scale(as.data.frame(lapply(Q600[samples,],as.numeric))))

mre  <- as.data.frame(lapply(mre[samples,],as.numeric))
cov <- as.data.frame(lapply(cov[samples,],as.numeric))
cov_imputed <- cov
for (i in 1:ncol(cov)) {
  cov_imputed[is.na(cov[, i]), i] <- median(cov[, i], na.rm = TRUE)
}


##修改菌属名
colnames(otu)  <- colnames_otu
colnames(Q300) <- colnames_Q300
colnames(Q600) <- colnames_Q600
colnames(mre)  <- colnames_mre

colnames(otu) <- paste0("Genus_", colnames(otu))

# 初始化空矩阵存储相关性系数和P值
cor_matrix <- matrix(NA, nrow = ncol(mre), ncol = ncol(otu))
p_values <- matrix(NA, nrow = ncol(mre), ncol = ncol(otu))

cor_result <- data.frame(matrix(nrow = ncol(mre), ncol = ncol(otu)))
rownames(cor_result) <- colnames(mre)
colnames(cor_result) <- colnames(otu)

p_result <- data.frame(matrix(nrow = ncol(mre), ncol = ncol(otu)))
rownames(p_result) <- colnames(mre)
colnames(p_result) <- colnames(otu)


# 循环遍历 mre 的每一列和 otu 的每一列进行偏相关性分析
for (i in 1:ncol(mre)) {
  for (j in 1:ncol(otu)) {
    # 创建联合数据框，包括 mre 列, otu 列和混杂因素
    joint_data <- data.frame(mre_column = mre[, i], otu_column = otu[, j], cov_imputed)
    # 计算偏相关性和P值
    cor_ <- pcor.test(joint_data$mre_column, joint_data$otu_column, joint_data[,3:(2+length(COV_fea))], 
                      method="spearman")
    cor_matrix[i, j] <- cor_$estimate
    p_values[i, j] <- cor_$p.value
    # 存储偏相关性系数和P值
    cor_result[i, j] <- cor_$estimate
    p_result[i, j] <- cor_$p.value
  }
}

# 对这个向量进行 FDR 校正
fdr_values_vector <- p.adjust(p_values, method = "fdr")
fdr_result_matrix <- matrix(fdr_values_vector, nrow = nrow(p_values), ncol = ncol(p_values))
fdr_result_df <- as.data.frame(fdr_result_matrix)
rownames(fdr_result_df) <- colnames(mre)
colnames(fdr_result_df) <- colnames(otu)
write.csv(p_result,file="20240508_BD2_Spearman_pcor_p_result_7mre_12otu.csv")
write.csv(fdr_result_df,file="20240508_BD2_Spearman_pcor_fdr_result_7mre_12otu.csv")

#筛选p<0.05 并且 |cor|>0.25的数据
# 初始化一个空的数据框来存储结果
select_result <- data.frame(rowName = character(), colName = character(), stringsAsFactors = FALSE)

# 遍历每个元素
for (i in rownames(p_result)) {
  for (j in colnames(p_result)) {
    # 检查每个元素是否满足条件
    if (p_result[i,j] < 0.05 && abs(cor_result[i,j]) > 0.25) {
      # 如果满足条件，那么将行名和列名添加到结果数据框中
      select_result <- rbind(select_result, data.frame(otu = i, mre = j, p_value = p_result[i,j], cor_value = cor_result[i,j], fdr = fdr_result_df[i,j], stringsAsFactors = FALSE))
    }
  }
}
write.csv(select_result,file="20240508_偏相关BD2_Spearman_otu_7mre_select_result.csv")

otu_select_names <- unique(select_result$otu)
mre_select_names <- unique(select_result$mre)
cor_result_select <- cor_result[otu_select_names,mre_select_names]
p_result_select <- p_result[otu_select_names,mre_select_names]



# 标有显著性的热图
getPSig <- function(dc) {
  sc <- ''
  if (dc < 0.001) sc <- '***'
  else if (dc < 0.01) sc <- '**'
  else if (dc < 0.05) sc <- '*'
  sc
}



sig_mat <- matrix(sapply(as.matrix(p_result_select), getPSig), nrow = nrow(as.matrix(p_result_select)))
p_values <- format(as.matrix(p_result_select), scientific = FALSE, digits = 3)

pdf("20240508_偏相关heatmap_OTU_7MRE_p_heng.pdf",width = 12,height = 8)
sig_mat <- matrix(sapply(as.matrix(p_result_select), getPSig), nrow=nrow(as.matrix(p_result_select)))
heatmap_pic_new <- pheatmap(
  as.matrix(cor_result_select),
  cellwidth = 30,
  cellheight = 30,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle = 90,
  display_numbers = sig_mat,
  fontsize_row = 12,  # 设置横坐标标签的字体大小
  fontsize_col = 12,  # 设置纵坐标标签的字体大小
  fontface_row = "bold",  # 设置横坐标标签为粗体
  fontface_col = "bold",  # 设置纵坐标标签为粗体
  color = colors <- colorRampPalette(c("#4E81BE", "lightyellow", "#E64B35CC"))(11)
)
while (!is.null(dev.list()))  dev.off()


########mre q300

cor_matrix <- matrix(NA, nrow = ncol(mre), ncol = ncol(Q300))
p_values <- matrix(NA, nrow = ncol(mre), ncol = ncol(Q300))

cor_result <- data.frame(matrix(nrow = ncol(mre), ncol = ncol(Q300)))
rownames(cor_result) <- colnames(mre)
colnames(cor_result) <- colnames(Q300)

p_result <- data.frame(matrix(nrow = ncol(mre), ncol = ncol(Q300)))
rownames(p_result) <- colnames(mre)
colnames(p_result) <- colnames(Q300)


# 循环遍历 mre 的每一列和 Q300 的每一列进行偏相关性分析
for (i in 1:ncol(mre)) {
  for (j in 1:ncol(Q300)) {
    # 创建联合数据框，包括 mre 列, otu 列和混杂因素
    joint_data <- data.frame(mre_column = mre[, i], Q300_column = Q300[, j], cov_imputed)
    # 计算偏相关性和P值
    cor_ <- pcor.test(joint_data$mre_column, joint_data$Q300_column, joint_data[,3:(2+length(COV_fea))], 
                      method="spearman")
    cor_matrix[i, j] <- cor_$estimate
    p_values[i, j] <- cor_$p.value
    # 存储偏相关性系数和P值
    cor_result[i, j] <- cor_$estimate
    p_result[i, j] <- cor_$p.value
  }
}

# 对这个向量进行 FDR 校正
fdr_values_vector <- p.adjust(p_values, method = "fdr")
fdr_result_matrix <- matrix(fdr_values_vector, nrow = nrow(p_values), ncol = ncol(p_values))
fdr_result_df <- as.data.frame(fdr_result_matrix)
rownames(fdr_result_df) <- colnames(mre)
colnames(fdr_result_df) <- colnames(Q300)
write.csv(p_result,file="20240508_BD2_Spearman_pcor_p_result_7mre_18Q300.csv")
write.csv(fdr_result_df,file="20240508_BD2_Spearman_pcor_fdr_result_7mre_18Q300.csv")
#筛选p<0.05 并且 |cor|>0.25的数据


select_result <- data.frame(rowName = character(), colName = character(), stringsAsFactors = FALSE)

# 遍历每个元素
for (i in rownames(p_result)) {
  for (j in colnames(p_result)) {
    # 检查每个元素是否满足条件
    if (p_result[i,j] < 0.05 && abs(cor_result[i,j]) > 0.25) {
      # 如果满足条件，那么将行名和列名添加到结果数据框中
      select_result <- rbind(select_result, data.frame(Q300 = i, mre = j, p_value = p_result[i,j], cor_value = cor_result[i,j], fdr = fdr_result_df[i,j], stringsAsFactors = FALSE))
    }
  }
}
write.csv(select_result,file="20240508_偏相关BD2_Spearman_Q300_7mre_select_result.csv")

Q300_select_names <- unique(select_result$Q300)
mre_select_names <- unique(select_result$mre)
cor_result_select <- cor_result[Q300_select_names,mre_select_names]
p_result_select <- p_result[Q300_select_names,mre_select_names]



# 标有显著性的热图
getPSig <- function(dc) {
  sc <- ''
  if (dc < 0.001) sc <- '***'
  else if (dc < 0.01) sc <- '**'
  else if (dc < 0.05) sc <- '*'
  sc
}



sig_mat <- matrix(sapply(as.matrix(p_result_select), getPSig), nrow = nrow(as.matrix(p_result_select)))
p_values <- format(as.matrix(p_result_select), scientific = FALSE, digits = 3)

pdf("20240508_偏相关heatmap_Q300_7MRE_p_heng.pdf",width = 12,height = 8)
sig_mat <- matrix(sapply(as.matrix(p_result_select), getPSig), nrow=nrow(as.matrix(p_result_select)))
heatmap_pic_new <- pheatmap(
  as.matrix(cor_result_select),
  cellwidth = 30,
  cellheight = 30,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle = 90,
  display_numbers = sig_mat,
  fontsize_row = 12,  # 设置横坐标标签的字体大小
  fontsize_col = 12,  # 设置纵坐标标签的字体大小
  fontface_row = "bold",  # 设置横坐标标签为粗体
  fontface_col = "bold",  # 设置纵坐标标签为粗体
  color = colors <- colorRampPalette(c("#4E81BE", "lightyellow", "#E64B35CC"))(11)
)
while (!is.null(dev.list()))  dev.off()




######mre q600
# 初始化空矩阵存储相关性系数和P值
cor_matrix <- matrix(NA, nrow = ncol(mre), ncol = ncol(Q600))
p_values <- matrix(NA, nrow = ncol(mre), ncol = ncol(Q600))

cor_result <- data.frame(matrix(nrow = ncol(mre), ncol = ncol(Q600)))
rownames(cor_result) <- colnames(mre)
colnames(cor_result) <- colnames(Q600)

p_result <- data.frame(matrix(nrow = ncol(mre), ncol = ncol(Q600)))
rownames(p_result) <- colnames(mre)
colnames(p_result) <- colnames(Q600)

# 循环遍历 mre 的每一列和 Q600 的每一列进行偏相关性分析
for (i in 1:ncol(mre)) {
  for (j in 1:ncol(Q600)) {
    # 创建联合数据框，包括 mre 列, otu 列和混杂因素
    joint_data <- data.frame(mre_column = mre[, i], Q600_column = Q600[, j], cov_imputed)
    # 计算偏相关性和P值
    cor_ <- pcor.test(joint_data$mre_column, joint_data$Q600_column, joint_data[,3:(2+length(COV_fea))], 
                      method="spearman")
    cor_matrix[i, j] <- cor_$estimate
    p_values[i, j] <- cor_$p.value
    # 存储偏相关性系数和P值
    cor_result[i, j] <- cor_$estimate
    p_result[i, j] <- cor_$p.value
  }
}

# 对这个向量进行 FDR 校正
fdr_values_vector <- p.adjust(p_values, method = "fdr")
fdr_result_matrix <- matrix(fdr_values_vector, nrow = nrow(p_values), ncol = ncol(p_values))
fdr_result_df <- as.data.frame(fdr_result_matrix)
rownames(fdr_result_df) <- colnames(mre)
colnames(fdr_result_df) <- colnames(Q600)
write.csv(p_result,file="BD2_Spearman_pcor_p_result_7mre_67Q600.csv")
write.csv(fdr_result_df,file="BD2_Spearman_pcor_fdr_result_7mre_67Q600.csv")




#筛选p<0.05 并且 |cor|>0.25的数据
select_result <- data.frame(rowName = character(), colName = character(), stringsAsFactors = FALSE)

# 遍历每个元素
for (i in rownames(p_result)) {
  for (j in colnames(p_result)) {
    # 检查每个元素是否满足条件
    if (p_result[i,j] < 0.05 && abs(cor_result[i,j]) > 0.25) {
      # 如果满足条件，那么将行名和列名添加到结果数据框中
      select_result <- rbind(select_result, data.frame(Q600 = i, mre = j, p_value = p_result[i,j], cor_value = cor_result[i,j], fdr = fdr_result_df[i,j], stringsAsFactors = FALSE))
    }
  }
}
write.csv(select_result,file="20240508_偏相关BD2_Spearman_Q600_7mre_select_result.csv")
Q600_select_names <- unique(select_result$Q600)
mre_select_names <- unique(select_result$mre)
cor_result_select <- cor_result[Q600_select_names,mre_select_names]
p_result_select <- p_result[Q600_select_names,mre_select_names]



# 标有显著性的热图
getPSig <- function(dc) {
  sc <- ''
  if (dc < 0.001) sc <- '***'
  else if (dc < 0.01) sc <- '**'
  else if (dc < 0.05) sc <- '*'
  sc
}



sig_mat <- matrix(sapply(as.matrix(p_result_select), getPSig), nrow = nrow(as.matrix(p_result_select)))
p_values <- format(as.matrix(p_result_select), scientific = FALSE, digits = 3)

pdf("20240508_偏相关heatmap_Q600_7MRE_p_heng.pdf",width = 12,height = 8)
sig_mat <- matrix(sapply(as.matrix(p_result_select), getPSig), nrow=nrow(as.matrix(p_result_select)))
heatmap_pic_new <- pheatmap(
  as.matrix(cor_result_select),
  cellwidth = 30,
  cellheight = 30,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle = 90,
  display_numbers = sig_mat,
  fontsize_row = 12,  # 设置横坐标标签的字体大小
  fontsize_col = 12,  # 设置纵坐标标签的字体大小
  fontface_row = "bold",  # 设置横坐标标签为粗体
  fontface_col = "bold",  # 设置纵坐标标签为粗体
  color = colors <- colorRampPalette(c("#4E81BE", "lightyellow", "#E64B35CC"))(11)
)
while (!is.null(dev.list()))  dev.off()



