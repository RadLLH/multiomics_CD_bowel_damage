options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
BiocManager::install("microbiomeMarker")
BiocManager::install("DescTools")

library(tidyverse)
library(magrittr)
library(DescTools)
library(microbiomeMarker)
library(phyloseq)


####### HC VS CD

#otu <- as.data.frame(t(read.csv('otu_table_genus.csv',row.names = 1)))
otu <- as.data.frame(t(read.csv('./genus_select_FDR/73otu_result.csv',row.names = 1)))

col_names <- colnames(otu)

# 分割每个列名
split_names <- strsplit(col_names, ";")
# 将list中的vector转化为character
id_chars <- paste("otu_ID", 1:ncol(otu), sep='_')
# 将结果转换为dataframe
tax <- data.frame(
  rbind(do.call(rbind, split_names))
)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax) <- id_chars
colnames(otu) <- id_chars


cov <- as.data.frame((read.csv('./lefse/Covariables.csv',row.names = 1)))
samples <- intersect(rownames(cov),rownames(otu))
cov <- cov[samples,]
otu  <- otu[samples,]
# Replace NAs with median in numeric columns
numeric_columns <- sapply(cov, is.numeric)
cov[, numeric_columns] <- lapply(cov[, numeric_columns], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

otu <- t(otu)

#cov <- cov[(cov$Group != "BD2"),]
otu <- otu[,colnames(otu) %in% rownames(cov)]
cov$Group <- ifelse(cov$Group %in% c("BD1", "BD2"), "CD", cov$Group)


otu %<>% as.matrix() 
tax %<>% as.matrix() 


# 1.2.1 构建phyloseq对象

physeq <- phyloseq(
  otu_table(otu,taxa_are_rows = TRUE), 
  tax_table(tax), # 数据格式必须是矩阵，否则会改变行列名
  sample_data(cov)
)


# 1.2.2 数据标准化-抽平
##?normalize # 查看数据标准化方法对应字符缩写。
set.seed(12345)
data <- physeq
data

# 2.1.1 lefse分析-LDA阈值设置为4

lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "Group",
  taxa_rank = "none",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.05, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.05, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 30, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。

# 2.1.2 提取biomarker鉴定结果
##没有注释名的分类单元，将会按照前几个分类单元的注释结果进行命名。
lefse %>% marker_table() -> res.diff



dim(res.diff) # 61个biomarker。
res.diff %>% head()

result = data.frame(feature=res.diff$feature,enrich_group=res.diff$enrich_group,ef_lda=res.diff$ef_lda,pvalue=res.diff$pvalue,padj=res.diff$padj)
last_char <- sapply(strsplit(result$feature, "_"), function(x) x[length(x)])
last_num <- as.numeric(last_char)
result$feature <- col_names[last_num]
write.csv(result,file="./lefse/HC_CD_p0.05.csv")


####### HC VS NBD

#otu <- as.data.frame(t(read.csv('otu_table_genus.csv',row.names = 1)))
otu <- as.data.frame(t(read.csv('./genus_select_FDR/73otu_result.csv',row.names = 1)))

col_names <- colnames(otu)

# 分割每个列名
split_names <- strsplit(col_names, ";")
# 将list中的vector转化为character
id_chars <- paste("otu_ID", 1:ncol(otu), sep='_')
# 将结果转换为dataframe
tax <- data.frame(
  rbind(do.call(rbind, split_names))
)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax) <- id_chars
colnames(otu) <- id_chars


cov <- as.data.frame((read.csv('./lefse/Covariables.csv',row.names = 1)))
samples <- intersect(rownames(cov),rownames(otu))
cov <- cov[samples,]
otu  <- otu[samples,]
# Replace NAs with median in numeric columns
numeric_columns <- sapply(cov, is.numeric)
cov[, numeric_columns] <- lapply(cov[, numeric_columns], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

otu <- t(otu)

cov <- cov[(cov$Group != "BD2"),]
otu <- otu[,colnames(otu) %in% rownames(cov)]
cov$Group <- ifelse(cov$Group %in% c("BD1"), "NBD", cov$Group)


otu %<>% as.matrix() 
tax %<>% as.matrix() 


# 1.2.1 构建phyloseq对象

physeq <- phyloseq(
  otu_table(otu,taxa_are_rows = TRUE), 
  tax_table(tax), # 数据格式必须是矩阵，否则会改变行列名
  sample_data(cov)
)


# 1.2.2 数据标准化-抽平
##?normalize # 查看数据标准化方法对应字符缩写。
set.seed(12345)
data <- physeq
data

# 2.1.1 lefse分析-LDA阈值设置为4

lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "Group",
  taxa_rank = "none",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.05, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.05, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 30, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。

# 2.1.2 提取biomarker鉴定结果
##没有注释名的分类单元，将会按照前几个分类单元的注释结果进行命名。
lefse %>% marker_table() -> res.diff



dim(res.diff) # 61个biomarker。
res.diff %>% head()

result = data.frame(feature=res.diff$feature,enrich_group=res.diff$enrich_group,ef_lda=res.diff$ef_lda,pvalue=res.diff$pvalue,padj=res.diff$padj)
last_char <- sapply(strsplit(result$feature, "_"), function(x) x[length(x)])
last_num <- as.numeric(last_char)
result$feature <- col_names[last_num]
write.csv(result,file="./lefse/HC_NBD_p0.05.csv")





####### HC VS BD

#otu <- as.data.frame(t(read.csv('otu_table_genus.csv',row.names = 1)))
otu <- as.data.frame(t(read.csv('./genus_select_FDR/73otu_result.csv',row.names = 1)))

col_names <- colnames(otu)

# 分割每个列名
split_names <- strsplit(col_names, ";")
# 将list中的vector转化为character
id_chars <- paste("otu_ID", 1:ncol(otu), sep='_')
# 将结果转换为dataframe
tax <- data.frame(
  rbind(do.call(rbind, split_names))
)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax) <- id_chars
colnames(otu) <- id_chars


cov <- as.data.frame((read.csv('./lefse/Covariables.csv',row.names = 1)))
samples <- intersect(rownames(cov),rownames(otu))
cov <- cov[samples,]
otu  <- otu[samples,]
# Replace NAs with median in numeric columns
numeric_columns <- sapply(cov, is.numeric)
cov[, numeric_columns] <- lapply(cov[, numeric_columns], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

otu <- t(otu)

cov <- cov[(cov$Group != "BD1"),]
otu <- otu[,colnames(otu) %in% rownames(cov)]
cov$Group <- ifelse(cov$Group %in% c("BD2"), "BD", cov$Group)


otu %<>% as.matrix() 
tax %<>% as.matrix() 


# 1.2.1 构建phyloseq对象

physeq <- phyloseq(
  otu_table(otu,taxa_are_rows = TRUE), 
  tax_table(tax), # 数据格式必须是矩阵，否则会改变行列名
  sample_data(cov)
)


# 1.2.2 数据标准化-抽平
##?normalize # 查看数据标准化方法对应字符缩写。
set.seed(12345)
data <- physeq
data

# 2.1.1 lefse分析-LDA阈值设置为4

lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "Group",
  taxa_rank = "none",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.05, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.05, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 30, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。

# 2.1.2 提取biomarker鉴定结果
##没有注释名的分类单元，将会按照前几个分类单元的注释结果进行命名。
lefse %>% marker_table() -> res.diff



dim(res.diff) # 61个biomarker。
res.diff %>% head()

result = data.frame(feature=res.diff$feature,enrich_group=res.diff$enrich_group,ef_lda=res.diff$ef_lda,pvalue=res.diff$pvalue,padj=res.diff$padj)
last_char <- sapply(strsplit(result$feature, "_"), function(x) x[length(x)])
last_num <- as.numeric(last_char)
result$feature <- col_names[last_num]
write.csv(result,file="./lefse/HC_BD_p0.05.csv")




####### BD VS NBD

#otu <- as.data.frame(t(read.csv('otu_table_genus.csv',row.names = 1)))
otu <- as.data.frame(t(read.csv('./genus_select_FDR/73otu_result.csv',row.names = 1)))

col_names <- colnames(otu)

# 分割每个列名
split_names <- strsplit(col_names, ";")
# 将list中的vector转化为character
id_chars <- paste("otu_ID", 1:ncol(otu), sep='_')
# 将结果转换为dataframe
tax <- data.frame(
  rbind(do.call(rbind, split_names))
)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax) <- id_chars
colnames(otu) <- id_chars


cov <- as.data.frame((read.csv('./lefse/Covariables.csv',row.names = 1)))
samples <- intersect(rownames(cov),rownames(otu))
cov <- cov[samples,]
otu  <- otu[samples,]
# Replace NAs with median in numeric columns
numeric_columns <- sapply(cov, is.numeric)
cov[, numeric_columns] <- lapply(cov[, numeric_columns], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

otu <- t(otu)

cov <- cov[(cov$Group != "HC"),]
otu <- otu[,colnames(otu) %in% rownames(cov)]
cov$Group <- ifelse(cov$Group %in% c("BD2"), "BD", cov$Group)
cov$Group <- ifelse(cov$Group %in% c("BD1"), "NBD", cov$Group)


otu %<>% as.matrix() 
tax %<>% as.matrix() 


# 1.2.1 构建phyloseq对象

physeq <- phyloseq(
  otu_table(otu,taxa_are_rows = TRUE), 
  tax_table(tax), # 数据格式必须是矩阵，否则会改变行列名
  sample_data(cov)
)


# 1.2.2 数据标准化-抽平
##?normalize # 查看数据标准化方法对应字符缩写。
set.seed(12345)
data <- physeq
data

# 2.1.1 lefse分析-LDA阈值设置为4

lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "Group",
  taxa_rank = "none",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.05, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.05, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 30, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。

# 2.1.2 提取biomarker鉴定结果
##没有注释名的分类单元，将会按照前几个分类单元的注释结果进行命名。
lefse %>% marker_table() -> res.diff



dim(res.diff) # 61个biomarker。
res.diff %>% head()

result = data.frame(feature=res.diff$feature,enrich_group=res.diff$enrich_group,ef_lda=res.diff$ef_lda,pvalue=res.diff$pvalue,padj=res.diff$padj)
last_char <- sapply(strsplit(result$feature, "_"), function(x) x[length(x)])
last_num <- as.numeric(last_char)
result$feature <- col_names[last_num]
write.csv(result,file="./lefse/BD_NBD_p0.05.csv")






############ p=0.1

####### HC VS CD

#otu <- as.data.frame(t(read.csv('otu_table_genus.csv',row.names = 1)))
otu <- as.data.frame(t(read.csv('./genus_select_FDR/73otu_result.csv',row.names = 1)))

col_names <- colnames(otu)

# 分割每个列名
split_names <- strsplit(col_names, ";")
# 将list中的vector转化为character
id_chars <- paste("otu_ID", 1:ncol(otu), sep='_')
# 将结果转换为dataframe
tax <- data.frame(
  rbind(do.call(rbind, split_names))
)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax) <- id_chars
colnames(otu) <- id_chars


cov <- as.data.frame((read.csv('./lefse/Covariables.csv',row.names = 1)))
samples <- intersect(rownames(cov),rownames(otu))
cov <- cov[samples,]
otu  <- otu[samples,]
# Replace NAs with median in numeric columns
numeric_columns <- sapply(cov, is.numeric)
cov[, numeric_columns] <- lapply(cov[, numeric_columns], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

otu <- t(otu)

#cov <- cov[(cov$Group != "BD2"),]
otu <- otu[,colnames(otu) %in% rownames(cov)]
cov$Group <- ifelse(cov$Group %in% c("BD1", "BD2"), "CD", cov$Group)


otu %<>% as.matrix() 
tax %<>% as.matrix() 


# 1.2.1 构建phyloseq对象

physeq <- phyloseq(
  otu_table(otu,taxa_are_rows = TRUE), 
  tax_table(tax), # 数据格式必须是矩阵，否则会改变行列名
  sample_data(cov)
)


# 1.2.2 数据标准化-抽平
##?normalize # 查看数据标准化方法对应字符缩写。
set.seed(12345)
data <- physeq
data

# 2.1.1 lefse分析-LDA阈值设置为4

lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "Group",
  taxa_rank = "none",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.1, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.1, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 30, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。

# 2.1.2 提取biomarker鉴定结果
##没有注释名的分类单元，将会按照前几个分类单元的注释结果进行命名。
lefse %>% marker_table() -> res.diff



dim(res.diff) # 61个biomarker。
res.diff %>% head()

result = data.frame(feature=res.diff$feature,enrich_group=res.diff$enrich_group,ef_lda=res.diff$ef_lda,pvalue=res.diff$pvalue,padj=res.diff$padj)
last_char <- sapply(strsplit(result$feature, "_"), function(x) x[length(x)])
last_num <- as.numeric(last_char)
result$feature <- col_names[last_num]
write.csv(result,file="./lefse/HC_CD_p0.1.csv")
result$name <- sapply(strsplit(result$feature, ";"), function(x) x[length(x)])

####### HC VS NBD

#otu <- as.data.frame(t(read.csv('otu_table_genus.csv',row.names = 1)))
otu <- as.data.frame(t(read.csv('./genus_select_FDR/73otu_result.csv',row.names = 1)))

col_names <- colnames(otu)

# 分割每个列名
split_names <- strsplit(col_names, ";")
# 将list中的vector转化为character
id_chars <- paste("otu_ID", 1:ncol(otu), sep='_')
# 将结果转换为dataframe
tax <- data.frame(
  rbind(do.call(rbind, split_names))
)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax) <- id_chars
colnames(otu) <- id_chars


cov <- as.data.frame((read.csv('./lefse/Covariables.csv',row.names = 1)))
samples <- intersect(rownames(cov),rownames(otu))
cov <- cov[samples,]
otu  <- otu[samples,]
# Replace NAs with median in numeric columns
numeric_columns <- sapply(cov, is.numeric)
cov[, numeric_columns] <- lapply(cov[, numeric_columns], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

otu <- t(otu)

cov <- cov[(cov$Group != "BD2"),]
otu <- otu[,colnames(otu) %in% rownames(cov)]
cov$Group <- ifelse(cov$Group %in% c("BD1"), "NBD", cov$Group)


otu %<>% as.matrix() 
tax %<>% as.matrix() 


# 1.2.1 构建phyloseq对象

physeq <- phyloseq(
  otu_table(otu,taxa_are_rows = TRUE), 
  tax_table(tax), # 数据格式必须是矩阵，否则会改变行列名
  sample_data(cov)
)


# 1.2.2 数据标准化-抽平
##?normalize # 查看数据标准化方法对应字符缩写。
set.seed(12345)
data <- physeq
data

# 2.1.1 lefse分析-LDA阈值设置为4

lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "Group",
  taxa_rank = "none",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.1, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.1, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 30, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。

# 2.1.2 提取biomarker鉴定结果
##没有注释名的分类单元，将会按照前几个分类单元的注释结果进行命名。
lefse %>% marker_table() -> res.diff



dim(res.diff) # 61个biomarker。
res.diff %>% head()

result = data.frame(feature=res.diff$feature,enrich_group=res.diff$enrich_group,ef_lda=res.diff$ef_lda,pvalue=res.diff$pvalue,padj=res.diff$padj)
last_char <- sapply(strsplit(result$feature, "_"), function(x) x[length(x)])
last_num <- as.numeric(last_char)
result$feature <- col_names[last_num]
write.csv(result,file="./lefse/HC_NBD_p0.1.csv")

result$name <- sapply(strsplit(result$feature, ";"), function(x) x[length(x)])



####### HC VS BD

#otu <- as.data.frame(t(read.csv('otu_table_genus.csv',row.names = 1)))
otu <- as.data.frame(t(read.csv('./genus_select_FDR/73otu_result.csv',row.names = 1)))

col_names <- colnames(otu)

# 分割每个列名
split_names <- strsplit(col_names, ";")
# 将list中的vector转化为character
id_chars <- paste("otu_ID", 1:ncol(otu), sep='_')
# 将结果转换为dataframe
tax <- data.frame(
  rbind(do.call(rbind, split_names))
)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax) <- id_chars
colnames(otu) <- id_chars


cov <- as.data.frame((read.csv('./lefse/Covariables.csv',row.names = 1)))
samples <- intersect(rownames(cov),rownames(otu))
cov <- cov[samples,]
otu  <- otu[samples,]
# Replace NAs with median in numeric columns
numeric_columns <- sapply(cov, is.numeric)
cov[, numeric_columns] <- lapply(cov[, numeric_columns], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

otu <- t(otu)

cov <- cov[(cov$Group != "BD1"),]
otu <- otu[,colnames(otu) %in% rownames(cov)]
cov$Group <- ifelse(cov$Group %in% c("BD2"), "BD", cov$Group)


otu %<>% as.matrix() 
tax %<>% as.matrix() 


# 1.2.1 构建phyloseq对象

physeq <- phyloseq(
  otu_table(otu,taxa_are_rows = TRUE), 
  tax_table(tax), # 数据格式必须是矩阵，否则会改变行列名
  sample_data(cov)
)


# 1.2.2 数据标准化-抽平
##?normalize # 查看数据标准化方法对应字符缩写。
set.seed(12345)
data <- physeq
data

# 2.1.1 lefse分析-LDA阈值设置为4

lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "Group",
  taxa_rank = "none",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.1, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.1, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 30, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。

# 2.1.2 提取biomarker鉴定结果
##没有注释名的分类单元，将会按照前几个分类单元的注释结果进行命名。
lefse %>% marker_table() -> res.diff



dim(res.diff) # 61个biomarker。
res.diff %>% head()

result = data.frame(feature=res.diff$feature,enrich_group=res.diff$enrich_group,ef_lda=res.diff$ef_lda,pvalue=res.diff$pvalue,padj=res.diff$padj)
last_char <- sapply(strsplit(result$feature, "_"), function(x) x[length(x)])
last_num <- as.numeric(last_char)
result$feature <- col_names[last_num]
write.csv(result,file="./lefse/HC_BD_p0.1.csv")
result$name <- sapply(strsplit(result$feature, ";"), function(x) x[length(x)])



####### BD VS NBD

#otu <- as.data.frame(t(read.csv('otu_table_genus.csv',row.names = 1)))
otu <- as.data.frame(t(read.csv('./genus_select_FDR/73otu_result.csv',row.names = 1)))

col_names <- colnames(otu)

# 分割每个列名
split_names <- strsplit(col_names, ";")
# 将list中的vector转化为character
id_chars <- paste("otu_ID", 1:ncol(otu), sep='_')
# 将结果转换为dataframe
tax <- data.frame(
  rbind(do.call(rbind, split_names))
)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax) <- id_chars
colnames(otu) <- id_chars


cov <- as.data.frame((read.csv('./lefse/Covariables.csv',row.names = 1)))
samples <- intersect(rownames(cov),rownames(otu))
cov <- cov[samples,]
otu  <- otu[samples,]
# Replace NAs with median in numeric columns
numeric_columns <- sapply(cov, is.numeric)
cov[, numeric_columns] <- lapply(cov[, numeric_columns], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

otu <- t(otu)

cov <- cov[(cov$Group != "HC"),]
otu <- otu[,colnames(otu) %in% rownames(cov)]
cov$Group <- ifelse(cov$Group %in% c("BD2"), "BD", cov$Group)
cov$Group <- ifelse(cov$Group %in% c("BD1"), "NBD", cov$Group)


otu %<>% as.matrix() 
tax %<>% as.matrix() 


# 1.2.1 构建phyloseq对象

physeq <- phyloseq(
  otu_table(otu,taxa_are_rows = TRUE), 
  tax_table(tax), # 数据格式必须是矩阵，否则会改变行列名
  sample_data(cov)
)


# 1.2.2 数据标准化-抽平
##?normalize # 查看数据标准化方法对应字符缩写。
set.seed(12345)
data <- physeq
data

# 2.1.1 lefse分析-LDA阈值设置为4

lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "Group",
  taxa_rank = "none",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.1, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.1, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 30, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。

# 2.1.2 提取biomarker鉴定结果
##没有注释名的分类单元，将会按照前几个分类单元的注释结果进行命名。
lefse %>% marker_table() -> res.diff



dim(res.diff) # 61个biomarker。
res.diff %>% head()

result = data.frame(feature=res.diff$feature,enrich_group=res.diff$enrich_group,ef_lda=res.diff$ef_lda,pvalue=res.diff$pvalue,padj=res.diff$padj)
last_char <- sapply(strsplit(result$feature, "_"), function(x) x[length(x)])
last_num <- as.numeric(last_char)
result$feature <- col_names[last_num]
write.csv(result,file="./lefse/BD_NBD_p0.1.csv")
result$name <- sapply(strsplit(result$feature, ";"), function(x) x[length(x)])



#画图

#物种分类层级图 需去掉taxa_rank，或指定 taxa_rank="all"
lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "Group",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.1, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.1, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 30, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。


cols <- RColorBrewer::brewer.pal(8, "Dark2")
pdf("./lefse/cladogram_lda_marke_only.pdf",
    width = unit(16,"cm"),
    height = unit(16,"cm"),
    family="Times")
plot_cladogram(
  lefse, 
  # 颜色数需要与enrich_group分类水平数保持一致。#
  color = cols[seq_along(res.diff$enrich_group %>% unique)],
  only_marker = TRUE,# FALSE图中展示所有数据，TRUE表示只展示差异显著数据。
  branch_size = 0.2, # 分支粗细
  alpha = 0.2,# 阴影颜色透明度
  node_size_scale = 1, # 节点大小
  node_size_offset = 1.1,
  clade_label_level = 7, # 添加分支标签的最大分类单元，其余分类单元注释作为图例，数值越大代表的分类单元越高。
  clade_label_font_size = 2,# 分支标签大小
  annotation_shape = 22,
  annotation_shape_size = 2,
  marker_legend_param = list(
    ncol = 2, # 图例排2列。
    direction = "horizontal" # 水平放置。
  )
) +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "bottom"
  )
dev.off()




####ggplot画条形图 lda score

files <- list.files(path = "./lefse", pattern = "p0\\.05\\.csv$", full.names = TRUE)
# 创建颜色映射列表
color_mapping <- c("HC" = "#00c16e", "BD" = "#037ef3", 
                   "NBD" = "#f85a40", "CD" = "#7552cc")
for(i in seq_along(files)) {
  # 使用read.csv读取文件，然后添加到data_list
  result <- as.data.frame(read.csv(files[i],row.names = 1))
  result$name <- sapply(strsplit(result$feature, ";"), function(x) x[length(x)])
  
  # 获取group列的第二个值
  second_group_value <- unique(result$enrich_group)[2]
  # 使用ifelse()将group列为第二个值的value列取负
  result$ef_lda <- ifelse(result$enrich_group == second_group_value, -result$ef_lda, result$ef_lda)
  # 对name进行排序
  result <- result[order(result$ef_lda, decreasing = FALSE), ]
  # 将name列改为因子，并按照ef_lda的排序顺序设置因子水平
  result$name <- factor(result$name, levels = unique(result$name))
  plot <- ggplot(result,aes(name,ef_lda))+
    geom_col(aes(fill=enrich_group))+
    scale_fill_manual(values=color_mapping)+#自定义颜色
    coord_flip() +
    labs(x = "Genus", y = "LDA Score", fill = "Enrich Group")
  plot_name <- sub(pattern = "\\.csv$", replacement = "_lda.pdf", x = files[i])
  ggsave(filename = plot_name, plot = plot, device = "pdf", width = 10, height = 7)
}

files <- list.files(path = "./lefse", pattern = "p0\\.1\\.csv$", full.names = TRUE)
for(i in seq_along(files)) {
  # 使用read.csv读取文件，然后添加到data_list
  result <- as.data.frame(read.csv(files[i],row.names = 1))
  result$name <- sapply(strsplit(result$feature, ";"), function(x) x[length(x)])
  # 获取group列的第二个值
  second_group_value <- unique(result$enrich_group)[2]
  # 使用ifelse()将group列为第二个值的value列取负
  result$ef_lda <- ifelse(result$enrich_group == second_group_value, -result$ef_lda, result$ef_lda)
  # 对name进行排序
  result <- result[order(result$ef_lda, decreasing = FALSE), ]
  # 将name列改为因子，并按照ef_lda的排序顺序设置因子水平
  result$name <- factor(result$name, levels = unique(result$name))
  plot <- ggplot(result,aes(name,ef_lda))+
    geom_col(aes(fill=enrich_group))+
    scale_fill_manual(values=color_mapping)+#自定义颜色
    coord_flip() +
    labs(x = "Genus", y = "LDA Score", fill = "Enrich Group")
  plot_name <- sub(pattern = "\\.csv$", replacement = "_lda.pdf", x = files[i])
  ggsave(filename = plot_name, plot = plot, device = "pdf", width = 10, height = 7)
}

