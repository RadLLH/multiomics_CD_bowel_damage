setwd('G:\\BD_Genues_metabolites_MRE\\BD_code\\lefse')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiomeMarker")

library(tidyverse)
library(magrittr)
library(DescTools)
library(microbiomeMarker)


otu <- as.data.frame(t(read.csv('otu_table_genus.csv',row.names = 1)))

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


cov <- as.data.frame((read.csv('Covariables.csv',row.names = 1)))
samples <- intersect(rownames(cov),rownames(otu))
cov <- cov[samples,]
otu  <- otu[samples,]
# Replace NAs with median in numeric columns
numeric_columns <- sapply(cov, is.numeric)
cov[, numeric_columns] <- lapply(cov[, numeric_columns], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

otu <- t(otu)

otu %<>% as.matrix() 
tax %<>% as.matrix() 


# 1.2.1 构建phyloseq对象
library(phyloseq)
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
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.05, # 筛选阈值根据自己数据设置。
  kw_cutoff = 0.05, # 筛选阈值根据自己数据设置,此处进行bonferroni校正。
  bootstrap_n = 50, # LDA的bootstrap迭代次数。
  lda_cutoff = 2 # 筛选阈值根据自己数据设置。
) # 此函数没有提供差异检验p值校正参数，可以自己校正。

# 2.1.2 提取biomarker鉴定结果
##没有注释名的分类单元，将会按照前几个分类单元的注释结果进行命名。
res.diff <- lefse %>% marker_table() 
write.csv(res.diff,"lda4.diff.csv",quote = FALSE)
dim(res.diff) # 61个biomarker。
res.diff %>% head()


cols <- RColorBrewer::brewer.pal(8, "Dark2")
pdf("cladogram_lda_marke_only.pdf",
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




#intersection_HCBD_BDNBD
# 读取第一个CSV文件
data1 <- read.csv("BD1_BD2_diff73.csv")

# 读取第二个CSV文件
data2 <- read.csv("HC_BD2_diff73.csv")

# 读取第三个CSV文件
data3 <- read.csv("HC_CD_diff73.csv")


# 提取A列的值
values1 <- data1$feature
values2 <- data2$feature
values3 <- data3$feature
# 计算交集
# 计算交集
intersection <- Reduce(intersect, list(values1, values2, values3))


# 打印交集
print(intersection)

# 读取第三个CSV文件
data4 <- read.csv("otu_table_genus.csv")

# 筛选第三个文件中C列在交集之中的行
matching_rows <- data4[grep(paste(intersection, collapse="|"), data4$Taxon), ]

# 保存到新的CSV文件
write.csv(matching_rows, "select_73OTU_2_lefse_0.05_P.csv", row.names = FALSE)