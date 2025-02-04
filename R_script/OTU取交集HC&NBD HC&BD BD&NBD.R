##################OTU取交集HC&NBD HC&BD BD&NBD [P < 0.05,][LDA > 2,]
setwd('G:\\BD_Genues_metabolites_MRE\\BD_code\\lefse')
# 读取第一个CSV文件
data1 <- read.csv("HC_BD2_diff73.csv")

# 读取第二个CSV文件
data2 <- read.csv("BD1_BD2_diff73.csv")

# 读取第三个CSV文件
data3 <- read.csv("HC_BD1_diff73.csv")


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

data4 <- data4[sapply(data4$Taxon, function(x) any(sapply(intersection, function(y) grepl(y, x)))), ]
print(data4$Taxon)

# 筛选第三个文件中C列在交集之中的行
data4 <- data4[data4$Taxon %in% intersection, ]

# 保存到新的CSV文件
write.csv(data4, "filtered_genus_6_select_lefse_HC_BD1_BD2.csv", row.names = FALSE)
