# MAPFA
Metagenomic Analysis Pipeline for Food Animals

## 流程

### 数据处理部分

原始数据的质控，过滤，过滤后质控，组装，组装结果评估，分箱（三个分箱算法），分箱结果的优化，分箱结果评估。

### 数据分析部分

1. 物种注释（kraken2）
2. ORF注释（metagenemark，cd-hit）
3. 基因丰度估计salmon
4. 毒力因子注释（usearch比对至VFDB数据库）
5. 耐药基因注释（diamond比对至可选的多个耐药数据库（CARD，Deep-ARG））
6. 插入序列注释（usearch比对）
7. 原噬菌体、质粒注释（usearch比对）
8. 共现关系（耐药基因与毒力因子，可移动遗传原件）
9. 功能注释（GO，KEGG）
10. 对高质量分箱结果分析，单菌注释，prokka预测基因

