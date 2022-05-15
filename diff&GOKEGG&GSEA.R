DEG_analysis <- function(countdata, clinicaldata, target_gene, log2FC = 1, 
                         signif_DEG = 'p.adj',DEG_pvaluecutoff = 0.05,DEG_qvaluecutoff = 0.05,
                         signif_GO = 'p.adj&q.value',GO_pvaluecutoff = 0.05,GO_qvaluecutoff = 0.05,
                         GSEA_padjcutoff = 0.25
){
  if(!require(dplyr))install.packages('dplyr')
  if(!require(BiocManager))install.packages('BiocManager')
  if(!require(ggplot2))install.packages('ggplot2')
  if(!require(ggsci))install.packages('ggsci')
  if(!require(ggplotify))install.packages("ggplotify")
  if(!require(patchwork))install.packages("patchwork")
  if(!require(cowplot))install.packages("cowplot")
  if(!require(DESeq2))BiocManager::install('DESeq2')
  if(!require(edgeR))BiocManager::install('edgeR')
  if(!require(limma))BiocManager::install('limma')
  if(!require(clusterProfiler))BiocManager::install('clusterProfiler')
  target <- countdata[target_gene,] #选取目的基因，countdata为已标准化的数据
  target<-t(target)
  countdata_min <- countdata[-which(rownames(countdata)==target_gene),]
  clinicaldata <- merge(x = clinicaldata, y = target, by = 'row.names', all = F)
  row.names(clinicaldata) <- clinicaldata[, 1] #第一行设置为行名
  clinicaldata <- clinicaldata[, -1]
  a <- data.frame(matrix(NA,lengths(rownames(clinicaldata)),1))
  colnames(a) <- paste(target_gene,'_status')
  clinicaldata <- cbind(clinicaldata,a)
  clinicaldata[,paste(target_gene,'_status')][clinicaldata[,target_gene] < median(clinicaldata[,target_gene])]<- 'low'
  clinicaldata[,paste(target_gene,'_status')][clinicaldata[,target_gene] >= median(clinicaldata[,target_gene])] <- 'high'
  clinicaldata[,paste(target_gene,'_status')] <- as.factor(clinicaldata[,paste(target_gene,'_status')])
  # 构建limma分组
  design <- model.matrix(~0+factor(clinicaldata[,paste(target_gene,'_status')]))
  colnames(design)=levels(factor(clinicaldata[,paste(target_gene,'_status')]))
  rownames(design)=colnames(countdata_min)
  # contrast.matrix<-makeContrasts(paste0(levels(clinicaldata[,paste(target_gene,'_status')]),
  #   collapse = "-"),levels = design) 这一行我怎么运行都报错，没辙了换成下面这一行，固定了high-low
  contrast.matrix<-makeContrasts('high-low',levels = design)
  fit1 <- lmFit(countdata_min,design)
  fit2 <- contrasts.fit(fit1, contrast.matrix)
  fit2 <- eBayes(fit2)
  all_diff <- topTable(fit2, adjust.method = 'fdr',coef = 1, p.value = 1, lfc = log(1,2), number = 30000, sort.by = 'logFC')
  if (signif_DEG == 'p.adj'){
    DEG_signif <- 'adj.P.Val'
    DEG_cutoff <- DEG_pvaluecutoff
  } else if (signif_DEG=='p.value') { 
    DEG_signif <- 'P.Value'
    DEG_cutoff <- DEG_qvaluecutoff
  }
  all_diff$change = as.factor(
    ifelse(all_diff[,DEG_signif] < DEG_cutoff & abs(all_diff$logFC) > log2FC,
           ifelse(all_diff$logFC > log2FC ,'UP','DOWN'),'NOT'))
  # 火山图
  this_tile <- paste0('Cutoff for logFC is ',round(log2FC,3),
                      '\nThe number of up gene is ',nrow(all_diff[all_diff$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(all_diff[all_diff$change =='DOWN',])
  )
  
  vol_plot <- ggplot(data=all_diff, aes(x=logFC, y=-log10(all_diff[,DEG_signif]),color=change)) + geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+ xlab("log2 fold change") + ylab(paste("-log10 ",DEG_signif)) + 
    ggtitle( this_tile ) +   theme(plot.title = element_text(size=15,hjust = 0.5)) + 
    scale_colour_manual(values = c('blue','black','red')) 	
  
  # GO/KEGG富集分析
  DEG_diff <- all_diff[-which(all_diff$change=='NOT'),]
  library(org.Hs.eg.db)
  if (signif_GO == 'p.value'){
    go <- enrichGO(rownames(DEG_diff),OrgDb = org.Hs.eg.db, 
                   ont='ALL',pAdjustMethod = 'none',pvalueCutoff = GO_pvaluecutoff, 
                   qvalueCutoff = 1,keyType = 'SYMBOL')
  } else if (signif_GO == 'p.adj'){
    go <- enrichGO(rownames(DEG_diff),OrgDb = org.Hs.eg.db, 
                   ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = GO_pvaluecutoff, 
                   qvalueCutoff = 1,keyType = 'SYMBOL')
  } else if (signif_GO == 'q.value') {
    go <- enrichGO(rownames(DEG_diff),OrgDb = org.Hs.eg.db, 
                   ont='ALL',pAdjustMethod = 'none',pvalueCutoff = 1, 
                   qvalueCutoff = GO_qvaluecutoff,keyType = 'SYMBOL')
  } else if (signif_GO == 'p.adj&q.value'){
    go <- enrichGO(rownames(DEG_diff),OrgDb = org.Hs.eg.db, 
                   ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = GO_pvaluecutoff, 
                   qvalueCutoff = GO_qvaluecutoff,keyType = 'SYMBOL')
  }
  
  # KEGG
  id_ENTREZID <- bitr(rownames(DEG_diff), fromType = 'SYMBOL', toType = 'ENTREZID',
                      OrgDb = 'org.Hs.eg.db')
  if (signif_GO == 'p.value'){
    kegg <- enrichKEGG(id_ENTREZID$ENTREZID, organism = 'hsa', keyType = 'kegg', 
                       pvalueCutoff = GO_pvaluecutoff, pAdjustMethod = 'none', minGSSize = 3, maxGSSize = 500, 
                       qvalueCutoff = 1, use_internal_data = FALSE)
    func <- kegg@result %>% dplyr::filter(pvalue <=GO_pvaluecutoff)
    if (nrow(func)==0){
      go <- go
    } else {
      func$ONTOLOGY <- 'KEGG'
      go@result <- merge(x = go@result, y = func, all = T)
    }
  } else if (signif_GO == 'p.adj'){
    kegg <- enrichKEGG(id_ENTREZID$ENTREZID, organism = 'hsa', keyType = 'kegg', 
                       pvalueCutoff = GO_pvaluecutoff, pAdjustMethod = 'BH', minGSSize = 3, maxGSSize = 500, 
                       qvalueCutoff = 1, use_internal_data = FALSE)
    func <- kegg@result %>% dplyr::filter(p.adjust <= GO_pvaluecutoff)
    if (nrow(func)==0){
      go <- go
    } else {
      func$ONTOLOGY <- 'KEGG'
      go@result <- merge(x = go@result, y = func, all = T)
    }
  } else if (signif_GO == 'q.value') {
    kegg <- enrichKEGG(id_ENTREZID$ENTREZID, organism = 'hsa', keyType = 'kegg', 
                       pvalueCutoff = 1, pAdjustMethod = 'none', minGSSize = 3, maxGSSize = 500, 
                       qvalueCutoff = GO_qvaluecutoff, use_internal_data = FALSE)
    func <- kegg@result %>% dplyr::filter(qvalue <= GO_qvaluecutoff)
    if (nrow(func)==0){
      go <- go
    } else {
      func$ONTOLOGY <- 'KEGG'
      go@result <- merge(x = go@result, y = func, all = T)
    }
  } else if (signif_GO == 'p.adj&q.value'){
    kegg <- enrichKEGG(id_ENTREZID$ENTREZID, organism = 'hsa', keyType = 'kegg', 
                       pvalueCutoff = GO_pvaluecutoff, pAdjustMethod = 'BH', minGSSize = 3, maxGSSize = 500, 
                       qvalueCutoff = GO_qvaluecutoff, use_internal_data = FALSE)
    func <- kegg@result %>% dplyr::filter(p.adjust <= GO_pvaluecutoff & qvalue <= GO_qvaluecutoff)
    if (nrow(func)==0){
      go <- go
    } else {
      func$ONTOLOGY <- 'KEGG'
      go@result <- merge(x = go@result, y = func, all = T)
    }
  }
  
  GO_plot <- barplot(go, showCategory = nrow(go@result), split="ONTOLOGY",color = "qvalue")+ facet_grid(ONTOLOGY~.,scale="free",space = "free_y")
  
  # 获取相关性表
  correlation <- data.frame()
  genelist <- rownames(countdata)
  gene <- target_gene
  genedata <- as.numeric(countdata[gene,])
  for(i in 1:length(genelist)){
    dd = cor.test(genedata,as.numeric(countdata[i,]),method="spearman")
    correlation[i,1] = gene
    correlation[i,2] = genelist[i]
    correlation[i,3] = dd$estimate
    correlation[i,4] = dd$p.value
  }
  
  colnames(correlation) <- c("gene1","gene2","cor","p.value")
  rownames(correlation) <- correlation[,2]
  correlation <- correlation[,-(1:2)]
  correlation <- correlation[-which(rownames(correlation)==target_gene),]
  
  # 差异基因GSEA
  # 因为在线GSEA没有hallmark和biocarta，所以我用的是读取线下文件的方式
  geneList <- all_diff$logFC
  names(geneList) = rownames(all_diff)
  geneList = sort(geneList, decreasing = TRUE)
  hallmarks <- read.gmt("G:\\R\\RData\\TCGA GENE\\h.all.v7.5.1.symbols.gmt")
  biocarta <- read.gmt("G:\\R\\RData\\TCGA GENE\\c2.cp.biocarta.v7.5.1.symbols.gmt")
  cp.kegg <- read.gmt("G:\\R\\RData\\TCGA GENE\\c2.cp.kegg.v7.5.1.symbols.gmt")
  go.bp <- read.gmt("G:\\R\\RData\\TCGA GENE\\c5.go.bp.v7.5.1.symbols.gmt")
  go.cc <- read.gmt("G:\\R\\RData\\TCGA GENE\\c5.go.cc.v7.5.1.symbols.gmt")
  go.mf <- read.gmt("G:\\R\\RData\\TCGA GENE\\c5.go.mf.v7.5.1.symbols.gmt")
  
  set.seed(1)
  y1 <- GSEA(geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = hallmarks)
  y2 <- GSEA(geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = biocarta)
  y3 <- GSEA(geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = cp.kegg)
  y4 <- GSEA(geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = go.bp)
  y5 <- GSEA(geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = go.cc)
  y6 <- GSEA(geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = go.mf)
  gsea_res <- y1
  if (nrow(y2@result > 0)){
    gsea_res@result <- merge(x = gsea_res@result, y = y2@result, all = T)
  } else {
    gsea_res <- gsea_res
  }
  if (nrow(y3@result > 0)){
    gsea_res@result <- merge(x = gsea_res@result, y = y3@result, all = T)
  } else {
    gsea_res <- gsea_res
  }
  if (nrow(y4@result > 0)){
    gsea_res@result <- merge(x = gsea_res@result, y = y4@result, all = T)
  } else {
    gsea_res <- gsea_res
  }
  if (nrow(y5@result > 0)){
    gsea_res@result <- merge(x = gsea_res@result, y = y5@result, all = T)
  } else {
    gsea_res <- gsea_res
  }
  if (nrow(y6@result > 0)){
    gsea_res@result <- merge(x = gsea_res@result, y = y6@result, all = T)
  } else {
    gsea_res <- gsea_res
  }
  
  deg_gseaplot <- dotplot(gsea_res,showCategory=nrow(gsea_res@result),split=".sign")+facet_grid(~.sign)
  
  # 基因相关性GSEA
  cor_geneList <- correlation$cor
  names(cor_geneList) = rownames(correlation)
  cor_geneList = sort(cor_geneList, decreasing = TRUE)
  set.seed(1)
  x1 <- GSEA(cor_geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = hallmarks)
  x2 <- GSEA(cor_geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = biocarta)
  x3 <- GSEA(cor_geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = cp.kegg)
  x4 <- GSEA(cor_geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = go.bp)
  x5 <- GSEA(cor_geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = go.cc)
  x6 <- GSEA(cor_geneList,minGSSize = 15,maxGSSize = 200,pvalueCutoff = GSEA_padjcutoff,TERM2GENE = go.mf)
  cor_gsea <- x1
  if (nrow(x2@result > 0)){
    cor_gsea@result <- merge(x = cor_gsea@result, y = x2@result, all = T)
  } else {
    cor_gsea <- cor_gsea
  }
  if (nrow(x3@result > 0)){
    cor_gsea@result <- merge(x = cor_gsea@result, y = x3@result, all = T)
  } else {
    cor_gsea <- cor_gsea
  }
  if (nrow(x4@result > 0)){
    cor_gsea@result <- merge(x = cor_gsea@result, y = x4@result, all = T)
  } else {
    cor_gsea <- cor_gsea
  }
  if (nrow(x5@result > 0)){
    cor_gsea@result <- merge(x = cor_gsea@result, y = x5@result, all = T)
  } else {
    cor_gsea <- cor_gsea
  }
  if (nrow(x6@result > 0)){
    cor_gsea@result <- merge(x = cor_gsea@result, y = x6@result, all = T)
  } else {
    cor_gsea <- cor_gsea
  }
  
  corgene_gseaplot <- dotplot(cor_gsea,showCategory=nrow(cor_gsea@result),split=".sign")+facet_grid(~.sign)
  
  return(list(all_diff,vol_plot,GO_plot,correlation,deg_gseaplot,corgene_gseaplot))
}