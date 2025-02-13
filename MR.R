
library(TwoSampleMR)
library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)
vcf_outcome <- readVcf("t2d.vcf.gz")
outcome <- gwasvcf_to_TwoSampleMR(vcf = vcf_outcome)
head(outcome)


library(data.table)
eqtl=fread("eqtlgen_p_5e-08.csv.gz")
eqtl=data.frame(eqtl)
head(eqtl)


library(readxl)
gene <- read_excel("gene416.xlsx")
gene=gene$gene
gene=gene[!duplicated(gene)]


library(clusterProfiler) 
library(org.Hs.eg.db)
df=bitr(geneID = gene,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
head(df)
df <- read.table("df.csv",header=T, sep=",",row.names = 1) 
list_Probe=df$ENSEMBL



library(TwoSampleMR)
library(dplyr)
library(vroom)
library(VariantAnnotation)
library(gwasvcf)
library(gwasglue)
library(ieugwasr)
library(plinkbinr)
library(MRPRESSO)

result = NULL
heterogeneity = NULL
pleiotropy = NULL
exposure_datsave = NULL
outcome_datsave = NULL
harmonise_datsave = NULL

for (i in 1:length(list_Probe)) {
  print(i)
  exposure = subset(eqtl, Probe == list_Probe[i])
  tryCatch({
    exposure = format_data(exposure,
                           type ="exposure",
                           phenotype_col ="Probe",
                           snp_col ="SNP",
                           beta_col = "b",
                           se_col = "SE",
                           eaf_col = "Freq",
                           effect_allele_col = "A1",
                           other_allele_col = "A2",
                           pval_col = "p",
                           chr_col = "Chr",
                           pos_col = "BP")
    exposure$id.exposure = exposure$exposure
    exposure_clump = ieugwasr::ld_clump(dplyr::tibble(rsid = exposure$SNP,
                                                      pval = exposure$pval.exposure),
                                        clump_kb = 10000,
                                        clump_r2 = 0.001,
                                        clump_p = 1,
                                        bfile = "/home/data/t090502/MR/EUR",
                                        plink_bin = plinkbinr::get_plink_exe(),
                                        pop = "EUR")
    exposure_clump = exposure[exposure$SNP %in% exposure_clump$rsid,]
    exposure_datsave = rbind(exposure_datsave, exposure_clump)
    outcome_dat = subset( outcome, outcome$SNP %in% exposure_clump$SNP)
    outcome_dat = format_data(dat = outcome_dat,
                              type = "outcome",
                              snp_col = "SNP",
                              beta_col = "beta.exposure",
                              pval_col = "pval.exposure",
                              se_col = "se.exposure",
                              eaf_col = "eaf.exposure",
                              effect_allele_col = "effect_allele.exposure",
                              other_allele_col = "other_allele.exposure")
    outcome_dat$id.outcome = "T2D"
    outcome_datsave = rbind(outcome_datsave, outcome_dat)
    mr_data = harmonise_data(exposure_dat =exposure_clump,
                             outcome_dat = outcome_dat,
                             action= 2)  
    harmonise_datsave = rbind(harmonise_datsave, mr_data)
    mr_res = mr(mr_data)
    mr_res = generate_odds_ratios(mr_res)
    result = rbind(result, mr_res)
    het= mr_heterogeneity(mr_data)
    heterogeneity = rbind(heterogeneity, het)
    pleio = mr_pleiotropy_test(mr_data)
    pleiotropy = rbind(pleiotropy, pleio)
  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
  })
}

head(result)

result$ENSEMBL =result$id.exposure
result2= left_join(result,df, by = "ENSEMBL")
write.csv(result2, file = "MR2eqtlgen_mr_res_ebi-a-GCST005180-416.csv")
write.csv(heterogeneity, file = "MR2eqtlgen_heterogeneity_ebi-a-GCST005180-416.csv")
write.csv(pleiotropy, file = "MR2eqtlgen_pleiotropy_ebi-a-GCST005180-416.csv")

res = result2[result2$method %in% c("Wald ratio","Inverse variance weighted"),]
res = res[order(res$pval),]
res$FDR = p.adjust(res$pval, method ="BH")
res = res[res$pval <0.05,]
het = heterogeneity[heterogeneity$id.exposure %in% res$id.exposure,]
pleio = pleiotropy[pleiotropy$id.exposure %in% res$id.exposure,]
id_1 = het$id.exposure[het$Q_pval <0.05]
id_2 = pleio$id.exposure[pleio$Q_pval <0.05]
id = unique(c(id_1,id_2))
res =res[!res$id.exposure %in% id,]  #!是去掉的意思
het =het[!het$id.exposure %in% id,]
pleio =pleio[!pleio$id.exposure %in% id,]
write.csv(res, file = "MR2eqtlgen_mr_res_ebi-a-GCST005180-416—0.05.csv")
write.csv(pleio, file = "MR2eqtlgen_pleiotropy_ebi-a-GCST005180-416—0.05.csv")
write.csv(het, file = "MR2eqtlgen_heterogeneity_ebi-a-GCST005180-416—0.05.csv")


head(exposure_datsave)
exposure_datsave$ENSEMBL =exposure_datsave$exposure
exposure_dat= left_join(exposure_datsave,df, by = "ENSEMBL")
exposure_dat$samplesize.exposure = 31684 
R2a=2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2b=2*exposure_dat$se.exposure*exposure_dat$se.exposure*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2=R2a/(R2a+R2b)
exposure_dat$F_statistics=R2*(exposure_dat$samplesize.exposure-2)/(1-R2)
write.csv(exposure_dat, file = "MR2eqtlgen_gene-instruments-416.csv")
exposure_dat <- read.table("MR2eqtlgen_gene-instruments-416.csv",header=T, sep=",",row.names = 1) 
save(exposure_dat,outcome_datsave,result2,harmonise_datsave,file ='MR2eqtlgen_mr_input_res.Rdata')
load('MR2eqtlgen_mr_input_res.Rdata')


table1 <- res %>% 
  left_join(exposure_dat, by = "exposure")

table1 <- table1 %>% 
  generate_odds_ratios()%>% 
  mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",or,or_lci95, or_uci95),
         `P value` = scales::scientific(pval),
         `F statistics` = sprintf("%.2f",F_statistics)) %>% 
  dplyr::select(Gene = SYMBOL.y, `ENSEMBL ID` = ENSEMBL.y,
                SNP, `Effect allele` = effect_allele.exposure, 
                `OR (95% CI)`, `P value`, 
                 `F statistics`)

save(table1,file ='MR2eqtlgen_table1.Rdata')
load('MR2eqtlgen_table1.Rdata')


library(grid)
library(forestploter)

mydata=read.table("MR2eqtlgen_mr_res_ebi-a-GCST005180-416—0.05.csv",header=T, sep=",",row.names = 1)
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.3f (%.3f - %.3f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
mydata$pval=scales::scientific(mydata$pval)
mydata$exposure=mydata$SYMBOL
head(mydata)
p=forest(mydata[,c(4:6,9,18,19)],
         est = mydata$or,
         lower =mydata$or_lci95, 
         upper = mydata$or_uci95,
         sizes =0.3,
         ci_column =5 ,
         ref_line = 1,
         xlim = c(0.05, 3),
)
ggsave("MR2eqtlgen_forest416.pdf",p,width = 10, height = 12)



library(TwoSampleMR)

bimr_T2D<- extract_instruments('ebi-a-GCST006867', p1=5e-08, clump=TRUE)
write.csv(bimr_T2D,file ='bimr_T2D.csv',quote = F)
bimr_T2D <- read.table("bimr_T2D.csv",header=T, sep=",",row.names = 1)
head(bimr_T2D)

library(TwoSampleMR)
library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)
vcf_outcome <- readVcf("eqtl-a-ENSG00000146535.vcf.gz")
outcome <- gwasvcf_to_TwoSampleMR(vcf = vcf_outcome)
outcome_dat = subset( outcome, outcome$SNP %in% bimr_T2D$SNP)
outcome_dat = format_data(dat = outcome_dat,
                          type = "outcome",
                          snp_col = "SNP",
                          beta_col = "beta.exposure",
                          pval_col = "pval.exposure",
                          se_col = "se.exposure",
                          eaf_col = "eaf.exposure",
                          effect_allele_col = "effect_allele.exposure",
                          other_allele_col = "other_allele.exposure")
outcome_dat$id.outcome = "GNA12"
mr_data = harmonise_data(exposure_dat =bimr_T2D,
                         outcome_dat = outcome_dat,
                         action= 2)  #等于2标记回文snp，等于1不标记，默认等于2
mr_res = mr(mr_data)
mr_res = generate_odds_ratios(mr_res)
biMR_result = rbind(biMR_result, mr_res)
write.csv(biMR_result, file = "MR2eqtlgen_biMR_result13.csv")




library(TwoSampleMR)
library(dplyr)
load('MR2eqtlgen_mr_input_res.Rdata')
load('MR2eqtlgen_table1.Rdata')
harmonise_datsave$samplesize.outcome = 655666
harmonise_datsave$samplesize.exposure = 31684
library(readxl)
outcome_gene <- read_excel("bi_MR_gene.xlsx")
outcome_gene$Gene =outcome_gene$bi_gene
head(outcome_gene)
outcome_gene= left_join(outcome_gene, table1, by = "Gene")


steiger_res= harmonise_datsave %>%
  filter(SNP %in% outcome_gene$SNP) %>% 
  steiger_filtering()
steiger_res$steiger_pval




library(TwoSampleMR)
library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)

outcome <- read.table("outcome.csv",header=T, sep=",",row.names = 1)
load('MR2eqtlgen_mr_input_res.Rdata')

data <- vcfR::read.vcfR("eqtl-a-ENSG00000112742.vcf")  
gt=data@gt
gt=as.data.frame(gt)

colnames(gt)
gt$FORMAT[1]

library(tidyverse)

gt$`eqtl-a-ENSG00000112742`[1] 

##!!!
gt=separate(gt,col='eqtl-a-ENSG00000112742',into = c('ES', 'SE',
                                                     'LP','AF','SS',
                                                     'ID'),sep = '\\:')

gc()



gt=na.omit(gt)
colnames(gt)=c('format','beta','se','logpvalue','eaf','samplesize','snp')
gt$beta=as.numeric(gt$beta)
gt$se=as.numeric(gt$se)
gt$logpvalue=as.numeric(gt$logpvalue)
gt$eaf=as.numeric(gt$eaf)
gt$samplesize=as.numeric(gt$samplesize)

gc()
gt$format=NULL


fix=data@fix
fix=as.data.frame(fix)
colnames(fix)
colnames(fix)=c('chr','pos','snp','ref','alt')
fix=fix[,1:5]



eqtl=left_join(fix,gt,by='snp')
eqtl=na.omit(eqtl)
eqtl$maf = ifelse(eqtl$eaf < 0.5, 
                  eqtl$eaf,
                  1 - eqtl$eaf)
eqtl$eaf=NULL
eqtl=eqtl[eqtl$chr==6,]  
eqtl$logpvalue=as.numeric(eqtl$logpvalue)
eqtl$p_value=10^(-eqtl$logpvalue)
eqtl$pos=as.numeric(eqtl$pos)
eqtl=eqtl[eqtl$pos >80760668 -1000000 ,]
eqtl=eqtl[eqtl$pos <80760668 +1000000 ,]

my_eqtl=eqtl[,c('snp','p_value','maf')]

colnames(my_eqtl)=c('snp','pvalues','MAF')
my_eqtl=na.omit(my_eqtl)

my_eqtl=my_eqtl[my_eqtl$MAF>0 ,]


library(TwoSampleMR)

outcome_dat = subset( outcome, outcome$SNP %in% my_eqtl$snp)
outcome_dat = format_data(dat = outcome_dat,
                          type = "outcome",
                          snp_col = "SNP",
                          pos_col = "pos.exposure",
                          chr_col = "chr.exposure",
                          beta_col = "beta.exposure",
                          pval_col = "pval.exposure",
                          se_col = "se.exposure",
                          eaf_col = "eaf.exposure",
                          effect_allele_col = "effect_allele.exposure",
                          other_allele_col = "other_allele.exposure")
coloc_SLE_dat <- outcome_dat %>% 
  mutate(chr.outcome = as.numeric(chr.outcome),
         pos.outcome = as.numeric(pos.outcome),
         outcome = "Type 2 diabetes", 
         id.outcome = "ebi-a-GCST006867")


gwas=coloc_SLE_dat
gwas$beta=as.numeric(gwas$beta.outcome)
gwas$se=as.numeric(gwas$se.outcome)
gwas$varbeta=(gwas$se)^2
gwas=gwas[,c('SNP','pval.outcome',"beta",'varbeta')]
colnames(gwas)=c('snp','pvalues','beta','varbeta')
gwas=na.omit(gwas)
library(coloc)
input <- merge(my_eqtl, gwas, by="snp", all=FALSE, suffixes=c("_eqtl","_gwas"))

library(coloc)
result<-coloc.abf(dataset1=list(pvalues=input$pvalues_gwas, type="cc", s=0.09, N=655666,snp=input$snp),
                  dataset2=list(pvalues=input$pvalues_eqtl, type="quant", N=31684,snp=input$snp),
                  MAF=input$MAF)
need_result<-result$results %>% filter(SNP.PP.H4>0.8)



library(locuscomparer)
snp_overlap <- intersect(my_eqtl[['snp']],
                         gwas[['snp']])
snp_overlap <- unique(snp_overlap)

exposure_dat <-my_eqtl[my_eqtl[['snp']] %in% snp_overlap,]
outcome_dat <- gwas[gwas[['snp']] %in% snp_overlap,]

exposure_dat <- exposure_dat[order(exposure_dat[['snp']]),]
outcome_dat <- outcome_dat[order(outcome_dat[['snp']]),]

exposure_dat <- data.frame(
  rsid = exposure_dat[['snp']],
  pval = exposure_dat[['pvalues']])

outcome_dat <- data.frame(
  rsid = outcome_dat[['snp']],
  pval = outcome_dat[['pvalues']])


library(locuscomparer)
p <- locuscompare(in_fn1 ="coloc_outcome_TTK.tsv", 
                  in_fn2 ="coloc_exposure_TTK.tsv", 
                  title1 = paste0('T2D'," GWAS"),  #！！！ 改改
                  title2 = paste0('TTK'," eQTL"), # !!!!改改
                  combine =F,
                  legend = F
)$locuscompare
ggsave("TTK.png",p,width = 4, height = 4)










