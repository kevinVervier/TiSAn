library(GenomicRanges)

# create eQTL db
x = read.table('/sdata/GTex_eQTl/GTex_eQTL/V6/GTEx_Analysis_V6_eQTLs/Brain_snps.db',header=TRUE)

brain.eqtl = GRanges(seqnames = Rle(x$snp_chrom),
                        ranges = IRanges(start= x$snp_pos, end=x$snp_pos),
                        rsid = x$rs_id_dbSNP142_GRCh37p13,
                        gene = Rle(x$gene)
)
#newStyle <- mapSeqlevels(seqlevels(brain.eqtl), "UCSC")
#brain.eqtl <- renameSeqlevels(brain.eqtl, newStyle)
brain.eqtl <- sortSeqlevels(brain.eqtl)
brain.eqtl <- sort(brain.eqtl)
# remove duplicated entries found in more than 1 tissue
brain.eqtl <- brain.eqtl[-which(duplicated(brain.eqtl))] # 380,059
save(brain.eqtl,file='brain_eqtl_gr.Rdata')
#saveRDS(brain.eqtl,file='brain_eqtl_gr.rds')

# create eQTL db
x = read.table('/sdata/GTex_eQTl/GTex_eQTL/V6/GTEx_Analysis_V6_eQTLs/non_Brain_snps.db',header=TRUE)

other.eqtl = GRanges(seqnames = Rle(x$snp_chrom),
                     ranges = IRanges(start= x$snp_pos, end=x$snp_pos),
                        rsid = x$rs_id_dbSNP142_GRCh37p13,
                       gene = Rle(x$gene)
)

other.eqtl <- sortSeqlevels(other.eqtl)
other.eqtl <- sort(other.eqtl)
# remove duplicated entries found in more than 1 tissue
other.eqtl <- other.eqtl[-which(duplicated(other.eqtl))] # 380,059
save(other.eqtl,file='non_brain_eqtl_gr.Rdata')
#saveRDS(other.eqtl,file='non_brain_eqtl_gr.rds')


################
# brain genes db
x = read.table('~/varann/aim1/output/brain-gene-analysis/brain_gene.db',header = TRUE)
brain.gene = GRanges(seqnames = Rle(x$chromosome_name),
                     ranges = IRanges(start= x$start_position, end=x$end_position),
                     symbol = x$hgnc_symbol
                     #  gene = Rle(x$gene)
)
#newStyle <- mapSeqlevels(seqlevels(brain.eqtl), "UCSC")
#brain.eqtl <- renameSeqlevels(brain.eqtl, newStyle)
brain.gene <- sortSeqlevels(brain.gene)
brain.gene <- sort(brain.gene)
save(brain.gene,file='brain_gene_gr.Rdata')

################
# non-brain genes db
x = read.table('~/varann/aim1/output/brain-gene-analysis/non-brain_gene.db',header = TRUE)
other.gene = GRanges(seqnames = Rle(x$chromosome_name),
                     ranges = IRanges(start= x$start_position, end=x$end_position),
                     symbol = x$hgnc_symbol
                     #  gene = Rle(x$gene)
)
#newStyle <- mapSeqlevels(seqlevels(brain.eqtl), "UCSC")
#brain.eqtl <- renameSeqlevels(brain.eqtl, newStyle)
other.gene <- sortSeqlevels(other.gene)
other.gene <- sort(other.gene)
save(other.gene,file='non_brain_gene_gr.Rdata')

#######
# RME

x = read.table('/sdata/ChIP-seq/RMEpigenomics/WGBS_DMRs_v2_brain_average.txt',header=TRUE)
y = read.table('/sdata/ChIP-seq/RMEpigenomics/WGBS_DMRs_v2_other_average.txt',header=TRUE)
  
rme = GRanges(seqnames = Rle(x$chr),
                     ranges = IRanges(start= x$start, end=x$end),
                     brain.methyl = round(x$average_methylation,3),
                     other.methyl = round(y$average_methylation,3)
)
rme <- sortSeqlevels(rme)
rme <- sort(rme)
save(rme,file='rme_gr.Rdata')


##################################
# Fetal brain DDMRs

x = read.csv('~/varann/aim1/data/Bray2015.csv',header=TRUE)
ddmr = GRanges(seqnames = Rle(x$chrom),
              ranges = IRanges(start= x$start, end=x$end)
)
ddmr <- sortSeqlevels(ddmr)
ddmr <- sort(ddmr)
save(ddmr,file='ddmr_gr.Rdata')



######################
# Mongo DB version
# 
# library(rjson)
# x = read.table('/sdata/GTex_eQTl/GTex_eQTL/V6/GTEx_Analysis_V6_eQTLs/Brain_snps.db',header=TRUE)
# x = data.frame(x)
# exportJson <- rjson::toJSON(x[1:10,])
# exportJson
# 
# write(exportJson, "/sdata/GTex_eQTl/GTex_eQTL/V6/GTEx_Analysis_V6_eQTLs/Brain_snps.json")
# 
# library(mongolite)
# dmd <- mongo("diamonds")
# dmd$import(url("http://jeroen.github.io/data/diamonds.json"))
# dmd$count()
# dmd$drop()
# 
# brain.eqtl.db <- mongo("brain.eqtl")
# brain.eqtl.db$import(file('/sdata/GTex_eQTl/GTex_eQTL/V6/GTEx_Analysis_V6_eQTLs/Brain_snps.json'))
