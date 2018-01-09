require('data.table')
require(GenomicRanges)

#steps=1000000

#i=1
#dat = GRanges()


tmp = fread("zcat /sdata/vcfannotations/TiSAn_Brain_light.bed.gz") 
#tmp$V4 = round(tmp$V4,digits=1)
tmp = tmp[tmp$V4 > 0,]


#tmp2 = rep(0,nrow(tmp))

#for(i in 1:100000){
#  tmp$V4[(1+(i-1)*21200):(i*21200)] = round(tmp$V4[(1+(i-1)*2120):(i*2120)],digits=1)
#}

#tmp2$V4[(1+i*21200):nrow(tmp)] = round(tmp$V4[(1+i*21200):nrow(tmp)],digits=1)

#tmp$V4 = tmp2$V4
#rm(tmp2)
#gc()
#tmp = makeGRangesFromDataFrame(tmp,start.field = "V2",end.field = "V3",seqnames.field = "V1")


# while(TRUE){
#   tmp = fread("TiSAn_Brain.bed",nrows = steps,skip = steps*(i-1)+1) 
#   tmp$V4 = round(tmp$V4,digits=1)
tmp2 = makeGRangesFromDataFrame(tmp,keep.extra.columns = TRUE,start.field = "V2",end.field = "V3",seqnames.field = "V1")
#   tmp2 = tmp2[-which(elementMetadata(tmp2)[,1]==0)]
#   dat = c(dat,tmp2)
#   save(dat,file='rounded_brain_score.Rdata')
#   i = i+1
#   cat(steps*(i-1)+1,'\n')
# }

#tmp3 = unlist(reduce(split(tmp2, elementMetadata(tmp2)$V4)))
#tmp3@elementMetadata$score = names(tmp3)
#names(tmp3) = NULL

#tmp3 = sort(tmp3)
brain.gr = tmp2
brain.gr@elementMetadata$score = brain.gr@elementMetadata$V4
brain.gr@elementMetadata$V4 = NULL
  
save(brain.gr,file='rounded_brain_score.Rdata')


#################################################
# Heart scores

tmp = fread("zcat /sdata/vcfannotations/TiSAn_Heart_light.bed.gz") 
tmp = tmp[tmp$V4 > 0,]
tmp2 = makeGRangesFromDataFrame(tmp,keep.extra.columns = TRUE,start.field = "V2",end.field = "V3",seqnames.field = "V1")
heart.gr = tmp2
heart.gr@elementMetadata$score = heart.gr@elementMetadata$V4
heart.gr@elementMetadata$V4 = NULL

save(heart.gr,file='rounded_heart_score.Rdata')
