###################################
#Meta-Analysis Code
#
###################################

###Combine all the results
final<-mapply(function(x,y,z){
  #####Testing the function
  #x<-gtex[[5]]
  #y<-cm[[5]]
  #z<-loci[[5]]
  
  
  #print(z)
  
  temp.gtex<-x[!is.na(x$Pvalue),]
  temp.cm<-y[!is.na(y$Pvalue),]
  
  
  temp.gtex$eqtl<-with(temp.gtex,paste(SNP,Gene,sep="-"))
  temp.gtex<-temp.gtex[!duplicated(temp.gtex$eqtl),]
  temp.cm$eqtl<-with(temp.cm,paste(SNP,Gene,sep="-"))
  temp.cm<-temp.cm[!duplicated(temp.cm$eqtl),]
  
  full<-inner_join(temp.gtex,temp.cm,by="eqtl") 
  print(sum(duplicated(full$eqtl)))
  #full<-full %>% arrange(Gene.y,Pvalue_FE)
  #full$newFDR_GTEx<-unlist(lapply(split(full$Pvalue_FE,full$Gene.y),function(y1)p.adjust(y1,method="fdr")))
  
  
  #full<-full %>% arrange(Gene.y,Pvalue)
  #full$newFDR_CM<-unlist(lapply(split(full$Pvalue,full$Gene.y),function(y1)p.adjust(y1,method="fdr")))
  
  
  temp<-full %>% select(Gene=Gene.x,GeneName=GeneName.x,SNP=SNP.x,
                        Beta_GTEx=Beta.x,Stat_GTEx=Statistic.x,PVALUE_GTEx=Pvalue.x,
                        Beta_CM=Beta.y,Stat_CM=Statistic.y,PVALUE_CM=Pvalue.y)
  temp
  
},gtex,cm,loci,SIMPLIFY=F)
names(final)<-loci





#y1=1
######UPDATED

final.tables<-mapply(function(x,z){
  #x<-final[[5]]
  #z<-loci[[5]]
  
  print(z)
  
  gwassnp<-sapply(strsplit(z,"-"),'[',3)
  type<-sapply(strsplit(z,"-"),'[',1)
  
  
  ###ID THE LD for SNPS
  system(sprintf("plink --vcf /home/lai-00/rkl1/nousome/COLOC_targeted/GENO/1000G.vcf.gz --r2 --ld-snp %s --ld-window-kb 100000 --ld-window 99999 --ld-window-r2 0",
                 gwassnp))
  plink.ld<-fread("plink.ld")
  plink.ld$CHR_B<-as.character(plink.ld$CHR_B)
  
  ld<-ldsnp[[which(names(ldsnp) %in% gwassnp)[1]]]
  temp1<-strsplit(ld$Coord,":")
  temp2<-strsplit(ld$Alleles,"[(/)]")
  ld.snp.pos<-data.table(chr=sapply(temp1,function(x)gsub("chr","",x[[1]])),
                         start=as.numeric(sapply(temp1,'[',2)))
  ld.snp.pos$end<-as.numeric(ld.snp.pos$start)+sapply(strsplit(ld$Alleles,"[(/)]"),function(y)sum(nchar(y))-1)
  ld.snp.pos$A1<-sapply(strsplit(ld$Alleles,"[(/)]"),'[',2)
  ld.snp.pos$A2<-sapply(strsplit(ld$Alleles,"[(/)]"),'[',3)
  
  ld.snp.pos<-cbind(ld.snp.pos,R2=ld$R2)
  
  tempsnp1<-sapply(strsplit(x$SNP,";"),function(x1)x1[[grep("rs",x1)[[1]]]])
  #rs<-sapply(tempsnp1,function(x)paste(x,collapse=","))
  rs<-tempsnp1
  
  tempsnp1<-sapply(strsplit(x$SNP,";"),function(x)unlist(sapply(x,function(y)y[grep("[0-9]",substr(y,1,1))])))
  tempsnp2<-data.table(do.call(rbind,strsplit(tempsnp1,"_")))[,1:4]
  tempsnp2$V2<-as.numeric(tempsnp2$V2)
  
  
  
  x1<-bind_cols(data.frame(rs),x,tempsnp2)
  x2<-left_join(x1,ld.snp.pos,by=c("V1"="chr","V2"="start")) 
  x3<-plink.ld %>% unique() %>% left_join(x2,.,by=c("V1"="CHR_B","V2"="BP_B"))
  x3$R2<-ifelse(is.na(x3$R2.x),ifelse(is.na(x3$R2.y),NA,x3$R2.y),x3$R2.x)
  x3<-x3 %>% filter(!is.na(R2))
  x3<-x3 %>% arrange(Gene,PVALUE_CM)
  x3$QVALUE_CM<-unlist(lapply(split(x3$PVALUE_CM,x3$Gene),function(y1)p.adjust(y1,method="fdr")))
  
  
  x3<-x3 %>% arrange(Gene,PVALUE_GTEx)
  x3$QVALUE_GTEx<-unlist(lapply(split(x3$PVALUE_GTEx,x3$Gene),function(y1)p.adjust(y1,method="fdr")))
  
  check.o<-with(x3,GRanges(seqnames=paste0("chr",CHR_A),IRanges(start=V2,width=1)))
  func.over<-do.call(cbind,lapply(functional,function(x){
    check.o %over% x
  }))
  
  
  func.over1<-apply(
    sapply(1:ncol(func.over),function(y1){
      s<-ifelse(func.over[,y1]==T,y1,NA)
    }),1,function(z){
      z1<-z[!is.na(z)]
      paste(z1,collapse=",")
    })
  
  x3$func_over<-func.over1
  x4<-x3 %>% mutate(SE_CM_INV=(1/(Beta_CM/Stat_CM)^2),SE_GTEx_INV=(1/(Beta_GTEx/Stat_GTEx)^2),
                    Bmeta=((Beta_CM*SE_CM_INV)+(Beta_GTEx*SE_GTEx_INV))/(SE_CM_INV+SE_GTEx_INV),
                    PoolSE=sqrt(1/(SE_CM_INV+SE_GTEx_INV)),
                    Z=Bmeta/PoolSE,Pmeta=dnorm(Z)) %>%
    arrange(Gene,Pmeta)
  #Vmeta=1/(SE_CM_INV^2+SE_GTEx_INV^2),
  #X2=(Bmeta^2)/Vmeta,
  #Pmeta=pchisq(X2,1,lower.tail=F)) %>% arrange(Gene,Pmeta)
  x4$FDRmeta<-unlist(lapply(split(x4$Pmeta,x4$Gene),function(y)p.adjust(y,method="fdr",n=length(y))))
  
  x5<-x4 %>% mutate(rs1=rs,TargetGene=GeneName,
                    CM_B_SE=paste0(round(Beta_CM,2)," (",round(Beta_CM/Stat_CM,2),")"),
                    GTEx_B_SE=paste0(round(Beta_GTEx,2)," (",round(Beta_GTEx/Stat_GTEx,2),")")) %>%
    #filter(FDRmeta<0.1) %>%
    select(rsid=SNP,SNP=rs1,LD=R2,TargetGene,GeneName=Gene,
           CM_Beta=CM_B_SE,CM_Pvalue=PVALUE_CM,
           GTEx_Beta=GTEx_B_SE,GTEx_Pvalue=PVALUE_GTEx,
           Meta_Beta=Bmeta,SE=PoolSE,Z,Pmeta,FDRmeta,Func_Over=func_over) %>%
    #Meta_Beta=Bmeta,Var=Vmeta,Pmeta,FDRmeta,Func_Over=func_over) %>%
    arrange(TargetGene,Pmeta)
  x5
  
},final,names(final),SIMPLIFY = F)