
###Library install###
librarylist_cran=c("shinyBS","shinyjs","shinydashboard","R.utils","data.table","png",
                   "fields","rmarkdown","DT","ggpubr","sparklyr","flexdashboard","shinycssloaders","highcharter","fastICA","rgl","pals","circlize")

librarylist_bc=c("ComplexHeatmap","org.Hs.eg.db","org.Mm.eg.db")

is.exist_cran <- librarylist_cran %in% rownames(installed.packages())
is.exist_bc <- librarylist_bc   %in% rownames(installed.packages())
if(length(librarylist_cran[!(is.exist_cran)])>0){install.packages(librarylist_cran[!(is.exist_cran)])}
if(length(librarylist_bc[!(is.exist_bc)])>0){
  install.packages("BiocManager")
  BiocManager::install(librarylist_bc[!(is.exist_bc)])
}
###Library install###

###Library loading###
library(R.utils)
library(fastICA)
library(rmarkdown)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(rgl)
library(pals)
library(circlize)
library(ComplexHeatmap)
###Library loading###



##user defined IC
run_ICA<-function(dm,IC){
  gc()
  gc()
  
  dm_loci=paste(dm[,1],dm[,2],dm[,3],dm[,4],sep="__")
  dm_data=dm[,c(5:ncol(dm))]
  rownames(dm_data)=dm_loci  
  
  set.seed(1974414)
  
  resICA  <- fastICA(X=dm_data, n.comp=IC,  maxit=100000, alg.typ="parallel",tol = 10^-6)
  colnames(resICA$A)=colnames(dm_data)
  resICA$A=t(resICA$A)
  
  return(resICA)
  
}


##PC >=80
run_ICA_PC80<-function(dm){
  gc()
  gc()
  
  dm_loci=paste(dm[,1],dm[,2],dm[,3],dm[,4],sep="__")
  dm_data=dm[,c(5:ncol(dm))]
  rownames(dm_data)=dm_loci
  
  set.seed(1974414)
  
  ##num PC >= 80
  pcaobj <- prcomp(t(dm_data), scale = TRUE)
  pCs=pcaobj$sdev ^ 2
  kiyo=pCs / sum(pCs)
  IC=sum(cumsum(kiyo)*100<80)+1
  
  resICA  <- fastICA(X=dm_data, n.comp=IC,  maxit=100000, alg.typ="parallel",tol = 10^-6)
  colnames(resICA$A)=colnames(dm_data)
  resICA$A=t(resICA$A)
  
  return(resICA)
  
}



###tabName=ICs
make_ICA_heatmap_shiny<-function(resICA,sampleList,GenomicFeature){
  
  #collist=c("blue","cyan","gray","yellow","red","pink","purple","lightgreen","magenta","green","orange")
  #collist=c("magenta","green","red","blue")
  set.seed(0414) 
  collist=grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][sample(433,100,replace=F)]
  
  if (GenomicFeature==1){Loci="CGI"}
  if (GenomicFeature==2){Loci="Gene body"}
  if (GenomicFeature==3){Loci="1st intron"}
  if (GenomicFeature==4){Loci="Promoter"}

  ###all data
  A=resICA$A

  ###remove common ICs
  #A=t(A0[(-1)*order(abs(rowMeans(A0)),decreasing=T)[1],])
  ######A0=t(resICA$A)
  #A=t(A0[(-1)*order(rowMax(abs(A0)),decreasing=T)[1],])
  ######A=t(A0[(-1)*order(as.numeric(rowMax(abs(A0))),decreasing=T)[1],])

  colnames(A)=paste("IC",seq(1: ncol(A)),sep="")
  Dist<-'ward.D2';
  coLlegend=min(abs(min(A)),abs(max(A)))
  
  coL = colorRamp2(c(-1*(coLlegend), 0, coLlegend), c("blue", "white", "red"))
  
  ColList=NULL
  for (i in 1:ncol(sampleList)){
    
    #classnum=length(levels(sampleList[,i]))
    classnum=length(unique(as.vector(sampleList[,i])))
    
    ColList.tmp=NULL
    
    for (j in 1:classnum){
      ColList.tmp=c(ColList.tmp, collist[j])
    }
    
    #names(ColList.tmp)=levels(sampleList[,i])
    #names(ColList.tmp)=sort(as.character(unique(as.vector(sampleList[,i]))))
    names(ColList.tmp)=sort(unique(as.vector(sampleList[,i])))
    ColList=c(ColList,list(ColList.tmp))
    
  }
  
  names(ColList)=colnames(sampleList)
  
  
  ha = HeatmapAnnotation(df = sampleList,
                         col = ColList,
                         show_annotation_name = TRUE,
                         annotation_name_side="left")
  
  xx=Heatmap(t(A[rownames(sampleList),]), 
             top_annotation = ha,
             col = coL,
             clustering_method_rows="ward.D2",
             clustering_method_columns="ward.D2",
             heatmap_legend_param = list(color_bar = "continuous"),
             row_dend_width = unit(50, "mm"),
             column_dend_heigh = unit(50, "mm"), 
             name="A",
             column_title=Loci,
             show_row_names = TRUE,
             show_column_names = TRUE,
             
             row_names_gp = gpar(fontsize = 9),
             column_names_gp = gpar(fontsize = 9)
  )
  
  return(xx)
}


##########make_highLF_heatmap_target_shiny<-function(resICA,sampleList,datamat,Zscore,Target,GenomicFeature){
make_highLF_heatmap_target_shiny<-function(resICA,sampleList,datamat,Zscore,Target,GenomicFeature,Platform){
  library(circlize)
  library(ComplexHeatmap)
  
  dm_loci=paste(datamat[,1],datamat[,2],datamat[,3],datamat[,4],sep="__")
  dm_data=datamat[,c(5:ncol(datamat))]
  rownames(dm_data)=dm_loci  
  
  
  #widthunit=3.5*NROW(sampleList)
  widthunit=140
  

  ###remove common ICs
  ######A0=t(resICA$A)
  #resICA$S=resICA$S[,(-1)*order(rowMax(abs(A0)),decreasing=T)[1]]
  ######resICA$S=resICA$S[,(-1)*order(as.numeric(rowMax(abs(A0))),decreasing=T)[1]]
  
    
  
  #collist=c("blue","cyan","gray","yellow","red","pink","purple","lightgreen","magenta","green","orange")
  #collist=c("magenta","green","red","blue")
  set.seed(0414) 
  collist=grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][sample(433,100,replace=F)]
  
  
  if (GenomicFeature==1){Loci="CGI"}
  if (GenomicFeature==2){Loci="Gene body"}
  if (GenomicFeature==3){Loci="1st intron"}
  if (GenomicFeature==4){Loci="Promoter"}
  
  ColList=NULL
  for (i in 1:ncol(sampleList)){
    
    #classnum=length(levels(sampleList[,i]))
    classnum=length(unique(as.vector(sampleList[,i])))
    
    ColList.tmp=NULL
    
    for (j in 1:classnum){
      ColList.tmp=c(ColList.tmp, collist[j])
    }
    
    #names(ColList.tmp)=levels(sampleList[,i])
    names(ColList.tmp)=sort(unique(as.vector(sampleList[,i])))
    ColList=c(ColList,list(ColList.tmp))
    
  }
  
  names(ColList)=colnames(sampleList)
  
  i=Target
  target_dm=subset(resICA$S[,i],abs(resICA$S[,i])>=Zscore)
  ZhighNum=NROW(subset(resICA$S[,i],resICA$S[,i]>=Zscore))
  ZlowNum =NROW(subset(resICA$S[,i],resICA$S[,i]<=(-1)*Zscore))       
  if (length(target_dm)>0){
    sort_names=names(sort(target_dm,decreasing=T))
    
    datamat2=data.frame(dm_data[sort_names,rownames(sampleList)],IC=target_dm[sort_names])
    datamat0=dm_data[sort_names,rownames(sampleList)]
    
    
    value = datamat2$IC      
    ha = HeatmapAnnotation(df = sampleList,
                           col = ColList,
                           show_annotation_name = TRUE,
                           annotation_name_side="left"
    )
    
    row_ha = rowAnnotation(Loadings = row_anno_barplot(value, axis = TRUE, axis_side = "bottom",baseline=0,ylim=c(-15,15)),
                           width = unit(3, "cm"),
                           show_annotation_name = TRUE,
                           annotation_name_offset = unit(1, "cm"),
                           annotation_name_rot = c(0, 0, 20)
                           
    )
    
    coL=colorRampPalette(c("darkblue","yellow"))(10)
    
    if (Platform==1){meth_legend="meth level"}
    if (Platform==2){meth_legend="beta value"}
    
    xx<-row_ha+Heatmap(datamat0,col=coL,
                       width = unit(widthunit, "mm"), cluster_rows = FALSE,
                       top_annotation = ha, 
                       cluster_columns = TRUE,
                       clustering_distance_columns = "euclidean",
                       clustering_method_columns="ward.D2",
                       show_row_names = FALSE,
                       ##########name="meth level",
                       name=meth_legend,
                       column_title = paste(Loci,"\nIC",i,"(Loadings>",Zscore,":",ZhighNum, ", Loadings<-", Zscore , ":", ZlowNum,")",sep=""),
                       bottom_annotation_height = unit(2, "cm"))
    
    
    draw(xx,heatmap_legend_side = "left",annotation_legend_side = "left")
    
  }else{
    next;   
  }
  
}



make_highLF_meth_target_shiny<-function(resICA,sampleList,datamat,Zscore,Target){
  
  dm_loci=paste(datamat[,1],datamat[,2],datamat[,3],datamat[,4],sep="__")
  dm_data=datamat[,c(5:ncol(datamat))]
  rownames(dm_data)=dm_loci  
  
  
  ###remove common ICs
  ######A0=t(resICA$A)
  #resICA$S=resICA$S[,(-1)*order(rowMax(abs(A0)),decreasing=T)[1]]
  ######resICA$S=resICA$S[,(-1)*order(as.numeric(rowMax(abs(A0))),decreasing=T)[1]]

  
  
  #collist=c("blue","cyan","gray","yellow","red","pink","purple","lightgreen","magenta","green","orange")
  #collist=c("magenta","green","red","blue")
  set.seed(0414) 
  collist=grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][sample(433,100,replace=F)]
  
  
  i=Target
  target_dm=subset(resICA$S[,i],abs(resICA$S[,i])>=Zscore)
  
  if (length(target_dm)>0){
    sort_names=names(sort(target_dm,decreasing=T))
    datamat2=data.frame(dm_data[sort_names,rownames(sampleList)],IC=target_dm[sort_names])
    datamat0=dm_data[sort_names,rownames(sampleList)]
    return(datamat0)
  }else{
    next;   
  }
  
}



generate_region2loadings<-function(resICA,Genome,Feature){
  
  ###remove common ICs
  ######A0=t(resICA$A)
  #resICA$S=resICA$S[,(-1)*order(rowMax(abs(A0)),decreasing=T)[1]]
  ######resICA$S=resICA$S[,(-1)*order(as.numeric(rowMax(abs(A0))),decreasing=T)[1]]
  
  
  ###all data
  targetregions=NULL
  
  if (Feature>=2){
    targetloci=matrix(unlist(strsplit(rownames(resICA$S), "__")),ncol=4,byrow=T)
    
    if (Genome<=2){geneSymbols2Id <- mapIds(org.Hs.eg.db, keys=targetloci[,4], column="ENTREZID", keytype="SYMBOL", multiVals="first")}
    if (Genome>=3){geneSymbols2Id <- mapIds(org.Mm.eg.db, keys=targetloci[,4], column="ENTREZID", keytype="SYMBOL", multiVals="first")}
    
    for(i in 1:ncol(resICA$S)){
      targetregions.tmp=data.frame(chr=targetloci[,1],
                                   start=as.numeric(targetloci[,2]),
                                   end=as.numeric(targetloci[,3]),
                                   GeneSymbol=paste0("<a href=\"https://www.ncbi.nlm.nih.gov/gene/",geneSymbols2Id,"\" target=\"blank_\">",targetloci[,4],"</a>"),
                                   genesymbol=targetloci[,4],
                                   Loadings=as.numeric(round(resICA$S[,i],2)),
                                   IC=paste("IC",i,sep=""))
      targetregions=rbind(targetregions,targetregions.tmp)
    }
  }
  
  if (Feature==1){
    targetloci=matrix(unlist(strsplit(rownames(resICA$S), "__")),ncol=4,byrow=T)
    targetloci0=gsub("__CGI", "", rownames(resICA$S))
    #targetloci=matrix(unlist(strsplit(rownames(resICA$S), "__")),ncol=4,byrow=T)
    
    if (Genome==1){CGI_closestTSS=read.table("~/annotation/hg38_cgi2closestTSS_nr",header=T,row.names=1)}
    if (Genome==2){CGI_closestTSS=read.table("~/annotation/hg19_cgi2closestTSS_nr",header=T,row.names=1)}
    if (Genome==3){CGI_closestTSS=read.table("~/annotation/mm10_cgi2closestTSS_nr",header=T,row.names=1)}
    if (Genome==4){CGI_closestTSS=read.table("~/annotation/mm9_cgi2closestTSS_nr",header=T,row.names=1)}
    

    for(i in 1:ncol(resICA$S)){
      targetregions.tmp=data.frame(chr=targetloci[,1],
                                   start=as.numeric(targetloci[,2]),
                                   end=as.numeric(targetloci[,3]),
                                   Loadings=as.numeric(round(resICA$S[,i],2)),
                                   #ClosestTSS=as.vector(CGI_closestTSS[rownames(resICA$S),"symbol"]),
                                   ClosestTSS=as.vector(CGI_closestTSS[as.matrix(targetloci0),"symbol"]),
                                   #DistToClosestTSS=as.numeric(CGI_closestTSS[rownames(resICA$S),"dist"]),
                                   DistToClosestTSS=as.numeric(CGI_closestTSS[as.matrix(targetloci0),"dist"]),
                                   IC=paste("IC",i,sep=""))
      targetregions=rbind(targetregions,targetregions.tmp)
    }
  }
  
  
  targetregions=subset(targetregions, abs(Loadings)>=2)
  
  return(targetregions)
  
}




generate_highLF_target_shiny<-function(resICA,TargetIC,Zscore,Genome,Feature){
  
  ###remove common ICs
  ######A0=t(resICA$A)
  #resICA$S=resICA$S[,(-1)*order(rowMax(abs(A0)),decreasing=T)[1]]
  ######resICA$S=resICA$S[,(-1)*order(as.numeric(rowMax(abs(A0))),decreasing=T)[1]]
  
  
  if (Genome==1){genomeref="hg38"}
  if (Genome==2){genomeref="hg19"}
  if (Genome==3){genomeref="mm10"}
  if (Genome==4){genomeref="mm9"}
  
  
  ####all data
  target_dm=subset(resICA$S[,TargetIC],abs(resICA$S[,TargetIC])>=Zscore)
  
  
  
  if (length(target_dm)>1){
    sort_names=names(sort(target_dm,decreasing=T))
    
    #targetloci=matrix(unlist(strsplit(names(target_dm), "_")),ncol=4,byrow=T)
    
    if (Feature>=2){
      targetloci=matrix(unlist(strsplit(sort_names, "__")),ncol=4,byrow=T)
      
      if (Genome <=2){
        geneSymbols2Id <- mapIds(org.Hs.eg.db, keys=targetloci[,4], column="ENTREZID", keytype="SYMBOL", multiVals="first")
      }
      
      if (Genome >=3){
        geneSymbols2Id <- mapIds(org.Mm.eg.db, keys=targetloci[,4], column="ENTREZID", keytype="SYMBOL", multiVals="first")  
      }
      
      targetgenes=data.frame(chr=targetloci[,1],
                             start=as.numeric(targetloci[,2]),
                             end=as.numeric(targetloci[,3]),
                             GeneSymbol=paste0("<a href=\"https://www.ncbi.nlm.nih.gov/gene/",geneSymbols2Id,"\" target=\"blank_\">",targetloci[,4],"</a>"),
                             #genesymbol="CGI",
                             genesymbol=targetloci[,4],
                             symbol=targetloci[,4],
                             Loadings=as.numeric(round(target_dm[sort_names],2)))
    }
    
    if (Feature==1){
      
      targetloci=matrix(unlist(strsplit(sort_names, "__")),ncol=4,byrow=T)
      sort_names0=gsub("__CGI", "", sort_names)
      
      if (Genome==1){CGI_closestTSS=read.table("~/annotation/hg38_cgi2closestTSS_nr",header=T,row.names=1)}
      if (Genome==2){CGI_closestTSS=read.table("~/annotation/hg19_cgi2closestTSS_nr",header=T,row.names=1)}
      if (Genome==3){CGI_closestTSS=read.table("~/annotation/mm10_cgi2closestTSS_nr",header=T,row.names=1)}
      if (Genome==4){CGI_closestTSS=read.table("~/annotation/mm9_cgi2closestTSS_nr",header=T,row.names=1)}
      
      #ClosestTSS=as.vector(CGI_closestTSS[as.matrix(targetloci0),"symbol"]),
      
      #target_CGI_closestTSS=CGI_closestTSS[sort_names,]
      target_CGI_closestTSS=CGI_closestTSS[sort_names0,]
      
      if (Genome <=2){
        #  #https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr10%3A82116202-82117120&showTracks=1
        #  #CpGi=paste0("<a href=\"https://genome.ucsc.edu/cgi-bin/",Genome,"Tracks?db=",version, "&position=",chr,"%3A",start,"-",end,"&showTracks=1","\" target=\"blank_\">UCSC Genome Browser</a>")
        
        CGI_URL=paste0("<a href=\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=",genomeref, "&position=",targetloci[,1],"%3A",as.numeric(targetloci[,2]),"-",as.numeric(targetloci[,3]),"&showTracks=1","\" target=\"blank_\">GenomeBrowser</a>")

        #GeneSymbol=paste0("<a href=\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=",chr,"%3A",start,"-",end,"&showTracks=1","\" target=\"blank_\">GenomeBrowser</a>"),
      }
      
      if (Genome >=3){
        CGI_URL=paste0("<a href=\"https://genome.ucsc.edu/cgi-bin/mmTracks?db=",genomeref, "&position=",targetloci[,1],"%3A",as.numeric(targetloci[,2]),"-",as.numeric(targetloci[,3]),"&showTracks=1","\" target=\"blank_\">GenomeBrowser</a>")
      }
      
      targetgenes=data.frame(chr=targetloci[,1],
                             start=as.numeric(targetloci[,2]),
                             end=as.numeric(targetloci[,3]),
                             
                             #GeneSymbol=paste0("<a href=\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=",Genome, "&position=",targetloci[,1],"%3A",as.numeric(targetloci[,2]),"-",as.numeric(targetloci[,3]),"&showTracks=1","\" target=\"blank_\">GenomeBrowser</a>"),
                             GenomeBrowser=CGI_URL,
                             
                             genesymbol="CGI",
                             Loadings=as.numeric(round(target_dm[sort_names],2)),
                             ClosestTSS=as.vector(target_CGI_closestTSS[,"symbol"]),
                             DistToClosestTSS=as.numeric(target_CGI_closestTSS[,"dist"])

                             
      )
      
    }
    
  } 
  return(targetgenes)
} 



ggplot2boxplot_IC_shiny<-function(resICA,sampleList,Subgroup){
  
  #collist=c("blue","cyan","gray","yellow","red","pink","purple","lightgreen","magenta","green","orange")
  #collist=c("magenta","green","red","blue")
  set.seed(0414) 
  collist=grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][sample(433,100,replace=F)]
  
  
  ###all data
  resICA_A=resICA$A
  
  
  ###remove common ICs
  #A0=t(resICA$A)
  #resICA_A=t(A0[(-1)*order(abs(rowMeans(A0)),decreasing=T)[1],])
  
  ###remove common ICs
  ######A0=t(resICA$A)
  #A=t(A0[(-1)*order(abs(rowMeans(A0)),decreasing=T)[1],])
  #resICA_A=t(A0[(-1)*order(rowMax(abs(A0)),decreasing=T)[1],])
  ######resICA_A=t(A0[(-1)*order(as.numeric(rowMax(abs(A0))),decreasing=T)[1],])
  
  
  
  colnames(resICA_A)=paste("IC",c(1:ncol(resICA_A)),sep="")
  #df.tmp=data.frame(Labels=sampleList[rownames(resICA_A),Subgroup],resICA_A)
  df.tmp=data.frame(Labels=as.character(sampleList[rownames(resICA_A),Subgroup]),resICA_A)
  df.m <- as.data.table(melt(df.tmp,id.var= "Labels"))
  
  unique_group=NROW(unique(df.m$Labels))
  
  df.m[, ymin:=min(value)-0.02, by = .(variable)]
  df.m[, ymax:=min(value)+0.02, by = .(variable)]
  
  if (NROW(levels(sampleList[,Subgroup]))==2){
    p<-ggplot(df.m, aes(x=Labels, y = value)) + 
      geom_boxplot(aes(fill=Labels)) +
      scale_fill_manual(values = collist[1:unique_group]) +
      theme(axis.text.x = element_text()) + 
      geom_blank(aes(y=ymin))+
      geom_blank(aes(y=ymax))+
      facet_wrap(~variable, scales = 'free') + 
      stat_compare_means(aes(label = paste0("p = ", ..p.format..)),label.x.npc="left",label.y.npc="bottom",method="t.test") +
      ylab("A") 
  }else{
    p<-ggplot(df.m, aes(x=Labels, y = value)) + 
      geom_boxplot(aes(fill=Labels)) +
      scale_fill_manual(values = collist[1:unique_group]) +
      theme(axis.text.x = element_blank()) +
      geom_blank(aes(y=ymin))+
      geom_blank(aes(y=ymax))+
      facet_wrap(~variable, scales = 'free') + 
      stat_compare_means(aes(label = paste0("p = ", ..p.format..)),label.x.npc="left",label.y.npc="bottom",method="anova") +
      ylab("A") 
  }
  
  #p <- p + guides(fill=guide_legend(title="Legend_Title"))
  return(p)
}


