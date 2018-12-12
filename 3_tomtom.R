library(foreach)
library(doParallel)
require(XML)
library(methods)
library(gplots)

setwd("./motif_supp_code/meme_out/")

# experiment names 
types = c("Asthma","RA","SLE","T1D","T2D")


# this script uses tomtom utility and meme motif databeses, make sure you have them installed 
# on your system and provide the corresponding paths here
meme_bin__path = "/home/<user>/meme/bin/"
motif_dbs="/home/megisc/<user>/motif_databases/"

#tf_db="HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
tf_db="JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme"

#setup parallel backend to use 8 processors
#if you don't wanto to run the analisys on multiple cores comment the following twolines
cl<-makeCluster(8)
registerDoParallel(cl)


res_all = foreach(type_can = types,.combine = rbind) %dopar% {
  require(XML)
  library("methods")
  
  tomtom_comm=paste0(meme_bin__path,"/tomtom -no-ssc -oc tmp -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 ",
                     type_can,"/dreme.txt ",motif_dbs,"/",tf_db, " -oc ",type_can,"/tomtom")
  system(tomtom_comm)
  
  res=read.table(paste0(type_can,"/tomtom/tomtom.txt"),stringsAsFactors = F)
  colnames(res) = c("Query_ID",	"Target_ID"	,"Optimal_offset",	"p.value"	,"E.value",	"q.value",	"Overlap"	,"Query_consensus"	,"Target_consensus",	"Orientation")
  res$type=type_can
  res$tfbs = sapply(res$Target_ID,function(x) strsplit(x,"_")[[1]][1])
  
  #parse dreme xml output 
  xml_out = xmlParse(paste0(type_can,"/tomtom/tomtom.xml"))
  rootnode = xmlRoot(xml_out)
  motifs = xmlChildren(rootnode[["targets"]])
  tab_names = do.call(rbind,
                      lapply(as.list(motifs), 
                             function(x) { 
                                c(as.character(xmlAttrs(x)["id"]),as.character(xmlAttrs(x)["alt"]))
                              } 
                             )
                      )
  matches = match(res$Target_ID,tab_names[,1])
  res$tfbs2 = tab_names[matches,2]
  
  return(res)
}

#comment this if you dont use multi-core
stopCluster(cl)



#annoatate motifs with results from tomtom

load("summary.Rdata")

ann_tfbs <- types_mot %>% dplyr::select(type,motif,main_dir) %>% dplyr::inner_join(res_all, by = c( "type"="type","motif"="Query_ID"))

rows <- unique(ann_tfbs$tfbs2)
columns <- unique(ann_tfbs$type)

#create a binary matrix where we mark the presence of the TFs in the experimnts 
mat = matrix(0,nrow = length(rows), ncol = length(columns),dimnames = list(rows,columns)) 
for (i in 1:nrow(ann_tfbs)){
  mat[ann_tfbs$tfbs2[i],ann_tfbs$type[i]]=1
}

heatmap.2(mat,Colv=T,Rowv = F,dendrogram = "col",tracecol = "black",colsep=1:ncol(mat),sepwidth=c(0.001,0.001),
          rowsep=1:nrow(mat), sepcolor="black", col=c("white","gray","red"),cexRow = 1,cexCol = 1, margins = c(5,10), key=F)
 

save.image("tomtom.Rdata")
