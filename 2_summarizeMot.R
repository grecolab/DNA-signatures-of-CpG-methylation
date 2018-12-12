library(Biostrings)
library(tibble)
library(gplots)
library(matrixStats)
library(factoextra)
library(dendextend)

setwd("./motif_supp_code/meme_out/")
load("./motif_supp_code/iDMCs.Rdata")


# motif_summary -----------------------------------------------------------

#set minumum percentage of CpGs to call a support set hyper- or hypo-methylated 
unbalance_tresh <-  0.7

#experiment names
types = c("Asthma","RA","SLE","T1D","T2D")


#retrives methylation information from the support set of the motif
get_meth_sign <- function(x,idm_type) {
  
  idm_type <- idm_type[!is.na(idm_type$hyper),]
  
  motif_file <- paste0(x["type"],"/",x["motif"],".txt")
  
  mot_cpg = read.table(motif_file , quote="\"", comment.char="",header = T, stringsAsFactors=FALSE)
  
  
  hyper_type=rownames(idm_type)[idm_type$hyper]
  n_hyper = length(intersect(hyper_type,mot_cpg$id))
  
  hypo_type=rownames(idm_type)[!idm_type$hyper]
  n_hypo = length(intersect(hypo_type,mot_cpg$id))
  
  tot_hyper = length(hyper_type)
  tot_hypo = length(hypo_type)
  
  #test for enrichment of hyper- or hypo-methylated CpGs in the support set
  contingency_mat <- matrix(c(n_hyper, n_hypo, tot_hyper-n_hyper, tot_hypo-n_hypo) ,nrow=2,ncol=2)
  fisher = fisher.test(contingency_mat)
  
  #save collected info in a data.frame
  motif_out_data <- data.frame(recovered_pos= nrow(mot_cpg), 
                               n_hyper=n_hyper,n_hypo=n_hypo, 
                               tot_hyper=tot_hyper, 
                               tot_hypo = tot_hypo, 
                               fisher.p = fisher$p.value)
  return(motif_out_data)
}


types_mot = tibble::tibble(type = character(), motif = character(),length = numeric(),
                           positives = numeric(), tot_pos = numeric(),
                           negatives = numeric(), tot_neg = numeric(),
                           pos_ratio = numeric(), neg_ratio = numeric(),
                           pvalue = numeric(), evalue = numeric(),
                           recovered_pos=numeric(),
                           n_hyper= numeric(), n_hypo = numeric(),
                           tot_hyper= numeric(),tot_hypo = numeric(),
                           fisher.p = numeric())

for (type in types){
  type_cpg_of_int <- iDMCs[[type]]
  
  types_tab = read.table(paste0(type,"/summary.txt"), header = T,stringsAsFactors = F)
  types_tab=cbind(types_tab, 
                  do.call(rbind, 
                          apply(types_tab, 1 ,
                                function(x,idmcs) get_meth_sign(x,idmcs), idmcs=type_cpg_of_int)
                          )
                  )
  types_mot = dplyr::bind_rows(types_mot, types_tab)
}


#code methylation status of the support set as 1 if > unbalance_tresh % are hyper, 
# -1 if > unbalance_tresh % are hypo, 0 otherwise
types_mot$main_dir = 0
types_mot[types_mot$n_hyper/(types_mot$recovered_pos) > unbalance_tresh, ]$main_dir = 1  
types_mot[types_mot$n_hypo/(types_mot$recovered_pos) > unbalance_tresh, ]$main_dir = -1  

write.table(types_mot,file = "summary.txt",quote = F, row.names = F, col.names = T)



# motif_similiarity -------------------------------------------------------

types_mot$name  <-  types_mot$type

cat(paste0(">",types_mot$motif,"-->",types_mot$type,"\n",types_mot$motif,"\n",collapse = ""))
names= paste0(types_mot$name,"-->",types_mot$motif,"-->",types_mot$main_dir)

#create similiarity matrix
mat_aln_all = matrix(NA,nrow = nrow(types_mot), ncol = nrow(types_mot),dimnames = list(names,names))

for(i in 1:nrow(mat_aln_all))
  for(j in 1:nrow(mat_aln_all))
    #here you can provide another function to compute the similarity between i-th and j-th motif
    #for example you can provide a function that akes two PWMs in input
    mat_aln_all[i,j] <- pairwiseAlignment(types_mot$motif[i], types_mot$motif[j], 
                                        substitutionMatrix=nucleotideSubstitutionMatrix(),
                                        type="global",scoreOnly=T,gapOpening=1, gapExtension=0.1)

# min max normalization in the range 0-1
range01 <- function(x){
  (x - min(x))/(max(x) - min(x))
}

mat_aln_all_0_1 = apply(mat_aln_all, 2, range01)
#mat_aln_d=apply(mat_aln_all,2,function(x) 1 - x/max(x))
mat_aln_d_all = as.dist(1 - mat_aln_all_0_1)


#represent distances as heatmap
heatmap(as.matrix(mat_aln_d_all),scale = "none")

#represent distances as dendrogram

par(mar=c(15,5,2,2))

hc_all=hclust(mat_aln_d_all)

col_dend  <- hc_all %>% as.dendrogram
dir= as.factor(sapply(strsplit(colnames(mat_aln_all_0_1),"-->"), function(x){x[3]}))
plot(col_dend)
#add methylation info as a colored bar
colored_bars(as.character(c("gray","red","green")[dir]), col_dend)

save.image("summary.Rdata")
