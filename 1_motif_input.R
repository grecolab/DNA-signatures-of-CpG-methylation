library(readr)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

flankSize=20

dir.create("beds",showWarnings = F)
dir.create("fasta",showWarnings = F)
dir.create("meme_out",showWarnings = F)



#get hg19 locations of 450k probes
data(Locations)
hg19_loc = as.data.frame(Locations)

#save background locations in a tab delimited file
background  = tibble::tibble(id=rownames(hg19_loc), chr = hg19_loc$chr, pos = hg19_loc$pos)
write.table(background,paste0("beds/450k_hg19_background.bed"),sep="\t",row.names = F,quote = F)

#target CpG analises names
targs = c("Asthma","RA","SLE","T1D","T2D")

# store data from multiple experiments in a list
iDMCs = list()

for (targ_name in targs){
  #### read the  file with differential methylation output for CpGs of interest 
  file = paste0("./testData/",targ_name,".tab")
  
  iDMC = read_delim(file,delim = "\t")
  iDMC <- iDMC[,c("ProbeID","status")]
  
  
  #get coordinates for the Cpgs of interest
  matches = match(iDMC$ProbeID, rownames(hg19_loc)) 
  iDMC = cbind(iDMC, hg19_loc[matches,])
  
  #define hyper- of hypo- methylation status 
  iDMC$hyper = (iDMC$status == "hyper")
  
  # store position and methylation status only
  iDMCs[[targ_name]] = iDMC
  
  #we need the genomic positions and ids only 
  iDMC_pos = iDMC[,c("ProbeID","chr","pos")]
  colnames(iDMC_pos) = c("id", "chr", "pos")
  
  #save the position file of CpGs in the beds folder 
  write.table(iDMC_pos,paste0("beds/",targ_name,".bed"),sep="\t",row.names = F,quote = F)
  
  #call the motifs_lib.R script passing the id of the analsyis, the id of the backgrounf file, the folder with position files
  # the folder where to store fasta files, the folder where to sore dreme output
  
  system(paste("Rscript ../motifs_lib.R", targ_name,"450k_hg19_background","beds","fasta","meme_out",flankSize), wait=FALSE)
}

#the above for loop runs all the motif search analyses in parallel, if you want to serialize set wait to FALSE

save(iDMCs,file = "iDMCs.Rdata")

#Wait for all the motif search processes to finish before executing the next script

