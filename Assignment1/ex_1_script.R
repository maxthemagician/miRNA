# script to compute master table
if(!require(data.table)){
  install.packages("data.table")
}
if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(data.table)
library(ggplot2)
outputfile = ""

# functions ---------------------------------------------------------------

species_mapping <- function(short){
  return(organisms$`#name`[organisms$`#organism`==short])
}
versionmapping <- function(acc){
  s = as.integer(substr(acc,3,9))
  if(s<=217)return(1)
  if(s<=541)return(2)
  if(s<=750)return(3)
  if(s<=1286)return(4)
  if(s<=1432)return(5)
  if(s<=1719)return(6)
  if(s<=3200)return(7)
  if(s<=3724)return(8)
  if(s<=5378)return(9)
  if(s<=6236)return(10)
  if(s<=8204)return(11)
  if(s<=10321)return(12)
  if(s<=10694)return(13)
  if(s<=12988)return(14)
  if(s<=15994)return(15)
  if(s<=16745)return(16)
  if(s<=18482)return(17)
  if(s<=19661)return(18)
  if(s<=23566)return(19)
  if(s<=26421)return(20)
  if(s<=31825)return(21)
  return(22)
}


# main --------------------------------------------------------------------


cel = fread("cel.gff3")
hsa = fread("hsa.gff3")
mmu = fread("mmu.gff3")

organisms = fread("organisms.txt")
mirna.xls = fread("miRNA.csv")


names = c("species","miRNA_prec_name","miRNA_prec_ID","miRNA_prec_seq","5p_mat_name","5p_mat_ID","5p_mat_seq","3p_mat_name","3p_mat_ID","3p_mat_seq","version")
mastertable = matrix(nrow = length(rownames(mirna.xls)),ncol = 11)
colnames(mastertable) <- names

for(i in 1:length(rownames(mirna.xls))){
  mastertable[i,1] = species_mapping(strsplit(mirna.xls$ID[i],split="-")[[1]][1])
  mastertable[i,2] = mirna.xls$ID[i]
  mastertable[i,3] = mirna.xls$Accession[i]
  mastertable[i,4] = mirna.xls$Sequence[i]
  mastertable[i,5] = mirna.xls$Mature1_ID[i]
  mastertable[i,6] = mirna.xls$Mature1_Acc[i]
  mastertable[i,7] = mirna.xls$Mature1_Seq[i]
  mastertable[i,8] = mirna.xls$Mature2_ID[i]
  mastertable[i,9] = mirna.xls$Mature2_Acc[i]
  mastertable[i,10] = mirna.xls$Mature2_Seq[i]
  mastertable[i,11] = versionmapping(mirna.xls$Accession[i])
}

colnames(mastertable) <- names ### comes before write.csv
write.csv(mastertable, paste0(outputfile,"mastertable.csv"))


types = c(species_mapping("cel"),species_mapping("mmu"),species_mapping("hsa"))
ids = which(mastertable[,1] %in% types)


new_table = matrix(mastertable[ids,],ncol = 11)
colnames(new_table) <- names
write.csv(new_table, paste0(outputfile,"mastertable_cleared.csv"))


# task 3
# a)
# number of miRNAs for each species in version 22
number = table(new_table[,1])

# contains Versionnumber, Species, added miRNA for this version and species
version_number = as.data.frame(table(new_table[,11],new_table[,1]))

# version_number sorted by version 1-22
c = version_number[order(as.numeric(as.character(version_number$Var1))),] 

# sum up count of miRNAs per species per version
for(i in 4:length(rownames(c))){ 
    c$Freq[i] = c$Freq[i]+c$Freq[i-3]
}

c = as.data.frame(c)
jpeg("stat_3a.jpg", width = 1024, height = 720)
ggplot(c, aes (y = Freq, x=as.numeric(as.character(c$Var1)), colour = Var2))+geom_point()+geom_line()+
  xlab("miRBase Version") +
  ylab("miRNA entrys") +
  ggtitle("Number of miRNA precursers per species and version") +
  theme(legend.title=element_blank())
dev.off()

# b)
mastertable = as.data.frame(new_table)
ids_mature = which(mastertable["3p_mat_ID"]!="")
filtered = mastertable[ids_mature,c("species", "version")]

t = as.data.frame(table(filtered),stringsAsFactors = FALSE)
t = rbind(t , c("Caenorhabditis elegans", 5, 0))
t = rbind(t , c("Mus musculus", 5, 0))
t = rbind(t , c("Homo sapiens", 5, 0))
t = t[order(as.numeric(as.character(t$version))),] 
t = as.data.frame(t,stringsAsFactors = T)
t$version = as.numeric(t$version)
t$Freq = as.numeric(t$Freq)
t$Freq[1] =t$Freq[1]*2
t$Freq[2] =t$Freq[2]*2
t$Freq[3] =t$Freq[3]*2
for (org in t$species[1:3]) {
  print(org)
  index = which(t$species==org)
  for(i in 2:22){
    p = which(t$version[index]==i)
    p2 = which(t$version[index]==i-1)
    p2 = t$Freq[index[p2]]
    t$Freq[index[p]] = as.numeric(t$Freq[index[p]]) + as.numeric(p2)
  }
}

mature = t[order(as.numeric(as.character(t$version))),] 
mature = as.data.frame(mature,stringsAsFactors = T )
mature$version = as.numeric(mature$version)
mature$Freq = as.numeric(mature$Freq)

# print
jpeg("stat_3b.jpg", width = 1024, height = 720)
ggplot(mature, aes (y = mature$Freq, x=as.numeric(as.character(mature$version)), colour = mature$species))+geom_point()+geom_line()+
  xlab("miRBase Version") +
  ylab("miRNA entrys") +
  ggtitle("Number of miRNA matures per species and version") +
  theme(legend.title=element_blank())
dev.off()





ids_pre = which(mastertable["3p_mat_ID"]=="")
filtered = mastertable[ids_pre,c("species", "version")]

t = as.data.frame(table(filtered),stringsAsFactors = FALSE)
t = rbind(t , c("Caenorhabditis elegans", 5, 0))
t = rbind(t , c("Mus musculus", 5, 0))
t = rbind(t , c("Homo sapiens", 5, 0))
t = t[order(as.numeric(as.character(t$version))),] 
t = as.data.frame(t,stringsAsFactors = T)
t$version = as.numeric(t$version)
t$Freq = as.numeric(t$Freq)
t$Freq[1] =t$Freq[1]*2
t$Freq[2] =t$Freq[2]*2
t$Freq[3] =t$Freq[3]*2
for (org in t$species[1:3]) {
  print(org)
  index = which(t$species==org)
  for(i in 2:22){
    p = which(t$version[index]==i)
    p2 = which(t$version[index]==i-1)
    p2 = t$Freq[index[p2]]
    t$Freq[index[p]] = as.numeric(t$Freq[index[p]]) + as.numeric(p2)
  }
}

pre = t[order(as.numeric(as.character(t$version))),] 
pre = as.data.frame(pre,stringsAsFactors = T )
pre$version = as.numeric(pre$version)
pre$Freq = as.numeric(pre$Freq)

#jpeg("stat_3b2.jpg", width = 1024, height = 720)
#ggplot(pre, aes (y = pre$Freq, x=as.numeric(as.character(pre$version)), colour = pre$species))+geom_point()+geom_line()+
#  xlab("miRBase Version") +
#  ylab("miRNA entrys") +
#  ggtitle("Number of miRNA precurser per species and version") +
#  theme(legend.title=element_blank())
#dev.off()


# task 4

#mastertable["precursor_A"] <- NA
#mastertable["precursor_C"] <- NA
#mastertable["precursor_G"] <- NA
#mastertable["precursor_U"] <- NA

#mastertable["mature_5p_A"] <- NA
#mastertable["mature_5p_C"] <- NA
#mastertable["mature_5p_G"] <- NA
#mastertable["mature_5p_U"] <- NA

#mastertable["mature_3p_A"] <- NA
#mastertable["mature_3p_C"] <- NA
#mastertable["mature_3p_G"] <- NA
#mastertable["mature_3p_U"] <- NA


#count A, C, G, U
library(stringr)
library(stringi)

mastertable$precursor_A <- str_count(mastertable$miRNA_prec_seq, "A")/stri_length(mastertable$miRNA_prec_seq)
mastertable$precursor_C <- str_count(mastertable$miRNA_prec_seq, "C")/stri_length(mastertable$miRNA_prec_seq)
mastertable$precursor_G <- str_count(mastertable$miRNA_prec_seq, "G")/stri_length(mastertable$miRNA_prec_seq)
mastertable$precursor_U <- str_count(mastertable$miRNA_prec_seq, "U")/stri_length(mastertable$miRNA_prec_seq)

mastertable$mature_5p_A <- str_count(mastertable$`5p_mat_seq`, "A")/stri_length(mastertable$`5p_mat_seq`)
mastertable$mature_5p_C <- str_count(mastertable$`5p_mat_seq`, "C")/stri_length(mastertable$`5p_mat_seq`)
mastertable$mature_5p_G <- str_count(mastertable$`5p_mat_seq`, "G")/stri_length(mastertable$`5p_mat_seq`)
mastertable$mature_5p_U <- str_count(mastertable$`5p_mat_seq`, "U")/stri_length(mastertable$`5p_mat_seq`)

mastertable$mature_3p_A <- str_count(mastertable$`3p_mat_seq`, "A")/stri_length(mastertable$`3p_mat_seq`)
mastertable$mature_3p_C <- str_count(mastertable$`3p_mat_seq`, "C")/stri_length(mastertable$`3p_mat_seq`)
mastertable$mature_3p_G <- str_count(mastertable$`3p_mat_seq`, "G")/stri_length(mastertable$`3p_mat_seq`)
mastertable$mature_3p_U <- str_count(mastertable$`3p_mat_seq`, "U")/stri_length(mastertable$`3p_mat_seq`)

write.csv(mastertable, paste0(outputfile,"mastertable_ACGU.csv"))
dfm <- melt(mastertable, id.vars=c("version", "species"),measure.vars=colnames(mastertable)[16:19])
dfm['DBversions'] = cut(as.numeric(as.character(dfm$version)),breaks = c(0,7,16,19,22))
colnames(dfm)[3] = "Nucleotide"
colnames(dfm)[4] = "Percentage"

for (s in unique(dfm$species)){
  p = ggplot(subset(dfm,species == s)) + geom_boxplot(aes(y=Percentage,fill=Nucleotide,x=DBversions))+ggtitle(s)
  ggsave(p,file = paste(s,"nulecotidePercentages.pdf",sep='_'))
}

