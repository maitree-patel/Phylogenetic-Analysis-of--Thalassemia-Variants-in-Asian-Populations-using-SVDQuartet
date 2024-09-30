library(vcfR)

#reading in vcf file using vcfR package function
chr11.HBB <- vcfR::read.vcfR("/Users/maitreepatel/Desktop/CME/Chr11HBB.vcf")

install.packages("remotes")
remotes::install_github("shankarkshakya/mypackage")

vcf2nexus <- function(vcf, file = "file.nex") {
  
  vcf <- extract.indels(vcf)
  vcf <- vcf[is.biallelic(vcf),]
  
  gt.filtered <- extract.gt(vcf, element = "GT", as.numeric = T, convertNA = T)
                            #, return.alleles=T)
  
  gt.filtered[is.na(gt.filtered)] <- "?"
  
  gt.filtered <- t(gt.filtered)
  
  ape::write.nexus.data(gt.filtered, file)
  
  nex.file <- scan(file, what = "character", sep = "\n",
                   quiet = TRUE)
  
  # bgn <- grep("BEGIN", snapp.file)
  # snapp.file[bgn] <- "BEGIN CHARACTERS;"
  fmt <- grep("FORMAT", nex.file)
  nex.file[fmt] <- "  FORMAT DATATYPE=INTEGER MISSING=? GAP=- SYMBOLS=\"012\" LABELS=LEFT TRANSPOSE=NO INTERLEAVE=NO;"
  
  #return(snapp.file)
  #write(snapp.file, file)
  
}

vcf2nexus(chr11.HBB,
          file = "Chr11_HBB.nex")



library(phangorn)
library(phytools)
#reading nexus file
hbb.nexus <- read.nexus.data("/Users/maitreepatel/Desktop/CME/Chr11_HBB.nex")
hbb.df <- as.data.frame(hbb.nexus[1:10])
class(hbb.df)
hbb.nj <- nj(dist.gene(hbb.df))
plot(hbb.nj)

BiocManager::install("adegenet")
library(adegenet)
library(vcfR)
hbb.vcf <- read.vcfR("/Users/maitreepatel/Downloads/11.5177714-5243592.ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf",
                     verbose = FALSE)
hbb.genlight <- vcfR2genlight(hbb.vcf)
hbb.genlight@gen
dist.matrix <- dist.gene(as.matrix(hbb.genlight),
                         method = "percentage")
hbb.nj <- nj(dist.matrix)
plot(hbb.nj)

beta.vcf <- read.vcfR("/Users/maitreepatel/Desktop/CME/betasubunit.vcf")
beta.genlight <- vcfR2genlight(beta.vcf)
dist.beta <- dist.gene(as.matrix(beta.genlight),
                       method = "percentage")
beta.nj <- nj(dist.beta)
plot(beta.nj)


#svdq tree
library(phangorn)
library(phytools)
svd.hbb.tree <- read.nexus("/Users/maitreepatel/Desktop/CME/paup_hbb_tree.nex")
labels <- read.csv("/Users/maitreepatel/Downloads/CME - Sheet1.csv")
svd.hbb.tree$tip.label <- paste0(labels$Population, "-" , labels$Sample)
#labels$Sample, "-", 

plotTree(svd.hbb.tree,
         edge.width = 2,
         font = 1)
nodelabels(svd.hbb.tree$node.label,
           cex=0.8)
save.plot(file="Rtree_hbb.png")
write.tree(svd.hbb.tree,
           file = "svd_hbb_tree.tre")
#FigTree
#colour code by region - characterise regions
  #plot it on a map
#colour code by snp?
  #are these variants significant
#radial tree
#plot a map (same colours on the tree)

#also try getting a tree from regions around the hbb (10,000 or a million bp on each side)


#radial tree using phytools
library(phytools)
hbb.tre <- read.tree(file="/Users/maitreepatel/Desktop/CME/svd_hbb_tree.tre")

fanCladogram<-function(tree,use.edge.length=FALSE,...){
  if(use.edge.length==FALSE){
    if(hasArg(power)) power<-list(...)$power
    else power<-0.6
    tree<-compute.brlen(tree,power=power)
  }
  args<-list(...)
  args$power<-NULL
  args$tree<-tree
  args$type<-"fan"
  args$color<-"transparent"
  do.call(plotTree,args)
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  X<-cbind(obj$xx[1:Ntip(tree)],
           obj$yy[1:Ntip(tree)])
  rownames(X)<-tree$tip.label
  A<-cbind(obj$xx[1:tree$Nnode+Ntip(tree)],
           obj$yy[1:tree$Nnode+Ntip(tree)])
  rownames(A)<-1:tree$Nnode+Ntip(tree)
  pp.args<-list(...)
  pp.args$power<-NULL
  pp.args$tree<-tree
  pp.args$X<-X
  pp.args$A<-A
  pp.args$ftype<-"off"
  pp.args$node.size<-c(0,0)
  pp.args$xlim<-obj$x.lim
  pp.args$ylim<-obj$y.lim
  pp.args$xaxt<-"n"
  pp.args$type<-NULL
  pp.args$add<-TRUE
  do.call(phylomorphospace,pp.args)
}

fanCladogram(hbb.tre,lwd=1,ftype="i",fsize=0.3)

install.packages("ggtree")
library(ggtree)

# Assuming you have a phylogenetic tree object called 'tree' and SNP data called 'snp_data'
# Annotate the tree with SNP data
library(vcfR)
vcf.file <- read.vcfR("chr11_hbb.vcf")
snp_data <- extract.gt(vcf.file, return.alleles = T)
names <- rownames(snp_data)
names[1500:5000]
nrow(snp_data)
p <- ggtree(svd.hbb.tree) + geom_tippoint(aes(color = snp_data))
p



snps <- c("11_5225477_2038", "11_5225483_2039", "11_5225492_2040", "11_5225505_2041", "11_5225519_2042", 
          "11_5225558_2043", "11_5225635_2044", "11_5225660_2045", "11_5225698_2046", "11_5225714_2047", 
          "11_5225757_2048", "11_5225765_2049", "11_5225766_2050", "11_5225900_2051", "11_5225901_2052",
          "11_5225917_2053", "11_5225941_2054", "11_5225943_2055", "11_5225966_2056", "11_5225989_2057",
          "11_5226006_2058", "11_5226025_2059", "11_5226096_2060", "11_5226121_2061", "11_5226186_2062", 
          "11_5226251_2063", "11_5226258_2064", "11_5226274_2065", "11_5226279_2066", "11_5226314_2067", 
          "11_5226316_2068", "11_5226318_2069", "11_5226341_2070", "11_5226359_2071", "11_5226403_2072", 
          "11_5226469_2073", "11_5226471_2074", "11_5226530_2075", "11_5226534_2076", "11_5226578_2077",
          "11_5226579_2078", "11_5226583_2079", "11_5226586_2080", "11_5226655_2081", "11_5226725_2082",
          "11_5226745_2083", "11_5226777_2084", "11_5226865_2085", "11_5226890_2086", "11_5226939_2087", 
          "11_5226965_2088", "11_5227041_2089")
length(snps)

positions <- match(snps, rownames(snp_data))
filtered <- snp_data[positions,]
mutants <- c()
for (snp in 1:nrow(filtered)) {
  for (sample in 1:ncol(filtered)) {
    if (filtered[snp, sample]=="1|1") {
      mutants <- c(mutants, colnames(filtered)[sample])
    } else if (filtered[snp, sample]=="1|0") {
      mutants <- c(mutants, colnames(filtered)[sample])
    } else if (filtered[snp, sample]=="0|1") {
      mutants <- c(mutants, colnames(filtered)[sample])
    } else {
      mutants <- mutants
    }
  }
  print(mutants)
}
samples <- unique(mutants)
labels[match(samples, labels$Sample),1:2]

beta <- read.vcfR("betasubunit.vcf")
snp_beta <- extract.gt(beta)
nrow(snp_beta)
