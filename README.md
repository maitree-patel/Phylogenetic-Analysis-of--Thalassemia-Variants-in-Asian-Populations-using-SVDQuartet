# Phylogenetic-Analysis-of-Beta-Thalassemia-Variants-in-Asian-Populations-using-SVDQuartet

Given a high prevalence in South-Asian populations and some information about common SNPs found in these regions I was interested in unraveling the phylogeny of β-thalassemia variants. In this study, I conducted a phylogenetic analysis of SNP data, i.e. variants from chromosome 11, as made available by the 1000 Genomes Project to understand the evolution of  variants found in the beta-globin gene region. The phylogenetic analysis was performed using SVDQuartets as a part of the PAUP package run through command line interface. SVDQuartet is a phylogenetic analysis method used for multi-locus or SNP data using a coalescent-based model.

## Data Availability
As mentioned above, variant data was accessed from The 1000 Genomes Phase 3 release. 1000 Genomes is a database with reference genomic variations from a diverse group of populations. Asian Populations of Punjabi in Lahore, Pakistan; Sri Lankan Tamil in UK; Japanese in Tokyo, Japan; Han Chinese South and Han Chinese in Beijing; Kinh in Ho Chi Minh City, Vietnam; Bengali in Bangladesh; Gujarati Indians in Texas and Indian Telugus in UK (Indian immigrant populations) were used in this study. The Asian populations were chosen based on data availability from Phase 3 release. The populations included are given below in Table 1 and graphically visualized in Figure 1.  Five samples from each population were selected based on its availability through the tool described below, i.e. a total of 50 samples were chosen for this analysis.

<img width="400" alt="Screenshot 2024-09-30 at 5 44 26 PM" src="https://github.com/user-attachments/assets/436f2d05-e241-4742-88e4-a4b727bfab74">

<img width="400" alt="Screenshot 2024-09-30 at 5 44 34 PM" src="https://github.com/user-attachments/assets/07a75b5f-4d36-4674-ad6f-a7477b2d20de">

The data was accessed through the 1000 Genomes International Genome Sample Resource (IGSR) data portal (10) and downloaded through the Ensembl Data Slicer Web Interface (11). The data slicer uses VCFTools to get specific regions of the genome using genomic coordinates that are described by the user. Hemoglobin β-subunit is found on chromosome 11 with location coordinates 11:5225464-5227071. Additionally, the coordinates of chromosome 11, i.e. 11:5177714-5243592 (that contains approximately 1,313 genes) as described by NCBI (12) were both input as a way to exploratorily analyze relationships between samples.

## Data Processing
In order to perform a phylogenetic analysis using SVDQuartets we require an input file in NEXUS format. To convert the vcf file downloaded, vcf2phylip was used to create a NEXUS file from the examples given as a part of the README file documentation on GitHub (13). Additionally during tree construction using FigTree (as will be described in detail below), The SVDQuartet tree was converted to .tre format (the input format required for FigTree) using the write.tree function from the ape package in R.

## Tree Construction using SVDQuartets
To perform the phylogenetic analysis using the SVDQuartet method, the program was run through command line interface with the eval parameter set to all (eval=all) to evaluate all quartets with bootstrap set to default (bootstrap=standard) which performs 100 bootstrap iterations. The resulting best tree is saved using savetrees with parameters brlens=yes to save the bootstrap support values and SupportValues=Both to also save branch length information and node labels. The tree was then plotted in R by reading, plotting and editing using functions from the phangorn and phytools packages. 
	
The total number of SNPs found in chromosome 11 was 2,807 and in the β-globin region was 55. SNPs from the β-globin region alone were first analyzed using SVDQuartets in the PAUP package. This resulted in a tree with a large polytomy. This means the information from the variants from the β-globin region does not give us enough information to determine the evolutionary relationship between the samples. This demanded for an increase in the sequence length for analysis, i.e. the use of the entire chromosome 11 region. 

As will be described below, the SVDQuartet tree was further visualized in FigTree v1.4.4 (15). It is a software used to view, visualize and annotate trees graphically. Two types of trees were visualized in the radial topology: i) by region and ii) by the presence of SNP in the β-globin region. This was done by manually color-coding sample branches in the tree based on the above two criteria. For the former criterion the populations from the same country or region were given the same color. Countries characterized as the South-East Asia region were visualized using warm colors whereas those that form South Asia were given cooler tones for better and broader construction of the relationships on the tree. The information about the presence or absence of SNP was analyzed from the genotype information in the VCF files in R using the vcfR package functions prior to annotating the tree based on the latter criteria.

## Results
Figure 2 shows the phylogenetic tree constructed using SNP data from chromosome 11. We can see a major clade consisting of South-East Asian populations (which includes the Chinese, Japanese and Vietnamese populations) and two major clades formed from South-Asian populations (Including Indian, Pakistani, Sri Lankan and Bangladeshi). This tree in itself does not provide enough information to interpret relationships between these samples. This is also true because rooting of the tree is not possible since they are all human samples and using another species as the root would translate to difference in the location of chromosome 11 and the β-globin genes along with the inconsistency in the characterisation of beta-globin gene variants that cause the disorder in that specie. Therefore, the tree was plotted in FigTree application to analyze trends within these major clades as described by the workflow above in section 2.4 Tree Construction using SVDQuartet. 

*Figure 2*

<img src="https://github.com/user-attachments/assets/6f4aeeae-e421-49cf-89b1-002f1e93ac48" width="550">

Figure 3 visualizes the phylogenetic tree based on regions of the populations of the samples. We can interpret the same as described from Figure 2 in better clarity and confidence using a radial tree. We can broadly see that most individuals from the South-East Asian populations that include Chinese, Japanese and Vietnamese samples cluster together on the top of the tree. And the South Asian Populations of Indian, Pakistani, Sri Lankan and Bangladeshi samples mostly cluster towards the end. These clusters are expected to be genetically similar given their histories of partition and geographic proximity. This can be due to higher degree of immigration between the populations (Peter et al., 2019) and shared ancestry. A study published in Nature analyzed SNPs from European populations and performed a Principal Component Analysis to find that populations found together on the map clustered closer together on the PCA plot (Novembre et al., 2008) which might also be the reason for clustering of populations in certain clades in Figure 3.

*Figure 3*
<img width="700" alt="Screenshot 2024-09-30 at 5 56 02 PM" src="https://github.com/user-attachments/assets/f6384250-53c0-49e9-857a-8971532ee981">




