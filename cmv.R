
# Author : Osman Merdan
# Date : 
#......................................................................#
# HUMAN HERPES VIRUS 5 NGS DATA ANALYSIS SUMMARY STATISTICS
#......................................................................#

# ======================== #
## 1- Loading Libraries ####
# =======================  #
suppressMessages(library(tidyverse))
suppressMessages(library(ggthemes))
suppressMessages(library(RColorBrewer))
suppressMessages(library(wesanderson))
suppressMessages(library(ggpubr))

# ======================================= #
## 2- Trimmed Fastq Summary Statistics ####
# ======================================= #

# MultiQC created multiqc_fastqc.txt file. 
# This file contains summary statistics about trimmed raw sequence data.
# Loading the summary statistics 
data <- read_delim(file = "~/cmv-project/qc/multiqc_data/multiqc_fastqc.txt",
                   delim = "\t", 
                   col_names = T)
# 2nd,3rd, and 4th columns includes metadata.
data <- data[,-c(2,3,4)]
# Separation of Forward(1) and Reverse(2) reads into Pair column
data <- separate(data = data,
                 col = Sample,
                 into = c("Sample","Pair"),
                 sep = "_paired_")
data

# =========================================== #
## 3- Visualizing Raw Sequence Counts Plot ####
# =========================================== #

# MultiQC HTML creates fastqc_sequence_counts_plot. 
# Using raw data of this plot create rmarkdown friendly sequence count plot. 
# Visualizing important variables.
read_delim("fastqc_sequence_counts_plot.tsv",
           delim = "\t",
           col_names = T)%>%
  pivot_longer(cols = c("Unique Reads","Duplicate Reads"),
               names_to = "Filter",
               values_to = "Number of Reads")%>%
  ggplot(aes(x = Category, y = `Number of Reads`, fill= Filter))+
  geom_col (stat="identitiy")+
  scale_fill_brewer()+
  theme_minimal()+
  labs(x="Sample", title = "Number of Reads Per Sample")+
  theme(axis.text.x = element_text(angle = 90),
        legend.key.size = unit(x=4, units = "mm"),
        plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = c(100000,20000),
             linetype = "dashed",
             linewidth = 0.2,
             color=c("black","red"))

# ====================================== #
## 4- Cumulative Genome Coverage Plot ####
# ====================================== #

# Empirical cumulative density function visualization 
# For just low-coverage, low-depth samples 
# In present project those samples were:
# "ERR7018443","ERR7018474","ERR7021883","ERR7029111","ERR7039821"
# Those sample names were listed in low-meandepth-samples.txt file.
x<-read_delim(
  "~/cmv-project/deduplicated-bam/low-meandepth-samples/low-meandepth-samples.txt",
  delim = "\t",
  col_names = F)
x <- as_vector(x)
# Load data
data<-read_delim("~/cmv-project/coverage/piled-depth.txt",
                 delim = "\t",
                 col_names = F)%>%
  filter(X3 %in% c(x))
# Visualize
ggplot(data=data, aes(x=X2, color=X3))+
  geom_line(aes(y = 1-..y..), # y=1-..y.. creates complementary ecdf plot
            stat = "ecdf",
            size=0.7)+
  coord_cartesian(xlim=c(0,50))+
  labs(x="Depth (X)",
       y="Fraction of Reference",
       title = "Genome Fraction Covered By At Least X Reads")+
  theme_minimal()+
  scale_color_brewer("Samples",palette = "Paired")+
  theme(legend.key.size = unit(x=5, units = "mm"),
        plot.title = element_text(hjust = 0.5))

# ======================= #
## 5- Coverage Heatmap ####
# ======================= #

# Important Note: This plot requres lots of computing power and memory
# Reading piled up coverage data 
coverage <- read_delim("~/cmv-project/coverage/piled-depth.txt", 
                       col_names = FALSE)

# Adding proper column names
colnames(coverage)<- c("Position", "Depth", "SampleID")

# Checking for sample ID for any problems.
unique(coverage$SampleID)

# In order to simplify the coverage heatmap -->
# a) converting depth coverage values to log values 
# b) and rounding up to nearest integer
# Base is chosen according to min depth treshold. 
coverage$Depth <- sapply(X=coverage$Depth,
                         FUN = function(x) round(log(x+1, base = 10),1))
colnames(coverage)[2]<-"Depth(Log10)"
# Plot the coverage
plt<-ggplot(data=coverage,
       aes(x= Position, y= SampleID, fill=`Depth(Log10)`))+
  scale_fill_viridis_c(option = "D")+
  theme_minimal()+
  geom_tile()+
  labs(title = "Depth of Positions")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        title = element_text(size = 10),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 7))
# Save the image
ggsave(filename="coverage-heatmap.png",
       device = "png",
       dpi = "print",
       plot = plt,
       width = unit(8,"cm"),
       height = unit(4,"cm"))

# ===================================================== #
## 6- Coverage Plot For Specific Regions of Interest ####
# ===================================================== #

# Important Note: This plot requres lots of computing power and memory
# Reading piled up coverage data 
coverage <- read_delim("~/cmv-project/coverage/piled-depth.txt", 
                       col_names = FALSE,
                       delim = "\t")

# Adding proper column names
colnames(coverage)<- c("Position", "Depth", "SampleID")

# Filter data for position and samples
# Create plot
plt<-coverage%>%
  filter(SampleID %in% c("ERR7036363",
                         "ERR7032218",
                         "ERR7024258"))%>%
  filter(Position>5000 & Position<13000)%>%
  ggplot(aes(x=Position,y=Depth,color=SampleID))+
  geom_line()+
  theme_minimal()+
  labs(title = "Depth of Positions")+
  scale_color_brewer(palette = "Paired")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        title = element_text(size = 10),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
# Save the image
ggsave(filename="coverage-limited.png",
       device = "png",
       dpi = "print",
       plot = plt,
       width = unit(8,"cm"),
       height = unit(4,"cm"))


# ===================================================== #
## 7- Heat map For Number of Non-synonymous Variants ####
# ===================================================== #

# Load data created with bash script. 
# This data has all annotated variants for each sample.
data<-read_delim("~/cmv-project/extracted-vcf/merged-variants.txt",
                 delim = "\t",
                 col_names = FALSE)[,-1]

# Rename columns
colnames(data)<-c("Position", "ID", "REF", "ALT",
                  "FILTER", "EFFECT", "IMPACT", 
                  "GENE", "GENEID", "FEATURE",
                  "FEATUREID", "HGVS_C", "HGVS_P",
                  "SampleID")

# What kind of unique effects annotated with SnpEff. 
unique(data[,"EFFECT"])

# Plot 
data%>%
  filter(EFFECT!= "synonymous_variant")%>%
  count(SampleID, GENE)%>%
  complete(SampleID, GENE, fill = list(n=0))%>%
  ggplot(aes(x= GENE,y= SampleID))+
  geom_tile(aes(fill= n), color="black")+
  scale_fill_distiller(palette = "RdYlBu")+
  coord_fixed()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 3),
        axis.text.y = element_text(size = 5))+
  guides(fill=guide_colorbar(barwidth = 0.2, barheight = 5))

# ===================================================== #
## 8- Summarizing Number of Total Variations.        ####
# ===================================================== #

# Load data created with bash script. 
# This data has all annotated variants for each sample.
data<-read_delim("~/cmv-project/extracted-vcf/merged-variants.txt",
                 delim = "\t",
                 col_names = FALSE)[,-1]

# Rename columns
colnames(data)<-c("Position", "ID", "REF", "ALT",
                  "FILTER", "EFFECT", "IMPACT", 
                  "GENE", "GENEID", "FEATURE",
                  "FEATUREID", "HGVS_C", "HGVS_P",
                  "SampleID")

# Summarize variation counts
data%>%
  group_by(SampleID, EFFECT)%>%
  count()%>%rename("Number of Variations"=n)


# ===================================================== #
## 9- Plotting UL8 Gene Mutations ####
# ===================================================== #

# Get the data 
data<-read_delim("~/cmv-project/extracted-vcf/merged-variants.txt",
                 delim = "\t",
                 col_names = FALSE)[,-1]
colnames(data)<-c("Position", "ID", "REF", "ALT",
                  "FILTER", "EFFECT", "IMPACT", 
                  "GENE", "GENEID", "FEATURE",
                  "FEATUREID", "HGVS_C", "HGVS_P",
                  "SampleID")
# Plot how many samples have non-synonymous variants for each position in UL8.
# color, fill by EFFECT column.
# Save plot as plot1 for later combinind plots
plot1<- data%>%
filter(EFFECT!="synonymous_variant" & GENE == "UL8")%>%
  ggplot(aes(x=Position, fill=EFFECT, color=EFFECT))+
  geom_bar()+ # geom bar stat count 
  theme_minimal()+
  labs(y="Number of Samples")+
  theme(legend.position = "top",
        legend.key.size = unit(0.7,"lines"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette= "Set1")+
  # remove x axis
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Create a genes data frame to use for gggenes plot.
# For detailed instructions: https://github.com/wilkox/gggenes
# This step can be automatized for large number of samples using gff files.
genes <- data_frame(molecule = "AY446894", # molecule name (chromosome)
                    gene = "UL8", # gene names (as vector)
                    start = 15865, # start position for each gene (as vector)
                    end = 16941,  # end position for each gene (as a vector)
                    strand = "forward", # strans ("forward" or "reverse")
                    # subgenes or domains in the proteins (as vector)
                    subgene = c("Transmembrane-1",
                                "Ig-Domain",
                                "Transmembrane-2",
                                "Ig strand B (motif)",
                                "Ig strand E (motif)",
                                "Ig strand F (motif)"),
                    # from = where the subgenes start using amono acid counts
                    from = c(15865+(22*3),
                             15865+(41*3),
                             16552+(84*3),
                             15865+(50*3),
                             15865+(98*3),
                             15865+(112*3)),
                    # to = where the subgenes end
                    # Transmembrane-1 domain spans between 23-45th amino acids
                    to = c(15865+(45*3)-1,
                           15865+(114*3)-1,
                           16552+(105*3)-1,
                           15865+(55*3)-1,
                           15865+(103*3)-1,
                           15865+(118*3)-1),
                    orientation = 1)

# Using gggenes package create gene tracks. (https://github.com/wilkox/gggenes)
# Save it as plot2 for later use.
library(gggenes)
plot2 <- ggplot(genes, aes(xmin = start, # start positions for genes
                           xmax = end, # end positions for genes
                           y = molecule, # y axis label (chromosome name)
                           label = gene))+ # add gene name lebels
  # place gene arrows
  geom_gene_arrow(fill = "white",
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(1, "mm"))+
  # place subgene segments
  geom_subgene_arrow(data = genes,aes(xmin = start, # for genes
                                      xmax = end, # for genes
                                      y = molecule,
                                      fill = subgene, # fill by subgene
                                      xsubmin = from, # start for each subgene
                                      xsubmax = to),# end for each subgene
                     color="black", # contour color for subgenes
                     alpha=.7) +
  geom_gene_label(align = "centre")+ # set gene label position
  scale_fill_brewer(palette = "Set3")+ # set subgene fill palette
  theme_genes()+ # add pretty genes theme
  theme(legend.position = "bottom",
        legend.key.size = unit(0.7,"lines"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  labs(y=NULL, x="Position")

# Using ggpubr::ggarrange combine two plots and align them by x axis
ggpubr::ggarrange(plot1,
                  plot2,
                  ncol = 1,
                  nrow = 2,
                  align = "v", # align vertically by x axis
                  heights = c(2,1))  # upper plot will be 2x higher
detach("package:gggenes", unload = TRUE)

# ===================================================== #
## 10- Plotting Resistance Mutations ####
# ===================================================== #

# Using trackViewer package create lollipop plots for mutations in 
# UL97 and UL54.
# Below example shows it only for UL54.
library(trackViewer)
# Detailed use here --> 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/trackViewer/inst/doc/lollipopPlot.html#Lolliplot
# https://jianhong.github.io/trackViewer/articles/lollipopPlot.html#jitter-the-label-1
# Get data
data<-read_delim("~/cmv-project/extracted-vcf/merged-variants.txt",
                 delim = "\t",
                 col_names = FALSE)[,-1]
colnames(data)<-c("Position", "ID", "REF", "ALT",
                  "FILTER", "EFFECT", "IMPACT", 
                  "GENE", "GENEID", "FEATURE",
                  "FEATUREID", "HGVS_C", "HGVS_P",
                  "SampleID")
# Filter and create counts amino acid change and effect for each position
data <- data%>%
  filter(GENE == "UL54")%>%
  filter(EFFECT != "synonymous_variant")%>%
  group_by(Position,HGVS_P,EFFECT)%>%count()

# Create Number of Samples GRanges object.
# This object defines mutations.
# Reason to name it like this is --> 
# because this object name will be used to name y axis.
`Number of Samples` <- GRanges(seqnames = "AY446894", #chromosome 
                               # mutation positions
                               ranges = IRanges(data$Position,
                                                width = 1,
                                                names = data$HGVS_P),
                               # score = hight of lollipops
                               score = data$n)

# Add color metadata for Number of Samples object.
# Modify the code accordingly.
`Number of Samples`$color <- ifelse(test = data$HGVS_P %in% c("p.Gln330Aln",
                                                              "p.Thr1122Ile"),
                                    yes = "coral1",
                                    no = 'black')

# Create features.gr Granges object to display genes/features
features.gr <- GRanges(seqnames = "AY446894", # chromosome 
                       ranges = IRanges(start = 78194, # gene start positions
                                        end = 81922, # gene end positions
                                        names="UL54"), # gene names
                       # fill color for genes as a vector
                       # this vector must be equal size with gene names
                       fill = c("#DFA32D"),  
                       height = c(0.03),
                       strand = "-")

# Legend for EFFECT types.
# Legend must be a named vector.
m<-list(labels=c("Novel","Previously characterized"),
        fill=c("coral1",'black'))

lolliplot(`Number of Samples`, # First argument is mutations obj
          features = features.gr, # Second argument is features obj
          cex = 0.9, # controls the lollipop ball sizes
          legend=m) # displays the legend
detach("package:trackViewer", unload = TRUE)

# ============================================================= #
## 11- Finding Samples With Double AA Change in One Position ####
# ============================================================= #

# Variant analysis suggests that:
# some samples have two concequtive SNPs in UL54
# Those two SNPs were pictured as two different amino acid changes 
# (Thr1122Ile and Thr1122Ala) in the annotated vcf files.
# Finding samples which has two amino acid change in same position
# Load the data

data<-read_delim("~/cmv-project/extracted-vcf/merged-variants.txt",
                 delim = "\t",
                 col_names = FALSE)[,-1]
colnames(data)<-c("Position", "ID", "REF", "ALT",
                  "FILTER", "EFFECT", "IMPACT", 
                  "GENE", "GENEID", "FEATURE",
                  "FEATUREID", "HGVS_C", "HGVS_P",
                  "SampleID")
data <- data%>%
  filter(GENE == "UL54")%>% # filter for UL54 gene
  filter(EFFECT != "synonymous_variant") # filter synonymous variants

Ile<-(data%>% # Samples which have "p.Thr1122Ile" mutation
        filter(HGVS_P=="p.Thr1122Ile"))$SampleID
Ala<-(data%>% # Samples which have "p.Thr1122Ala" mutation
        filter(HGVS_P=="p.Thr1122Ala"))$SampleID
# Intersections
print("== Samples p.Thr1122Ile and p.Thr1122Ala mutations were detected ==")
Reduce(f = intersect, x = list(Ile,Ala))
# Difference
print("== Samples only p.Thr1122Ala mutation were detected ==")
setdiff(Ala,Ile)
# Samples with neither mutation
print("Samples neither mutation is detected")
setdiff(data$SampleID, c(Ala,Ile))


# ===================================================== #
## 12- Plotting Genotype  ####
# ===================================================== #


 
