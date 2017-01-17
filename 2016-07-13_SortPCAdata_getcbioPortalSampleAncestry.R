library(cgdsr)

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[28,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[6,1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[7,1]

# Get data slices for a specified list of genes, genetic profile and case list
Mutations <- getProfileData(mycgds,c('APC','TP53'), mygeneticprofile, mycaselist)

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)


#Import 189 used TCGA samples
TCGA_189_SAMPLES <- read.delim("~/Documents/DataFromOthers/Rosa/TCGA_189_SAMPLES.txt", header=FALSE, stringsAsFactors=FALSE)

#Import sample map
Filename_SampleID_Map_Simplified <- read.delim("E:/TCGA_ClinicalData/Filename_SampleID_Map_Simplified.txt", stringsAsFactors=FALSE)

#Get list of samples that are in both.  Do all of the ones we're interested in have Affy 6.0 data?
exist <- Filename_SampleID_Map_Simplified[Filename_SampleID_Map_Simplified$Comment..TCGA.Barcode. %in% TCGA_189_SAMPLES$V3,]

#Import PCA data from all TCGA samples
PCA_TCGASamples <- read.delim("E:/COMPUTER_BACKUP-last-2016-02-09/_2015_CNV/2015-10-06_PCA-output_TCGA_Whites/2016-07-13_PCA_Output_AllSamples.txt", stringsAsFactors=FALSE)


#Add Sample Barcode to PCA_TCGASamples
SampleIDs <- vector(mode="character",length = 1411)
for(i in 1:nrow(PCA_TCGASamples)){
     
     ID <- Filename_SampleID_Map_Simplified$Comment..TCGA.Barcode.[ paste0(Filename_SampleID_Map_Simplified$Hybridization.Name, ".txt") == PCA_TCGASamples[i,1] ]
     if(length(ID) == 0){
          ID <- NA
     }
     SampleIDs[i] <- ID
     rm(ID)
}

PCA_TCGASamples$Barcode <- SampleIDs
PCA_TCGASamples$The189 <- NA

for(i in 1:nrow(PCA_TCGASamples)){
     if( PCA_TCGASamples$Barcode[i] %in% TCGA_189_SAMPLES$V3 ){
          PCA_TCGASamples$The189[i] <- 1
     }
     else{ PCA_TCGASamples$The189[i] <- 0}
}

library(ggplot2)
library(dplyr)

PCA_TCGASamples2 <- PCA_TCGASamples %>% arrange(Barcode) %>% distinct(Barcode)

ggplot(PCA_TCGASamples2, aes(x=EV...26.9183, y=EV...4.22101, col=factor(The189))) + geom_jitter() + scale_colour_manual(values=c("gray70","red"), name=c("Dataset"), labels=c("Other","cbioPortal")) + ggtitle("PCA of all TCGA Samples with Affy 6.0 data") + xlab("PC 1 (African Ancestry)") + ylab("PC 2 (Asian Ancestry)") + geom_text(data=PCA_TCGASamples2[PCA_TCGASamples2$Barcode %in% TCGA_189_SAMPLES$V3 & (PCA_TCGASamples2$EV...26.9183 < -0.02 | PCA_TCGASamples2$EV...4.22101 > 0.05),], aes(label=Barcode), hjust=0.5, vjust=-0.8)

NonWhites <- c("TCGA-AA-3858","TCGA-AF-3913","TCGA-AY-4070")
