library(cgdsr)

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

# Get available case lists (collection of samples) for a given cancer study
#Choose coadread_tcga_pub
mycancerstudy = getCancerStudies(mycgds)[36,1]
#choose "coadread_tcga_pub_nonhypermut_all" 
mycaselist = getCaseLists(mycgds,mycancerstudy)[6,1]

# Get available genetic profiles
#Choose coadread_tcga_mutations
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[1,1]

GeneList.BCL9L <- c("BCL9L")


# Get data slices for a specified list of genes, genetic profile and case list
Mutations <- getProfileData(mycgds, GeneList.BCL9L, mygeneticprofile, mycaselist)
is.na(Mutations) <- Mutations == "NaN"

NonWhites <- c("TCGA.AA.3858.01","TCGA.AF.3913.01","TCGA.AY.4070.01")

Mutations.Whites <- Mutations[!(rownames(Mutations) %in% NonWhites),]
library(dplyr)
summary(Mutations.Whites)
which(!is.na(Mutations.Whites))
#apply(Mutations.Whites, 2, function(x) length(which(!is.na(x))))

SampleswithMutationsW <- names(which(!is.na(Mutations.Whites)))
SampleswithMutations <- names(which(!is.na(Mutations)))

#Get clinical data
Clinical <- getClinicalData(mycgds, mycaselist)
ClinicalforSamples <- Clinical[]


Mutations.Fig2a <- Mutations.Whites[,colnames(Mutations.Whites) %in% GeneList.Fig2a]
Mutation.Counts.2a <- as.data.frame(apply(Mutations.Fig2a, 2, function(x) length(which(!is.na(x)))))
colnames(Mutation.Counts.2a) <- c("TCGA_mut")
Mutation.Counts.2a$TCGA_norm <- 186 - Mutation.Counts.2a$TCGA_mut
Mutation.Counts.2a$CCCC_mut <- c(0, 0,27, 1,0,4,0,0,22,0,3,0,5,0,2,24,12)
Mutation.Counts.2a$CCCC_norm <- 43 - Mutation.Counts.2a$CCCC_mut

get_fisher <- function(df){
     mat <- matrix(as.numeric(unlist(df[c(1:4)])), ncol=2)
     f <- fisher.test(as.table(mat), alt="two.sided")
     return(f$p.value)
}

fishers <- apply(Mutation.Counts.2a, 1,  get_fisher)
Mutation.Counts.2a$FET_p <- round(fishers,5)
Mutation.Counts.2a$Bonf_p <- round(p.adjust(fishers, method="bonf"), 5)

Mutation.Counts.2a$TCGA_perc <- round(Mutation.Counts.2a$TCGA_mut / (Mutation.Counts.2a$TCGA_mut + Mutation.Counts.2a$TCGA_norm) * 100, 3)
Mutation.Counts.2a$CCCC_perc <- round(Mutation.Counts.2a$CCCC_mut / (Mutation.Counts.2a$CCCC_mut + Mutation.Counts.2a$CCCC_norm) * 100, 3)
write.table(Mutation.Counts.2a, "../../2016_Exome/2016-07-26_ExomePaper_Fig2a.txt", quote=FALSE, col.names = TRUE, row.names = TRUE, sep="\t")


library(tidyr)
Fig2a <- Mutation.Counts.2a %>% mutate(Gene = rownames(Mutation.Counts.2a)) %>% select(-c(TCGA_norm, CCCC_norm)) %>% gather(Sample, Percent.Mutated, TCGA_perc, CCCC_perc)

png("../../2016_Exome/2016-07-26_Figure2a.png", height=4, width=6, units = "in", res=300)
ggplot(Fig2a, aes(x=Gene, y=Percent.Mutated)) + geom_bar(aes(fill=Sample), position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_fill_manual(values=c("#1b9e77","#7570b3"), name="Cohort", labels=c("CRC","TCGA")) + ylab("Percent of Patients with Mutations")
graphics.off()



Mutations.Fig2b<- Mutations.Whites[,colnames(Mutations.Whites) %in% GeneList.Fig2b]
Mutation.Counts.2b <- as.data.frame(apply(Mutations.Fig2b, 2, function(x) length(which(!is.na(x)))))
colnames(Mutation.Counts.2b) <- c("TCGA_mut")
Mutation.Counts.2b$TCGA_norm <- 186 - Mutation.Counts.2b$TCGA_mut
Mutation.Counts.2b$CCCC_mut <- c(2,4,3,2,2,1,3,4,2,3,5,2,2)
Mutation.Counts.2b$CCCC_norm <- 43 - Mutation.Counts.2b$CCCC_mut

fishers <- apply(Mutation.Counts.2b, 1,  get_fisher)
Mutation.Counts.2b$FET_p <- round(fishers,5)
Mutation.Counts.2b$Bonf_p <- round(p.adjust(fishers, method="bonf"), 5)

Mutation.Counts.2b$TCGA_perc <- round(Mutation.Counts.2b$TCGA_mut / (Mutation.Counts.2b$TCGA_mut + Mutation.Counts.2b$TCGA_norm) * 100, 3)
Mutation.Counts.2b$CCCC_perc <- round(Mutation.Counts.2b$CCCC_mut / (Mutation.Counts.2b$CCCC_mut + Mutation.Counts.2b$CCCC_norm) * 100, 3)
write.table(Mutation.Counts.2b, "../../2016_Exome/2016-07-26_ExomePaper_Fig2b.txt", quote=FALSE, col.names = TRUE, row.names = TRUE, sep="\t")


Fig2b <- Mutation.Counts.2b %>% mutate(Gene = rownames(Mutation.Counts.2b)) %>% select(-c(TCGA_norm, CCCC_norm)) %>% gather(Sample, Percent.Mutated, TCGA_perc, CCCC_perc)

png("../../2016_Exome/2016-07-26_Figure2b.png", height=4, width=6, units = "in", res=300)
ggplot(Fig2b, aes(x=Gene, y=Percent.Mutated)) + geom_bar(aes(fill=Sample), position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_fill_manual(values=c("#1b9e77","#7570b3"), name="Cohort", labels=c("CRC","TCGA")) + ylab("Percent of Patients with Mutations")
graphics.off()


Mutations.Fig2c <- Mutations.Whites[,colnames(Mutations.Whites) %in% GeneList.Fig2c]
Mutation.Counts.2c <- as.data.frame(apply(Mutations.Fig2c, 2, function(x) length(which(!is.na(x)))))
colnames(Mutation.Counts.2c) <- c("TCGA_mut")
Mutation.Counts.2c$TCGA_norm <- 186 - Mutation.Counts.2c$TCGA_mut
Mutation.Counts.2c$CCCC_mut <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,3,0,0,1,0,0)
Mutation.Counts.2c$CCCC_norm <- 43 - Mutation.Counts.2c$CCCC_mut

fishers <- apply(Mutation.Counts.2c, 1,  get_fisher)
Mutation.Counts.2c$FET_p <- round(fishers,5)
Mutation.Counts.2c$Bonf_p <- round(p.adjust(fishers, method="bonf"), 5)

Mutation.Counts.2c$TCGA_perc <- round(Mutation.Counts.2c$TCGA_mut / (Mutation.Counts.2c$TCGA_mut + Mutation.Counts.2c$TCGA_norm) * 100, 3)
Mutation.Counts.2c$CCCC_perc <- round(Mutation.Counts.2c$CCCC_mut / (Mutation.Counts.2c$CCCC_mut + Mutation.Counts.2c$CCCC_norm) * 100, 3)
write.table(Mutation.Counts.2c, "../../2016_Exome/2016-07-26_ExomePaper_Fig2c.txt", quote=FALSE, col.names = TRUE, row.names = TRUE, sep="\t")


Fig2c <- Mutation.Counts.2c %>% mutate(Gene = rownames(Mutation.Counts.2c)) %>% select(-c(TCGA_norm, CCCC_norm)) %>% gather(Sample, Percent.Mutated, TCGA_perc, CCCC_perc)

png("../../2016_Exome/2016-07-26_Figure2c.png", height=4, width=6, units = "in", res=300)
ggplot(Fig2c, aes(x=Gene, y=Percent.Mutated)) + geom_bar(aes(fill=Sample), position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_fill_manual(values=c("#1b9e77","#7570b3"), name="Cohort", labels=c("CRC","TCGA")) + ylab("Percent of Patients with Mutations")
graphics.off()