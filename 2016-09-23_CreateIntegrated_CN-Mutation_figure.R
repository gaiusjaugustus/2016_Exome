# This file is to create the mutation portion of the INTEGRATED ANALYSIS figure for
# the Exome paper

#########
# Using 2016-09-13_GainLossAcrossGenome_Plot.R as my guide from the CNV folder
#########

# Import annotated mutation file
library(data.table)
Mutations <- fread("C:/Users/gaugustus/Documents/DataFromOthers/Zarko/2016-09-26_Results_of_all_available_samples.mutations.geneannotated.txt")

Genes <- Mutations %>% select(Hugo_Symbol, Chromosome, start_position, end_position) %>% arrange(Hugo_Symbol) %>% distinct()


#Some genes didn't come through for some reason
MissingGenes <- Summary2[is.na(Summary2$start_position) & Summary2$Count > 2,]
library(biomaRt)
ensembl <-  useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes2 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "external_gene_name"), filters = c("external_gene_name"), mart = ensembl, values = MissingGenes$Hugo_Symbol )

Genes[Genes$Hugo_Symbol == "AC090825.1",] <- c("AC090825.1", 15, 99792297, 99792379)
Genes[Genes$Hugo_Symbol == "RP11-32B5.7",] <- c("RP11-32B5.7", 15, 21298233, 21325241)
Genes[Genes$Hugo_Symbol == "RP11-337C18.8",] <- c("RP11-337C18.8", 1, 147172771, 147211568)
Genes[Genes$Hugo_Symbol == "RP11-114H24.7",] <- c("RP11-114H24.7", 15, 77916522, 77922019)
Genes[Genes$Hugo_Symbol == "RP11-403I13.8",] <- c("RP11-403I13.8", 1, 145164099, 145216058)




# Need to count up how many mutations there are per gene

Summary <- Mutations %>% group_by(Hugo_Symbol) %>% summarise(Count = length(Hugo_Symbol))
Summary2 <- merge(Summary, Genes, by="Hugo_Symbol", all.x=TRUE)
Summary2 <- Summary2 %>% filter(!Chromosome %in% c("X","Y"))
Summary2$start_position <- as.numeric(Summary2$start_position)
Summary2$end_position <- as.numeric(Summary2$end_position)
Summary2$Chromosome <- factor(Summary2$Chromosome, levels=c(1:22))
Summary2$Position <- rowMeans(Summary2[, c("start_position", "end_position")])


SummarybyPatient <- Mutations %>% group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarise(Count = length(Hugo_Symbol), Proportion = Count/43)
SummarybyPatient2 <- merge(SummarybyPatient, Genes, by="Hugo_Symbol", all.x=TRUE)
SummarybyPatient2 <- SummarybyPatient2 %>% filter(!Chromosome %in% c("X","Y"))
SummarybyPatient2$Chromosome <- factor(SummarybyPatient2$Chromosome, levels=c(1:22))
SummarybyPatient2$start_position <- as.numeric(SummarybyPatient2$start_position)
SummarybyPatient2$end_position <- as.numeric(SummarybyPatient2$end_position)
SummarybyPatient2$Position <- rowMeans(SummarybyPatient2[, c("start_position", "end_position")])


library(ggplot2)
library(scales)

cust_theme <- theme_bw() + theme(legend.position="none", 
                                 axis.title = element_blank(), axis.ticks.x = element_blank(), 
                                 axis.text.x = element_blank(),  
                                 strip.background = element_blank(), panel.margin = unit(0, "lines"), 
                                 #panel.border = element_rect(size = 0.25, color = "black"), 
                                 panel.grid = element_blank())


# Final figure
ggplot(Summary2, aes(x=Position, y=Count), fill="#694929", col="#694929") + geom_line(  ) + geom_ribbon( aes(ymin=0, ymax=Count), col="#694929", fill="#694929" ) + facet_grid( . ~ Chromosome, scales="free_x", space="free_x") + cust_theme + scale_y_continuous(breaks=seq(0, 45, 5)) + scale_x_continuous(expand = c(0,0)) + geom_text(data = Summary2[Summary2$Count > 8,], aes(x = Position, y=Count, label=Hugo_Symbol), size=2, vjust=-1 )

#ggplot(Summary2[ Summary2$Chromosome ==14,], aes(x=Position, y=Count), fill="#694929", col="#694929") + geom_line(  ) + geom_ribbon( aes(ymin=0, ymax=Count), col="#694929", fill="#694929" ) + facet_grid( . ~ Chromosome, scales="free_x", space="free_x") + cust_theme + scale_y_continuous(breaks=seq(0, 45, 5)) + scale_x_continuous(expand = c(0,0)) + geom_text(data = Summary2[Summary2$Count > 8,], aes(x = Position, y=Count, label=Hugo_Symbol), size=2, vjust=-1 )


# Final figure for per Patient data
ggplot(SummarybyPatient2, aes(x=Position, y=Proportion), fill="#694929", col="#694929") + geom_line(  ) + geom_ribbon( aes(ymin=0, ymax=Proportion), col="#694929", fill="#694929" ) + facet_grid( . ~ Chromosome, scales="free_x", space="free_x") + cust_theme + scale_y_continuous(breaks=seq(0, 1, 0.1)) + scale_x_continuous(expand = c(0,0)) + geom_hline(yintercept=3.5/43) + geom_text(data = SummarybyPatient2[SummarybyPatient2$Proportion > 0.07,], aes(x = Position, y=Proportion, label=Hugo_Symbol), size=2, vjust=-1 )
