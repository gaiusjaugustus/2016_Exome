############### Get genomic locations for genes that are mutated

# Import list of mutations

library(data.table)
mutations <- fread("C:/Users/gaugustus/Documents/DataFromOthers/Zarko/2016-06-01_Results_of_all_available_samples.mutations.txt")

# Query biomaRt for their start and end positions
library(dplyr)
genes <- mutations %>% select(Hugo_Symbol) %>% arrange(Hugo_Symbol) %>% distinct(Hugo_Symbol)

library(biomaRt)
ensembl <-  useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes2 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "external_gene_name"), filters = c("hgnc_symbol"), mart = ensembl, values = list(genes) )
colnames(genes2)[5] <- "Hugo_Symbol"


# There were some duplicates in biomart.  I don't know which to choose so I just chose at random, since duplicates are generally close to each other.

genes3 <- genes2 %>% arrange(Hugo_Symbol)
genes3 <- genes3[ !duplicated(genes3$Hugo_Symbol),]

# Merged dataframes by Hugo_Symbol (gene name)

matching <- match(mutations$Hugo_Symbol, genes3$Hugo_Symbol)

matching <- merge(mutations, genes3, by="Hugo_Symbol", all.x=TRUE)

# Exported to a new file for reuse.
write.table(matching, "C:/Users/gaugustus/Documents/DataFromOthers/Zarko/2016-09-26_Results_of_all_available_samples.mutations.geneannotated.txt", quote=FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

