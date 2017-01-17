# Needed to get mutation information for certain genes in our dataset to compare to CN

library(data.table)
AllMutations <- fread("~/Documents/DataFromOthers/Zarko/2016-06-01_Results_of_all_available_samples.mutations.txt")
setDT(AllMutations)

SimplifiedMutations <- AllMutations[, AllMuts := lapply(HGVSc, paste0, collapse=", "), by = c("Tumor_Sample_Barcode", "Hugo_Symbol")]

SimpMutations <- aggregate(AllMutations$HGVSc, by=list(AllMutations$Tumor_Sample_Barcode, AllMutations$Hugo_Symbol), paste)
setDT(SimpMutations)
colnames(SimpMutations) <- c("Tumor_Sample_Barcode","Hugo_Symbol", "Mutations")

MutsNeeded <- SimpMutations[Hugo_Symbol %in% c("SMAD4","SYCP2","APC","TP53","KRAS")]
MutsNeeded <- MutsNeeded[!Tumor_Sample_Barcode %in% c("UICE_0001","UICE_0004","UICE_0005","UICE_0006","UICE_0007","UICE_0008","UICE_0012","UICE_0013","UICE_0015","UICE_0018","UICE_0025","UICE_0030","UICE_0031","UICE_0032","UICE_0037","UICE_0039","UICE_0040","UICE_0043","UICE_0044","UICE_0046","UICE_0047","UICE_0049","UICE_0050","UICE_0051","UICE_0055","UICE_0056")]
