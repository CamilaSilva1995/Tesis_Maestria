#ultimos datos
library("phyloseq") 
setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data3")
crop_vs_native_rawCounts_S <- import_biom("crop_vs_native_rawCounts_S.biom")
class(crop_vs_native_rawCounts_S) # objeto phyloseq

metadata_crop_vs_native_rawCounts_S <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data3/crop_vs_native_rawCounts_S_metadata.csv",header =  FALSE, row.names = 1, sep = ",")



