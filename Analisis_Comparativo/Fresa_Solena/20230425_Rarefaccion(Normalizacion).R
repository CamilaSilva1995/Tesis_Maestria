###https://cran.r-project.org/web/packages/Rarefy/vignettes/Rarefy_basics.html
###

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("stringi")
library("dplyr")
library("plyr")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
### Cargado de datos originales 
fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data),1,6)
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1/metadata.csv",header =  FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
fresa_kraken@sam_data$Sample<-row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data)<-c('Treatment','Samples')
samples_to_remove <- c("MP2079","MP2080","MP2088","MP2109","MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)


#### rarefy in phyloseq

#fijamos la semilla
fresa_kraken_rarefy <- rarefy_even_depth(fresa_kraken_fil, rngseed = 123, sample.size = min(sample_sums(fresa_kraken_fil)))

## queremos ver laprofundidad delas muestras
dprof <- data.frame(Samples = colnames(fresa_kraken_rarefy@otu_table@.Data),
                    Reads = sample_sums(fresa_kraken_rarefy),
                    Treatment = fresa_kraken_rarefy@sam_data@.Data[[1]])

mu_Samples <- mean(dprof$Reads) 

ggplot(data = dprof, mapping = aes(x = Samples, y = Reads))+
  geom_bar(stat = "identity", aes( fill = Treatment)) +
  geom_hline(yintercept=mu_Samples, color = "black", linetype="dashed", size=1.5) +
  scale_fill_manual(values = c("cyan3","darkmagenta")) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


## luego dela rarefacción sacar los indices de diversidad alfa


index = estimate_richness(fresa_kraken_rarefy)
index

p<-plot_richness(physeq = fresa_kraken_rarefy, measures = c("Observed","Chao1","Shannon","simpson"),x = "Treatment", color = "Treatment") 
p + geom_point(size=3, alpha=0.5) +
  theme(axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))

# y diversidad beta

percentages_fil <- transform_sample_counts(fresa_kraken_rarefy, function(x) x*100 / sum(x) )
head(percentages_fil@otu_table@.Data)
meta_ord_fil <- ordinate(physeq = percentages_fil, method = "NMDS", distance = "bray") 

plot_ordination(physeq = percentages_fil, ordination = meta_ord_fil, color = "Treatment") +
  geom_text(mapping = aes(label = colnames(fresa_kraken_rarefy@otu_table@.Data)), size = 3, vjust = 1.5)


####################################################################
# media de Shannon y varianza de Chao1 para prueba de hipotesis 

## prueba dehipotesis con medias con indice Shannon

fresa_kraken_rarefy_df <- psmelt(fresa_kraken_rarefy)

## solo tomamos las columnas Sample, OTU y Abundance
fresa_kraken_rarefy_df2 <- fresa_kraken_rarefy_df[c(1,2,3)]
head(fresa_kraken_rarefy_df2)

# queremos pasar de dataframe a table, para calcular el alfa diversidad
rarefy_df <- reshape(fresa_kraken_rarefy_df2, idvar = "Sample",v.names= c("Abundance"), timevar = "OTU", direction = "wide")
rownames(rarefy_df) <- rarefy_df$Sample
rarefy_df <- select(rarefy_df, -Sample)
head(rarefy_df)

Shannon_OTU <- diversity(rarefy_df, "shannon")
Shannon_OTU_df <- data.frame(Shannon=Shannon_OTU)

total <-cbind(fresa_kraken_rarefy@sam_data,Shannon_OTU_df)
#
mu_Shannon <- ddply(total, "Treatment", summarise, grp.mean=mean(Shannon))

p <- ggplot(total, aes(x=Shannon))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
q <- p + geom_vline(data=mu_Shannon, aes(xintercept=grp.mean, color="red"),linetype="dashed")




# Prueba de Hipotesis con varianzas para indice Chao1












