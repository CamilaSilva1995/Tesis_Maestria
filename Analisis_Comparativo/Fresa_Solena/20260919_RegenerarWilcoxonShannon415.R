## Figura 4.15 - Prueba de Mann-Whitney sobre el indice de Shannon (Img/cap3/Wilcoxon_Shannon.png)
## Regenera la figura en alta resolucion (300 dpi) con etiquetas en espanol.
## La version anterior estaba a 960x540 px y sin informacion en el eje Y.

library("phyloseq")
library("ggplot2")
library("vegan")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"

fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data, 4, 100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data), 1, 6)
metadata_fresa <- read.csv2("metadata.csv", header = FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
fresa_kraken@sam_data$Sample <- row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data) <- c('Treatment', 'Samples')

samples_to_remove <- c("MP2079", "MP2080", "MP2088", "MP2109", "MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)

## Indice de Shannon por muestra y su rango (igual que 20230510_PruebaWilcoxon.R, lineas 45-50)
OTU <- t(fresa_kraken_fil@otu_table@.Data)
SAM <- fresa_kraken_fil@sam_data
Shannon_OTU <- diversity(OTU, "shannon")
Shannon_OTU_df <- data.frame(sample = names(Shannon_OTU), value = Shannon_OTU,
                             measure = rep("Shannon", length(Shannon_OTU)))
total <- cbind(Shannon_OTU_df, SAM)
total$Rank <- rank(total$value)

total$Treatment <- factor(total$Treatment, levels = c("healthy", "wilted"),
                          labels = c("Saludable", "No saludable"))

figura <- ggplot(data = total, aes(x = Rank, y = value)) +
  geom_point(aes(colour = Treatment), size = 3) +
  ## El script original usaba ylab("") y axis.text.y = element_blank(), lo que dejaba
  ## el eje Y sin informacion; aqui se rotula con la magnitud que realmente representa.
  labs(x = "Rango", y = "Índice de Shannon", colour = "Tratamiento") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        text = element_text(size = 12),
        axis.text = element_text(size = 10))

ggsave("Wilcoxon_Shannon.png", plot = figura, path = outpath,
       width = 30, height = 15, dpi = 300, units = "cm", bg = "white")

cat("Listo 4.14\n")

ggsave("Wilcoxon_Shannon.png", plot = figura, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img",
       width = 30, height = 15, dpi = 300, units = "cm", bg = "white")
