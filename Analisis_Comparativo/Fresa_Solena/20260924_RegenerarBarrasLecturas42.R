## Figura 4.2 - Numero de lecturas por muestra segun tratamiento (Img/cap3/Barras.png)
## Regenera la figura en alta resolucion (300 dpi) con etiquetas en espanol.
##
## La version anterior estaba a 960x540 px y, sobre todo, salia practicamente
## negra: se generaba con plot_bar(fresa_kraken_fil, fill = "Treatment")
## (20230227_Funciones&Graficas.R, linea 120), que dibuja un rectangulo con
## borde por cada OTU. Con miles de OTUs apilados los bordes negros tapaban por
## completo el color del tratamiento y la figura no mostraba lo que afirma su
## pie. Aqui se suman las lecturas por muestra (sample_sums), de modo que cada
## muestra es una sola barra sin bordes internos y el color del tratamiento si
## se distingue. Se agregan ademas las medias por grupo, que son justamente lo
## que el pie de figura compara.

library("phyloseq")
library("ggplot2")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath_analisis <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
outpath_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"

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

## Total de lecturas por muestra
lecturas <- data.frame(
  Muestra = sample_names(fresa_kraken_fil),
  Lecturas = as.numeric(sample_sums(fresa_kraken_fil)),
  Tratamiento = factor(as.character(sample_data(fresa_kraken_fil)$Treatment),
                       levels = c("healthy", "wilted"),
                       labels = c("Saludable", "No saludable"))
)
lecturas$Muestra <- factor(lecturas$Muestra, levels = lecturas$Muestra)

medias <- aggregate(Lecturas ~ Tratamiento, data = lecturas, FUN = mean)
n_grupo <- table(lecturas$Tratamiento)
cat("n por grupo:\n"); print(n_grupo)
cat("media de lecturas por grupo:\n"); print(medias)

figura <- ggplot(lecturas, aes(x = Muestra, y = Lecturas, fill = Tratamiento)) +
  geom_col(width = 0.8) +
  geom_hline(data = medias, aes(yintercept = Lecturas, colour = Tratamiento),
             linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
  scale_y_continuous(labels = function(x) x / 1e6) +
  labs(x = "Muestra", y = "Número de lecturas (millones)", fill = "Tratamiento") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 7, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 9))

ggsave("Barras.png", plot = figura, path = outpath_analisis,
       width = 30, height = 15, dpi = 300, units = "cm", bg = "white")
ggsave("Barras.png", plot = figura, path = outpath_tesis,
       width = 30, height = 15, dpi = 300, units = "cm", bg = "white")

cat("Listo 4.2\n")
