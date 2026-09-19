## Figura 4.13 - Prueba t de Student sobre el indice de Shannon en Fusarium
## (Img/cap3/tTest_Shannon_Fusarium.png). Regenera en alta resolucion (300 dpi),
## en espanol. La version anterior estaba a 960x512 px y era un montaje manual de
## tres imagenes sueltas de Results_img.

library("phyloseq")
library("ggplot2")
library("vegan")
library("dplyr")
library("plyr")
library("gginference")
library("cowplot")

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

## Subconjunto de Eukaryota aglomerado a nivel de Genero
## (igual que 20230410_PruebasdeHipotesisMedias.R, lineas 282-297)
merge_Eukaryota <- subset_taxa(fresa_kraken_fil, Kingdom == "Eukaryota")
glom <- tax_glom(merge_Eukaryota, taxrank = "Genus")
glom_df <- psmelt(glom)

glom_df2 <- glom_df[c(1, 2, 3)]
df <- reshape(glom_df2, idvar = "Sample", v.names = c("Abundance"),
              timevar = "OTU", direction = "wide")
rownames(df) <- df$Sample
df <- dplyr::select(df, -Sample)

Shannon_OTU <- diversity(df, "shannon")
Shannon_OTU_df <- data.frame(sample = names(Shannon_OTU), value = Shannon_OTU,
                             measure = rep("Shannon", length(Shannon_OTU)))
total_Shannon <- cbind(glom@sam_data, Shannon_OTU_df)

## Tratamiento en espanol
total_Shannon$Treatment <- factor(total_Shannon$Treatment,
                                  levels = c("healthy", "wilted"),
                                  labels = c("Saludable", "No Saludable"))
mu_Shannon <- ddply(total_Shannon, "Treatment", summarise, grp.mean = mean(value))

## ---- Panel A: histograma con la media por grupo ----
panelA <- ggplot(total_Shannon, aes(x = value)) +
  geom_histogram(color = "pink", fill = "black", bins = 30) +
  facet_grid(Treatment ~ .) +
  geom_vline(data = mu_Shannon, aes(xintercept = grp.mean, color = "Media"),
             linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("Media" = "red")) +
  labs(x = "Índice de diversidad de Shannon", y = "Número de muestras", color = NULL) +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        strip.text = element_text(size = 10))

## ---- Paneles B y C: prueba t con varianzas iguales y diferentes ----
total_H <- total_Shannon[total_Shannon$Treatment == "Saludable", ]
total_W <- total_Shannon[total_Shannon$Treatment == "No Saludable", ]

pruebat  <- t.test(total_H$value, total_W$value, var.equal = TRUE,  alternative = "two.sided")
pruebat2 <- t.test(total_H$value, total_W$value, var.equal = FALSE, alternative = "two.sided")

gl1 <- round(as.numeric(pruebat$parameter))
gl2 <- round(as.numeric(pruebat2$parameter))
cat("grados de libertad:", gl1, "y", gl2, "\n")

etiqueta_t <- function(p, gl, titulo) {
  g <- ggttest(p)
  ## ggttest rotula en ingles. El texto del estadistico vive en el mapping de la
  ## capa de texto (GeomText); se reemplaza por una etiqueta fija en espanol.
  for (i in seq_along(g$layers)) {
    if (inherits(g$layers[[i]]$geom, "GeomText") &&
        !is.null(g$layers[[i]]$mapping$label)) {
      g$layers[[i]]$mapping$label <- NULL
      g$layers[[i]]$aes_params$label <- paste("Estadístico =",
                                              round(as.numeric(p$statistic), 4))
    }
  }
  g +
    labs(title = titulo,
         subtitle = "Hipótesis alternativa: bilateral",
         caption = "Alfa = 0.05",
         x = paste("Distribución t con", gl, "grados de libertad")) +
    theme(plot.title = element_text(hjust = 0, size = 12),
          plot.subtitle = element_text(hjust = 0, size = 10),
          text = element_text(size = 11))
}

panelB <- etiqueta_t(pruebat,  gl1, "Asumiendo varianzas iguales")
panelC <- etiqueta_t(pruebat2, gl2, "Asumiendo varianzas diferentes")

derecha <- plot_grid(panelB, panelC, ncol = 1, labels = c("B)", "C)"),
                     label_size = 13, hjust = 0, label_x = 0.01)
figura <- plot_grid(panelA, derecha, ncol = 2, labels = c("A)", ""),
                    label_size = 13, hjust = 0, label_x = 0.01,
                    rel_widths = c(1, 1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("tTest_Shannon_Fusarium.png", plot = figura, path = outpath,
       width = 30, height = 16, dpi = 300, units = "cm", bg = "white")

cat("Listo 4.13\n")

ggsave("tTest_Shannon_Fusarium.png", plot = figura, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img",
       width = 30, height = 16, dpi = 300, units = "cm", bg = "white")
