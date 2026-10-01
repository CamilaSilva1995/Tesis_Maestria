## Figura 4.4 - PERMANOVA y PERMDISP sobre la composicion del microbioma rizosferico de fresa (Img/cap3/PERMANOVA_PERMDISP.png)
## (Analisis de varianza multivariado por permutaciones; Anderson 2001, 2006).
##
## Pregunta: la composicion de la comunidad difiere entre plantas saludables y no saludables?
##   - PERMANOVA (adonis2): contrasta la igualdad de centroides entre grupos.
##   - PERMDISP (betadisper): contrasta la igualdad de dispersion alrededor del centroide.
##     Es el diagnostico obligado del PERMANOVA, porque este tambien reacciona a diferencias
##     de dispersion, sobre todo con grupos desbalanceados (35 frente a 18).
##
## Entrada : Data/fresa_solena/Data1/fresa_kraken.biom + metadata.csv (mismo filtro que la tesis)
## Salida  : Results_img/PERMANOVA_PERMDISP.png, latex/Img/cap3/PERMANOVA_PERMDISP.png
##           Results_img/PERMANOVA_PERMDISP_resultados.csv

suppressMessages({library("phyloseq"); library("vegan"); library("ggplot2"); library("cowplot")})

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
out_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"
out_res   <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"

## ---- Carga y filtro de calidad (identico al resto de scripts de la tesis) ----
fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data, 4, 100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data), 1, 6)
metadata_fresa <- read.csv2("metadata.csv", header = FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
colnames(fresa_kraken@sam_data) <- "Treatment"
samples_to_remove <- c("MP2079", "MP2080", "MP2088", "MP2109", "MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)

## Abundancias relativas y matriz de disimilitud de Bray-Curtis (53 x 53)
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x / sum(x))
D_bray <- phyloseq::distance(percentages_fil, method = "bray")
D_jacc <- phyloseq::distance(percentages_fil, method = "jaccard")   # version cuantitativa (vegan)
meta <- data.frame(sample_data(fresa_kraken_fil))
meta$Treatment <- factor(meta$Treatment, levels = c("healthy", "wilted"),
                         labels = c("Saludable", "No saludable"))
## Comprobacion explicita del emparejamiento muestra-tratamiento
stopifnot(identical(rownames(meta), attr(D_bray, "Labels")))

## ---- PERMANOVA ----
## Semilla unica de la tesis para todo procedimiento aleatorio: 2026
set.seed(2026)
perm_bray <- adonis2(D_bray ~ Treatment, data = meta, permutations = 999)
set.seed(2026)
perm_jacc <- adonis2(D_jacc ~ Treatment, data = meta, permutations = 999)
cat("== PERMANOVA Bray-Curtis ==\n"); print(perm_bray)
cat("== PERMANOVA Jaccard ==\n");     print(perm_jacc)

## ---- PERMDISP ----
disp_bray <- betadisper(D_bray, meta$Treatment, type = "centroid")
set.seed(2026)
disp_test <- permutest(disp_bray, permutations = 999)
cat("== PERMDISP Bray-Curtis ==\n"); print(disp_test)
cat("Distancia media al centroide por grupo:\n"); print(round(tapply(disp_bray$distances, meta$Treatment, mean), 4))

## ---- Tabla de resultados ----
res <- data.frame(
  prueba       = c("PERMANOVA", "PERMANOVA", "PERMDISP"),
  disimilitud  = c("Bray-Curtis", "Jaccard cuantitativo", "Bray-Curtis"),
  estadistico  = c(perm_bray$F[1], perm_jacc$F[1], disp_test$tab$F[1]),
  R2           = c(perm_bray$R2[1], perm_jacc$R2[1], NA),
  gl           = c("1, 51", "1, 51", "1, 51"),
  permutaciones = 999,
  p            = c(perm_bray$`Pr(>F)`[1], perm_jacc$`Pr(>F)`[1], disp_test$tab$`Pr(>F)`[1]))
res$estadistico <- round(res$estadistico, 3); res$R2 <- round(res$R2, 4)
print(res)
write.csv(res, file.path(out_res, "PERMANOVA_PERMDISP_resultados.csv"), row.names = FALSE)

## ---- Figura: PCoA con centroides (A) y distancia al centroide por grupo (B) ----
colores <- c("Saludable" = "#F8766D", "No saludable" = "#00BFC4")
pc <- data.frame(disp_bray$vectors[, 1:2], Grupo = meta$Treatment)
colnames(pc)[1:2] <- c("PCoA1", "PCoA2")
cen <- data.frame(disp_bray$centroids[, 1:2], Grupo = rownames(disp_bray$centroids))
colnames(cen)[1:2] <- c("PCoA1", "PCoA2")
seg <- merge(pc, cen, by = "Grupo", suffixes = c("", "_c"))
var_exp <- round(100 * disp_bray$eig[1:2] / sum(disp_bray$eig[disp_bray$eig > 0]), 1)

pA <- ggplot() +
  geom_segment(data = seg, aes(x = PCoA1_c, y = PCoA2_c, xend = PCoA1, yend = PCoA2, colour = Grupo), alpha = 0.4) +
  geom_point(data = pc, aes(PCoA1, PCoA2, colour = Grupo), size = 2.5) +
  geom_point(data = cen, aes(PCoA1, PCoA2, fill = Grupo), shape = 23, size = 5, colour = "black") +
  scale_colour_manual(values = colores) + scale_fill_manual(values = colores, guide = "none") +
  labs(x = paste0("PCoA 1 (", var_exp[1], " %)"), y = paste0("PCoA 2 (", var_exp[2], " %)"),
       colour = "Tratamiento",
       title = sprintf("PERMANOVA: pseudo-F = %.2f, R² = %.3f, p = %.3f",
                       perm_bray$F[1], perm_bray$R2[1], perm_bray$`Pr(>F)`[1])) +
  theme_bw() + theme(legend.position = "bottom", text = element_text(size = 12))

dd <- data.frame(dist = disp_bray$distances, Grupo = meta$Treatment)
pB <- ggplot(dd, aes(Grupo, dist, colour = Grupo)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = colores, guide = "none") +
  labs(x = NULL, y = "Distancia al centroide del grupo",
       title = sprintf("PERMDISP: F = %.2f, p = %.3f", disp_test$tab$F[1], disp_test$tab$`Pr(>F)`[1])) +
  theme_bw() + theme(text = element_text(size = 12))

fig <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1.35, 1), labels = c("A)", "B)"), label_size = 13)
ggsave("PERMANOVA_PERMDISP.png", plot = fig, path = out_tesis, width = 30, height = 14, dpi = 300, units = "cm", bg = "white")
ggsave("PERMANOVA_PERMDISP.png", plot = fig, path = out_res,   width = 30, height = 14, dpi = 300, units = "cm", bg = "white")
cat("Figura guardada en", out_tesis, "y", out_res, "\n")
