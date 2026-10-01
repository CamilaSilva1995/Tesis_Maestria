## Figura 4.17 - Proporcion de Fusarium y de los generos candidatos por grupo
##               (Img/cap3/Proporciones_Fusarium.png)
##
## Para cada muestra se calcula el porcentaje de lecturas clasificadas asignadas a un genero:
##   (a) respecto de todas las lecturas clasificadas de la muestra (denominador: comunidad completa)
##   (b) respecto de las lecturas asignadas a eucariotas (denominador: fraccion eucariota)
## y se compara entre grupos con la prueba de Mann-Whitney (exacta, bilateral).
## Generos: Fusarium (eje del estudio) y los candidatos identificados por inspeccion de barras:
##   benéficos Pseudomonas, Bacillus, Streptomyces, Paenibacillus; patogenos Ralstonia, Phytophthora.
## Semilla unica de la tesis: 2026 (solo afecta al jitter de la figura).

suppressMessages({library("phyloseq"); library("ggplot2"); library("cowplot")})
setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
out_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"
out_res   <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
colores <- c("Saludable" = "#F8766D", "No saludable" = "#00BFC4")
set.seed(2026)

## ---- Carga y filtro de calidad (identico al resto de scripts) ----
fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data, 4, 100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data), 1, 6)
metadata_fresa <- read.csv2("metadata.csv", header = FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
colnames(fresa_kraken@sam_data) <- "Treatment"
samples_to_remove <- c("MP2079", "MP2080", "MP2088", "MP2109", "MP2137")
fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)
grupo <- factor(data.frame(sample_data(fil))$Treatment, levels = c("healthy", "wilted"),
                labels = c("Saludable", "No saludable"))
names(grupo) <- sample_names(fil)

X   <- as(otu_table(fil), "matrix")                   # taxones x muestras
tax <- as(tax_table(fil), "matrix")
stopifnot(identical(colnames(X), names(grupo)))      # emparejamiento explicito muestra-grupo
total   <- colSums(X)
total_e <- colSums(X[tax[, "Kingdom"] == "Eukaryota", , drop = FALSE])

pct <- function(genero, denominador = total) {
  100 * colSums(X[tax[, "Genus"] == genero, , drop = FALSE]) / denominador
}
prueba <- function(v) {
  x <- v[grupo == "Saludable"]; y <- v[grupo == "No saludable"]
  mw <- suppressWarnings(wilcox.test(x, y, exact = TRUE))
  c(media_sal = mean(x), mediana_sal = median(x), media_nosal = mean(y), mediana_nosal = median(y),
    W = unname(mw$statistic), p_MW = mw$p.value)
}

## ---- Fusarium con los dos denominadores ----
fus_total <- pct("Fusarium"); fus_euk <- pct("Fusarium", total_e)
cat("Lecturas de Fusarium por muestra: min =", min(colSums(X[tax[, "Genus"] == "Fusarium", , drop = FALSE])),
    " max =", max(colSums(X[tax[, "Genus"] == "Fusarium", , drop = FALSE])), "\n")
cat("Fraccion eucariota de las lecturas clasificadas (%): media =", round(mean(100 * total_e / total), 3), "\n")

## ---- Generos candidatos ----
candidatos <- c("Fusarium", "Phytophthora", "Ralstonia", "Pseudomonas", "Bacillus", "Streptomyces", "Paenibacillus")
tabla <- t(sapply(candidatos, function(g) prueba(pct(g))))
tabla <- rbind(tabla, "Fusarium (fraccion eucariota)" = prueba(fus_euk))
tabla <- cbind(tabla, p_BH = p.adjust(tabla[, "p_MW"], method = "BH"))
cat("\n== Porcentaje de lecturas por genero y grupo, prueba de Mann-Whitney (p_BH: Benjamini-Hochberg) ==\n")
print(round(tabla, 4))
write.csv(round(tabla, 5), file.path(out_res, "Proporciones_generos_candidatos.csv"))

## ---- Figura ----
df_fus <- rbind(data.frame(Grupo = grupo, Valor = fus_total, Denominador = "Comunidad completa"),
                data.frame(Grupo = grupo, Valor = fus_euk,   Denominador = "Fracción eucariota"))
etq <- sprintf("p = %.3f", c(tabla["Fusarium", "p_MW"], tabla["Fusarium (fraccion eucariota)", "p_MW"]))
pA <- ggplot(subset(df_fus, Denominador == "Comunidad completa"), aes(Grupo, Valor, colour = Grupo)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.12, size = 2.3, alpha = 0.8) +
  scale_colour_manual(values = colores, guide = "none") +
  labs(x = NULL, y = "% de las lecturas clasificadas", title = "Fusarium, comunidad completa", subtitle = paste("Mann-Whitney", etq[1])) +
  theme_bw() + theme(text = element_text(size = 12))
pB <- ggplot(subset(df_fus, Denominador == "Fracción eucariota"), aes(Grupo, Valor, colour = Grupo)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.12, size = 2.3, alpha = 0.8) +
  scale_colour_manual(values = colores, guide = "none") +
  labs(x = NULL, y = "% de las lecturas de eucariotas", title = "Fusarium, fracción eucariota", subtitle = paste("Mann-Whitney", etq[2])) +
  theme_bw() + theme(text = element_text(size = 12))

otros <- setdiff(candidatos, "Fusarium")
df_c <- do.call(rbind, lapply(otros, function(g) data.frame(Grupo = grupo, Valor = pct(g), Genero = g)))
df_c$Genero <- factor(df_c$Genero, levels = otros,
                      labels = sprintf("%s (p = %.3f)", otros, tabla[otros, "p_MW"]))
pC <- ggplot(df_c, aes(Grupo, Valor, colour = Grupo)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.12, size = 1.6, alpha = 0.8) +
  facet_wrap(~ Genero, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = colores) +
  labs(x = NULL, y = "% de las lecturas clasificadas", colour = "Tratamiento", title = "Géneros candidatos") +
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom",
                     axis.text.x = element_blank(), axis.ticks.x = element_blank())

arriba <- plot_grid(pA, pB, ncol = 2, labels = c("A)", "B)"), label_size = 13)
fig <- plot_grid(arriba, pC, ncol = 1, rel_heights = c(1, 1.05), labels = c("", "C)"), label_size = 13)
ggsave("Proporciones_Fusarium.png", plot = fig, path = out_tesis, width = 30, height = 24, dpi = 300, units = "cm", bg = "white")
ggsave("Proporciones_Fusarium.png", plot = fig, path = out_res,   width = 30, height = 24, dpi = 300, units = "cm", bg = "white")
cat("\nFigura guardada en", out_tesis, "\n")
