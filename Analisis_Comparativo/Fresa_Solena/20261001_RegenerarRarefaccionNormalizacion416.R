## Figura 4.16 - Rarefaccion y normalizacion: curvas de rarefaccion, diversidad alfa tras rarificar
##               a profundidad comun y sensibilidad del PERMANOVA a la normalizacion TMM
##               (Img/cap3/Rarefaccion_Normalizacion.png)
##
## Pregunta: la tendencia hacia mayor diversidad en las plantas saludables se debe a que tienen
## mas lecturas (esfuerzo de muestreo) o se mantiene cuando todas las muestras se llevan a la
## misma profundidad?
##
## Pasos:
##   1. Curvas de rarefaccion (vegan::rarecurve): numero esperado de taxones en funcion de la
##      profundidad. Si las curvas se aplanan, el muestreo esta saturado.
##   2. Rarefaccion a profundidad comun (phyloseq::rarefy_even_depth) con la profundidad minima
##      del conjunto filtrado, semilla 2026. Diversidad alfa y pruebas de hipotesis sobre los
##      datos rarificados.
##   3. Normalizacion TMM (edgeR::calcNormFactors, Robinson & Oshlack 2010). Como Shannon y Simpson
##      dependen solo de las proporciones dentro de cada muestra, un factor de escala por muestra
##      no los cambia; TMM se evalua donde si importa: la disimilitud entre muestras (PERMANOVA).
##
## Semilla unica de la tesis para todo procedimiento aleatorio: 2026

suppressMessages({library("phyloseq"); library("vegan"); library("ggplot2"); library("cowplot"); library("reshape2")})
setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
out_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"
out_res   <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
colores <- c("Saludable" = "#F8766D", "No saludable" = "#00BFC4")

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
X <- t(as(otu_table(fil), "matrix"))            # muestras x taxones
stopifnot(identical(rownames(X), names(grupo)))  # emparejamiento explicito muestra-grupo

## ---- 1. Curvas de rarefaccion ----
profundidad <- rowSums(X)
cat("Profundidad (lecturas clasificadas): min =", min(profundidad), " max =", max(profundidad), "\n")
paso <- 200000
curvas <- rarecurve(X, step = paso, label = FALSE, tidy = TRUE)   # columnas: Site, Sample (profundidad), Species
curvas$Grupo <- grupo[as.character(curvas$Site)]
## saturacion: fraccion de la riqueza final alcanzada a la profundidad minima
sat <- sapply(rownames(X), function(s) rarefy(X[s, , drop = FALSE], sample = min(profundidad)) / sum(X[s, ] > 0))
cat("Fraccion de la riqueza observada alcanzada a la profundidad minima: media =", round(mean(sat), 3),
    " min =", round(min(sat), 3), "\n")

## ---- 2. Rarefaccion a profundidad comun ----
set.seed(2026)
rar <- rarefy_even_depth(fil, sample.size = min(profundidad), rngseed = 2026, replace = FALSE, verbose = FALSE)
cat("Taxones eliminados por quedar sin lecturas tras rarificar:", ntaxa(fil) - ntaxa(rar), "\n")
alfa_orig <- estimate_richness(fil, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
alfa_rar  <- estimate_richness(rar, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
stopifnot(identical(rownames(alfa_orig), rownames(alfa_rar)), identical(rownames(alfa_rar), names(grupo)))
alfa_rar$Grupo <- grupo; alfa_orig$Grupo <- grupo

medias <- function(d) { a <- aggregate(. ~ Grupo, d[, c("Observed","Chao1","Shannon","Simpson","Grupo")], mean); a[,-1] <- round(a[,-1], 3); a }
cat("\n== Medias por grupo, datos originales ==\n"); print(medias(alfa_orig))
cat("== Medias por grupo, datos rarificados a", min(profundidad), "lecturas ==\n"); print(medias(alfa_rar))

prueba <- function(v, g) {
  x <- v[g == "Saludable"]; y <- v[g == "No saludable"]
  Ft <- var.test(x, y); tw <- t.test(x, y); mw <- wilcox.test(x, y, exact = TRUE)
  c(media_sal = mean(x), media_nosal = mean(y), F = unname(Ft$statistic), p_F = Ft$p.value,
    t_Welch = unname(tw$statistic), gl_Welch = unname(tw$parameter), p_Welch = tw$p.value,
    W = unname(mw$statistic), p_MW = mw$p.value)
}
tabla <- rbind(
  "Observados (original)"   = prueba(alfa_orig$Observed, grupo),
  "Observados (rarificado)" = prueba(alfa_rar$Observed, grupo),
  "Chao1 (original)"        = prueba(alfa_orig$Chao1, grupo),
  "Chao1 (rarificado)"      = prueba(alfa_rar$Chao1, grupo),
  "Shannon (original)"      = prueba(alfa_orig$Shannon, grupo),
  "Shannon (rarificado)"    = prueba(alfa_rar$Shannon, grupo),
  "Simpson (original)"      = prueba(alfa_orig$Simpson, grupo),
  "Simpson (rarificado)"    = prueba(alfa_rar$Simpson, grupo))
cat("\n== Pruebas sobre los indices, originales y rarificados ==\n"); print(round(tabla, 4))
write.csv(round(tabla, 5), file.path(out_res, "Rarefaccion_pruebas.csv"))

## ---- 3. Normalizacion TMM y sensibilidad del PERMANOVA ----
tmm_manual <- function(counts, logratioTrim = 0.3, sumTrim = 0.05, Acutoff = -1e10) {
  ## Implementacion directa de Robinson & Oshlack (2010): se usa solo si edgeR no esta instalado.
  lib <- colSums(counts); f75 <- apply(counts, 2, function(x) quantile(x / sum(x), 0.75))
  ref <- which.min(abs(f75 - mean(f75)))
  factores <- sapply(seq_len(ncol(counts)), function(i) {
    obs <- counts[, i] / lib[i]; refp <- counts[, ref] / lib[ref]
    ok <- obs > 0 & refp > 0
    M <- log2(obs[ok] / refp[ok]); A <- 0.5 * log2(obs[ok] * refp[ok])
    v <- (lib[i] - counts[ok, i]) / lib[i] / counts[ok, i] + (lib[ref] - counts[ok, ref]) / lib[ref] / counts[ok, ref]
    n <- length(M); loM <- floor(n * logratioTrim) + 1; hiM <- n + 1 - loM
    loA <- floor(n * sumTrim) + 1; hiA <- n + 1 - loA
    keep <- (rank(M) >= loM & rank(M) <= hiM) & (rank(A) >= loA & rank(A) <= hiA) & A > Acutoff
    2^(sum(M[keep] / v[keep]) / sum(1 / v[keep]))
  })
  factores / exp(mean(log(factores)))
}
counts <- t(X)                                   # taxones x muestras (formato edgeR)
if (requireNamespace("edgeR", quietly = TRUE)) {
  y <- edgeR::DGEList(counts = counts, remove.zeros = TRUE)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  nf <- y$samples$norm.factors; metodo_tmm <- "edgeR::calcNormFactors(method = 'TMM')"
} else {
  nf <- tmm_manual(counts); metodo_tmm <- "implementacion propia de TMM (edgeR no instalado)"
}
names(nf) <- colnames(counts)
cat("\nTMM calculado con:", metodo_tmm, "\n")
cat("Factores de normalizacion TMM: min =", round(min(nf), 3), " max =", round(max(nf), 3), "\n")
cat("Factor TMM medio por grupo:\n"); print(round(tapply(nf, grupo, mean), 3))
## conteos normalizados: dividir cada muestra por su tamano de biblioteca efectivo (lib * factor)
tmm <- sweep(counts, 2, colSums(counts) * nf, "/")        # ahora cada columna suma 1/nf
tss <- sweep(counts, 2, colSums(counts), "/")             # abundancias relativas (TSS)

set.seed(2026); perm_tss <- adonis2(vegdist(t(tss), "bray") ~ grupo, permutations = 999)
set.seed(2026); perm_tmm <- adonis2(vegdist(t(tmm), "bray") ~ grupo, permutations = 999)
set.seed(2026); perm_rar <- adonis2(vegdist(t(as(otu_table(rar), "matrix")), "bray") ~ grupo, permutations = 999)
sens <- data.frame(normalizacion = c("Abundancia relativa (TSS)", "TMM", "Rarefaccion a profundidad comun"),
                   pseudo_F = c(perm_tss$F[1], perm_tmm$F[1], perm_rar$F[1]),
                   R2 = c(perm_tss$R2[1], perm_tmm$R2[1], perm_rar$R2[1]),
                   p = c(perm_tss$`Pr(>F)`[1], perm_tmm$`Pr(>F)`[1], perm_rar$`Pr(>F)`[1]))
cat("\n== PERMANOVA (Bray-Curtis) segun la normalizacion ==\n"); print(sens, digits = 3)
write.csv(sens, file.path(out_res, "Normalizacion_PERMANOVA.csv"), row.names = FALSE)

## ---- Figura ----
pA <- ggplot(curvas, aes(Sample / 1e6, Species, group = Site, colour = Grupo)) +
  geom_line(alpha = 0.8) +
  geom_vline(xintercept = min(profundidad) / 1e6, linetype = "dashed") +
  scale_colour_manual(values = colores) +
  labs(x = "Lecturas submuestreadas (millones)", y = "Taxones esperados", colour = "Tratamiento",
       title = "Curvas de rarefacción") +
  theme_bw() + theme(legend.position = "bottom", text = element_text(size = 12))

larga <- melt(alfa_rar, id.vars = "Grupo", measure.vars = c("Observed", "Chao1", "Shannon", "Simpson"),
              variable.name = "Indice", value.name = "Valor")
levels(larga$Indice) <- c("Observado", "Chao1", "Shannon", "Simpson")
pB <- ggplot(larga, aes(Grupo, Valor, colour = Grupo)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.12, size = 1.8, alpha = 0.8) +
  facet_wrap(~ Indice, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = colores, guide = "none") +
  labs(x = NULL, y = "Medida de diversidad alfa",
       title = sprintf("Diversidad alfa tras rarificar a %.2f millones de lecturas", min(profundidad) / 1e6)) +
  theme_bw() + theme(text = element_text(size = 12), axis.text.x = element_text(angle = 30, hjust = 1))

fig <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1, 1.6), labels = c("A)", "B)"), label_size = 13)
ggsave("Rarefaccion_Normalizacion.png", plot = fig, path = out_tesis, width = 32, height = 13, dpi = 300, units = "cm", bg = "white")
ggsave("Rarefaccion_Normalizacion.png", plot = fig, path = out_res,   width = 32, height = 13, dpi = 300, units = "cm", bg = "white")
cat("\nFigura guardada en", out_tesis, "\n")
