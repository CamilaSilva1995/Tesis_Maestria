## Figura 3.6 - Tercer conjunto complementario: cultivo frente a suelo nativo (Img/cap2/Complementarios_nativo.png)
## Regenera crop_vs_native_rawCounts_S_Alfa/Beta_Diversidad.png de 20230502_NuevosDatos.R con buen diseno,
## y agrega Mann-Whitney por indice y PERMANOVA (4 frente a 4 muestras).

source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_cargar_complementarios.R")

## ---- Panel A: diversidad alfa ----
alfa <- estimate_richness(nat, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
## estimate_richness convierte los guiones de los nombres en puntos (make.names); se restauran los nombres originales
stopifnot(identical(rownames(alfa), make.names(names(cat2)))); rownames(alfa) <- names(cat2); alfa$Categoria <- cat2
larga <- melt(alfa, id.vars = "Categoria", measure.vars = c("Observed", "Chao1", "Shannon", "Simpson"), variable.name = "Indice", value.name = "Valor")
levels(larga$Indice) <- c("Observado", "Chao1", "Shannon", "Simpson")
mw <- sapply(c("Observed", "Chao1", "Shannon", "Simpson"), function(v) wilcox.test(alfa[[v]] ~ alfa$Categoria, exact = TRUE)$p.value)
cat("Medias por categoria:\n"); m <- aggregate(. ~ Categoria, alfa[, c("Observed","Chao1","Shannon","Simpson","Categoria")], mean); m[,-1] <- round(m[,-1], 3); print(m)
cat("Mann-Whitney exacto (4 frente a 4; el valor p minimo posible es 0.029):\n"); print(round(mw, 3))
cat("Lecturas clasificadas (millones):\n"); print(round(tapply(sample_sums(nat) / 1e6, cat2, mean), 1))
etq <- data.frame(Indice = factor(c("Observado", "Chao1", "Shannon", "Simpson"), levels = levels(larga$Indice)), texto = sprintf("p = %.3f", mw))
pA <- ggplot(larga, aes(Categoria, Valor, colour = Categoria)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.1, size = 2.4, alpha = 0.85) +
  facet_wrap(~ Indice, scales = "free_y", nrow = 1) +
  geom_text(data = etq, aes(x = 1.5, y = Inf, label = texto), inherit.aes = FALSE, vjust = 1.4, size = 3.3) +
  scale_colour_manual(values = colores2, guide = "none") +
  labs(x = NULL, y = "Medida de diversidad alfa", title = "A) Diversidad alfa, cultivo frente a suelo nativo") +
  theme_bw() + theme(text = element_text(size = 12), plot.title = element_text(face = "bold"))

## ---- Panel B: NMDS Bray-Curtis + PERMANOVA ----
D <- phyloseq::distance(rel_nat, method = "bray")
## Con 8 muestras el NMDS no es confiable (estres cercano a cero); se usa PCoA, que es determinista
ord <- ordinate(rel_nat, method = "PCoA", distance = D)
sc <- data.frame(ord$vectors[, 1:2]); colnames(sc) <- c("PCoA1", "PCoA2"); sc$Categoria <- cat2[rownames(sc)]
var_exp <- round(100 * ord$values$Relative_eig[1:2], 1)
meta <- data.frame(sample_data(nat)); set.seed(2026); perm <- adonis2(D ~ Categoria, data = meta, permutations = 999)
cat("== PERMANOVA cultivo frente a nativo ==\n"); print(perm)
pB <- ggplot(sc, aes(PCoA1, PCoA2, colour = Categoria)) + geom_point(size = 4, alpha = 0.9) +
  scale_colour_manual(values = colores2) +
  labs(x = sprintf("PCoA 1 (%.1f %%)", var_exp[1]), y = sprintf("PCoA 2 (%.1f %%)", var_exp[2]),
       title = "B) PCoA, Bray-Curtis",
       subtitle = sprintf("PERMANOVA: R² = %.2f, p = %.3f", perm$R2[1], perm$`Pr(>F)`[1]), colour = "Categoría") +
  theme_bw() + theme(legend.position = "bottom", text = element_text(size = 12), plot.title = element_text(face = "bold"))

fig <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1.7, 1))
ggsave("Complementarios_nativo.png", plot = fig, path = out_tesis, width = 32, height = 13, dpi = 300, units = "cm", bg = "white")
ggsave("Complementarios_nativo.png", plot = fig, path = out_res,   width = 32, height = 13, dpi = 300, units = "cm", bg = "white")
cat("Figura guardada en", out_tesis, "\n")
