## Figura 3.5 - Conjunto complementario de cinco categorias: composicion y diversidad beta
##               (Img/cap2/Complementarios_beta.png)
## Regenera AllData_Beta_diversidad.png y el StackBar de 20230321/20230327 con buen diseno, y agrega
## PERMANOVA y PERMDISP con la categoria como factor (cinco niveles).

source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_cargar_complementarios.R")

## ---- Panel A: abundancia relativa por filo, muestras agrupadas por categoria ----
filo <- tax_glom(rel_todo, taxrank = "Phylum")
df <- psmelt(filo)
top <- names(sort(tapply(df$Abundance, df$Phylum, mean), decreasing = TRUE))[1:8]
df$Filo <- ifelse(df$Phylum %in% top, df$Phylum, "Otros")
df$Filo <- factor(df$Filo, levels = c(top, "Otros"))
orden <- unique(df[order(df$Categoria, df$Sample), "Sample"])
df$Sample <- factor(df$Sample, levels = orden)
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_paleta_tesis.R")
pal <- paleta_para_filos(c(top, "Otros"))   # mismo color por filo en toda la tesis
pA <- ggplot(df, aes(Sample, 100 * Abundance, fill = Filo)) + geom_col(width = 1) +
  facet_grid(~ Categoria, scales = "free_x", space = "free_x", labeller = label_wrap_gen(12)) +
  scale_fill_manual(values = pal) +
  labs(x = "Muestras", y = "Abundancia relativa (%)", fill = "Filo", title = "A) Composición por filo") +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 12),
                     strip.text = element_text(size = 9), legend.position = "right", plot.title = element_text(face = "bold"))
cat("Abundancia relativa media (%) de los dos filos dominantes por categoria:\n")
print(round(100 * with(df[df$Phylum %in% top[1:2], ], tapply(Abundance, list(Categoria, Phylum), mean)), 1))

## ---- Panel B: NMDS Bray-Curtis ----
D <- phyloseq::distance(rel_todo, method = "bray")
set.seed(2026); ord <- ordinate(rel_todo, method = "NMDS", distance = D, trymax = 50, trace = 0)
sc <- data.frame(scores(ord, display = "sites")); sc$Categoria <- cat5[rownames(sc)]
cen <- aggregate(cbind(NMDS1, NMDS2) ~ Categoria, sc, mean)
pB <- ggplot(sc, aes(NMDS1, NMDS2, colour = Categoria)) +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_point(data = cen, aes(fill = Categoria), shape = 23, size = 5, colour = "black") +
  scale_colour_manual(values = colores5) + scale_fill_manual(values = colores5, guide = "none") +
  labs(title = sprintf("B) NMDS, Bray-Curtis (estrés = %.3f)", ord$stress), colour = "Categoría") +
  theme_bw() + theme(legend.position = "bottom", text = element_text(size = 12), plot.title = element_text(face = "bold")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3)))

## ---- PERMANOVA y PERMDISP con cinco categorias ----
meta <- data.frame(sample_data(todo_fil)); stopifnot(identical(rownames(meta), attr(D, "Labels")))
set.seed(2026); perm <- adonis2(D ~ Categoria, data = meta, permutations = 999)
disp <- betadisper(D, meta$Categoria, type = "centroid"); set.seed(2026); pdisp <- permutest(disp, permutations = 999)
cat("== PERMANOVA (5 categorias) ==\n"); print(perm)
cat("== PERMDISP ==\n"); print(pdisp$tab); cat("Distancia media al centroide por categoria:\n"); print(round(tapply(disp$distances, meta$Categoria, mean), 3))
## PERMANOVA por pares entre categorias (BH)
pares <- combn(levels(meta$Categoria), 2)
pp <- apply(pares, 2, function(p) { s <- meta$Categoria %in% p; set.seed(2026)
  a <- adonis2(as.dist(as.matrix(D)[s, s]) ~ droplevels(meta$Categoria[s]), permutations = 999); c(R2 = a$R2[1], p = a$`Pr(>F)`[1]) })
pares_df <- data.frame(grupo1 = pares[1, ], grupo2 = pares[2, ], R2 = round(pp["R2", ], 3), p = pp["p", ], p_BH = round(p.adjust(pp["p", ], "BH"), 3))
cat("== PERMANOVA por pares ==\n"); print(pares_df)
write.csv(pares_df, file.path(out_res, "Complementarios_permanova_pares.csv"), row.names = FALSE)
pdisp_df <- data.frame(prueba = c("PERMANOVA", "PERMDISP"), F = c(perm$F[1], pdisp$tab$F[1]), R2 = c(perm$R2[1], NA),
                       p = c(perm$`Pr(>F)`[1], pdisp$tab$`Pr(>F)`[1]))
write.csv(pdisp_df, file.path(out_res, "Complementarios_permanova.csv"), row.names = FALSE)

fig <- plot_grid(pA, pB, ncol = 1, rel_heights = c(1, 1.1))
ggsave("Complementarios_beta.png", plot = fig, path = out_tesis, width = 30, height = 26, dpi = 300, units = "cm", bg = "white")
ggsave("Complementarios_beta.png", plot = fig, path = out_res,   width = 30, height = 26, dpi = 300, units = "cm", bg = "white")
cat("Figura guardada en", out_tesis, "\n")
