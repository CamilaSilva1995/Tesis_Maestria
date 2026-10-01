## Figura 3.4 - Conjunto complementario de cinco categorias: profundidad y diversidad alfa
##               (Img/cap2/Complementarios_alfa.png)
## Regenera, con buen diseno y 300 dpi, la figura AllData_Alfa_diversidad.png de 20230321_NuevosDatosAll.R
## y la acompana de una prueba de Kruskal-Wallis entre las cinco categorias para cada indice.

source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_cargar_complementarios.R")

## ---- Panel A: lecturas clasificadas por muestra ----
prof <- data.frame(muestra = sample_names(todo_fil), lecturas = sample_sums(todo_fil) / 1e6, Categoria = cat5)
prof <- prof[order(prof$Categoria, -prof$lecturas), ]; prof$muestra <- factor(prof$muestra, levels = prof$muestra)
pA <- ggplot(prof, aes(muestra, lecturas, fill = Categoria)) + geom_col() +
  scale_fill_manual(values = colores5) +
  labs(x = "Muestra", y = "Millones de lecturas clasificadas", fill = "Categoría",
       title = "A) Profundidad por muestra") +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     legend.position = "bottom", text = element_text(size = 12), plot.title = element_text(face = "bold")) +
  guides(fill = guide_legend(nrow = 2))
cat("Lecturas clasificadas (millones), media por categoria:\n"); print(round(tapply(prof$lecturas, prof$Categoria, mean), 1))

## ---- Panel B: diversidad alfa por categoria ----
alfa <- estimate_richness(todo_fil, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
stopifnot(identical(rownames(alfa), names(cat5)))
alfa$Categoria <- cat5
larga <- melt(alfa, id.vars = "Categoria", measure.vars = c("Observed", "Chao1", "Shannon", "Simpson"),
              variable.name = "Indice", value.name = "Valor")
levels(larga$Indice) <- c("Observado", "Chao1", "Shannon", "Simpson")
## Kruskal-Wallis por indice (cinco grupos) y prueba post hoc de Dunn simplificada con Mann-Whitney + BH
kw <- sapply(c("Observed", "Chao1", "Shannon", "Simpson"), function(v) kruskal.test(alfa[[v]] ~ alfa$Categoria)$p.value)
cat("Kruskal-Wallis (5 categorias), valor p por indice:\n"); print(signif(kw, 3))
cat("Medias por categoria:\n"); m <- aggregate(. ~ Categoria, alfa[, c("Observed","Chao1","Shannon","Simpson","Categoria")], mean); m[,-1] <- round(m[,-1], 3); print(m)
pw <- pairwise.wilcox.test(alfa$Shannon, alfa$Categoria, p.adjust.method = "BH", exact = FALSE)
cat("Shannon, Mann-Whitney por pares con ajuste BH:\n"); print(round(pw$p.value, 3))
write.csv(data.frame(indice = names(kw), p_kruskal = signif(kw, 4)), file.path(out_res, "Complementarios_kruskal.csv"), row.names = FALSE)
etq <- data.frame(Indice = factor(c("Observado", "Chao1", "Shannon", "Simpson"), levels = levels(larga$Indice)),
                  texto = sprintf("Kruskal-Wallis p = %.3g", kw))
pB <- ggplot(larga, aes(Categoria, Valor, colour = Categoria)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) + geom_jitter(width = 0.12, size = 1.6, alpha = 0.8) +
  facet_wrap(~ Indice, scales = "free_y", nrow = 1) +
  geom_text(data = etq, aes(x = 3, y = Inf, label = texto), inherit.aes = FALSE, vjust = 1.4, size = 3.2) +
  scale_colour_manual(values = colores5, guide = "none") +
  labs(x = NULL, y = "Medida de diversidad alfa", title = "B) Diversidad alfa por categoría") +
  theme_bw() + theme(text = element_text(size = 12), axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
                     plot.title = element_text(face = "bold"))

fig <- plot_grid(pA, pB, ncol = 1, rel_heights = c(0.9, 1.3))
ggsave("Complementarios_alfa.png", plot = fig, path = out_tesis, width = 30, height = 24, dpi = 300, units = "cm", bg = "white")
ggsave("Complementarios_alfa.png", plot = fig, path = out_res,   width = 30, height = 24, dpi = 300, units = "cm", bg = "white")
cat("Figura guardada en", out_tesis, "\n")
