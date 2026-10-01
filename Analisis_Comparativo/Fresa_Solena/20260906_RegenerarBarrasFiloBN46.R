## Figura 4.6 - Composicion media por filo en cada grupo (Img/cap3/BarrasFilo_blackandwhite.png)
## Sustituye la antigua grafica de barras con contorno blanco/negro por tratamiento (funcion
## Abundance_barras() de 20230227_Funciones&Graficas.R), que repetia la Figura 4.5 sin anadir
## informacion. Ahora:
##   A) abundancia relativa media de cada filo en el grupo saludable y en el no saludable
##   B) los mismos filos, muestra por muestra y en escala logaritmica, para que los filos
##      minoritarios (menos del 3 %) sean visibles. Valor p de Mann-Whitney por filo (ajuste BH).
## El nombre del archivo se conserva para no romper las referencias de la tesis.

source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_paleta_tesis.R")
suppressMessages(library("reshape2"))
datos <- cargar_fresa()

n_top <- 10
df <- psmelt(tax_glom(datos$rel, taxrank = "Phylum"))
top <- names(sort(tapply(df$Abundance, df$Phylum, mean), decreasing = TRUE))[1:n_top]
df$Filo <- factor(ifelse(df$Phylum %in% top, df$Phylum, "Otros"), levels = c(top, "Otros"))
pal <- paleta_para_filos(levels(df$Filo))

## ---- Panel A: media por grupo ----
medias <- aggregate(Abundance ~ Filo + Tratamiento, df, mean)
cat("Abundancia relativa media (%) por filo y grupo:\n")
tab_medias <- dcast(medias, Filo ~ Tratamiento, value.var = "Abundance"); tab_medias[, -1] <- round(tab_medias[, -1], 3); print(tab_medias)
pA <- ggplot(medias, aes(x = Tratamiento, y = Abundance, fill = Filo)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Abundancia relativa media (%)", fill = "Filo", title = "A) Composición media por grupo") +
  tema_tesis() + theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))

## ---- Panel B: cada filo, muestra por muestra, escala logaritmica ----
df_top <- df[df$Filo != "Otros", ]
df_top <- aggregate(Abundance ~ Sample + Tratamiento + Filo, df_top, sum)
mw <- sapply(top, function(f) { s <- df_top[df_top$Filo == f, ]; wilcox.test(Abundance ~ Tratamiento, s, exact = FALSE)$p.value })
mw_bh <- p.adjust(mw, "BH")
cat("Mann-Whitney por filo (p y p ajustado BH):\n"); print(round(cbind(p = mw, p_BH = mw_bh), 3))
etq <- data.frame(Filo = factor(top, levels = levels(df$Filo)), texto = sprintf("p = %.2f", mw_bh))
pB <- ggplot(df_top, aes(x = Tratamiento, y = Abundance, colour = Tratamiento)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  geom_text(data = etq, aes(x = 1.5, y = Inf, label = texto), inherit.aes = FALSE, vjust = 1.5, size = 3) +
  facet_wrap(~ Filo, nrow = 2, scales = "free_y") +
  scale_y_log10() +
  scale_colour_manual(values = colores_trat) +
  labs(x = NULL, y = "Abundancia relativa (%, escala logarítmica)", colour = "Tratamiento",
       title = "B) Filos principales por muestra; valor p de Mann-Whitney ajustado (BH)") +
  tema_tesis() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

fig <- plot_grid(pA, pB, ncol = 1, rel_heights = c(0.9, 1.3))
guardar_figura("BarrasFilo_blackandwhite.png", fig, out_cap3, ancho = 30, alto = 26)
