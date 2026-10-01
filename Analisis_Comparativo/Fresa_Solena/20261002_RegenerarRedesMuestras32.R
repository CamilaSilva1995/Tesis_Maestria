## Figura 3.2 - Redes de muestras con phyloseq (Img/cap2/Redes_muestras.png)
##
## Regenera, con buen diseno y 300 dpi, las dos redes simples del script 20230314_Redes.R:
##   A) make_network + plot_network: distancia de Jaccard (cuantitativa, sobre abundancias
##      relativas), arista si d <= 0.16 (aprox. el 25 % de los pares mas parecidos)
##   B) plot_net: misma distancia, arista si d <= 0.13 (aprox. el 10 % de los pares mas parecidos)
## Los umbrales originales de 2023 (0.8 y 0.5) producen el grafo completo, porque la distancia de
## Jaccard cuantitativa entre estas muestras va de 0.08 a 0.47; por eso se fijan a partir de los
## cuantiles de la distribucion de distancias. Cada nodo es una muestra (53, tras el filtro de
## calidad); dos muestras se conectan si su composicion es suficientemente parecida. Se busca ver
## si las muestras saludables y no saludables forman grupos separados, y se cuantifica con la
## fraccion de aristas dentro de un mismo grupo frente a la esperada por azar (prueba por permutacion).
## Semilla unica de la tesis: 2026 (afecta al trazado fruchterman-reingold).

suppressMessages({library("phyloseq"); library("ggplot2"); library("igraph"); library("cowplot")})
setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
out_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap2"
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
sample_data(fil)$Tratamiento <- factor(sample_data(fil)$Treatment, levels = c("healthy", "wilted"),
                                       labels = c("Saludable", "No saludable"))
rel <- transform_sample_counts(fil, function(x) x / sum(x))

## ---- Red A: Jaccard, d <= 0.8 (make_network + plot_network) ----
set.seed(2026)
igA <- make_network(rel, type = "samples", distance = "jaccard", max.dist = 0.16)
cat("Red A (Jaccard <= 0.16): nodos =", vcount(igA), " aristas =", ecount(igA),
    " componentes =", components(igA)$no, "\n")
set.seed(2026)
pA <- plot_network(igA, rel, type = "samples", color = "Tratamiento", shape = "Tratamiento",
                   label = NULL, point_size = 4, line_alpha = 0.25, line_weight = 0.4) +
  scale_colour_manual(values = colores) +
  labs(colour = "Tratamiento", shape = "Tratamiento",
       title = "A) Jaccard, arista si d ≤ 0.16",
       subtitle = sprintf("25 %% de los pares más parecidos: %d muestras, %d aristas", vcount(igA), ecount(igA))) +
  guides(size = "none", linewidth = "none", alpha = "none", colour = guide_legend(override.aes = list(size = 4))) +
  theme_void() + theme(legend.position = "bottom", text = element_text(size = 12),
                       plot.title = element_text(face = "bold"))

## ---- Red B: Jaccard, d <= 0.5 (plot_net) ----
set.seed(2026)
igB <- make_network(rel, type = "samples", distance = "jaccard", max.dist = 0.13)
cat("Red B (Jaccard <= 0.13): nodos =", vcount(igB), " aristas =", ecount(igB),
    " componentes =", components(igB)$no, " aislados =", sum(degree(igB) == 0), "\n")
set.seed(2026)
pB <- plot_net(rel, distance = "jaccard", type = "samples", maxdist = 0.13,
               color = "Tratamiento", shape = "Tratamiento", point_size = 4,
               laymeth = "fruchterman.reingold") +
  scale_colour_manual(values = colores) +
  labs(colour = "Tratamiento", shape = "Tratamiento",
       title = "B) Jaccard, arista si d ≤ 0.13",
       subtitle = sprintf("10 %% de los pares más parecidos: %d muestras, %d aristas", sum(degree(igB) > 0), ecount(igB))) +
  guides(size = "none", linewidth = "none", alpha = "none", colour = guide_legend(override.aes = list(size = 4))) +
  theme_void() + theme(legend.position = "bottom", text = element_text(size = 12),
                       plot.title = element_text(face = "bold"))

## ---- Metricas por grupo: proporcion de aristas dentro y entre grupos ----
grupo <- sample_data(rel)$Tratamiento; names(grupo) <- sample_names(rel)
for (nm in c("igA", "igB")) {
  g <- get(nm); e <- as_edgelist(g)
  mismo <- grupo[e[, 1]] == grupo[e[, 2]]
  cat(sprintf("%s: aristas dentro del mismo grupo = %d (%.1f%%), entre grupos = %d (%.1f%%)\n",
              nm, sum(mismo), 100 * mean(mismo), sum(!mismo), 100 * mean(!mismo)))
}
## Esperado bajo independencia: fraccion de pares posibles que son del mismo grupo
n1 <- sum(grupo == "Saludable"); n2 <- sum(grupo == "No saludable")
cat(sprintf("Fraccion esperada de pares del mismo grupo si el grupo no importara: %.1f%%\n",
            100 * (choose(n1, 2) + choose(n2, 2)) / choose(n1 + n2, 2)))
## Prueba por permutacion: se barajan las etiquetas de grupo 999 veces y se recalcula la fraccion
## de aristas dentro de un mismo grupo (asortatividad por grupo)
set.seed(2026)
for (nm in c("igA", "igB")) {
  g <- get(nm); e <- as_edgelist(g)
  obs <- mean(grupo[e[, 1]] == grupo[e[, 2]])
  perm <- replicate(999, { gp <- sample(grupo); names(gp) <- names(grupo); mean(gp[e[, 1]] == gp[e[, 2]]) })
  cat(sprintf("%s: fraccion observada dentro de grupo = %.3f, media permutada = %.3f, p = %.3f, asortatividad = %.3f\n",
              nm, obs, mean(perm), (sum(perm >= obs) + 1) / 1000, assortativity_nominal(g, as.integer(grupo[V(g)$name]))))
}
## Distancia media dentro de cada grupo
D <- as.matrix(phyloseq::distance(rel, "jaccard"))
cat(sprintf("Distancia de Jaccard media: saludable-saludable = %.3f, no saludable-no saludable = %.3f, entre grupos = %.3f\n",
            mean(D[grupo == "Saludable", grupo == "Saludable"][lower.tri(diag(n1))]),
            mean(D[grupo == "No saludable", grupo == "No saludable"][lower.tri(diag(n2))]),
            mean(D[grupo == "Saludable", grupo == "No saludable"])))

fig <- plot_grid(pA, pB, ncol = 2)
ggsave("Redes_muestras.png", plot = fig, path = out_tesis, width = 30, height = 15, dpi = 300, units = "cm", bg = "white")
ggsave("Redes_muestras.png", plot = fig, path = out_res,   width = 30, height = 15, dpi = 300, units = "cm", bg = "white")
cat("Figura guardada en", out_tesis, "\n")
