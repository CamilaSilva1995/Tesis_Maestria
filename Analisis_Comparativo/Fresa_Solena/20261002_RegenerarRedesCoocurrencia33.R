## Figura 3.3 - Redes de coocurrencia entre generos (Img/cap2/Redes_coocurrencia.png)
##
## Usa las salidas originales de MicNet y Alnitak guardadas en Data1/Redes (marzo de 2023):
##   - SparCC_Output.csv: matriz de correlacion SparCC 1795 x 1795 entre generos (MicNet)
##   - Output_UMAP_HDBSCAN .csv: cluster de cada genero (UMAP + HDBSCAN, MicNet)
##   - 20230314-1_raw_network.csv: vecinos de Fusarium segun Alnitak (Spearman >= 0.7 y Bray-Curtis <= 0.3)
## Paneles:
##   A) Red de coocurrencia entre generos con |r| >= 0.8 (componente gigante), color por filo,
##      tamano por grado. Fusarium no tiene vecinos a este umbral.
##   B) Vecindario de Fusarium con |r| >= 0.5: sus 40 vecinos mas correlacionados.
##   C) Vecinos de Fusarium segun Alnitak sobre la tabla completa.
## Tambien calcula metricas de la red a varios umbrales (Results_img/Redes_metricas.csv).
## Los identificadores de la matriz son los del BIOM de generos; se traducen a nombres con la
## taxonomia de fresa_kraken.biom. Semilla unica de la tesis: 2026.

suppressMessages({library("jsonlite"); library("igraph"); library("ggplot2"); library("cowplot"); library("RColorBrewer")})
setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
out_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap2"
out_res   <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
set.seed(2026)

## ---- Taxonomia por identificador ----
bk <- fromJSON("fresa_kraken.biom", simplifyVector = FALSE)
tax <- do.call(rbind, lapply(bk$rows, function(r) c(id = r$id, substr(unlist(r$metadata$taxonomy), 4, 100))))
colnames(tax) <- c("id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax <- as.data.frame(tax, stringsAsFactors = FALSE); rownames(tax) <- tax$id
nombre <- function(ids) sapply(ids, function(i) { v <- unlist(tax[i, 2:7]); v <- v[v != ""]; if (length(v)) unname(tail(v, 1)) else i })

## ---- Matriz SparCC ----
sp <- read.csv("Redes/SparCC_Output.csv", row.names = 1, check.names = FALSE)
M <- as.matrix(sp); rownames(M) <- colnames(M) <- as.character(rownames(sp)); diag(M) <- 0
stopifnot(all(rownames(M) %in% tax$id))
cl <- read.csv("Redes/Output_UMAP_HDBSCAN .csv", row.names = 1)
stopifnot(nrow(cl) == nrow(M))
filo <- tax[rownames(M), "Phylum"]; genero <- tax[rownames(M), "Genus"]
id_fus <- rownames(M)[genero == "Fusarium"]; stopifnot(length(id_fus) == 1)
cat("Generos:", nrow(M), "| Fusarium =", id_fus, "| cluster HDBSCAN de Fusarium:", cl$Cluster[match(id_fus, rownames(M))], "\n")

## ---- Metricas de la red a varios umbrales ----
metricas <- do.call(rbind, lapply(c(0.5, 0.6, 0.7, 0.8), function(u) {
  A <- abs(M) >= u; g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  g <- delete_vertices(g, degree(g) == 0)
  comp <- components(g); set.seed(2026); mod <- modularity(cluster_louvain(g))
  data.frame(umbral = u, nodos = vcount(g), aristas = ecount(g), negativas = sum(M[upper.tri(M)] <= -u),
             densidad = round(edge_density(g), 4), componentes = comp$no, gigante = max(comp$csize),
             coef_agrupamiento = round(transitivity(g, type = "global"), 3), modularidad = round(mod, 3),
             grado_medio = round(mean(degree(g)), 1), grado_Fusarium = sum(abs(M[id_fus, ]) >= u))
}))
print(metricas); write.csv(metricas, file.path(out_res, "Redes_metricas.csv"), row.names = FALSE)

## ---- Panel A: red |r| >= 0.8, componente gigante ----
u <- 0.8
gA <- graph_from_adjacency_matrix(abs(M) >= u, mode = "undirected", diag = FALSE)
gA <- delete_vertices(gA, degree(gA) == 0)
comp <- components(gA); gA <- induced_subgraph(gA, which(comp$membership == which.max(comp$csize)))
filos_top <- names(sort(table(filo), decreasing = TRUE))[1:7]
V(gA)$filo <- ifelse(tax[V(gA)$name, "Phylum"] %in% filos_top, tax[V(gA)$name, "Phylum"], "Otros")
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_paleta_tesis.R")
pal <- paleta_para_filos(c(filos_top, "Otros"))   # mismo color por filo en toda la tesis
set.seed(2026); layA <- layout_with_fr(gA)
dfA <- data.frame(x = layA[, 1], y = layA[, 2], filo = V(gA)$filo, grado = degree(gA))
eA <- as_edgelist(gA, names = FALSE)
segA <- data.frame(x = layA[eA[, 1], 1], y = layA[eA[, 1], 2], xend = layA[eA[, 2], 1], yend = layA[eA[, 2], 2])
pA <- ggplot() +
  geom_segment(data = segA, aes(x, y, xend = xend, yend = yend), colour = "grey80", linewidth = 0.15, alpha = 0.5) +
  geom_point(data = dfA, aes(x, y, colour = filo, size = grado), alpha = 0.9) +
  scale_colour_manual(values = pal, name = "Filo") + scale_size(range = c(0.8, 4), name = "Grado") +
  labs(title = sprintf("A) Coocurrencia entre géneros, |r| ≥ %.1f", u),
       subtitle = sprintf("Componente gigante: %d géneros, %d aristas", vcount(gA), ecount(gA))) +
  theme_void() + theme(legend.position = "right", text = element_text(size = 11), plot.title = element_text(face = "bold")) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

## ---- Panel B: vecindario de Fusarium, |r| >= 0.5 ----
uB <- 0.5
vec <- names(which(abs(M[id_fus, ]) >= uB))
gB <- graph_from_adjacency_matrix(abs(M[c(id_fus, vec), c(id_fus, vec)]) >= uB, mode = "undirected", diag = FALSE)
E(gB)$r <- M[c(id_fus, vec), c(id_fus, vec)][as_edgelist(gB, names = FALSE)]
set.seed(2026); layB <- layout_with_fr(gB)
dfB <- data.frame(x = layB[, 1], y = layB[, 2], nombre = nombre(V(gB)$name),
                  filo = ifelse(tax[V(gB)$name, "Phylum"] %in% filos_top, tax[V(gB)$name, "Phylum"], "Otros"),
                  es_fus = V(gB)$name == id_fus, r_fus = M[id_fus, V(gB)$name])
eB <- as_edgelist(gB, names = FALSE)
segB <- data.frame(x = layB[eB[, 1], 1], y = layB[eB[, 1], 2], xend = layB[eB[, 2], 1], yend = layB[eB[, 2], 2],
                   con_fus = (V(gB)$name[eB[, 1]] == id_fus | V(gB)$name[eB[, 2]] == id_fus), r = E(gB)$r)
top5 <- names(sort(M[id_fus, vec], decreasing = TRUE))[1:6]
dfB$etiqueta <- ifelse(dfB$es_fus | V(gB)$name %in% top5, dfB$nombre, "")
pB <- ggplot() +
  geom_segment(data = segB, aes(x, y, xend = xend, yend = yend, alpha = con_fus, linewidth = abs(r)), colour = "grey50") +
  scale_alpha_manual(values = c(`FALSE` = 0.15, `TRUE` = 0.8), guide = "none") + scale_linewidth(range = c(0.2, 1.2), guide = "none") +
  geom_point(data = dfB, aes(x, y, colour = filo, size = ifelse(es_fus, 7, 3.2)), alpha = 0.95) +
  scale_colour_manual(values = pal, guide = "none") + scale_size_identity() +
  geom_text(data = dfB, aes(x, y, label = etiqueta), size = 3, fontface = "italic", vjust = -1.1) +
  labs(title = sprintf("B) Vecindario de Fusarium, |r| ≥ %.1f", uB),
       subtitle = sprintf("%d vecinos; correlación máxima r = %.2f", length(vec), max(M[id_fus, vec]))) +
  theme_void() + theme(text = element_text(size = 11), plot.title = element_text(face = "bold")) +
  coord_cartesian(clip = "off")

## ---- Panel C: vecinos de Fusarium segun Alnitak ----
al <- read.csv("Redes/20230314-1_raw_network.csv", check.names = FALSE)
al$taxon1 <- as.character(al$taxon1); al$taxon2 <- as.character(al$taxon2)
etq <- function(i) { v <- unlist(tax[i, 2:7]); v <- v[v != ""]; n <- unname(tail(v, 1))
  if (tax[i, "Species"] != "") paste(tax[i, "Genus"], tax[i, "Species"]) else if (tax[i, "Genus"] == "") paste0(n, " (", c(Phylum="filo",Class="clase",Order="orden",Family="familia")[names(v)[length(v)]], ")") else n }
nodosC <- unique(c(al$taxon1, al$taxon2))
gC <- graph_from_data_frame(al[, c("taxon1", "taxon2")], directed = FALSE, vertices = nodosC)
set.seed(2026); layC <- layout_as_star(gC, center = V(gC)[name == "5506"])
dfC <- data.frame(x = layC[, 1], y = layC[, 2], nombre = sapply(V(gC)$name, etq), es_fus = V(gC)$name == "5506")
eC <- as_edgelist(gC, names = FALSE)
segC <- data.frame(x = layC[eC[, 1], 1], y = layC[eC[, 1], 2], xend = layC[eC[, 2], 1], yend = layC[eC[, 2], 2],
                   rho = al$`Spearman Correlation`, bc = al$`Bray Curtis Dissimilarity`)
pC <- ggplot() +
  geom_segment(data = segC, aes(x, y, xend = xend, yend = yend, linewidth = rho), colour = "grey40") +
  scale_linewidth(range = c(0.6, 2), guide = "none") +
  geom_label(data = segC, aes((x + xend) / 2, (y + yend) / 2, label = sprintf("ρ = %.2f\nd = %.2f", rho, bc)), size = 2.8, label.size = 0) +
  geom_point(data = dfC, aes(x, y, size = ifelse(es_fus, 9, 6)), colour = ifelse(dfC$es_fus, "#00BFC4", "#F8766D")) +
  scale_size_identity() +
  geom_text(data = dfC, aes(x, y, label = nombre), size = 3.1, fontface = "italic", vjust = -1.6) +
  labs(title = "C) Vecinos de Fusarium según Alnitak",
       subtitle = "Spearman ρ ≥ 0.7 y Bray-Curtis d ≤ 0.3, tabla completa") +
  theme_void() + theme(text = element_text(size = 11), plot.title = element_text(face = "bold")) +
  coord_cartesian(clip = "off", xlim = range(layC[, 1]) * 1.5, ylim = range(layC[, 2]) * 1.5)

abajo <- plot_grid(pB, pC, ncol = 2)
fig <- plot_grid(pA, abajo, ncol = 1, rel_heights = c(1.15, 1))
ggsave("Redes_coocurrencia.png", plot = fig, path = out_tesis, width = 30, height = 30, dpi = 300, units = "cm", bg = "white")
ggsave("Redes_coocurrencia.png", plot = fig, path = out_res,   width = 30, height = 30, dpi = 300, units = "cm", bg = "white")
cat("Figura guardada en", out_tesis, "\n")
cat("Vecinos de Fusarium (SparCC, |r| >= 0.5), los seis mayores:\n")
print(round(sort(M[id_fus, vec], decreasing = TRUE)[1:6], 3)); cat(paste(nombre(top5), collapse = ", "), "\n")
