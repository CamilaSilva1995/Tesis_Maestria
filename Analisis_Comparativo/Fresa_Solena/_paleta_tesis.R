## Paleta, tema y cargador comunes para todas las figuras de la tesis.
## Lo usan por source() los scripts 2026*_Regenerar*.R. Garantiza que:
##   - el grupo saludable es siempre "#F8766D" (coral) y el no saludable "#00BFC4" (turquesa),
##   - cada filo tiene siempre el mismo color en todas las figuras (barras y redes),
##   - los taxones de familia, genero y especie usan una misma paleta cualitativa,
##   - todas las figuras comparten tipografia, fondo blanco y 300 dpi.
## Semilla unica de la tesis para todo procedimiento aleatorio: 2026.

suppressMessages({library("phyloseq"); library("ggplot2"); library("cowplot"); library("RColorBrewer")})
set.seed(2026)

ruta_datos  <- "/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1"
out_res     <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
out_cap3    <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"   # figuras del capitulo 4 (Resultados)
out_cap2    <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap2"   # figuras del capitulo 3 (Metodologia)
out_cap1    <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap1"   # figuras del capitulo 1 (Preliminares)

## ---- Colores ----
colores_trat <- c("Saludable" = "#F8766D", "No saludable" = "#00BFC4")

## Filos principales del conjunto (en orden de abundancia media) con un color fijo cada uno
paleta_filos <- c(
  "Actinobacteria"   = "#8DA0CB", "Proteobacteria"  = "#66C2A5", "Planctomycetes"  = "#FFD92F",
  "Firmicutes"       = "#FC8D62", "Bacteroidetes"   = "#E78AC3", "Acidobacteria"   = "#1B9E77",
  "Ascomycota"       = "#E5C494", "Cyanobacteria"   = "#A6D854", "Chloroflexi"     = "#7570B3",
  "Verrucomicrobia"  = "#D95F02", "Gemmatimonadetes" = "#B3B3B3", "Oomycota"        = "#E7298A",
  "Basidiomycota"    = "#A6761D", "Nitrospirae"     = "#66A61E", "Deinococcus-Thermus" = "#1F78B4",
  "Otros"            = "grey75")
## Reserva para filos que no esten en la lista anterior
reserva_filos <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Pastel2"))

## Devuelve una paleta con nombres para el vector de filos dado (siempre el mismo color por filo)
paleta_para_filos <- function(nombres) {
  nombres <- unique(as.character(nombres))
  pal <- paleta_filos[nombres]; names(pal) <- nombres
  faltan <- is.na(pal)
  if (any(faltan)) pal[faltan] <- reserva_filos[seq_len(sum(faltan))]
  pal
}

## Paleta cualitativa para familias, generos y especies ("Otros" siempre gris)
paleta_taxones <- function(nombres) {
  nombres <- unique(as.character(nombres)); otros <- nombres == "Otros"
  base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF",
            "#1B9E77", "#E6AB02", "#7570B3", "#66A61E", "#E7298A", "#1F78B4", "#B2DF8A",
            "#FDBF6F", "#CAB2D6", "#FFFF99", "#B15928")
  pal <- setNames(rep("grey75", length(nombres)), nombres)
  pal[!otros] <- base[seq_len(sum(!otros))]
  pal
}

## ---- Tema ----
tema_tesis <- function(base = 12) {
  theme_bw(base_size = base) +
    theme(legend.position = "bottom", legend.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = base - 2),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
}

## Guarda una figura en la carpeta de la tesis y en Results_img
guardar_figura <- function(nombre, fig, carpeta_tesis, ancho = 30, alto = 16) {
  ggsave(nombre, plot = fig, path = carpeta_tesis, width = ancho, height = alto, dpi = 300, units = "cm", bg = "white")
  ggsave(nombre, plot = fig, path = out_res,       width = ancho, height = alto, dpi = 300, units = "cm", bg = "white")
  cat("Figura guardada:", file.path(carpeta_tesis, nombre), "\n")
}

## ---- Datos principales (Data1: 58 muestras, 53 tras el filtro de calidad) ----
muestras_eliminadas <- c("MP2079", "MP2080", "MP2088", "MP2109", "MP2137")

cargar_fresa <- function() {
  ps <- import_biom(file.path(ruta_datos, "fresa_kraken.biom"))
  colnames(ps@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  ps@tax_table@.Data <- substr(ps@tax_table@.Data, 4, 100)
  colnames(ps@otu_table@.Data) <- substr(colnames(ps@otu_table@.Data), 1, 6)
  md <- read.csv2(file.path(ruta_datos, "metadata.csv"), header = FALSE, row.names = 1, sep = ",")
  ps@sam_data <- sample_data(md)
  colnames(ps@sam_data) <- "Treatment"
  sample_data(ps)$Tratamiento <- factor(sample_data(ps)$Treatment, levels = c("healthy", "wilted"),
                                        labels = names(colores_trat))
  fil <- prune_samples(!(sample_names(ps) %in% muestras_eliminadas), ps)
  list(crudo = ps, fil = fil, rel = transform_sample_counts(fil, function(x) 100 * x / sum(x)))
}

## Vector con nombres: tratamiento de cada muestra
tratamiento_de <- function(ps) { t <- sample_data(ps)$Tratamiento; names(t) <- sample_names(ps); t }
