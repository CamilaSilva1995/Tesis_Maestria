## Cargador comun para los conjuntos de datos complementarios (Seccion 3.4 de la tesis).
## Lo usan los scripts 20261002_RegenerarComplementarios34.R, ..35.R y ..36.R mediante source().
##
## Conjunto "Data_all": 85 muestras = Data1 (58, cultivo saludable / marchito) + Data2 (27, tres
##   tipos de sitio). Se aplica el mismo filtro de calidad de la tesis (25 millones de lecturas tras
##   fastp), que elimina las cinco muestras ya conocidas; ninguna muestra de Data2 queda por debajo.
## Conjunto "Data3": 8 muestras, cultivo (4) frente a suelo nativo (4).
## Semilla unica de la tesis: 2026.

suppressMessages({library("phyloseq"); library("vegan"); library("ggplot2"); library("cowplot"); library("reshape2"); library("RColorBrewer")})
base <- "/home/camila/GIT/Tesis_Maestria/Data/fresa_solena"
out_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap2"
out_res   <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
set.seed(2026)

cargar_biom <- function(archivo) {
  ps <- import_biom(archivo)
  colnames(ps@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  ps@tax_table@.Data <- substr(ps@tax_table@.Data, 4, 100)
  sample_names(ps) <- sub("\\.kraken2\\.report$", "", sample_names(ps))
  ps
}

## ---- Data_all: cinco categorias ----
etiquetas5 <- c(cropHealthy = "Cultivo saludable", cropWilted = "Cultivo marchito", nearCrop = "Cerca del cultivo",
                uncultivatedLand = "Suelo no cultivado", forestDegraded = "Bosque degradado")
colores5 <- c("Cultivo saludable" = "#F8766D", "Cultivo marchito" = "#00BFC4", "Cerca del cultivo" = "#7CAE00",
              "Suelo no cultivado" = "#C77CFF", "Bosque degradado" = "#E69F00")
todo <- cargar_biom(file.path(base, "Data_all/fresa_kraken_all.biom"))
md_todo <- read.csv(file.path(base, "Data_all/metadata.csv"), row.names = 1)
stopifnot(all(sample_names(todo) %in% rownames(md_todo)))
md_todo <- md_todo[sample_names(todo), , drop = FALSE]
md_todo$Categoria <- factor(etiquetas5[md_todo$Category], levels = etiquetas5)
sample_data(todo) <- sample_data(md_todo)
samples_to_remove <- c("MP2079", "MP2080", "MP2088", "MP2109", "MP2137")
todo_fil <- prune_samples(!(sample_names(todo) %in% samples_to_remove), todo)
cat5 <- sample_data(todo_fil)$Categoria; names(cat5) <- sample_names(todo_fil)
stopifnot(identical(names(cat5), sample_names(todo_fil)))
rel_todo <- transform_sample_counts(todo_fil, function(x) x / sum(x))

## ---- Data3: cultivo frente a suelo nativo ----
etiquetas2 <- c(crop = "Cultivo", native = "Suelo nativo")
colores2 <- c("Cultivo" = "#F8766D", "Suelo nativo" = "#2E8B57")
nat <- cargar_biom(file.path(base, "Data3/crop_vs_native_kraken.biom"))
md_nat <- read.csv(file.path(base, "Data3/metadata.csv"), row.names = 1)
stopifnot(all(sample_names(nat) %in% rownames(md_nat)))
md_nat <- md_nat[sample_names(nat), , drop = FALSE]
md_nat$Categoria <- factor(etiquetas2[md_nat$category], levels = etiquetas2)
sample_data(nat) <- sample_data(md_nat)
cat2 <- sample_data(nat)$Categoria; names(cat2) <- sample_names(nat)
rel_nat <- transform_sample_counts(nat, function(x) x / sum(x))

cat("Data_all filtrado:", nsamples(todo_fil), "muestras |", paste(names(table(cat5)), table(cat5), sep = " = ", collapse = ", "), "\n")
cat("Data3:", nsamples(nat), "muestras |", paste(names(table(cat2)), table(cat2), sep = " = ", collapse = ", "), "\n")
