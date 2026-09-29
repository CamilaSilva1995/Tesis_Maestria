## Tamano del efecto, potencia estadistica y PERMANOVA para la comparacion de Shannon
## entre plantas saludables y no saludables. Respaldo cuantitativo de la seccion de
## conclusiones del articulo (potencia reducida con 35 y 18 muestras).
suppressMessages({library(phyloseq); library(vegan)})
setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
fk <- import_biom("fresa_kraken.biom")
colnames(fk@otu_table@.Data) <- substr(colnames(fk@otu_table@.Data), 1, 6)
md <- read.csv2("metadata.csv", header = FALSE, row.names = 1, sep = ",")
fk@sam_data <- sample_data(md); colnames(fk@sam_data) <- "Treatment"
fk <- prune_samples(!(sample_names(fk) %in% c("MP2079","MP2080","MP2088","MP2109","MP2137")), fk)

H <- diversity(t(fk@otu_table@.Data), "shannon")
g <- fk@sam_data$Treatment
print(table(g)); print(tapply(H, g, mean)); print(tapply(H, g, sd))

m <- tapply(H, g, mean); s <- tapply(H, g, sd); n <- table(g)
sp <- sqrt(((n[1]-1)*s[1]^2 + (n[2]-1)*s[2]^2) / (sum(n) - 2))
cat("d de Cohen:", (m[1] - m[2]) / sp, "\n")
nh <- 2 / (1/n[1] + 1/n[2])   # media armonica de los tamanos de grupo
print(power.t.test(n = nh, delta = m[1] - m[2], sd = sp, sig.level = 0.05))
print(power.t.test(power = 0.8, delta = m[1] - m[2], sd = sp, sig.level = 0.05))

print(var.test(H ~ g))
set.seed(1)
rel <- t(fk@otu_table@.Data) / colSums(fk@otu_table@.Data)
print(adonis2(vegdist(rel, "bray") ~ g, permutations = 999))
