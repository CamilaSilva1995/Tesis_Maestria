## Verificacion numerica de las formulas reportadas en el articulo (ecuaciones 1 a 3 y
## medidas de diversidad beta) contra lo que calculan phyloseq y vegan sobre los datos de Solena.
## Tambien comprueba los estadisticos de las pruebas de hipotesis (t, Welch, F, Mann-Whitney).
suppressMessages({library(phyloseq); library(vegan)})
setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
fk <- import_biom("fresa_kraken.biom")
colnames(fk@otu_table@.Data) <- substr(colnames(fk@otu_table@.Data), 1, 6)
md <- read.csv2("metadata.csv", header = FALSE, row.names = 1, sep = ",")
fk@sam_data <- sample_data(md); colnames(fk@sam_data) <- "Treatment"
fk <- prune_samples(!(sample_names(fk) %in% c("MP2079","MP2080","MP2088","MP2109","MP2137")), fk)

X <- t(fk@otu_table@.Data)          # muestras x taxones
x <- X[1, ]; y <- X[2, ]
p <- x / sum(x); q <- y / sum(y)

chk <- function(nombre, a, b) cat(sprintf("%-30s phyloseq/vegan = %.6f   formula = %.6f   %s\n",
                                          nombre, a, b, ifelse(abs(a - b) < 1e-6, "OK", "DIFIEREN")))

## ---- Diversidad alfa (ecuaciones 1 a 3 del articulo) ----
er <- estimate_richness(fk, measures = c("Observed","Chao1","Shannon","Simpson","InvSimpson"))
S_obs <- sum(x > 0); F1 <- sum(x == 1); F2 <- sum(x == 2)
chk("Chao1 corregido (ec. 1)", er$Chao1[1], S_obs + F1 * (F1 - 1) / (2 * (F2 + 1)))
pp <- p[p > 0]
chk("Shannon (ec. 2)", er$Shannon[1], -sum(pp * log(pp)))
chk("Simpson 1 - sum(p^2) (ec. 3)", er$Simpson[1], 1 - sum(p^2))
chk("InvSimpson (no usado)", er$InvSimpson[1], 1 / sum(p^2))

## ---- Diversidad beta sobre abundancias relativas ----
P <- transform_sample_counts(fk, function(v) v * 100 / sum(v))
D <- function(m) as.matrix(phyloseq::distance(P, method = m))[1, 2]
a <- 100 * p; b <- 100 * q
bc <- sum(abs(a - b)) / sum(a + b)
chk("Bray-Curtis", D("bray"), bc)
chk("Jaccard cuantitativo 2B/(1+B)", D("jaccard"), 2 * bc / (1 + bc))
chk("Euclidiana", D("euclidean"), sqrt(sum((a - b)^2)))
chk("Manhattan", D("manhattan"), sum(abs(a - b)))
M <- (p + q) / 2
KL <- function(u, v) { i <- u > 0; sum(u[i] * log(u[i] / v[i])) }
chk("Jensen-Shannon (divergencia)", D("jsd"), 0.5 * KL(p, M) + 0.5 * KL(q, M))

## ---- Pruebas de hipotesis sobre Shannon, comunidad completa ----
H <- er$Shannon; g <- fk@sam_data$Treatment
n1 <- sum(g == "healthy"); n2 <- sum(g == "wilted")
s1 <- var(H[g == "healthy"]); s2 <- var(H[g == "wilted"])
tw <- t.test(H ~ g); tp <- t.test(H ~ g, var.equal = TRUE)
chk("gl de Welch", tw$parameter, (s1/n1 + s2/n2)^2 / ((s1/n1)^2/(n1-1) + (s2/n2)^2/(n2-1)))
chk("gl t varianzas iguales", tp$parameter, n1 + n2 - 2)
cat(sprintf("p Welch = %.4f (articulo: 0.1501)\n", tw$p.value))
ft <- var.test(H ~ g); chk("F = S1^2/S2^2", ft$statistic, s1 / s2)
w <- wilcox.test(H ~ g, conf.int = TRUE)
r <- rank(H); chk("Mann-Whitney W", w$statistic, sum(r[g == "healthy"]) - n1 * (n1 + 1) / 2)
cat(sprintf("W = %g, p = %.4f, IC95 = (%.4f, %.4f), estimacion = %.4f (articulo: 376, 0.2586, (-0.0124, 0.0561), 0.017)\n",
            w$statistic, w$p.value, w$conf.int[1], w$conf.int[2], w$estimate))
