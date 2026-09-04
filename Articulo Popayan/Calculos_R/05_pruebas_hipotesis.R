# ============================================================================
# 05_pruebas_hipotesis.R
# Validación estadística de las diferencias de diversidad entre grupos.
#   - t de Student (varianzas iguales y de Welch) sobre el índice de Shannon
#   - t de Student sobre Shannon restringido al género Fusarium
#   - prueba de Wilcoxon / Mann-Whitney (no paramétrica) sobre Shannon
# Reproduce los resultados de las Figuras 5 y 6 del artículo.
# ============================================================================

if (!exists("fresa_kraken_fil")) source(file.path(PROJECT_ROOT, "Calculos_R", "01_preprocesamiento.R"))

## --- Índice de Shannon por muestra (paquete vegan) ------------------------
OTU <- t(fresa_kraken_fil@otu_table@.Data)          # muestras en filas
SAM <- fresa_kraken_fil@sam_data

Shannon_OTU    <- diversity(OTU, index = "shannon")
Shannon_OTU_df <- data.frame(
  sample  = names(Shannon_OTU),
  value   = Shannon_OTU,
  measure = "Shannon"
)
total <- cbind(Shannon_OTU_df, SAM)
total$Rank <- rank(total$value)

healthy <- total[total$Treatment == "healthy", ]$value
wilted  <- total[total$Treatment == "wilted",  ]$value

## --- (1) t de Student sobre Shannon: varianzas iguales y de Welch ---------
t_var_iguales <- t.test(healthy, wilted, var.equal = TRUE)   # varianzas iguales
t_welch       <- t.test(healthy, wilted, var.equal = FALSE)  # varianzas distintas (Welch)

cat("\n=== t de Student (Shannon, comunidad completa) ===\n")
cat("Varianzas iguales : t =", round(t_var_iguales$statistic, 4),
    "| p =", round(t_var_iguales$p.value, 4), "\n")
cat("Welch             : t =", round(t_welch$statistic, 4),
    "| p =", round(t_welch$p.value, 4), "\n")

## --- (2) t de Student sobre Shannon restringido a Fusarium ----------------
fusarium <- subset_taxa(fresa_kraken_fil, Genus == "Fusarium")
if (ntaxa(fusarium) > 0) {
  OTU_f <- t(fusarium@otu_table@.Data)
  sh_f  <- diversity(OTU_f, index = "shannon")
  df_f  <- cbind(data.frame(value = sh_f), fusarium@sam_data)
  h_f   <- df_f[df_f$Treatment == "healthy", ]$value
  w_f   <- df_f[df_f$Treatment == "wilted",  ]$value

  t_f_iguales <- t.test(h_f, w_f, var.equal = TRUE)
  t_f_welch   <- t.test(h_f, w_f, var.equal = FALSE)

  cat("\n=== t de Student (Shannon, género Fusarium) ===\n")
  cat("Varianzas iguales : t =", round(t_f_iguales$statistic, 4),
      "| p =", round(t_f_iguales$p.value, 4), "\n")
  cat("Welch             : t =", round(t_f_welch$statistic, 4),
      "| p =", round(t_f_welch$p.value, 4), "\n")
}

## --- (3) Prueba de Wilcoxon / Mann-Whitney (no paramétrica) ---------------
w_test <- wilcox.test(x = healthy, y = wilted, alternative = "two.sided",
                      mu = 0, paired = FALSE, conf.int = TRUE, conf.level = 0.95)
cat("\n=== Wilcoxon / Mann-Whitney (Shannon, comunidad completa) ===\n")
print(w_test)

## --- Figura de rangos (equivalente a la Figura 6) -------------------------
ggplot(data = total, aes(x = Rank, y = value)) +
  geom_point(aes(colour = Treatment), size = 3) +
  ylab("Índice de Shannon") + xlab("rango") +
  theme_bw() +
  theme(axis.text.y = element_blank())
save_fig("Wilcoxon_Shannon.png")

## --- Persistir resumen de resultados --------------------------------------
resumen <- data.frame(
  prueba = c("t Student (var. iguales)", "t Welch",
             "Wilcoxon/Mann-Whitney"),
  estadistico = c(unname(t_var_iguales$statistic),
                  unname(t_welch$statistic),
                  unname(w_test$statistic)),
  p_valor = c(t_var_iguales$p.value, t_welch$p.value, w_test$p.value)
)
write.csv(resumen, file.path(OUTPUT_DIR, "pruebas_hipotesis_resumen.csv"), row.names = FALSE)

message("05_pruebas_hipotesis.R completado.")
