## Figura 1.3 - Region de rechazo en una prueba de hipotesis (Img/cap1/pruebat_explicacion.png)
## Regenera la figura en alta resolucion (300 dpi) con etiquetas en espanol.
## La version anterior estaba a 960x540 px.

library("ggplot2")
library("cowplot")

outpath <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap1"

## Curva de la normal estandar (igual que 20230410_PruebasdeHipotesisMedias.R, lineas 163-193)
x <- seq(-4, 4, length = 1000)
y <- dnorm(x, mean = 0, sd = 1)
z_critico <- qnorm(0.95)
data <- data.frame(x, y)

curva <- ggplot(data, aes(x, y)) +
  geom_line(color = "#C6E2FF", linewidth = 1) +
  geom_area(data = subset(data, x < z_critico), aes(x = x, y = y),
            fill = "#C6E2FF", alpha = 0.8) +                      # region de aceptacion
  geom_area(data = subset(data, x > z_critico), aes(x = x, y = y),
            fill = "grey", alpha = 0.5) +                         # region de rechazo
  geom_vline(xintercept = z_critico, color = "royalblue4",
             linetype = "solid", linewidth = 1) +                 # dato observado
  annotate("text", x = z_critico + 0.22, y = 0.15, label = "Dato observado",
           color = "royalblue4", angle = 90, size = 4.5) +
  ## Esta anotacion estaba comentada en el script original (linea 186) y se habia
  ## agregado a mano sobre la imagen; aqui se genera desde el codigo.
  annotate("text", x = 2.45, y = 0.035, label = "Región de\nRechazo",
           color = "grey25", size = 4.5, lineheight = 0.95) +
  labs(title = "Región de rechazo y p-valor en una prueba de hipótesis",
       x = "Valor de la prueba", y = "Densidad") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 13),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

## Recuadros con las tres pruebas (estaban puestos a mano sobre la imagen)
pruebas <- data.frame(
  n = c("1)", "2)", "3)"),
  y = c(3, 2, 1),
  txt = c("Prueba t de Student\nsobre el Índice Shannon\n(Varianzas Iguales y\nDiferentes)",
          "Prueba t de Student\nsobre el Índice Shannon\nen Fusarium (Varianzas\nIguales y Diferentes)",
          "Prueba Mann-Whitney\nsobre las distribuciones\nde los valores de\nShannon.")
)

cajas <- ggplot(pruebas, aes(x = 1, y = y)) +
  geom_label(aes(label = txt), fill = "grey92", colour = "grey15", label.size = 0,
             size = 3.8, lineheight = 1.05, label.padding = unit(0.6, "lines")) +
  geom_text(aes(x = 0.55, label = n), size = 4.5) +
  scale_x_continuous(limits = c(0.4, 1.6)) +
  scale_y_continuous(limits = c(0.5, 3.5)) +
  theme_void()

figura <- plot_grid(curva, cajas, ncol = 2, rel_widths = c(1, 0.42)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("pruebat_explicacion.png", plot = figura, path = outpath,
       width = 30, height = 15, dpi = 300, units = "cm", bg = "white")

cat("Listo 1.3\n")

ggsave("pruebat_explicacion.png", plot = figura, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img",
       width = 30, height = 15, dpi = 300, units = "cm", bg = "white")
