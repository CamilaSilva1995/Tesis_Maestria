"""Construye el cuaderno Explorador_Interactivo_Pruebas_de_hipotesis_con_Plotly.ipynb.

Uso: python3 construir_explorador_interactivo.py   (luego ejecutarlo para incluir las salidas; ver README.md)

Es el tercer cuaderno de la serie. Los dos talleres (construir_cuaderno.py y construir_taller_rapido.py)
enseñan las pruebas; este toma los mismos datos y resultados, con las mismas semillas, y los convierte en
figuras interactivas de Plotly que reúne en una página HTML (explorador_pruebas_de_hipotesis.html).
Tiene la misma organización en episodios y la misma regla de tiempo que los talleres.
"""
import math
import re
import nbformat as nbf

NOMBRE = "Explorador_Interactivo_Pruebas_de_hipotesis_con_Plotly.ipynb"
REPO = "https://colab.research.google.com/github/CamilaSilva1995/Tesis_Maestria/blob/main/Taller_Practico/"
COMPLETO = REPO + "Taller_Practico_Analisis_estadistico_de_Datos_Metagenomicos.ipynb"
RAPIDO = REPO + "Taller_Rapido_Pruebas_de_hipotesis_con_Datos_Metagenomicos.ipynb"

# ------------------------------------------------------------------ regla de tiempo (la de los talleres)
# Explicación = texto que se explica en voz alta + celdas que se ejecutan + figuras que se leen + pausas.
PALABRAS_POR_MINUTO = 120    # ritmo de quien explica un texto técnico y lo va comentando
MINUTOS_POR_CELDA = 2.0      # leer los comentarios de una celda de código, ejecutarla y revisar su salida
MINUTOS_POR_FIGURA = 1.0     # explorar una figura entre todos
MINUTOS_POR_PAUSA = 1.0      # cada pregunta rápida «Para pensar»
REDONDEO = 5                 # la explicación se redondea hacia arriba al múltiplo de 5: deja margen para preguntas

EPISODIOS = []               # un diccionario por episodio, en orden
CIERRE = []                  # celdas que van después del último episodio (referencias y créditos)
n_ejercicios = 0             # los ejercicios se numeran de corrido en todo el cuaderno

def episodio(titulo, preguntas, objetivos):
    "Abre un episodio: lo que se agregue después con md(), code(), pausa()... queda dentro de él."
    EPISODIOS.append({"titulo": titulo, "preguntas": preguntas, "objetivos": objetivos,
                      "celdas": [], "palabras": 0, "codigo": 0, "figuras": 0, "pausas": 0, "ejercicios": []})

def palabras(texto):
    "Cuenta las palabras de un texto (sin números, símbolos ni fórmulas)."
    return len(re.findall(r"[A-Za-zÁÉÍÓÚÜÑáéíóúüñ]{2,}", texto))

def md(s):
    "Agrega una celda de texto (Markdown) al episodio abierto y suma sus palabras al tiempo de explicación."
    ep = EPISODIOS[-1]
    ep["celdas"].append(nbf.v4.new_markdown_cell(s.strip()))
    ep["palabras"] += palabras(s)

def code(s, figuras=0):
    "Agrega una celda de código; figuras = cuántas figuras dibuja (cuentan para el tiempo)."
    ep = EPISODIOS[-1]
    ep["celdas"].append(nbf.v4.new_code_cell(s.strip()))
    ep["codigo"] += 1
    ep["figuras"] += figuras

def pausa(texto):
    "Agrega una pregunta rápida para explorar la figura que se acaba de dibujar."
    ep = EPISODIOS[-1]
    ep["celdas"].append(nbf.v4.new_markdown_cell(f"> 💬 **Para pensar (1 min).** {texto.strip()}"))
    ep["pausas"] += 1

def ejercicio(titulo, minutos, enunciado, solucion_md, solucion_code=None):
    "Agrega un ejercicio: el enunciado y la solución (texto y, si hay, código) dentro de un desplegable."
    global n_ejercicios
    n_ejercicios += 1
    s = (f"### ✏️ Ejercicio {n_ejercicios}: {titulo}\n\n⏱️ *{minutos} min*\n\n{enunciado.strip()}\n\n<details>\n"
         f"<summary><b>👉 Ver solución</b> (haz clic para desplegar)</summary>\n\n{solucion_md.strip()}\n")
    if solucion_code:   # la solución en código es opcional
        s += f"\n```python\n{solucion_code.strip()}\n```\n"
    s += "\n</details>"
    ep = EPISODIOS[-1]
    ep["celdas"].append(nbf.v4.new_markdown_cell(s))
    ep["ejercicios"].append(minutos)

def puntos_clave(puntos):
    "Cierra el episodio abierto con sus puntos clave."
    ep = EPISODIOS[-1]
    ep["celdas"].append(nbf.v4.new_markdown_cell("> **🔑 Puntos clave**\n" + "\n".join(f"> - {p}" for p in puntos)))

def lista(elementos):
    "Lista con viñetas dentro de un bloque de cita."
    return "\n".join(f"> - {e}" for e in elementos)

# ================================================================== EPISODIO 1
episodio("Los datos y las pruebas, listos para dibujar",
         ["¿Qué necesito tener calculado antes de hacer una figura interactiva?"],
         ["Cargar los datos reales y generar los simulados, igual que en los talleres.",
          "Dejar en variables los resultados de las pruebas que se van a mostrar."])

md(r"""
Una figura interactiva no calcula nada: solo muestra lo que ya está calculado. Por eso este episodio deja listos los mismos datos y las mismas pruebas de los talleres. Son tres celdas largas que ya conoces: ejecútalas y sigue.

La primera prepara las herramientas. La única novedad es **Plotly**, la librería que dibuja las figuras interactivas.
""")
code(r"""
import html                                             # para incrustar la página dentro del cuaderno (Episodio 5)
import numpy as np                                      # arreglos numéricos y números aleatorios
import pandas as pd                                     # tablas de datos (DataFrame)
from scipy import stats                                 # pruebas estadísticas
from scipy.spatial.distance import pdist, squareform    # distancias entre muestras
import plotly                                           # figuras interactivas
import plotly.graph_objects as go                       # las piezas de una figura: trazas (go.Box, go.Scatter...) y diseño
import plotly.io as pio                                 # cómo se muestran y se guardan las figuras
from plotly.subplots import make_subplots               # varias gráficas dentro de una misma figura

try:
    import google.colab                                 # solo existe en Colab: allí Plotly ya sabe cómo mostrar las figuras
except ImportError:
    # fuera de Colab, la biblioteca de Plotly se carga desde internet para que el cuaderno guardado no pese varios megas
    pio.renderers.default = "plotly_mimetype+notebook_connected"

# Carpeta del repositorio de GitHub donde están los archivos CSV de la clase
URL = "https://raw.githubusercontent.com/CamilaSilva1995/Tesis_Maestria/main/Taller_Practico/datos/"

def cargar(nombre):
    "Lee un CSV desde GitHub; si no hay internet, desde la carpeta local datos/."
    try:
        return pd.read_csv(URL + nombre)        # primer intento: descargarlo de GitHub
    except Exception:
        return pd.read_csv("datos/" + nombre)   # si falla, leer la copia local

def formato_p(p):
    "Escribe un valor p con cuatro decimales; los muy pequeños, en notación científica."
    return f"{p:.4f}" if p >= 1e-4 else f"{p:.1e}"

def p_en_texto(p):
    "Escribe un valor p para un párrafo: tres decimales, o «p < 0.001» si es más pequeño."
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"

# Aspecto de las figuras y de la página: "noche" (fondo azul noche) o "claro". Cambia el tema y vuelve a ejecutar todo
TEMA = "noche"
PALETAS = {
    "noche": dict(fondo="#0B1026", panel="#111A3B", rejilla="#27335F", texto="#E8ECF8", suave="#9AA6CF", enlace="#7DD3FC",
                  neutro="#4A5A8F", acento="#FFD166", alerta="#FF4D8D", aviso="#FFB703", bien="#7CFFB2", saludable="#F8766D", no_saludable="#00BFC4",
                  control="#B9C4EC", control_texto="#0B1026"),
    "claro": dict(fondo="#F5F7FB", panel="#FFFFFF", rejilla="#E3E7F0", texto="#1F2937", suave="#6B7280", enlace="#0B6E99",
                  neutro="#B8BDCB", acento="#E09F00", alerta="#D62768", aviso="#F08C00", bien="#4C9F70", saludable="#F8766D", no_saludable="#00BFC4",
                  control="#E8ECF6", control_texto="#1F2937"),
}
P = PALETAS[TEMA]                                       # los colores del tema elegido
# Un color fijo por grupo, el mismo de los talleres y de la tesis (coral y turquesa) en los dos temas
COLORES = {"Saludable": P["saludable"], "No saludable": P["no_saludable"]}
GRUPOS = ["Saludable", "No saludable"]                  # el orden de los grupos en las figuras

def con_alfa(color, alfa):
    "Convierte un color #RRGGBB en uno con transparencia, para rellenos suaves."
    r, g, b = (int(color[i:i + 2], 16) for i in (1, 3, 5))
    return f"rgba({r},{g},{b},{alfa})"

def estilo(fig, titulo, alto=480):
    "Da a una figura el aspecto común: colores del tema, título centrado, altura y leyenda horizontal abajo."
    fig.update_layout(
        template="plotly_dark" if TEMA == "noche" else "plotly_white",
        paper_bgcolor=P["panel"], plot_bgcolor=P["panel"], font=dict(color=P["texto"], size=13),
        title=dict(text=titulo, x=0.5, xanchor="center", font=dict(size=17)), height=alto,
        margin=dict(l=70, r=30, t=80, b=60),
        legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.18, yanchor="top", bgcolor="rgba(0,0,0,0)"),
        hoverlabel=dict(bgcolor=P["fondo"], bordercolor=P["acento"], font=dict(color=P["texto"], size=13)))
    fig.update_xaxes(gridcolor=P["rejilla"], zerolinecolor=P["rejilla"], linecolor=P["rejilla"])
    fig.update_yaxes(gridcolor=P["rejilla"], zerolinecolor=P["rejilla"], linecolor=P["rejilla"])
    return fig

def paneles_relativos(fig):
    "Ancla el título de cada panel de make_subplots a su propio panel (y no a la figura): así los paneles se pueden apilar en pantallas angostas."
    for i, anotacion in enumerate(fig.layout.annotations[:2]):
        eje = "" if i == 0 else "2"
        anotacion.update(xref=f"x{eje} domain", x=0.5, yref=f"y{eje} domain", y=1.0, yanchor="bottom")

def violin(valores, g, unidad, visible=True, leyenda=True):
    "Medio violín con su caja y, al lado, todos los puntos; al pasar el cursor por un punto se lee qué muestra es."
    return go.Violin(
        y=valores, name=g, text=list(valores.index),         # el identificador de cada muestra
        side="positive", width=1.1,                          # solo la mitad derecha del violín: la forma de la distribución
        points="all", pointpos=-0.55, jitter=0.5,            # todos los puntos, a la izquierda del violín
        box_visible=True, meanline_visible=True,             # la caja (cuartiles) y la media, dentro del violín
        line_color=COLORES[g], fillcolor=con_alfa(COLORES[g], 0.30),
        marker=dict(size=8, opacity=0.9, line=dict(width=0.6, color=P["panel"])),
        hoveron="points", hovertemplate="%{text}<br>" + unidad + "<extra></extra>",
        legendgroup=g, showlegend=leyenda, visible=visible)

print("Listo. plotly", plotly.__version__, "| numpy", np.__version__, "| pandas", pd.__version__)
""")
code(r"""
# ---------- Datos reales: tres tablas unidas por el identificador de la muestra
meta = cargar("metadatos.csv").set_index("muestra")          # grupo y lecturas clasificadas de cada muestra
alfa = cargar("diversidad_alfa.csv").set_index("muestra")    # índices de diversidad alfa (Shannon, Chao1...)
conteos = cargar("conteos_por_genero.csv")                   # una fila por género y una columna por muestra

excluir = ["MP2079", "MP2080", "MP2088", "MP2109", "MP2137"] # filtro de calidad: menos de 25 millones de lecturas
datos = meta.join(alfa, how="inner").drop(index=excluir)     # une por identificador y quita esas cinco
muestras = list(datos.index)                                 # las 53 muestras reales
generos = conteos.set_index(["reino", "filo", "genero"])[muestras]
assert list(generos.columns) == muestras                     # comprobación: las tablas quedaron emparejadas

# ---------- Datos simulados: el caso ideal, con la misma semilla de los talleres
rng_sim = np.random.default_rng(2026)
n_sim = 30                                                   # muestras por grupo: diseño balanceado
media_sal, media_no, desviacion = 7.48, 7.42, 0.05           # Shannon normal, misma desviación en los dos grupos
ids = [f"S{i:02d}" for i in range(1, n_sim + 1)] + [f"N{i:02d}" for i in range(1, n_sim + 1)]
datos_sim = pd.DataFrame({
    "grupo": ["Saludable"] * n_sim + ["No saludable"] * n_sim,
    "shannon": np.concatenate([rng_sim.normal(media_sal, desviacion, n_sim),
                               rng_sim.normal(media_no, desviacion, n_sim)]),
}, index=pd.Index(ids, name="muestra"))

n_generos, lecturas_sim, variabilidad = 40, 100_000, 0.25
nombres = [f"G{i:02d}" for i in range(1, n_generos + 1)]        # G01, G02, ..., G40
patogeno = "G05"                                                 # el género que hacemos más abundante
comp_sal = pd.Series(1 / np.arange(1, n_generos + 1), index=nombres)
comp_sal = comp_sal / comp_sal.sum()                             # saludables: pocos géneros abundantes, muchos escasos
comp_no = comp_sal.copy()
comp_no[patogeno] = 4 * comp_no[patogeno]                        # no saludables: G05 cuatro veces más abundante
comp_no = comp_no / comp_no.sum()

def simular_muestras(composicion, n_muestras):
    "Simula n_muestras columnas de conteos alrededor de una composición promedio."
    columnas = []
    for _ in range(n_muestras):
        p = composicion.to_numpy() * np.exp(rng_sim.normal(0, variabilidad, len(composicion)))
        p = p / p.sum()
        columnas.append(rng_sim.multinomial(lecturas_sim, p))    # se "secuencian" 100 000 lecturas
    return np.array(columnas).T                                  # filas = géneros, columnas = muestras

generos_sim = pd.DataFrame(np.hstack([simular_muestras(comp_sal, n_sim), simular_muestras(comp_no, n_sim)]),
                           index=pd.Index(nombres, name="genero"), columns=datos_sim.index)
assert list(generos_sim.columns) == list(datos_sim.index)        # emparejadas por identificador

print("Reales:", len(datos), "muestras y", len(generos), "géneros | Simulados:", len(datos_sim), "muestras y", len(generos_sim), "géneros")
""")
code(r"""
# ---------- Las pruebas de los talleres, en funciones
def prueba_permutacion(x, y, n_perm=5000, semilla=2026):
    "Valor p por permutaciones: devuelve la diferencia de medias observada, las barajadas y el valor p."
    rng = np.random.default_rng(semilla)
    observada = x.mean() - y.mean()
    todos = np.concatenate([x, y])                       # todos los valores, sin etiqueta de grupo
    barajadas = []
    for _ in range(n_perm):
        rng.shuffle(todos)                               # barajar etiquetas
        barajadas.append(todos[:len(x)].mean() - todos[len(x):].mean())
    barajadas = np.array(barajadas)
    return observada, barajadas, np.mean(np.abs(barajadas) >= abs(observada))

def pseudo_F(D, etiquetas):
    "Pseudo-F y R2 del PERMANOVA a partir de la matriz de distancias y el grupo de cada muestra."
    n = len(etiquetas)
    grupos = np.unique(etiquetas)
    SS_total = (D ** 2).sum() / (2 * n)                  # variación total
    SS_dentro = sum((D[np.ix_(etiquetas == g, etiquetas == g)] ** 2).sum() / (2 * (etiquetas == g).sum())
                    for g in grupos)                     # variación dentro de los grupos
    SS_entre = SS_total - SS_dentro                      # variación entre grupos
    return (SS_entre / (len(grupos) - 1)) / (SS_dentro / (n - len(grupos))), SS_entre / SS_total

def permanova(D, etiquetas, n_perm=999, semilla=2026):
    "PERMANOVA: devuelve el pseudo-F observado, el R2 y el valor p por permutaciones."
    rng = np.random.default_rng(semilla)
    F_obs, R2 = pseudo_F(D, etiquetas)
    F_perm = np.array([pseudo_F(D, rng.permutation(etiquetas))[0] for _ in range(n_perm)])
    return F_obs, R2, (np.sum(F_perm >= F_obs) + 1) / (n_perm + 1)

def pcoa(D):
    "Coordenadas principales a partir de una matriz de distancias; las primeras columnas recogen más variación."
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n        # matriz de centrado
    B = J @ (-0.5 * D ** 2) @ J                # doble centrado de -1/2 por las distancias al cuadrado
    val, vec = np.linalg.eigh(B)               # valores y vectores propios, de menor a mayor
    val, vec = val[::-1], vec[:, ::-1]         # se invierte el orden: primero los ejes más importantes
    keep = val > 1e-10                         # se descartan los valores propios nulos o negativos
    return vec[:, keep] * np.sqrt(val[keep])

def permdisp(D, etiquetas, n_perm=999, semilla=2026):
    "PERMDISP: devuelve la distancia de cada muestra al centroide de su grupo, el estadístico F y el valor p."
    coords = pcoa(D)
    dist_centro = np.empty(len(etiquetas))
    for g in np.unique(etiquetas):
        m = etiquetas == g
        dist_centro[m] = np.sqrt(((coords[m] - coords[m].mean(axis=0)) ** 2).sum(axis=1))
    F_obs = stats.f_oneway(*[dist_centro[etiquetas == g] for g in np.unique(etiquetas)]).statistic
    rng = np.random.default_rng(semilla)
    F_perm = [stats.f_oneway(*[dist_centro[rng.permutation(etiquetas) == g] for g in np.unique(etiquetas)]).statistic
              for _ in range(n_perm)]
    return dist_centro, F_obs, (np.sum(np.array(F_perm) >= F_obs) + 1) / (n_perm + 1)

# ---------- Lo que van a mostrar las figuras
# Shannon de cada grupo (x = saludables, y = no saludables) en cada caso
x_sim = datos_sim.loc[datos_sim.grupo == "Saludable", "shannon"].to_numpy()
y_sim = datos_sim.loc[datos_sim.grupo == "No saludable", "shannon"].to_numpy()
x = datos.loc[datos.grupo == "Saludable", "shannon"].to_numpy()
y = datos.loc[datos.grupo == "No saludable", "shannon"].to_numpy()

# Porcentaje de lecturas de cada género en cada muestra
pct_sim = 100 * generos_sim / generos_sim.sum(axis=0)            # simulados: sobre el total de la muestra
por_genero = generos.groupby(level="genero").sum()               # reales: una fila por género (suma si está en varios filos)
pct = 100 * por_genero / datos.lecturas_clasificadas             # sobre TODAS las lecturas clasificadas de la muestra
# filo de cada género, para mostrarlo al pasar el cursor
filo_de = generos.reset_index().drop_duplicates("genero").set_index("genero").filo.fillna("sin filo asignado")

# Distancias de Bray-Curtis entre muestras y grupo de cada muestra
D_sim = squareform(pdist((generos_sim / generos_sim.sum(axis=0)).T.to_numpy(), metric="braycurtis"))
D = squareform(pdist((generos / generos.sum(axis=0)).T.to_numpy(), metric="braycurtis"))
etiquetas_sim, etiquetas = datos_sim.grupo.to_numpy(), datos.grupo.to_numpy()
print("Pruebas listas. Matrices de distancias:", D_sim.shape, "simulados |", D.shape, "reales")
""")
puntos_clave([
    "Una figura interactiva muestra resultados ya calculados: primero el análisis, después el dibujo.",
    "Los datos y las semillas son los de los talleres, así que los números coinciden.",
])

# ================================================================== EPISODIO 2
episodio("De una figura estática a una interactiva",
         ["¿Qué gano con una figura interactiva y cómo se construye una con Plotly?"],
         ["Construir una figura de Plotly a partir de trazas y de un diseño.",
          "Mostrar información de cada muestra al pasar el cursor.",
          "Agregar un deslizador que cambia lo que se ve."])

md(r"""
**Qué hacemos.** En los talleres las figuras eran imágenes: se miran, pero no se pueden interrogar. Una figura de Plotly responde: al pasar el cursor dice qué muestra es cada punto, se puede ampliar una zona y se pueden ocultar grupos desde la leyenda.

**Cómo se arma.** Toda figura de Plotly tiene dos partes:

- **Trazas** (`go.Box`, `go.Scatter`, `go.Bar`...): cada una es un conjunto de datos con su forma de dibujarse.
- **Diseño** (`layout`): títulos, ejes, leyenda y controles como menús y deslizadores.

La interactividad viaja dentro de la figura, sin servidor ni Python: por eso se puede guardar en un archivo HTML y abrirlo en cualquier navegador.

### 2.1 Una primera figura: la comunidad de un vistazo

Empezamos por la figura más conocida de un microbioma: la composición de cada muestra, en barras apiladas. Cada barra es una muestra y cada color, un filo. Sirve para ver cómo se construye una figura con una traza por serie.
""")
code(r"""
# % de lecturas de cada filo en cada muestra real (los géneros sin filo asignado no entran)
por_filo = generos.groupby(level="filo").sum()                   # una fila por filo
pct_filo = 100 * por_filo / por_filo.sum(axis=0)                 # cada muestra suma 100
principales = pct_filo.mean(axis=1).sort_values(ascending=False).index[:8]   # los ocho filos más abundantes
tabla_filos = pct_filo.loc[principales].copy()
tabla_filos.loc["Otros filos"] = (100 - tabla_filos.sum(axis=0)).clip(lower=0)   # el resto, en una sola categoría

# Muestras ordenadas: primero las saludables y después las no saludables
orden = list(datos.index[datos.grupo == "Saludable"]) + list(datos.index[datos.grupo == "No saludable"])
n_sal = int((datos.grupo == "Saludable").sum())
# colores vivos para los filos, distintos de los dos colores de los grupos
vivos = ["#4D96FF", "#B388FF", "#FFD166", "#7CFFB2", "#FF4D8D", "#C6F432", "#FDBA74", "#F9A8D4", "#7C88B5"]

fig_filos = go.Figure()
for filo, color in zip(tabla_filos.index, vivos):
    # una traza por filo: Plotly las apila porque el diseño dice barmode="stack"
    fig_filos.add_trace(go.Bar(x=orden, y=tabla_filos.loc[filo, orden], name=filo, marker_color=color,
                               hovertemplate="%{x}<br>" + filo + ": %{y:.1f} %<extra></extra>"))
estilo(fig_filos, "Datos reales: composición por filo de las 53 muestras", alto=520)
fig_filos.update_layout(barmode="stack", bargap=0.06, yaxis_title="% de las lecturas con filo asignado",
                        xaxis=dict(showticklabels=False, title="Cada barra es una muestra"), margin=dict(t=100),
                        legend_traceorder="normal")                  # la leyenda en el mismo orden en que se apilan las barras
# la muestra con más Firmicutes se sale del patrón; se señala porque volverá a aparecer en las demás figuras
atipica = pct_filo.loc["Firmicutes"].idxmax()
fig_filos.add_annotation(x=atipica, y=100, ax=-80, ay=55, text=f"{atipica}: {pct_filo.loc['Firmicutes', atipica]:.0f} % de Firmicutes",
                         showarrow=True, arrowcolor=P["texto"], font=dict(color=P["texto"], size=12),
                         bgcolor=P["fondo"], bordercolor=P["texto"], borderpad=3)
# una línea separa los dos grupos y un rótulo en cada extremo los nombra
fig_filos.add_vline(x=n_sal - 0.5, line_color=P["texto"], line_dash="dot")
for g, lado, x_rotulo in zip(GRUPOS, ["left", "right"], [0, 1]):
    fig_filos.add_annotation(x=x_rotulo, xref="paper", xanchor=lado, y=1.04, yref="paper", yanchor="bottom", showarrow=False,
                             font=dict(color=COLORES[g], size=14), text=f"<b>{g}s ({int((datos.grupo == g).sum())})</b>")
fig_filos.show()
""", figuras=1)
pausa(r"""
Haz clic en un filo de la leyenda para ocultarlo y doble clic para dejarlo solo. ¿Ves alguna diferencia evidente entre las plantas saludables y las no saludables?
""")
md(r"""
A simple vista las muestras se parecen mucho: los mismos filos dominan en todas. Una figura así describe, pero no decide: para saber si los grupos difieren hacen falta las pruebas.

### 2.2 Violines con puntos que se identifican

Un violín muestra la forma de la distribución; dentro lleva la caja con los cuartiles y, al lado, el punto de cada muestra.
""")
code(r"""
# Una fila y dos columnas: un panel por caso
fig_cajas = make_subplots(rows=1, cols=2, subplot_titles=["Datos simulados (caso ideal)", "Datos reales (fresa)"])
for columna, tabla in enumerate([datos_sim, datos], start=1):
    for g in GRUPOS:
        valores = tabla.loc[tabla.grupo == g, "shannon"]     # Shannon de las muestras de ese grupo
        # leyenda solo en el primer panel: una sola entrada por grupo
        fig_cajas.add_trace(violin(valores, g, "Shannon = %{y:.3f}", leyenda=(columna == 1)), row=1, col=columna)

# Debajo de cada panel: la desviación de cada grupo con la prueba F (¿varían lo mismo?), y el resultado de Welch
# (no supone varianzas iguales) y de Mann-Whitney (no supone normalidad)
for columna, (a, b) in enumerate([(x_sim, y_sim), (x, y)], start=1):
    F = a.var(ddof=1) / b.var(ddof=1)                                       # cociente de varianzas
    cola = stats.f.cdf(F, len(a) - 1, len(b) - 1)
    p_F = 2 * min(cola, 1 - cola)                                           # prueba F bilateral
    welch = stats.ttest_ind(a, b, equal_var=False).pvalue
    mw = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue
    fig_cajas.add_annotation(text=(f"DE = {a.std(ddof=1):.3f} y {b.std(ddof=1):.3f} (prueba F p = {formato_p(p_F)})<br>"
                                   f"Welch p = {formato_p(welch)} · Mann-Whitney p = {formato_p(mw)}"),
                             xref="x domain" if columna == 1 else "x2 domain", x=0.5,          # relativa a su panel
                             yref="y domain" if columna == 1 else "y2 domain", y=-0.13, yanchor="top",
                             showarrow=False, font=dict(color=P["acento"]))
fig_cajas.update_yaxes(title_text="Índice de Shannon", row=1, col=1)
estilo(fig_cajas, "Diversidad por grupo: cada punto es una muestra", alto=500)
fig_cajas.update_layout(showlegend=False, margin=dict(b=90))   # los ejes ya nombran los grupos; abajo van las anotaciones
paneles_relativos(fig_cajas)
fig_cajas.show()
""", figuras=1)
pausa(r"""
En el panel de los datos reales, busca la planta no saludable con el Shannon más bajo. ¿Cuál es su identificador? Es la muestra que hacía fallar los supuestos de la prueba t.
""")
md(r"""
Es la muestra `MD2086`. Volverá a aparecer en la ordenación del Episodio 4.

### 2.3 Un deslizador para ver cómo cambia el valor p

Un control no vuelve a calcular nada: solo decide **qué trazas se ven**. Por eso se calculan todas antes y el deslizador enciende unas y apaga otras. Aquí cada posición es una prueba de permutaciones con distinto número de muestras simuladas; la última posición son los datos reales. Las barras resaltadas son las barajadas tan extremas como lo observado: **su proporción es el valor p**.
""")
code(r"""
# Cada caso del deslizador: (etiqueta corta, saludables, no saludables, descripción)
casos = [(f"{n} + {n}", x_sim[:n], y_sim[:n], f"Datos simulados, {n} muestras por grupo") for n in (5, 10, 15, 20, 25, 30)]
casos.append(("reales", x, y, "Datos reales, 35 y 18 muestras"))
inicial = len(casos) - 2                                 # se abre en los datos simulados completos (30 + 30)

fig_perm = go.Figure()
titulos, limite = [], 0
for i, (etiqueta, a, b, descripcion) in enumerate(casos):
    observada, barajadas, p = prueba_permutacion(a, b)
    # el histograma se calcula aquí y se dibuja como barras: la página pesa mucho menos que con 5 000 valores
    alturas, bordes = np.histogram(barajadas, bins=40)
    centros = (bordes[:-1] + bordes[1:]) / 2
    extremas = np.abs(centros) >= abs(observada)         # las barras de las colas: tan extremas como lo observado
    # dos trazas de barras, para que cada una tenga su entrada en la leyenda: el centro y, en color de alerta, las colas
    for nombre, mascara, color in (("lo que produce el azar (H0)", ~extremas, P["neutro"]),
                                   ("tan extremas como lo observado: su proporción es el valor p", extremas, P["alerta"])):
        fig_perm.add_trace(go.Bar(x=centros[mascara], y=alturas[mascara], width=bordes[1] - bordes[0], marker_color=color,
                                  name=nombre, visible=(i == inicial),
                                  hovertemplate="diferencia ≈ %{x:.3f}<br>%{y} barajadas<extra></extra>"))
    # la diferencia observada y su opuesta (prueba bilateral), como líneas verticales
    for signo, raya in ((1, "solid"), (-1, "dash")):
        fig_perm.add_trace(go.Scatter(x=[signo * observada] * 2, y=[0, alturas.max()], mode="lines",
                                      line=dict(color=P["acento"], width=3, dash=raya), name="diferencia observada",
                                      visible=(i == inicial), showlegend=(signo == 1), hoverinfo="skip"))
    titulos.append(f"{descripcion}: diferencia = {observada:.3f}, p = {formato_p(p)}")
    limite = max(limite, np.abs(barajadas).max(), abs(observada))

# Cada paso del deslizador enciende las cuatro trazas de su caso (dos de barras y dos líneas) y cambia el título
pasos = [dict(method="update", label=etiqueta,
              args=[{"visible": [k // 4 == i for k in range(4 * len(casos))]}, {"title.text": titulos[i]}])
         for i, (etiqueta, *_) in enumerate(casos)]
estilo(fig_perm, titulos[inicial], alto=520)
fig_perm.update_layout(sliders=[dict(active=inicial, steps=pasos, currentvalue=dict(prefix="Muestras: "), pad=dict(t=60),
                                     bgcolor=P["rejilla"], activebgcolor=P["acento"], bordercolor=P["rejilla"], font=dict(color=P["texto"]))],
                       bargap=0, barmode="overlay", yaxis_title="Número de barajadas", legend=dict(y=1.0, yanchor="bottom"),
                       xaxis=dict(title="Diferencia de medias de Shannon con las etiquetas barajadas",
                                  range=[-1.1 * limite, 1.1 * limite]))   # mismo eje en todos los pasos, para poder comparar
fig_perm.show()
""", figuras=1)
md(r"""
Mueve el deslizador de izquierda a derecha. La diferencia verdadera es siempre la misma, pero con más muestras el histograma del azar se **estrecha**, la diferencia observada va quedando en la cola y las barras resaltadas son cada vez menos: el valor p baja. La última posición muestra los datos reales, con su histograma asimétrico y un p ≈ 0.04 que, como se ve en el taller, no es de fiar: barajar etiquetas supone que bajo $H_0$ los dos grupos tienen la misma distribución, varianza incluida.

### 2.4 El mismo deslizador, otra pregunta: ¿qué prueba elegir?

En los datos reales las pruebas no coinciden: permutación p ≈ 0.04, t de Student 0.058, Welch 0.150 y Mann-Whitney 0.259. ¿A cuál creer? A la que cumple sus supuestos. El Laboratorio 2 del taller completo lo muestra simulando miles de estudios **sin ninguna diferencia verdadera** y contando cuántas veces cada prueba «encuentra» una: lo correcto es cerca del 5 %. Aquí ese laboratorio se vuelve un deslizador con cuatro diseños del estudio; cada posición enciende las barras de su diseño.
""")
code(r"""
# Cuatro diseños: (etiqueta del deslizador, muestras saludables, no saludables, DE saludables, DE no saludables, descripción)
disenos = [("como los reales", 35, 18, 0.06, 0.17, "35 y 18 muestras, DE 0.06 y 0.17 (como los datos reales: el grupo pequeño es el más variable)"),
           ("DE iguales", 35, 18, 0.06, 0.06, "35 y 18 muestras, DE 0.06 y 0.06 (los dos grupos varían lo mismo)"),
           ("n iguales", 35, 35, 0.06, 0.17, "35 y 35 muestras, DE 0.06 y 0.17 (grupos del mismo tamaño)"),
           ("DE al revés", 35, 18, 0.17, 0.06, "35 y 18 muestras, DE 0.17 y 0.06 (el grupo grande es el más variable)")]
n_estudios = 10_000
pruebas = ["t de Student", "t de Welch", "Mann-Whitney"]

fig_alarmas = go.Figure()
titulos, alarmas = [], {}
for i, (etiqueta, n_a, n_b, de_a, de_b, descripcion) in enumerate(disenos):
    rng_lab = np.random.default_rng(2026)
    A = rng_lab.normal(7.45, de_a, (n_a, n_estudios))       # una columna por estudio; la media es la misma en los dos grupos
    B = rng_lab.normal(7.45, de_b, (n_b, n_estudios))
    # cada prueba se aplica a los 10 000 estudios a la vez (una columna por estudio) y se cuenta en cuántos da p < 0.05
    porcentajes = [100 * np.mean(stats.ttest_ind(A, B, equal_var=True).pvalue < 0.05),
                   100 * np.mean(stats.ttest_ind(A, B, equal_var=False).pvalue < 0.05),
                   100 * np.mean(stats.mannwhitneyu(A, B, alternative="two-sided").pvalue < 0.05)]
    alarmas[etiqueta] = porcentajes
    # una traza de barras por diseño; en color de alerta las pruebas que se alejan más de dos puntos del 5 %
    fig_alarmas.add_trace(go.Bar(x=pruebas, y=porcentajes, text=[f"{v:.1f} %" for v in porcentajes], textposition="outside",
                                 marker_color=[P["alerta"] if abs(v - 5) > 2 else P["bien"] for v in porcentajes],
                                 textfont=dict(color=P["texto"], size=14), visible=(i == 0), showlegend=False,
                                 hovertemplate="%{x}: ve una diferencia que no existe en el %{y:.1f} % de los estudios<extra></extra>"))
    titulos.append(f"Sin diferencia verdadera: {descripcion}")
# la línea del 5 % se dibuja como una traza más (siempre visible), para que tenga su entrada en la leyenda
fig_alarmas.add_trace(go.Scatter(x=[-0.5, len(pruebas) - 0.5], y=[5, 5], mode="lines", line=dict(color=P["acento"], dash="dot", width=2),
                                 name="5 %: lo esperado con α = 0.05", hoverinfo="skip"))

# Cada paso enciende la traza de barras de su diseño; la última traza (la línea del 5 %) queda siempre encendida
pasos = [dict(method="update", label=etiqueta,
              args=[{"visible": [k == i for k in range(len(disenos))] + [True]}, {"title.text": titulos[i]}])
         for i, (etiqueta, *_) in enumerate(disenos)]
estilo(fig_alarmas, titulos[0], alto=480)
fig_alarmas.update_layout(sliders=[dict(active=0, steps=pasos, currentvalue=dict(prefix="Diseño: "), pad=dict(t=50),
                                        bgcolor=P["rejilla"], activebgcolor=P["acento"], bordercolor=P["rejilla"], font=dict(color=P["texto"]))],
                          yaxis=dict(title="% de estudios en que la prueba ve una diferencia", range=[0, 1.3 * max(max(v) for v in alarmas.values())]),
                          xaxis=dict(title="Prueba (la verdad: los dos grupos tienen la misma media)", range=[-0.5, len(pruebas) - 0.5]),
                          legend=dict(y=1.0, yanchor="bottom"))
fig_alarmas.show()
for etiqueta, porcentajes in alarmas.items():
    print(f"{etiqueta:16s}", " | ".join(f"{prueba}: {v:4.1f} %" for prueba, v in zip(pruebas, porcentajes)))
""", figuras=1)
pausa(r"""
¿Qué prueba se queda cerca del 5 % en los cuatro diseños? ¿Con qué diseño se equivocan poco las tres?
""")
md(r"""
Con el diseño de los datos reales, la t de Student y Mann-Whitney ven diferencias que no existen bastante más del 5 % de las veces; Welch se mantiene cerca del 5 % en los cuatro diseños. Con grupos del mismo tamaño el problema casi desaparece: una buena razón para diseñar estudios balanceados. Por eso, en los datos reales, el resultado que vale es el de Welch: p = 0.150, no hay evidencia de que la diversidad promedio difiera.
""")
ejercicio("Un paso más en el deslizador", 4,
r"""
Agrega al deslizador un paso con solo **3 muestras por grupo** y vuelve a ejecutar la celda. ¿Qué valor p da? ¿Por qué el histograma tiene tan pocas barras?
""",
r"""
Basta con agregar el 3 a la lista de tamaños: `for n in (3, 5, 10, 15, 20, 25, 30)`. Con tres muestras por grupo solo hay 20 formas distintas de repartir las seis muestras, así que las 5 000 barajadas repiten una y otra vez las mismas pocas diferencias y el valor p no puede ser pequeño (aquí da cerca de 0.89). Con tan pocos datos ninguna prueba puede detectar la diferencia, aunque exista.
""",
r"""
# solo cambia esta línea de la celda del deslizador; lo demás queda igual
casos = [(f"{n} + {n}", x_sim[:n], y_sim[:n], f"Datos simulados, {n} muestras por grupo") for n in (3, 5, 10, 15, 20, 25, 30)]
""")
puntos_clave([
    "Una figura de Plotly se arma con trazas (los datos) y un diseño (títulos, ejes y controles).",
    "`hovertemplate` decide qué se lee al pasar el cursor; así cada punto se puede identificar.",
    "Un deslizador o un menú no recalcula: enciende y apaga trazas que ya están calculadas.",
])

# ================================================================== EPISODIO 3
episodio("Muchos géneros a la vez",
         ["¿Cómo exploro la abundancia de muchos géneros sin hacer una figura para cada uno?"],
         ["Agregar un menú desplegable que cambia el género que se muestra.",
          "Resumir cientos de pruebas en una sola figura.",
          "Ver el efecto de las comparaciones múltiples y de su corrección."])

md(r"""
**Qué hacemos.** En el taller comparamos *Fusarium* y, en un ejercicio, tres géneros más. Con una figura interactiva se pueden dejar muchos géneros en un menú, y hasta mirar los 1 795 a la vez.

**Por qué importa.** Cuantas más pruebas se hacen, más resultados «significativos» aparecen por puro azar. Verlo en una figura es la mejor forma de no olvidarlo.

### 3.1 Un menú para elegir el género
""")
code(r"""
# Géneros del menú: los del taller y, además, los doce con más lecturas en los datos reales
del_taller = ["Fusarium", "Phytophthora", "Ralstonia"]
mas_abundantes = por_genero.sum(axis=1).sort_values(ascending=False).index[:12]
candidatos = del_taller + [g for g in mas_abundantes if g not in del_taller]

fig_genero = go.Figure()
titulos, p_menu = [], {}                                 # p_menu: el valor p de cada género del menú
for i, genero in enumerate(candidatos):
    valores = pct.loc[genero]                            # % de lecturas de ese género en cada muestra
    a, b = valores[datos.grupo == "Saludable"], valores[datos.grupo == "No saludable"]
    p = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue
    p_menu[genero] = p
    for g, v in zip(GRUPOS, (a, b)):
        # al abrir solo se ve el primer género
        fig_genero.add_trace(violin(v, g, "%{y:.3f} % de las lecturas", visible=(i == 0), leyenda=False))
    titulos.append(f"<i>{genero}</i> ({filo_de[genero]}): {a.mean():.3f} % y {b.mean():.3f} % en promedio; Mann-Whitney p = {formato_p(p)}")

# Cada botón del menú enciende los dos violines de su género y cambia el título
botones = [dict(label=genero, method="update",
                args=[{"visible": [k // 2 == i for k in range(2 * len(candidatos))]}, {"title.text": titulos[i]}])
           for i, genero in enumerate(candidatos)]
estilo(fig_genero, titulos[0], alto=500)
fig_genero.update_layout(updatemenus=[dict(type="dropdown", buttons=botones, x=0, xanchor="left", y=1.02, yanchor="bottom",
                                           bgcolor=P["control"], bordercolor=P["acento"], font=dict(color=P["control_texto"]))],
                         yaxis_title="% de las lecturas clasificadas", margin=dict(t=110), showlegend=False)
fig_genero.show()
""", figuras=1)
pausa(r"""
Recorre el menú. ¿Cuántos de los quince géneros tienen p < 0.05? Si ninguno difiriera de verdad, ¿cuántos esperarías ver por azar entre quince?
""")
md(r"""
Solo uno de los quince, *Phytophthora*, tiene p < 0.05, y entre quince pruebas se espera cerca de uno por azar.

### 3.2 Todos los géneros en una figura

Ahora la misma comparación para **todos** los géneros a la vez. Cada burbuja es un género y su tamaño indica qué tan abundante es: a la derecha están los más abundantes en las plantas no saludables y a la izquierda los más abundantes en las saludables; cuanto más arriba, más pequeño es su valor p. Esta figura se conoce como *gráfico de volcán*.

Con tantas pruebas hay que **corregir los valores p**. El método de Benjamini-Hochberg los ajusta para controlar la proporción de falsos descubrimientos: ordena los valores p de menor a mayor y exige más a los primeros.
""")
code(r"""
def benjamini_hochberg(p):
    "Valores p ajustados por Benjamini-Hochberg, para controlar la proporción de falsos descubrimientos."
    p = np.asarray(p)
    n = len(p)
    orden = np.argsort(p)                                        # posiciones de menor a mayor valor p
    ajustados = p[orden] * n / np.arange(1, n + 1)               # p por (número de pruebas / puesto que ocupa)
    ajustados = np.minimum.accumulate(ajustados[::-1])[::-1]     # un p ajustado nunca supera al del puesto siguiente
    salida = np.empty(n)
    salida[orden] = np.minimum(ajustados, 1)                     # vuelve al orden original, sin pasar de 1
    return salida

def comparar_generos(porcentajes, grupo):
    "Compara cada género (fila) entre grupos: medias, cambio relativo, valor p de Mann-Whitney y valor p ajustado."
    a = porcentajes.loc[:, (grupo == "Saludable").to_numpy()].to_numpy()
    b = porcentajes.loc[:, (grupo == "No saludable").to_numpy()].to_numpy()
    # axis=1: una prueba por fila; con cientos de géneros y empates se usa la aproximación normal
    p = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic", axis=1).pvalue
    p = np.nan_to_num(p, nan=1.0)                                # géneros sin ninguna variación: no hay nada que comparar
    minimo = 1e-6                                                # evita dividir entre cero en géneros ausentes de un grupo
    return pd.DataFrame({"media_sal": a.mean(axis=1), "media_no": b.mean(axis=1),
                         "cambio": np.log2((b.mean(axis=1) + minimo) / (a.mean(axis=1) + minimo)),   # log2: +1 = el doble
                         "p": p, "p_ajustado": benjamini_hochberg(p)}, index=porcentajes.index)

comparacion_sim = comparar_generos(pct_sim, datos_sim.grupo)
comparacion = comparar_generos(pct, datos.grupo)
comparacion["filo"] = filo_de.reindex(comparacion.index)
comparacion_sim["filo"] = "género simulado"

# Tres clases de punto, cada una con su color
clases = [("p ≥ 0.05", P["neutro"], lambda t: t.p >= 0.05),
          ("p < 0.05, pero no tras corregir", P["aviso"], lambda t: (t.p < 0.05) & (t.p_ajustado >= 0.05)),
          ("p < 0.05 también tras corregir", P["alerta"], lambda t: t.p_ajustado < 0.05)]
fig_volcan = go.Figure()
titulos = []
for i, (caso, tabla) in enumerate([("Datos reales", comparacion), ("Datos simulados", comparacion_sim)]):
    for nombre, color, condicion in clases:
        sub = tabla[condicion(tabla)]
        # el tamaño de la burbuja crece con la abundancia del género (la mayor de sus dos medias)
        abundancia = sub[["media_sal", "media_no"]].max(axis=1)
        fig_volcan.add_trace(go.Scatter(
            x=sub.cambio, y=-np.log10(sub.p), mode="markers", name=nombre,
            marker=dict(color=color, size=np.clip(5 + 6 * np.sqrt(abundancia), 5, 28), opacity=0.75, line=dict(width=0.5, color=P["panel"])),
            text=[f"<b>{genero}</b> ({filo})" for genero, filo in zip(sub.index, sub.filo)],          # nombre y filo del género
            customdata=sub[["media_sal", "media_no", "p", "p_ajustado"]].to_numpy(),                  # los números que se leen al pasar el cursor
            hovertemplate=("%{text}<br>saludables: %{customdata[0]:.4f} %<br>no saludables: %{customdata[1]:.4f} %"
                           "<br>p = %{customdata[2]:.2g} · p ajustado = %{customdata[3]:.2g}<extra></extra>"),
            visible=(i == 0)))
    # una cuarta traza solo con texto: el nombre de los géneros que siguen siendo significativos tras corregir
    destacados = tabla[tabla.p_ajustado < 0.05]
    lado = ["middle right" if cambio > 0 else "middle left" for cambio in destacados.cambio]   # el nombre, hacia afuera
    fig_volcan.add_trace(go.Scatter(x=destacados.cambio, y=-np.log10(destacados.p), mode="text",
                                    text=[f"  {genero}  " for genero in destacados.index], textposition=lado,
                                    textfont=dict(color=P["texto"], size=12), hoverinfo="skip", showlegend=False, visible=(i == 0)))
    titulos.append(f"{caso}: {int((tabla.p < 0.05).sum())} de {len(tabla)} géneros con p < 0.05; "
                   f"{int((tabla.p_ajustado < 0.05).sum())} tras corregir por comparaciones múltiples")

# Cada caso tiene cuatro trazas (tres clases de punto y los nombres): el botón enciende las suyas
botones = [dict(label=caso, method="update",
                args=[{"visible": [k // 4 == i for k in range(8)]}, {"title.text": titulos[i]}])
           for i, caso in enumerate(["Datos reales", "Datos simulados"])]
estilo(fig_volcan, titulos[0], alto=560)
fig_volcan.add_hline(y=-np.log10(0.05), line_dash="dot", line_color=P["acento"], annotation_text="p = 0.05",
                     annotation_font_color=P["acento"])
fig_volcan.update_layout(updatemenus=[dict(type="buttons", direction="right", buttons=botones, x=0, xanchor="left", y=1.02, yanchor="bottom",
                                           bgcolor=P["control"], bordercolor=P["acento"], font=dict(color=P["control_texto"]))],
                         margin=dict(t=120), yaxis_title="−log10 del valor p (más arriba = p más pequeño)",
                         xaxis_title="Cambio en las no saludables (log2: +1 = el doble, −1 = la mitad)")
fig_volcan.show()
print(titulos[0]); print(titulos[1])
print("Por puro azar se esperaría cerca del 5 %:", round(0.05 * len(comparacion)), "géneros en los datos reales.")
""", figuras=1)
md(r"""
**Datos reales.** Muchos géneros superan la línea de p = 0.05, pero fíjate en dos cosas: son cambios muy pequeños (casi todos los puntos están pegados al cero del eje horizontal) y, tras la corrección, no queda **ninguno**. Con 1 795 pruebas, encontrar decenas de valores p pequeños no es un hallazgo: es lo que se espera.

**Datos simulados.** Cambia de caso con los botones. `G05`, el género que hicimos cuatro veces más abundante, queda muy arriba y a la derecha. Pero mira los demás: muchos aparecen a la **izquierda** y varios sobreviven a la corrección, aunque no los tocamos. Bajaron solo porque `G05` subió y las proporciones tienen que sumar 100. Es la **composicionalidad**: una diferencia real en un taxón produce diferencias aparentes en todos los demás.
""")
ejercicio("Otros géneros en el menú", 4,
r"""
Agrega al menú del apartado 3.1 los géneros *Bacillus* y *Paenibacillus*, dos bacterias muy estudiadas en la rizósfera, y vuelve a ejecutar la celda. ¿Alguno difiere entre grupos?
""",
r"""
Se agregan a la lista `del_taller`. Ninguno de los dos difiere de forma significativa (p ≈ 0.48 y 0.15). *Paenibacillus* es instructivo: su promedio es casi tres veces mayor en las plantas no saludables, pero se debe a una sola muestra con un 11.6 % de sus lecturas en ese género; pasa el cursor para ver cuál es (otra vez `MD2086`). Las medianas de los dos grupos son casi iguales. Mann-Whitney compara rangos, no promedios, y por eso una sola muestra no lo arrastra. Si al escribir un nombre aparece un error `KeyError`, ese género no está en la tabla: los nombres deben coincidir exactamente con los de la columna `genero`.
""",
r"""
# solo cambia esta línea de la celda del menú; lo demás queda igual
del_taller = ["Fusarium", "Phytophthora", "Ralstonia", "Bacillus", "Paenibacillus"]
""")
puntos_clave([
    "Un menú desplegable permite dejar muchas comparaciones en una sola figura.",
    "Con cientos de pruebas, decenas de valores p menores que 0.05 aparecen por azar: hay que corregir.",
    "En datos composicionales, un taxón que cambia de verdad arrastra a los demás.",
])

# ================================================================== EPISODIO 4
episodio("La comunidad completa y el resumen",
         ["¿Cómo muestro de forma interactiva la comparación de comunidades completas?"],
         ["Construir una ordenación en la que cada muestra se puede identificar.",
          "Presentar los resultados de todas las pruebas en una tabla."])

md(r"""
**Qué hacemos.** La ordenación del taller mostraba nubes de puntos. Ahora cada punto dice qué muestra es, se marca con una estrella el **centroide** de cada grupo y una línea une cada muestra con él. La longitud de esas líneas es exactamente lo que compara PERMDISP.

### 4.1 Una ordenación que se puede interrogar
""")
code(r"""
resultados = {}                                          # PERMANOVA y PERMDISP de cada caso, para los títulos
for caso, matriz, grupo in [("Simulados", D_sim, etiquetas_sim), ("Reales", D, etiquetas)]:
    F_caso, R2_caso, p_caso = permanova(matriz, grupo)
    _, F_disp, p_disp = permdisp(matriz, grupo)
    resultados[caso] = dict(F=F_caso, R2=R2_caso, p=p_caso, F_disp=F_disp, p_disp=p_disp)

fig_ord = make_subplots(rows=1, cols=2, subplot_titles=[
    f"{caso}<br>PERMANOVA R² = {r['R2']:.3f}, p = {r['p']:.3f} · PERMDISP p = {r['p_disp']:.3f}"
    for caso, r in zip(["Datos simulados (caso ideal)", "Datos reales (fresa)"], resultados.values())])
for columna, (tabla, matriz) in enumerate([(datos_sim, D_sim), (datos, D)], start=1):
    coords = pcoa(matriz)                                # columna 0 = eje 1, columna 1 = eje 2
    explicado = 100 * (coords ** 2).sum(axis=0) / (coords ** 2).sum()   # % de la variación que recoge cada eje
    lejanas = []                                         # la muestra más alejada de su centroide, por grupo
    for g in GRUPOS:
        m = (tabla.grupo == g).to_numpy()                # máscara: True en las muestras del grupo g
        cx, cy = coords[m, 0].mean(), coords[m, 1].mean()   # centroide del grupo: el promedio de sus coordenadas
        distancia = np.hypot(coords[m, 0] - cx, coords[m, 1] - cy)
        lejanas.append((distancia.max(), tabla.index[m][distancia.argmax()], coords[m, 0][distancia.argmax()], coords[m, 1][distancia.argmax()]))
        # una línea de cada muestra a su centroide: su longitud es lo que mide PERMDISP (None corta la línea entre muestras)
        rayos_x = [v for px in coords[m, 0] for v in (cx, px, None)]
        rayos_y = [v for py in coords[m, 1] for v in (cy, py, None)]
        fig_ord.add_trace(go.Scatter(x=rayos_x, y=rayos_y, mode="lines", line=dict(color=con_alfa(COLORES[g], 0.35), width=1),
                                     hoverinfo="skip", legendgroup=g, showlegend=False), row=1, col=columna)
        fig_ord.add_trace(go.Scatter(
            x=coords[m, 0], y=coords[m, 1], mode="markers", name=g,
            marker=dict(color=COLORES[g], size=10, opacity=0.9, line=dict(width=0.6, color=P["panel"])),
            text=list(tabla.index[m]), customdata=tabla.shannon[m].to_numpy(),
            hovertemplate="%{text}<br>Shannon = %{customdata:.3f}<extra>" + g + "</extra>",
            legendgroup=g, showlegend=(columna == 1)), row=1, col=columna)
        # el centroide, marcado con una estrella
        fig_ord.add_trace(go.Scatter(
            x=[cx], y=[cy], mode="markers", name="centroide",
            marker=dict(symbol="star", size=20, color=COLORES[g], line=dict(color=P["texto"], width=1.5)),
            hovertemplate="centroide de " + g.lower() + "<extra></extra>", legendgroup=g, showlegend=False), row=1, col=columna)
    if columna == 2:                                     # en los datos reales se nombra la muestra más alejada de su centroide
        _, lejana, px, py = max(lejanas)
        fig_ord.add_annotation(x=px, y=py, text=f"<b>{lejana}</b>", ax=-45, ay=-28, showarrow=True, arrowcolor=P["texto"],
                               font=dict(color=P["texto"], size=12), row=1, col=2)
    fig_ord.update_xaxes(title_text=f"Eje 1 ({explicado[0]:.0f} % de la variación)", row=1, col=columna)
    fig_ord.update_yaxes(title_text=f"Eje 2 ({explicado[1]:.0f} %)", row=1, col=columna)
estilo(fig_ord, "Ordenación (PCoA) sobre distancias de Bray-Curtis; ★ = centroide de cada grupo", alto=540)
fig_ord.update_layout(margin=dict(t=110))                # sitio para los títulos de dos líneas de los paneles
paneles_relativos(fig_ord)
fig_ord.show()
""", figuras=1)
pausa(r"""
En los datos reales hay una muestra muy alejada de las demás. ¿Cuál es? ¿Es la misma que tenía el Shannon más bajo en la primera figura?
""")
md(r"""
La muestra alejada es `MD2086`, la misma que tenía el Shannon más bajo: una sola planta explica buena parte de la mayor dispersión del grupo no saludable.

### 4.2 La tabla de resultados

Una tabla también puede ser una figura de Plotly. Se colorean las celdas con p < 0.05 para que los dos casos se comparen de un vistazo.
""")
code(r"""
def pruebas_de_un_caso(x, y, porcentaje_x, porcentaje_y, D, etiquetas):
    "Todas las pruebas de los talleres para un caso: lista de (variable, prueba, estadístico, valor p)."
    dif, _, p_perm = prueba_permutacion(x, y)
    F = x.var(ddof=1) / y.var(ddof=1)
    cola = stats.f.cdf(F, len(x) - 1, len(y) - 1)
    student, welch = stats.ttest_ind(x, y, equal_var=True), stats.ttest_ind(x, y, equal_var=False)
    mw = stats.mannwhitneyu(x, y, alternative="two-sided", method="exact")
    mw_genero = stats.mannwhitneyu(porcentaje_x, porcentaje_y, alternative="two-sided", method="exact")
    F_perm, R2, p_permanova = permanova(D, etiquetas)
    _, F_disp, p_disp = permdisp(D, etiquetas)
    return [("Índice de Shannon", "Permutación de medias", f"dif = {dif:.3f}", p_perm),
            ("Índice de Shannon", "F de varianzas", f"F = {F:.3f}", 2 * min(cola, 1 - cola)),
            ("Índice de Shannon", "t de Student", f"t = {student.statistic:.2f}", student.pvalue),
            ("Índice de Shannon", "t de Welch", f"t = {welch.statistic:.2f}", welch.pvalue),
            ("Índice de Shannon", "Mann-Whitney", f"U = {mw.statistic:.0f}", mw.pvalue),
            ("Género de interés (%): G05 · <i>Fusarium</i>", "Mann-Whitney", f"U = {mw_genero.statistic:.0f}", mw_genero.pvalue),
            ("Composición (Bray-Curtis)", "PERMANOVA", f"F = {F_perm:.2f}, R² = {R2:.3f}", p_permanova),
            ("Composición (Bray-Curtis)", "PERMDISP", f"F = {F_disp:.2f}", p_disp)]

sal_sim, no_sal_sim = (datos_sim.grupo == "Saludable").to_numpy(), (datos_sim.grupo == "No saludable").to_numpy()
sal, no_sal = (datos.grupo == "Saludable").to_numpy(), (datos.grupo == "No saludable").to_numpy()
filas_sim = pruebas_de_un_caso(x_sim, y_sim, pct_sim.loc[patogeno].to_numpy()[sal_sim], pct_sim.loc[patogeno].to_numpy()[no_sal_sim], D_sim, etiquetas_sim)
filas_real = pruebas_de_un_caso(x, y, pct.loc["Fusarium"].to_numpy()[sal], pct.loc["Fusarium"].to_numpy()[no_sal], D, etiquetas)

# Colores de la tabla, tomados del tema: fondo normal, celdas con p < 0.05 resaltadas y cabecera
fondo, resaltado, cabecera, texto = P["panel"], con_alfa(P["alerta"], 0.35), P["rejilla"], P["texto"]
fig_tabla = go.Figure(go.Table(
    columnwidth=[1.5, 1.1, 1.6, 1.6],
    header=dict(values=["<b>Variable</b>", "<b>Prueba</b>", "<b>Datos simulados (caso ideal)</b>", "<b>Datos reales (fresa)</b>"],
                fill_color=cabecera, line_color=P["rejilla"], font=dict(color=texto, size=13), align="left", height=34),
    cells=dict(values=[[f[0] for f in filas_sim], [f[1] for f in filas_sim],
                       [f"{f[2]}; p = {formato_p(f[3])}" for f in filas_sim],
                       [f"{f[2]}; p = {formato_p(f[3])}" for f in filas_real]],
               # un color por celda: las dos primeras columnas, fondo normal; las otras dos, según su valor p
               fill_color=[[fondo] * 8, [fondo] * 8,
                           [resaltado if f[3] < 0.05 else fondo for f in filas_sim],
                           [resaltado if f[3] < 0.05 else fondo for f in filas_real]],
               line_color=P["rejilla"], font=dict(color=texto, size=13), align="left", height=32)))
estilo(fig_tabla, "Todas las pruebas, lado a lado (resaltado: p < 0.05)", alto=430)
fig_tabla.update_layout(margin=dict(l=20, r=20, t=70, b=10))
fig_tabla.show()
""", figuras=1)
md(r"""
La tabla resume la lección de los talleres. En el **caso ideal** todas las pruebas sobre el índice de Shannon coinciden y la composición difiere sin que difiera la dispersión. En el **caso real** las pruebas no coinciden, porque los supuestos fallan, y lo único sólido es que las plantas no saludables varían más: en diversidad (prueba F) y en composición (PERMDISP).
""")
puntos_clave([
    "En una ordenación interactiva cada muestra se puede identificar, y eso ayuda a encontrar las atípicas.",
    "El centroide de cada grupo es la referencia de PERMDISP: lo que se compara es la distancia de las muestras a él.",
    "Una tabla con las celdas significativas resaltadas deja comparar los dos casos de un vistazo.",
])

# ================================================================== EPISODIO 5
episodio("Armar y guardar la página HTML",
         ["¿Cómo reúno las figuras en una página que cualquiera pueda abrir?"],
         ["Convertir cada figura en un bloque de HTML.",
          "Armar una página con títulos, texto y figuras, y guardarla en un archivo.",
          "Ver la página dentro del cuaderno y descargarla."])

md(r"""
**Qué hacemos.** Cada figura de Plotly sabe escribirse como un bloque de HTML con `fig.to_html()`. Una página es, entonces, una plantilla de texto en la que se pegan esos bloques entre títulos y párrafos.

**Tres decisiones.**

- **La biblioteca de Plotly.** La página la necesita para dibujar; se carga una sola vez en la cabecera y todas las figuras la reutilizan. Con `"cdn"` se descarga de internet al abrir la página (archivo pequeño, necesita conexión); con `True` va dentro del archivo (pesa unos 5 MB más, funciona sin conexión).
- **El texto.** Cada sección lleva un título, un **qué mirar** antes de la figura (con la instrucción de qué control usar) y un **qué concluir** después. Una figura interactiva sin esa guía se queda sin usar. Los números que citan los textos se toman de las variables calculadas arriba, para que texto y figura nunca se contradigan.
- **El orden.** Las secciones siguen el del taller: los datos, el valor p, elegir la prueba, un género de interés, todos los géneros, la composición completa y el resumen. Cada una indica a qué episodio del taller corresponde.

La página abre con unas **tarjetas de cifras**, que resumen el análisis antes de la primera figura, lleva un índice y toma sus colores del mismo tema que las figuras. Un detalle más: los botones de la barra de herramientas de Plotly vienen en inglés, así que la página registra una traducción (`Plotly.register`) y cada figura la pide con `locale: "es"`.
""")
code(r"""
PLOTLY_JS = "cdn"        # "cdn": la página carga Plotly de internet; True: Plotly va dentro del archivo (sin conexión)
TITULO = "Pruebas de hipótesis con datos metagenómicos: explorador interactivo"
REPOSITORIO = "https://github.com/CamilaSilva1995/Tesis_Maestria"

# ---------- Los números que citan los textos: salen de las pruebas de arriba, así los textos y las figuras no se contradicen
p_real = dict(zip(["permutación", "F", "Student", "Welch", "Mann-Whitney", "Fusarium", "PERMANOVA", "PERMDISP"], [f[3] for f in filas_real]))
p_sim = dict(zip(["permutación", "F", "Student", "Welch", "Mann-Whitney", "G05", "PERMANOVA", "PERMDISP"], [f[3] for f in filas_sim]))
r_sim, r_real = resultados["Simulados"], resultados["Reales"]
fus = comparacion.loc["Fusarium"]                                        # medias y valor p de Fusarium
alarma_real = dict(zip(pruebas, alarmas[disenos[0][0]]))                 # falsas alarmas con el diseño de los datos reales
menu_p05 = [g for g in candidatos if p_menu[g] < 0.05]                   # géneros del menú con p < 0.05
n_p05, n_bh = int((comparacion.p < 0.05).sum()), int((comparacion.p_ajustado < 0.05).sum())
n_bh_sim = int((comparacion_sim.p_ajustado < 0.05).sum())
miles = lambda n: f"{n:,}".replace(",", " ")                             # 1795 -> "1 795"
pct_txt = lambda v: f"{v:.0f} %"

# ---------- Las secciones, en el orden del taller: datos → valor p → elegir la prueba → un género → todos → composición → resumen
# Cada una: (ancla, nombre corto para el índice, título, episodio del taller, qué mirar, figura, qué concluir)
secciones = [
    ("comunidad", "Los datos: la comunidad", "1. Los datos: la comunidad de un vistazo", "Episodio 1 del taller",
     f"<p><b>Qué mirar.</b> Cada barra es una de las {len(datos)} muestras reales y cada color, un filo (los ocho más abundantes; el resto, en «Otros filos»). "
     "<b>Haz clic</b> en un filo de la leyenda para ocultarlo y doble clic para dejarlo solo; pasa el cursor para leer el porcentaje.</p>",
     fig_filos,
     f"<p class='conclusion'><b>Qué concluir.</b> Los mismos filos dominan en los dos grupos: a simple vista no hay diferencia, y por eso hacen falta las pruebas. "
     f"Una planta no saludable, {atipica}, se sale del patrón con un {pct_filo.loc['Firmicutes', atipica]:.0f} % de Firmicutes; vuelve a aparecer en las secciones 2 y 7.</p>"),
    ("diversidad", "Los datos: la diversidad", "2. Los datos: la diversidad que vamos a comparar", "Episodios 1 y 3 del taller",
     "<p><b>Qué mirar.</b> El índice de Shannon de cada muestra: el violín dibuja la distribución, la caja los cuartiles y la línea punteada, la media. "
     "<b>Pasa el cursor</b> sobre un punto para ver qué muestra es. Compara la escala de los dos paneles y cuánto se dispersan los puntos turquesa en los datos reales.</p>",
     fig_cajas,
     f"<p class='conclusion'><b>Qué concluir.</b> En el caso ideal los dos grupos varían lo mismo y las pruebas coinciden (Welch {p_en_texto(p_sim['Welch'])}). "
     f"En los datos reales la media es un poco mayor en las saludables, pero la diferencia no se distingue del azar (Welch {p_en_texto(p_real['Welch'])}; "
     f"Mann-Whitney {p_en_texto(p_real['Mann-Whitney'])}). Lo que sí es distinto es la dispersión: las plantas no saludables son más heterogéneas "
     f"(prueba F {p_en_texto(p_real['F'])}), y una de ellas, {atipica}, queda muy por debajo de las demás.</p>"),
    ("valor-p", "El valor p", "3. ¿Diferencia real o azar? El valor p", "Episodio 2 del taller",
     "<p><b>Qué mirar.</b> El histograma es lo que produce el azar cuando el grupo no importa (las etiquetas, barajadas 5 000 veces); las líneas amarillas, "
     "la diferencia observada y su opuesta. Las barras rosas son las barajadas tan extremas como lo observado: <b>su proporción es el valor p</b>. "
     "<b>Mueve el deslizador</b> de izquierda a derecha.</p>",
     fig_perm,
     f"<p class='conclusion'><b>Qué concluir.</b> La diferencia verdadera es siempre la misma, pero con más muestras el histograma del azar se estrecha y el valor p baja: "
     "el valor p depende de la diferencia, de cuánto varían los datos y de cuántas muestras hay. En la última posición, los datos reales dan "
     f"{p_en_texto(p_real['permutación'])} con un histograma asimétrico; guarda ese número: la siguiente sección muestra por qué no es de fiar.</p>"),
    ("elegir", "Elegir la prueba", "4. Elegir la prueba: ¿cuántas veces se equivoca cada una?", "Episodios 3 y 4 del taller",
     f"<p><b>Qué mirar.</b> Aquí no hay ninguna diferencia verdadera: se simulan {miles(n_estudios)} estudios con la misma media en los dos grupos y se cuenta "
     "en cuántos cada prueba dice que sí la hay. Lo correcto es cerca del 5 %, la línea amarilla. El primer diseño copia a los datos reales: 35 y 18 muestras "
     "y el grupo pequeño mucho más variable. <b>Mueve el deslizador</b> para cambiar el diseño del estudio.</p>",
     fig_alarmas,
     f"<p class='conclusion'><b>Qué concluir.</b> Con el diseño de los datos reales, la t de Student ve diferencias que no existen en el {alarma_real['t de Student']:.0f} % "
     f"de los estudios y Mann-Whitney en el {alarma_real['Mann-Whitney']:.0f} %; Welch se queda en el {alarma_real['t de Welch']:.0f} %, y lo hace en los cuatro diseños. "
     "Barajar etiquetas (sección 3) y la t de Student suponen que los dos grupos tienen la misma distribución, varianza incluida, y en los datos reales eso falla "
     f"(prueba F {p_en_texto(p_real['F'])}). Por eso el resultado que vale es el de Welch, {p_en_texto(p_real['Welch'])}: no hay evidencia de que la diversidad promedio difiera. "
     "Con grupos del mismo tamaño el problema casi desaparece: una razón para diseñar estudios balanceados.</p>"),
    ("genero", "Un género: Fusarium", "5. Un género de interés: ¿es más abundante en un grupo?", "Episodio 5 del taller",
     f"<p><b>Qué mirar.</b> La proporción de lecturas de un género en cada muestra. <b>Elige un género</b> en el menú: los del taller y los doce con más lecturas. "
     "El título da las medias de los dos grupos y el valor p de Mann-Whitney; pasa el cursor para ver cada muestra.</p>",
     fig_genero,
     f"<p class='conclusion'><b>Qué concluir.</b> <i>Fusarium</i>, el patógeno sospechoso, no es más abundante en las plantas no saludables: "
     f"{fus.media_sal:.3f} % frente a {fus.media_no:.3f} % de las lecturas, {p_en_texto(fus.p)}. De los {len(candidatos)} géneros del menú, "
     f"{'solo ' + ' y '.join(f'<i>{g}</i>' for g in menu_p05) + ' tiene' if len(menu_p05) == 1 else str(len(menu_p05)) + ' tienen'} p &lt; 0.05, "
     f"y entre {len(candidatos)} pruebas se espera cerca de uno por puro azar: un p pequeño aislado, sin corrección, no es un hallazgo.</p>"),
    ("volcan", "Todos los géneros", "6. Todos los géneros a la vez", "Episodio 5 del taller",
     f"<p><b>Qué mirar.</b> Cada burbuja es un género ({miles(len(comparacion))} en los datos reales) y su tamaño, su abundancia: a la derecha, más abundante "
     "en las plantas no saludables; más arriba, valor p más pequeño. La línea punteada es p = 0.05. <b>Cambia de caso</b> con los botones y pasa el cursor "
     "para leer el nombre del género, sus medias y su p ajustado por Benjamini-Hochberg.</p>",
     fig_volcan,
     f"<p class='conclusion'><b>Qué concluir.</b> Datos reales: {n_p05} géneros superan p = 0.05, cerca de los {round(0.05 * len(comparacion))} que se esperan por azar "
     f"entre {miles(len(comparacion))} pruebas, y tras corregir no queda {'ninguno' if n_bh == 0 else n_bh}: no hay evidencia de que algún género difiera entre grupos. "
     f"Datos simulados: G05, el género que hicimos cuatro veces más abundante, destaca arriba a la derecha, pero otros {n_bh_sim - 1} sobreviven a la corrección sin que los "
     "tocáramos: las proporciones suman 100 y un taxón que sube hace bajar a los demás (composicionalidad).</p>"),
    ("composicion", "La composición completa", "7. La composición completa: PERMANOVA y PERMDISP", "Episodio 6 del taller",
     "<p><b>Qué mirar.</b> Cada punto es una muestra y dos puntos cercanos tienen composiciones parecidas (ordenación PCoA sobre distancias de Bray-Curtis). "
     "La estrella es el centroide de cada grupo y cada línea une una muestra con el suyo: la separación entre estrellas es lo que mide PERMANOVA y la longitud "
     "de las líneas, lo que compara PERMDISP. <b>Pasa el cursor</b> sobre un punto para ver qué muestra es y su Shannon.</p>",
     fig_ord,
     f"<p class='conclusion'><b>Qué concluir.</b> En el caso ideal las nubes se separan (PERMANOVA R² = {r_sim['R2']:.2f}, p = {r_sim['p']:.3f}) y se dispersan igual "
     f"(PERMDISP p = {r_sim['p_disp']:.2f}). En los datos reales el grupo explica solo el {100 * r_real['R2']:.0f} % de la variación (p = {r_real['p']:.2f}), "
     f"pero la nube turquesa es más amplia (PERMDISP p = {r_real['p_disp']:.3f}): las plantas no saludables no tienen otra composición, sino composiciones más distintas "
     f"entre sí. {atipica} es, otra vez, la más alejada.</p>"),
    ("resumen", "Resumen", "8. Resumen: todas las pruebas, lado a lado", "Episodio 7 del taller",
     "<p><b>Qué mirar.</b> Cada fila es una prueba y cada columna, un caso; se resaltan las celdas con p &lt; 0.05.</p>",
     fig_tabla,
     f"<p class='conclusion'><b>Qué concluir.</b> Cuando los supuestos se cumplen (caso ideal) todas las pruebas coinciden. En el caso real no coinciden: permutación y "
     f"Student dan {p_en_texto(p_real['permutación'])} y {p_real['Student']:.3f}; Welch y Mann-Whitney, {p_real['Welch']:.3f} y "
     f"{p_real['Mann-Whitney']:.3f}, y la sección 4 explica a cuáles creer. Lo que queda en pie es que las plantas no saludables son más heterogéneas, "
     "en diversidad (prueba F) y en composición (PERMDISP); el estado de la planta explica el "
     f"{100 * r_real['R2']:.0f} % de la variación en la composición y <i>Fusarium</i> no difiere.</p>"),
]

# Las cifras de la cabecera: (número, qué significa)
cifras = [
    (str(len(datos)), f"muestras reales: {int(sal.sum())} saludables y {int(no_sal.sum())} no saludables"),
    (str(len(datos_sim)), "muestras simuladas: 30 y 30"),
    (miles(len(comparacion)), "géneros comparados en los datos reales"),
    (f"{n_p05} → {n_bh}", "géneros con p < 0.05, antes y después de corregir"),
    (f"{pct_txt(100 * r_sim['R2'])} · {pct_txt(100 * r_real['R2'])}", "de la variación en la composición la explica el grupo: caso ideal y caso real"),
]

# ---------- La biblioteca de Plotly se carga una sola vez, en la cabecera, y todas las figuras la reutilizan
from plotly.offline import get_plotlyjs, get_plotlyjs_version
if PLOTLY_JS == "cdn":
    plotly_js = f'<script charset="utf-8" src="https://cdn.plot.ly/plotly-{get_plotlyjs_version()}.min.js"></script>'
else:
    plotly_js = "<script>" + get_plotlyjs() + "</script>"
# Los botones de la barra de herramientas de Plotly vienen en inglés: se registra una traducción y cada figura la pide con locale="es"
TRADUCCION = '''<script>
Plotly.register({moduleType: "locale", name: "es", dictionary: {
  "Zoom": "Ampliar: arrastra un rectángulo", "Pan": "Desplazar", "Zoom in": "Acercar", "Zoom out": "Alejar",
  "Autoscale": "Ajustar los ejes", "Reset axes": "Volver a la vista inicial",
  "Download plot as a png": "Descargar la figura como PNG", "Download plot as a PNG": "Descargar la figura como PNG",
  "Download plot": "Descargar la figura", "Box Select": "Seleccionar con un rectángulo", "Lasso Select": "Seleccionar con lazo",
  "Toggle Spike Lines": "Líneas guía", "Show closest data on hover": "Mostrar el dato más cercano", "Compare data on hover": "Comparar los datos",
  "Double-click on legend to isolate one trace": "Doble clic en la leyenda para dejar una sola serie",
  "Double-click to zoom back out": "Doble clic para volver", "Produced with Plotly.js": "Hecho con Plotly.js",
  "Reset view": "Volver a la vista inicial", "Reset views": "Volver a la vista inicial",
  "Taking snapshot - this may take a few seconds": "Generando la imagen; puede tardar unos segundos", "Snapshot succeeded": "Imagen generada",
  "Sorry, there was a problem downloading your snapshot!": "No se pudo descargar la imagen",
  "Preparing image - this may take a few seconds": "Generando la imagen; puede tardar unos segundos", "Image download succeeded": "Imagen descargada",
  "Sorry, there was a problem downloading your image!": "No se pudo descargar la imagen", "Share chart...": "Compartir",
  "new text": "nuevo texto", "trace": "serie"
}});
</script>'''
# En pantallas angostas (celulares) Plotly no parte los títulos ni reacomoda los paneles: este pequeño script lo hace al cargar
# la página y al cambiar el ancho de la ventana. En una pantalla ancha no cambia nada.
ADAPTAR = '''<script>
(function () {
  const originales = new Map();                    // por figura: lo que hay que restaurar al volver a una pantalla ancha
  function partirLinea(texto, maximo) {            // inserta <br> en un espacio (fuera de etiquetas) cuando la línea pasa de «maximo» letras
    let salida = "", linea = 0, enEtiqueta = false;
    for (let i = 0; i < texto.length; i++) {
      const c = texto[i];
      if (c === "<" && /[A-Za-z\\/]/.test(texto[i + 1] || "")) enEtiqueta = true;   // «p < 0.05» no es una etiqueta
      if (c === " " && !enEtiqueta && linea >= maximo) { salida += "<br>"; linea = 0; continue; }
      salida += c;
      if (!enEtiqueta) linea += 1;
      if (c === ">") enEtiqueta = false;
    }
    return salida;
  }
  function partir(texto, maximo) {                 // respeta los saltos de línea que ya trae el texto y parte cada línea por separado
    if (!texto) return texto;
    return texto.split(/<br\\s*\\/?>/).map(function (linea) { return partirLinea(linea, maximo); }).join("<br>");
  }
  const EJES = ["xaxis", "xaxis2", "yaxis", "yaxis2"];
  function tituloEje(diseno, eje) { return diseno[eje] && diseno[eje].title && diseno[eje].title.text; }
  function adaptar() {
    document.querySelectorAll(".plotly-graph-div").forEach(function (figura) {
      const diseno = figura.layout;
      if (!diseno) return;
      if (!originales.has(figura)) {
        originales.set(figura, JSON.parse(JSON.stringify({
          alto: diseno.height, margen: diseno.margin && diseno.margin.t, letra: diseno.font && diseno.font.size,
          titulo: diseno.title && diseno.title.font && diseno.title.font.size,
          dominios: diseno.yaxis2 ? [diseno.xaxis.domain, diseno.xaxis2.domain, diseno.yaxis.domain, diseno.yaxis2.domain] : null,
          ejes: EJES.map(function (eje) { return tituloEje(diseno, eje) || null; }),
          anotaciones: (diseno.annotations || []).map(function (a) { return a.text || ""; }),
          letraAnotaciones: (diseno.annotations || []).map(function (a) { return (a.font && a.font.size) || null; }),
          leyendaY: diseno.legend && diseno.legend.y, margenB: diseno.margin && diseno.margin.b,
          deslizadores: diseno.sliders || null, menus: diseno.updatemenus || null })));
      }
      const o = originales.get(figura);
      const angosta = figura.clientWidth < 640;
      const maximo = Math.floor(figura.clientWidth / (angosta ? 10 : 9));      // letras que caben en una línea del título
      const cambios = { "font.size": angosta ? 11 : o.letra, "title.font.size": angosta ? 14 : o.titulo, "height": o.alto, "margin.t": o.margen };
      // los títulos largos se parten: el que se ve y los que cambian los deslizadores y los menús
      cambios["title.text"] = partir(diseno.title && diseno.title.text, maximo);
      const retitular = function (control) { control.args[1]["title.text"] = partir(control.args[1]["title.text"], maximo); };
      if (o.deslizadores) cambios["sliders"] = JSON.parse(JSON.stringify(o.deslizadores)).map(function (s) {
        s.steps.forEach(retitular);
        if (angosta) s.pad = Object.assign({}, s.pad, { t: ((s.pad && s.pad.t) || 0) + 30 });   // sitio para un título de eje de dos líneas
        return s;
      });
      if (o.menus) cambios["updatemenus"] = JSON.parse(JSON.stringify(o.menus)).map(function (m) { m.buttons.forEach(retitular); return m; });
      EJES.forEach(function (eje, i) {                                           // los títulos de los ejes, igual
        if (o.ejes[i]) cambios[eje + ".title.text"] = angosta ? partir(o.ejes[i], eje[0] === "x" ? maximo : 34) : o.ejes[i];
      });
      o.anotaciones.forEach(function (texto, i) {                                // y los títulos de los paneles y demás rótulos
        if (texto) cambios["annotations[" + i + "].text"] = angosta ? partir(texto, 52) : texto;
        cambios["annotations[" + i + "].font.size"] = angosta ? 11 : o.letraAnotaciones[i];
      });
      if (typeof o.leyendaY === "number" && o.leyendaY < 0) {                    // leyenda debajo: más abajo si el título del eje crece
        cambios["legend.y"] = angosta ? o.leyendaY - 0.1 : o.leyendaY;
        cambios["margin.b"] = angosta ? (o.margenB || 60) + 40 : o.margenB;
      }
      const lineas = (cambios["title.text"] || "").split("<br>").length;
      const leyendaArriba = diseno.legend && diseno.legend.y >= 1;
      if (angosta) cambios["margin.t"] = o.margen + 18 * (lineas - 1) + (leyendaArriba ? 60 : 0);
      if (o.dominios) {                                                          // dos paneles: lado a lado o apilados
        const d = angosta ? [[0, 1], [0, 1], [0.58, 1], [0, 0.42]] : o.dominios;
        cambios["xaxis.domain"] = d[0]; cambios["xaxis2.domain"] = d[1]; cambios["yaxis.domain"] = d[2]; cambios["yaxis2.domain"] = d[3];
        if (angosta) { cambios["height"] = Math.round(o.alto * 1.75); cambios["margin.t"] += 25; }   // sitio para el título del panel de arriba
      }
      Plotly.relayout(figura, cambios);
    });
  }
  let espera;
  window.addEventListener("resize", function () { clearTimeout(espera); espera = setTimeout(adaptar, 200); });
  window.addEventListener("load", adaptar);
})();
</script>'''
CONFIG = {"displaylogo": False, "responsive": True, "locale": "es",            # sin logo, se adapta al ancho, botones en español
          "modeBarButtonsToRemove": ["lasso2d", "select2d", "toggleSpikelines", "hoverCompareCartesian", "hoverClosestCartesian",
                                     "toggleHover", "sendChartToCloud"]}       # fuera los botones que no se necesitan (y el de subir a la nube)

PLANTILLA_HTML = '''<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITULO__</title>
__PLOTLY__
__TRADUCCION__
<style>
  * { box-sizing: border-box; }
  body { margin: 0; font-family: system-ui, "Segoe UI", Roboto, Arial, sans-serif; color: __TEXTO__; line-height: 1.55;
         background: __FONDO__; background-image: radial-gradient(circle at 15% 0%, __BRILLO__ 0%, __FONDO__ 55%); background-attachment: fixed; }
  main { max-width: 1150px; margin: 0 auto; padding: 20px 16px 40px; }
  h1 { font-size: clamp(1.6rem, 4vw, 2.5rem); line-height: 1.15; margin: 0.4em 0 0.3em;
       background: linear-gradient(90deg, __SALUDABLE__, __ACENTO__, __NO_SALUDABLE__); -webkit-background-clip: text; background-clip: text; color: transparent; }
  h2 { font-size: 1.3rem; margin: 0 0 0.3em; }
  p { margin: 0.4em 0 0.9em; }
  b { color: __ACENTO__; }
  a { color: __ENLACE__; }
  .docente { color: __SUAVE__; font-size: 0.95rem; }
  .cifras { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; margin: 1.2em 0; }
  .cifra { background: __PANEL__; border: 1px solid __REJILLA__; border-radius: 14px; padding: 14px 16px; }
  .cifra strong { display: block; font-size: 1.9rem; line-height: 1.1; color: __ACENTO__; }
  .cifra span { font-size: 0.88rem; color: __SUAVE__; }
  nav { display: flex; flex-wrap: wrap; gap: 8px; margin: 1.2em 0 1.6em; }
  nav a { text-decoration: none; padding: 6px 14px; border-radius: 999px; border: 1px solid __REJILLA__; background: __PANEL__; font-size: 0.92rem; }
  nav a:hover { border-color: __ACENTO__; color: __ACENTO__; }
  section { background: __PANEL__; border: 1px solid __REJILLA__; border-radius: 16px; padding: 18px 18px 8px; margin: 0 0 22px;
            box-shadow: 0 10px 30px rgba(0, 0, 0, 0.25); overflow: hidden; }
  .cabecera { display: flex; justify-content: space-between; align-items: baseline; flex-wrap: wrap; gap: 6px 12px; }
  .episodio { font-size: 0.8rem; color: __SUAVE__; border: 1px solid __REJILLA__; border-radius: 999px; padding: 2px 10px; white-space: nowrap; }
  .conclusion { border-left: 3px solid __ACENTO__; padding: 6px 12px; margin: 0.6em 0 0.8em; background: rgba(127, 127, 127, 0.08); border-radius: 0 8px 8px 0; }
  .nota { color: __SUAVE__; font-size: 0.9rem; text-align: center; line-height: 1.5; }
  @media (max-width: 600px) { main { padding: 12px 10px 30px; } section { padding: 12px 10px 4px; border-radius: 12px; } .cifra strong { font-size: 1.5rem; } }
</style>
</head>
<body>
<main>
<h1>__TITULO__</h1>
<p>Las mismas pruebas, aplicadas a dos casos: unos datos <b>simulados</b> (60 muestras, 30 por grupo), en los que conocemos la verdad y los supuestos
se cumplen, y unos datos <b>reales</b>: 53 metagenomas <i>shotgun</i> de la rizósfera de fresa, 35 de plantas saludables y 18 de plantas no saludables,
facilitados por Solena Ag. Todas las figuras son interactivas: pasa el cursor sobre un punto, amplía una zona arrastrando y haz doble clic para volver.</p>
<p class="docente">Para la clase: las ocho secciones siguen el orden del taller (los datos, el valor p, elegir la prueba, un género de interés, todos los géneros
y la composición completa). En cada una, pide una predicción antes de usar el control y lee la conclusión después.</p>
<div class="cifras">__CIFRAS__</div>
<nav>__MENU__</nav>
__SECCIONES__
<p class="nota">Fuente: Silva Gomez, P. C. (2026). <i>Análisis estadístico de la diversidad microbiana a distintos niveles taxonómicos en el microbioma
rizosférico de plantas de fresa saludables y no saludables</i>. Tesis de Maestría en Ciencias Matemáticas, Posgrado Conjunto UMSNH-UNAM.
Datos facilitados por Solena Ag (2023).<br>Código, datos y talleres: <a href="__REPOSITORIO__">__REPOSITORIO_CORTO__</a>.
Página generada con Python y Plotly; __CONEXION__.</p>
</main>
__ADAPTAR__
</body>
</html>'''

bloques = []
for ancla, corto, titulo, episodio_taller, mirar, fig, concluir in secciones:
    figura_html = fig.to_html(full_html=False, include_plotlyjs=False, config=CONFIG)   # la biblioteca ya está en la cabecera
    bloques.append(f'<section id="{ancla}"><div class="cabecera"><h2>{titulo}</h2><span class="episodio">{episodio_taller}</span></div>'
                   f'{mirar}{figura_html}{concluir}</section>')
menu = " ".join(f'<a href="#{ancla}">{i}. {corto}</a>' for i, (ancla, corto, *_) in enumerate(secciones, 1))
tarjetas = "".join(f'<div class="cifra"><strong>{numero}</strong><span>{significado}</span></div>' for numero, significado in cifras)

pagina = PLANTILLA_HTML
# Cada marca de la plantilla se cambia por su contenido o por un color del tema
marcas = {"__TITULO__": TITULO, "__MENU__": menu, "__CIFRAS__": tarjetas, "__SECCIONES__": "\n".join(bloques),
          "__PLOTLY__": plotly_js, "__TRADUCCION__": TRADUCCION, "__ADAPTAR__": ADAPTAR, "__REPOSITORIO__": REPOSITORIO, "__REPOSITORIO_CORTO__": REPOSITORIO.replace("https://", ""),
          "__CONEXION__": "solo necesita internet para cargar Plotly" if PLOTLY_JS == "cdn" else "funciona sin conexión",
          "__BRILLO__": "#1B2A6B" if TEMA == "noche" else "#FFFFFF"}
marcas.update({f"__{nombre.upper()}__": color for nombre, color in P.items()})   # __FONDO__, __PANEL__, __ACENTO__...
for marca, valor in marcas.items():
    pagina = pagina.replace(marca, valor)
assert "__" not in pagina.replace("__plotly", ""), "quedó una marca sin sustituir"

RUTA_HTML = "explorador_pruebas_de_hipotesis.html"
with open(RUTA_HTML, "w", encoding="utf-8") as archivo:
    archivo.write(pagina)
print(f"Página guardada en {RUTA_HTML}: {len(secciones)} figuras, {len(pagina.encode('utf-8')) / 1024:.0f} KB")
""")
md(r"""
La página ya está guardada. Para verla sin salir del cuaderno se incrusta en un marco (`iframe`), que la aísla del resto del cuaderno para que sus estilos no se mezclen.
""")
code(r"""
from IPython.display import HTML, display

# srcdoc recibe la página completa como texto; html.escape protege las comillas y los signos < >
display(HTML(f'<div><iframe srcdoc="{html.escape(pagina)}" width="100%" height="780" style="border:0; border-radius:12px"></iframe></div>'))
""")
md(r"""
Por último, el archivo. En Colab vive en el servidor de Google, así que hay que descargarlo; en tu computador ya quedó en la carpeta del cuaderno.
""")
code(r"""
try:
    from google.colab import files                       # solo existe en Colab
    files.download(RUTA_HTML)                            # el navegador descarga el archivo
except ImportError:
    print(f"Fuera de Colab no hace falta descargar: abre {RUTA_HTML} con doble clic desde la carpeta del cuaderno.")
""")
ejercicio("Tu propia versión de la página", 5,
r"""
1. En la primera celda del cuaderno cambia `TEMA = "noche"` por `TEMA = "claro"`.
2. En la celda que arma la página cambia `TITULO` por uno tuyo.
3. Ejecuta todo de nuevo (**Entorno de ejecución → Ejecutar todas**) y abre el archivo descargado.

¿Qué habría que cambiar para que la página funcione sin conexión a internet?
""",
r"""
Al ejecutar todo de nuevo, las figuras y la página pasan al fondo claro, porque todos sus colores salen de la paleta del tema elegido (`P = PALETAS[TEMA]`). Para que funcione sin conexión basta con `PLOTLY_JS = True`: la biblioteca de Plotly queda dentro del archivo, que pasa a pesar unos 5 MB más.
""")
puntos_clave([
    "`fig.to_html()` convierte una figura en un bloque de HTML; una página es una plantilla con esos bloques y la biblioteca de Plotly cargada una sola vez.",
    "La interactividad viaja en el archivo: no hace falta Python ni un servidor para abrirlo.",
    "Una figura interactiva necesita una frase que diga qué explorar.",
])

# ================================================================== CIERRE
CIERRE.append(nbf.v4.new_markdown_cell(r"""
---
## Referencias

**Figuras interactivas**

- Plotly Technologies Inc. (2015). *Collaborative data science*. Plotly Technologies Inc. https://plotly.com/python/
- Documentación de Plotly para Python: [diagramas de caja](https://plotly.com/python/box-plots/), [deslizadores](https://plotly.com/python/sliders/), [menús desplegables](https://plotly.com/python/dropdowns/), [tablas](https://plotly.com/python/table/) y [exportar a HTML](https://plotly.com/python/interactive-html-export/).

**Métodos estadísticos**

- Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26(1), 32–46. https://doi.org/10.1111/j.1442-9993.2001.01070.x
- Anderson, M. J. (2006). Distance-based tests for homogeneity of multivariate dispersions. *Biometrics*, 62(1), 245–253. https://doi.org/10.1111/j.1541-0420.2005.00440.x
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J. (2017). Microbiome datasets are compositional: and this is not optional. *Frontiers in Microbiology*, 8, 2224. https://doi.org/10.3389/fmicb.2017.02224
- Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. *The Annals of Mathematical Statistics*, 18(1), 50–60. https://doi.org/10.1214/aoms/1177730491

**Fuente de los datos y del análisis**

- Silva Gomez, P. C. (2026). *Análisis estadístico de la diversidad microbiana a distintos niveles taxonómicos en el microbioma rizosférico de plantas de fresa saludables y no saludables*. Tesis de Maestría en Ciencias Matemáticas, Posgrado Conjunto UMSNH-UNAM. Código y datos: https://github.com/CamilaSilva1995/Tesis_Maestria. Datos metagenómicos facilitados por Solena Ag.
- Solena Ag. (2023). [Metagenomas *shotgun* de la rizósfera de plantas de fresa saludables y no saludables] [Conjunto de datos no publicado]. https://www.solena.ag
""".strip()))

CIERRE.append(nbf.v4.new_markdown_cell(r"""
---
*Datos facilitados por Solena Ag. Código y datos en [GitHub](https://github.com/CamilaSilva1995/Tesis_Maestria/tree/main/Taller_Practico).*
""".strip()))

# ================================================================== TIEMPOS, PORTADA Y ARMADO
def minutos(ep):
    "Minutos de explicación (según la regla de tiempo, redondeados hacia arriba) y de ejercicios de un episodio."
    explicacion = (ep["palabras"] / PALABRAS_POR_MINUTO + ep["codigo"] * MINUTOS_POR_CELDA
                   + ep["figuras"] * MINUTOS_POR_FIGURA + ep["pausas"] * MINUTOS_POR_PAUSA)
    return REDONDEO * math.ceil(explicacion / REDONDEO), sum(ep["ejercicios"])

def reloj(m):
    "Minutos acumulados como h:mm."
    return f"{m // 60}:{m % 60:02d}"

def duracion(m):
    "Una duración en palabras: 75 -> '1 h 15 min'."
    return f"{m // 60} h {m % 60:02d} min" if m >= 60 else f"{m} min"

# Índice con la hora de inicio de cada episodio
filas, inicio = [], 0
for i, ep in enumerate(EPISODIOS, 1):
    exp, ej = minutos(ep)
    ep["exp"], ep["ej"] = exp, ej
    filas.append(f"| {reloj(inicio)} | **{i}. {ep['titulo']}** | {ep['preguntas'][0]} | {exp} min | {f'{ej} min' if ej else '—'} |")
    inicio += exp + ej
filas.append(f"| {reloj(inicio)} | *Fin* | | | |")
total_exp, total_ej = sum(ep["exp"] for ep in EPISODIOS), sum(ep["ej"] for ep in EPISODIOS)
INDICE = "| Inicio | Episodio | Pregunta que responde | Explicación | Ejercicio |\n|---|---|---|---|---|\n" + "\n".join(filas)

PORTADA = r"""
# 🍓 Explorador interactivo: pruebas de hipótesis con datos metagenómicos

## De las figuras del taller a una página HTML con Plotly

**Python · Plotly · Google Colab · «N_EP» episodios · «TOTAL»**

Los talleres terminan con figuras que se miran. Este cuaderno las convierte en figuras que **se exploran**, sobre un fondo azul noche: al pasar el cursor cada punto dice qué muestra es, un deslizador cambia el tamaño de la muestra, otro el diseño del estudio y un menú cambia el género. Al final, todas quedan reunidas en **una página HTML** que se abre en cualquier navegador, sin Python ni conexión a un servidor.

Es el tercer cuaderno de la serie y usa los mismos dos casos, datos **simulados** y datos **reales**, con los mismos resultados:

- [Taller completo](__COMPLETO__), de unas tres horas: las pruebas, paso a paso.
- [Taller rápido](__RAPIDO__), de una hora: la versión resumida.
- **Explorador interactivo** (este cuaderno): los resultados, en una página para explorar y compartir.

## ¿Para quién es?

Para quien ya trabajó alguno de los dos talleres, o conoce las pruebas, y quiere **presentar y explorar los resultados de forma interactiva**: para una clase, un informe o una reunión de laboratorio. No hace falta haber usado Plotly: el código está comentado paso a paso.

## Conocimientos previos

- **Las pruebas estadísticas:** cualquiera de los dos talleres de esta serie. Aquí no se vuelven a explicar: se dibujan.
- **Python y pandas:** [Plotting and Programming in Python](https://swcarpentry.github.io/python-novice-gapminder/), de Software Carpentry, y [Data Analysis and Visualization in Python for Ecologists](https://datacarpentry.github.io/python-ecology-lesson/), de Data Carpentry (en español: [Análisis y visualización de datos usando Python](https://datacarpentry.github.io/python-ecology-lesson-es/)).

## Objetivos generales

Al terminar podrás:

1. **Construir** una figura interactiva de Plotly a partir de trazas y de un diseño.
2. **Agregar** información al pasar el cursor, un deslizador y un menú desplegable.
3. **Resumir** cientos de pruebas en una figura y ver el efecto de corregir por comparaciones múltiples.
4. **Reunir** varias figuras en una página HTML, guardarla y descargarla.

## Índice

«INDICE»

Son «TOTAL_EXP» de explicación y «TOTAL_EJ» min de ejercicios. Hay **«N_EJ» ejercicios** cortos, con la solución escondida.

Cada episodio tiene la misma organización que los talleres: **preguntas y objetivos**, **explicación paso a paso** con el código comentado, una pregunta **para pensar** mientras exploras la figura, un **ejercicio** y los **puntos clave**.

## Cómo usar este cuaderno

1. En Google Colab: **Archivo → Subir notebook**, o ábrelo desde GitHub con **Archivo → Abrir notebook → GitHub**.
2. Ejecuta las celdas en orden con **Shift + Enter**. Los datos reales se descargan solos desde GitHub; no hay que subir nada.
3. Se usan `numpy`, `pandas`, `scipy` y `plotly`, que ya vienen instalados en Colab.
4. La última celda descarga la página `explorador_pruebas_de_hipotesis.html`.

Los datos reales fueron facilitados por la empresa [Solena Ag](https://www.solena.ag) (2023) y están publicados en
[GitHub](https://github.com/CamilaSilva1995/Tesis_Maestria/tree/main/Taller_Practico/datos).
"""
for marca, valor in {"«TOTAL»": duracion(inicio), "«N_EP»": str(len(EPISODIOS)), "«INDICE»": INDICE,
                     "«TOTAL_EXP»": duracion(total_exp), "«TOTAL_EJ»": str(total_ej), "«N_EJ»": str(n_ejercicios),
                     "__COMPLETO__": COMPLETO, "__RAPIDO__": RAPIDO}.items():
    PORTADA = PORTADA.replace(marca, valor)

# Las celdas del cuaderno, en orden: portada, episodios (encabezado + contenido) y cierre
C = [nbf.v4.new_markdown_cell(PORTADA.strip())]
for i, ep in enumerate(EPISODIOS, 1):
    tiempo = f"Explicación: {ep['exp']} min" + (f" · Ejercicio: {ep['ej']} min" if ep["ej"] else "")
    C.append(nbf.v4.new_markdown_cell(
        f"---\n## Episodio {i}. {ep['titulo']}\n\n⏱️ **{tiempo}**\n\n"
        f"> **❓ Preguntas**\n{lista(ep['preguntas'])}\n>\n> **🎯 Objetivos**\n{lista(ep['objetivos'])}"))
    C.extend(ep["celdas"])
C.extend(CIERRE)

nb = nbf.v4.new_notebook()
# Metadatos del cuaderno: kernel de Python 3 y opciones de Colab (nombre y tabla de contenido visible)
nb.metadata = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
               "language_info": {"name": "python"},
               "colab": {"name": NOMBRE, "toc_visible": True, "provenance": []}}
# El cuaderno se escribe sin salidas; para guardarlas hay que ejecutarlo después (ver README.md)
nb.cells = C
nbf.write(nb, NOMBRE)

# Resumen para quien prepara la clase: de dónde sale el tiempo de cada episodio
print("cuaderno escrito con", len(C), "celdas")
print("ep  palabras  celdas  figuras  pausas  explicación  ejercicios")
for i, ep in enumerate(EPISODIOS, 1):
    print(f"{i:2d}  {ep['palabras']:8d}  {ep['codigo']:6d}  {ep['figuras']:7d}  {ep['pausas']:6d}  {ep['exp']:8d} min  {ep['ej']:7d} min")
print(f"total: {duracion(total_exp)} de explicación + {total_ej} min de ejercicios = {duracion(inicio)}")
