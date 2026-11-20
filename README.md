# Anti-dsDNA: CLIA frente a radioinmunoanálisis (Farr)

**pMenor · Estadística clínica en R**

Código del análisis estadístico del artículo:

> Martín-Martínez MD, García-de la Rosa G, De Castro-Córdova J. Comparative
> evaluation of chemiluminescent immunoassay and radioimmunoassay for the
> detection of anti-double-stranded DNA antibodies. *J Rheum Dis.* 2025.
> [doi:10.4078/jrd.2025.0111](https://doi.org/10.4078/jrd.2025.0111)

## La pregunta

Muchos laboratorios están dejando el radioinmunoanálisis de Farr (RIA), el
método de referencia para los anticuerpos anti-dsDNA, por motivos logísticos y
de seguridad. ¿Puede un inmunoensayo quimioluminiscente (CLIA) sustituirlo sin
perder rendimiento diagnóstico?

## El análisis

807 muestras medidas con los dos métodos. El análisis tiene dos partes:

| Pregunta | Método |
|---|---|
| ¿Dan los dos métodos el mismo valor? | Spearman, Bland-Altman y regresión Tobit (por los valores censurados) |
| ¿Clasifican igual positivo/negativo? | Tabla 2×2, kappa de Cohen y PABAK |
| ¿Cómo rinde el CLIA frente al RIA? | Sensibilidad, especificidad, valores predictivos y cocientes de probabilidad, con IC 95 % |
| ¿Mejora con otro punto de corte? | Curva ROC e índice de Youden |
| ¿Depende el desacuerdo de la concentración? | Bland-Altman en porcentaje (análisis de sensibilidad, no publicado) |

Dos detalles que hacen interesante este caso:

- **Dos tercios de los RIA están en el límite de cuantificación** (< 2,5 UI/mL).
  Con tantos empates, la regresión de Passing-Bablok no es aplicable y hay que
  tratar la censura explícitamente.
- **Con pocos positivos, el kappa engaña.** Con un 7 % de positivos, los
  métodos coinciden en el 91 % de las muestras, pero el kappa es solo 0,37. El
  PABAK (0,81) muestra cuánto pesa la prevalencia.

### Resultados publicados

| Medida | Resultado |
|---|---|
| Correlación de Spearman | ρ = 0,108 (p = 0,002) |
| Bland-Altman | sesgo 14,2 UI/mL; límites −87,1 a 115,4 UI/mL |
| Kappa / PABAK | 0,365 / 0,812 |
| Sensibilidad / especificidad (corte 29 UI/mL) | 45 % / 94 % |
| AUC | 0,719 |
| Corte de Youden | 19,15 UI/mL (sensibilidad 57 %, especificidad 86 %) |

**Conclusión del artículo:** el CLIA no puede sustituir al RIA sin pérdida de
rendimiento diagnóstico.

## Datos

Este repositorio **no contiene datos de pacientes**. `R/simular_datos.R` genera
un dataset ficticio con la misma estructura (`data/antidsdna_simulado.csv`).
La tabla 2×2, y con ella el kappa, el PABAK y el rendimiento diagnóstico,
coinciden con el artículo por construcción. El resto de cifras se parecen pero
no son las publicadas.

## Cómo reproducirlo

1. Clona el repositorio y abre la carpeta como proyecto de RStudio.
2. Instala los paquetes:

   ```r
   install.packages(c("dplyr", "tidyr", "ggplot2", "pROC", "epiR",
                      "AER", "pwr", "rmarkdown", "scales", "svglite", "systemfonts"))
   ```

3. Ejecuta `rmarkdown::render("analisis.Rmd")`.

El resultado, con código, salidas y figuras, se puede leer directamente en
[`analisis.md`](analisis.md). Las figuras usan las fuentes IBM Plex Sans e IBM
Plex Mono; si no están instaladas, R usa otras y avisa.

## Estructura

```
├── analisis.Rmd          análisis completo
├── analisis.md           el mismo análisis ya ejecutado (se lee en GitHub)
├── R/
│   ├── simular_datos.R   genera los datos ficticios
│   └── theme_pmenor.R    estilo de los gráficos
├── data/                 datos simulados
└── output/               figuras
```

## Licencia

El código se publica con licencia [MIT](LICENSE): puedes usarlo y adaptarlo citando la fuente. Si lo usas en un trabajo, cita también el artículo.

---

**pMenor** · Estadística clínica en R · [pmenor.es](https://pmenor.es)
