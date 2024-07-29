library(dplyr)
library(fitdistrplus)

# Argumentos de entrada(conteos de referencia y muestra problema)
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Se requieren 3 argumentos: <archivo_conteos_referencia> <archivo_conteos_muestra_problema> <archivo_output> ")
} 

archivo_conteos_referencia = args[1]
archivo_conteos_muestra_problema = args[2]
output = args[3]

########################### Ordena datos (REFERENCIA) ##########################
df_conteos = read.csv(archivo_conteos_referencia, sep="\t", header=TRUE)
colnames(df_conteos)[-1] <- gsub("X", "", colnames(df_conteos)[-1])
total_anotaciones = df_conteos[, ncol(df_conteos)]
rownames(df_conteos) = df_conteos$Muestra 
df_conteos = df_conteos[-1]
df_conteos = df_conteos[, -ncol(df_conteos)]
df_rates = df_conteos/total_anotaciones * 1e6
distribuciones = c("normal", "gamma", "lognormal", "weibull")
ec_numbers = colnames(df_rates)
######################## Ordena datos (Muestra Problema) #######################
df_muestra_problema = read.csv(archivo_conteos_muestra_problema,
                               sep = "\t", header = TRUE)
colnames(df_muestra_problema)[-1] <- gsub("X", "", colnames(df_muestra_problema)[-1])
total_anotaciones_problema = df_muestra_problema[, ncol(df_muestra_problema)]
rownames(df_muestra_problema) = df_muestra_problema$Muestra
df_muestra_problema = df_muestra_problema[-1]
df_muestra_problema = df_muestra_problema[, -ncol(df_muestra_problema)]
df_rates_muestra_problema = df_muestra_problema/total_anotaciones_problema * 1e6
################### Ajusta distribuciones de todas las muestras ################
resultados_ajuste = list()
for (nombre_columna in ec_numbers) {
  col = df_rates[[nombre_columna]] 
  ajustes = lapply(distribuciones, function(dist) {
    tryCatch(
      fitdistr(col, dist),
      error = function(e) NULL
    )
  })
  resultados_ajuste[[nombre_columna]] = setNames(ajustes, distribuciones)
}
#resultados_ajuste
##################### Calcular los valores de AIC ##############################
lista_aic <- list()

for (variable in names(resultados_ajuste)) {
  ajustes <- resultados_ajuste[[variable]]
  aic_distribuciones <- list()
  
  # Verificar si hay ajustes disponibles para la variable actual
  if (length(ajustes) > 0) {
    # Iterar sobre cada distribución ajustada para la variable actual
    for (dist in names(ajustes)) {
      ajuste <- ajustes[[dist]]
      # Verificar si el ajuste no es nulo
      if (!is.null(ajuste)) {
        # Calcular el AIC y almacenarlo en la lista aic_distribuciones junto con el nombre de la distribución
        aic_value <- AIC(ajuste)
        aic_distribuciones[[dist]] <- aic_value
      }
    }
  }
  
  # Agregar los AIC de las distribuciones de la variable actual a la lista
  lista_aic[[variable]] <- aic_distribuciones
}
#lista_aic
###############################################################################
for (data in names(resultados_ajuste)) {
  datos_lista <- lista_aic[[data]]
  deslista = unlist(datos_lista)
  deslista = sort(deslista)
  #print(deslista[1])
}
###############################################################################
df_resultados = data.frame(variable = character(), AIC_min = numeric(),
                           distribucion_min = character(),
                           stringsAsFactors = FALSE)
for (data in names(resultados_ajuste)) {
  datos_lista <- lista_aic[[data]]
  if (length(datos_lista) > 0) {
    # Unlist y ordena los valores de AIC
    deslista <- sort(unlist(datos_lista))
    # Encuentra el valor de AIC más bajo
    aic_min <- deslista[1]
    # Encuentra la distribución correspondiente al valor de AIC más bajo
    for (dist in names(datos_lista)) {
      if (datos_lista[[dist]] == aic_min) {
        distribucion_min <- dist
        break
      }
    }
    df_resultados = rbind(df_resultados, data.frame(variable = data, 
                                                    AIC_min = aic_min,
                                                    distribucion_min = distribucion_min,
                                                    row.names = FALSE))
  } 
}
#df_resultados

rownames(df_resultados) <- df_resultados$variable
df_resultados$variable <- NULL
#print(df_resultados)

df_transpuesto <- t(df_resultados)
segunda_fila <- df_transpuesto[2, , drop = FALSE]
df_combinado <- rbind(df_rates, segunda_fila)

#write.csv(df_combinado, file = "/Users/brandonbuenohernandez/Documents/R_scipts/dist_ref3.csv",
#          sep = "\t")
##############################################################################

resultados_comparacion <- data.frame()

for (enzima in rownames(df_resultados)) {
  distribucion = df_resultados[enzima, "distribucion_min"]
  parametros = resultados_ajuste[[enzima]][[distribucion]]$estimate
  
  if (enzima %in% colnames(df_rates_muestra_problema)) {
    for (i in seq_along(df_rates_muestra_problema[[enzima]])) {
      valor = df_rates_muestra_problema[[enzima]][i]
      muestra = rownames(df_rates_muestra_problema)[i]
      
      if (!is.na(valor) && is.numeric(valor)) {
        if (distribucion == "normal") {
          mean = parametros["mean"]
          sd = parametros["sd"]
          p_value = pnorm(valor, mean, sd)
          z_score = qnorm(p_value)
        } else if (distribucion == "gamma") {
          shape = parametros["shape"]
          rate = parametros["rate"]
          p_value = pgamma(valor, shape, rate = rate)
          z_score = qnorm(p_value)
        } else if (distribucion == "lognormal") {
          meanlog = parametros["meanlog"]
          sdlog = parametros["sdlog"]
          p_value = plnorm(valor, meanlog, sdlog)
          z_score = qnorm(p_value)
        } else if (distribucion == "weibull") {
          shape = parametros["shape"]
          scale = parametros["scale"]
          p_value = pweibull(valor, shape, scale)
          z_score = qnorm(p_value)
        }
        
        resultados_comparacion <- rbind(resultados_comparacion, 
                                        data.frame(Muestra = muestra, 
                                                   Enzima = enzima, 
                                                   Valor = valor, 
                                                   Distribucion = distribucion, 
                                                   Z_Score = z_score))
        cat(muestra, enzima, "Rate:", valor, "Distribución:", distribucion, "Z-Score:", z_score, "\n")
      }
    }
  } else {
    for (muestra in rownames(df_rates_muestra_problema)) {
      resultados_comparacion <- rbind(resultados_comparacion, 
                                      data.frame(Muestra = muestra, 
                                                 Enzima = enzima, 
                                                 Valor = NA, 
                                                 Distribucion = distribucion, 
                                                 Z_Score = NA))
      cat(muestra, enzima, "Rate: NA Distribución:", distribucion, "Z-Score: NA\n")
    }
  }
}

resultados_comparacion
write.csv(resultados_comparacion,
          file = output,
          row.names = "Muestra")
