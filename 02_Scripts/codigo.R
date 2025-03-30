                        # CODIGO #

# Librerias a utilizar
library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(dplyr)

# Generar un objeto phyloseq a partir de datos
data("dietswap", package = "microbiome")
filo_sec <- dietswap
filo_sec

View(otu_table(filo_sec))

otu_table(filo_sec) -> abundancias

View(sample_data(filo_sec))

View(tax_table(filo_sec))
#-------------------------------------------------------------------#
                      # Curvas de rarefacción #



# Generar una curva de rarefacción de las abundancias
# Para que rarecurve funcione es necesario que la base de datos
# sea una matriz o un data frame, es por ello que converti la otu table
# a un data frame.
data_abundancias <- as.data.frame(abundancias)
str(data_abundancias)
View(data_abundancias)

# Dato que rare curve nececita que los sample esten en los renglones
# transpuse la tabla del data frame
t(data_abundancias) -> t_data_abundancias
View(t_data_abundancias)

# Le pedi un vector con colores obscuros a chat gpt para que cada linea se 
# vea distinta
colores_chat <- c("black", "darkblue", "darkcyan", "darkgoldenrod",
                  "darkgray", "darkgreen", "darkgrey", "darkkhaki",
                  "darkmagenta", "darkolivegreen", "darkorange", "darkorchid",
                  "darkred", "darksalmon", "darkseagreen", "darkslateblue",
                  "darkslategray", "darkslategrey", "darkturquoise", "darkviolet",
                  "midnightblue", "navy","brown","deepskyblue", "skyblue")


# Repetimos esta linea 9 veces, modificando de que renglon a que renglon va.
rarecurve(t_data_abundancias[1:25, ], 
          step = 25, 
          xlab = "Tamaño muestral", 
          ylab = "Especies",
          tidy = FALSE,
          col = colores_chat,  
          label = FALSE) 

rarecurve(t_data_abundancias[26:50, ], 
          step = 25, 
          xlab = "Tamaño muestral", 
          ylab = "Especies",
          tidy = FALSE,
          col = colores_chat,  
          label = FALSE) 

rarecurve(t_data_abundancias[51:75, ], 
          step = 25, 
          xlab = "Tamaño muestral", 
          ylab = "Especies",
          tidy = FALSE,
          col = colores_chat,  
          label = FALSE) 

rarecurve(t_data_abundancias[76:100, ], 
          step = 25, 
          xlab = "Tamaño muestral", 
          ylab = "Especies",
          tidy = FALSE,
          col = colores_chat,  
          label = FALSE) 

rarecurve(t_data_abundancias[100:125, ], 
          step = 25, 
          xlab = "Tamaño muestral", 
          ylab = "Especies",
          tidy = FALSE,
          col = colores_chat,  
          label = FALSE) 

rarecurve(t_data_abundancias[126:150, ], 
          step = 25, 
          xlab = "Tamaño muestral", 
          ylab = "Especies",
          tidy = FALSE,
          col = colores_chat,  
          label = FALSE) 

rarecurve(t_data_abundancias[151:175, ], 
          step = 25, 
          xlab = "Tamaño muestral", 
          ylab = "Especies",
          tidy = FALSE,
          col = colores_chat,  
          label = FALSE) 

rarecurve(t_data_abundancias[176:200, ], 
          step = 25, 
          xlab = "Tamaño muestral", 
          ylab = "Especies",
          tidy = FALSE,
          col = colores_chat,  
          label = FALSE) 

rarecurve(t_data_abundancias[201:222, ], 
          step = 25, 
          xlab = "Tamaño muestral", 
          ylab = "Especies",
          tidy = FALSE,
          col = colores_chat,  
          label = FALSE)


# Ciclo for para determinar la muestra que no llego a la asintota

t_data_abundancias[51:75,] -> rango 
View(rango)

sum(rango[1,])
sum(rango[1,1:5])
length(col(rango[[1]]))

x <- 50
for (i in 1:25) {
  sum(rango[i,]) -> tamaño_muestral
  x+i -> muestra
  print(paste0("El tamaño muestral es:",tamaño_muestral,
                " para la muestra ",muestra))
} # Se debe eliminar la mustra nnumero 56


#-------------------------------------------------------------------#
                     # Diversidad alfa (𝛼) #

Calcula y grafica los siguientes índices de diversidad alfa:

# Usa plot_richness() de phyloseq.  

  # Creamos un objeto que unicamente tome en cuenta las taxa mayores a 1   
filosec.prune <- prune_taxa(taxa_sums(filo_sec) > 1, filo_sec)

  
# Observed (Riqueza)
plot_richness(filosec.prune, x = "nationality", measures = "Observed", color = "sex")

# Shannon
  # Generamos una grafica separada por nacionalidad, bmi y sexo del indice de shannon
plot_richness(filosec.prune, x = "nationality", measures = "shannon")

plot_richness(filosec.prune, x = "bmi_group", measures = "shannon")

plot_richness(filosec.prune, x = "sex", measures = "shannon")


# separados por sexo y peso
plot_richness(filosec.prune, measures="Shannon", x="sex", color="bmi_group")

# separados por nacionalidad y sexo
plot_richness(filosec.prune, measures="Shannon", x="nationality", color="bmi_group")
  


# Simpson
# Generamos una grafica separada por nacionalidad del indice de simpson
plot_richness(filosec.prune, x = "nationality", measures = "simpson")









