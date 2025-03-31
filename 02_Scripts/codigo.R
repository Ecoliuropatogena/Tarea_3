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
tax_ta
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
estimate_richness(filosec.prune, measures = "Observed") 

# Shannon
  # Generamos una grafica separada por nacionalidad, bmi y sexo del indice de shannon
plot_richness(filosec.prune, x = "nationality", measures = "shannon")

  # Calculo de shannon
estimate_richness(filosec.prune, measures = "shannon")


# Shannon EXTRA
plot_richness(filosec.prune, x = "bmi_group", measures = "shannon")

plot_richness(filosec.prune, x = "sex", measures = "shannon")

# separados por sexo y peso
plot_richness(filosec.prune, measures="Shannon", x="sex", color="bmi_group")

# separados por nacionalidad y sexo
plot_richness(filosec.prune, measures="Shannon", x="nationality", color="bmi_group")
  


# Simpson
# Generamos una grafica separada por nacionalidad del indice de simpson
plot_richness(filosec.prune, x = "nationality", measures = "simpson")

# Calculo de simpson
estimate_richness(filosec.prune, measures = "simpson")


#--------------------------------------------------------------------#
                  # Filtrado y transformación
# Aplica un filtrado para quedarte solo con los géneros más 
# abundantes (por ejemplo, los que tienen más del 0.1% de 
# abundancia relativa en al menos 10% de las muestras).

# Tengo 222 muestras y 130 taxones
# Me estan pidiendo filtrar mis datos de acuerdo a los taxones más abundantes en 
# todas las muestras

# Primero generare un objeto que contenga la riqueza basal y haré comparaciones
filtro_basal <- prune_taxa(taxa_sums(filo_sec) > 1, filo_sec)
estimate_richness(filtro_basal, measures = "Observed") -> Riqueza_i # Riqueza inicial

filtro_abundancias <- prune_taxa(taxa_sums(filo_sec) > 300, filo_sec)
filtro_abundancias # Este objeto tiene las taxa que tienen una abundancia de al 
# menos 300 en todos los sample

estimate_richness(filtro_abundancias, measures = "Observed") -> Riqueza_f # Riqueza final

# Cálculo del cambio de la riqueza de acuerdo al filtrado por abundancias
Riqueza_i-Riqueza_f -> cambio # La riqueza de A debe ser más grande que la de B porque esta considerando más taxa

(cambio/Riqueza_i)*100 -> porcentaje_cambio
min(porcentaje_cambio) # Esto quiere decir que en Riqueza_f tego a los taxones que representan al menos
# un 4.9% de la riqueza en los sample
max(porcentaje_cambio) 


#----------------------------------------------------------------------------#

                            # Diversidad beta #
#Realiza una ordención PCoA utilizando distancia Bray-Curtis. Usa ordinate() y plot_ordination().
  #Responde:
#  ¿Los grupos se separan visiblemente?
#  ¿Qué podría estar causando esas diferencias?

BiocManager::install("MicrobiotaProcess")
library(MicrobiotaProcess)
library(vegan)
library(ggplot2)


# Para poder realizar el PCoA es necesario usar el objeto phylose
class(filo_sec)

# Primero usamos la funcion distance para obtener las distancias con el metodo
# Bray-curtis
ditancia_braycurtis <- distance(filo_sec, method = "bray")

# Ahora obtenemos las ordenadas
datos_ordenadas_pcoa <- ordinate(filo_sec, method = "PCoA", distance = ditancia_braycurtis)

# primero hacemos el grafico base generando los puntos y separandolos por color
# de acuerdo a las variables de los metadatos.
  # Aquí mismo se generan los ejes de x y y, seleccionado el punto en el
  # que intersectarán.
  # Ademas se añadiran los elipses
pcoa_plot<- plot_ordination(filo_sec, datos_ordenadas_pcoa, 
                             type = "sample", 
                             color = "group") +
  geom_point(size = 0.5) +  
  geom_hline(yintercept = 0, color = "darkgreen", linetype = 5) +
  geom_vline(xintercept = 0, color = "navy", linetype = 5) +
  stat_ellipse(aes(group = group), level = 0.8, color = "darkred", linetype = "dashed")

# Por ultimo se modificaran las etiquetas de los ejes y el titulo de la grafica
pcoa_plot +                         
  ggtitle("PCoA (Bray-Curtis)") +               
  labs(x = paste0("PCoA1 (", round(datos_ordenadas_pcoa$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(datos_ordenadas_pcoa$values$Relative_eig[2] * 100, 1), "%)")) +  # Etiquetas con porcentaje de varianza
  theme_bw() +                                 
  theme(plot.title = element_text(hjust = 0.5)) 




#----------------------------------------------------------------------------#
                            # Rank - abundance #
# calculo de de la abundancia de cad taxpn
abundancias_taxones <- taxa_sums(filo_sec)

# Abundancias ordenadas de mayor a menor
abundancias_ordenadas <- sort(abundancias_taxones, decreasing = TRUE)

# gráfica de barras rank-abundance

par(mgp = c(0, 1, -0.5))
barplot(abundancias_ordenadas, 
        main = "Rank-Abundance", 
        xlab = "Taxones (ordenados de mayor a menor)", 
        ylab = "Abundancia total", 
        col = "sienna", 
        las = 2,  
        cex.names = 0.5,
        cex.axis = 0.7)
#----------------------------------------------------------------------------#
                  # Gráficas apiladas de abundancia por taxón#

#Agrupa por phylum o género y grafica la composición de cada muestra como gráfica de barras apiladas.

# Podemos usar unicamente la funcion plot bar
plot_bar(filo_sec, fill = "Phylum")

# o complementarla
plot_bar(filo_sec, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

# Por grupo
plot_bar(filo_sec, fill = "group") + 
  geom_bar(aes(color=group), stat="identity", position="stack")



#----------------------------------------------------------------------------#
                        # Global patterns #

# Filtrar taxa con menos de 5 lecturas en al menos 20% 
# de las muestras 
# Aglomerar a nivel de Familia
# Transformar a abundancias relativas (%)
# Subset para incluir solo muestras de: Soil, Feces, Skin
# Entrega el código y muestra las dimensiones del objeto resultante

data("GlobalPatterns")


GP <- prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)

human <- get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")

sample_data(GP)$human <- factor(human)


cinco_lec_20<- filter_taxa(gp, function(x) sum(x > total*0.20) > 0, TRUE)


#----------------------------------------------------------------------------#
                        # Diversidad alfa #
# Calcular 3 índices de diversidad alfa (Shannon, Simpson, Observed)
indices_diversidad <-plot_richness(GP, x = "Primer", measures = c("Observed", "shannon", "simpson"), color = "SampleType")
indices_diversidad



# Crear boxplots comparativos de los índices entre tipos de muestra
alpha_meas = c("Observed","Shannon", "Simpson")
p <- plot_richness(GP, "human", "SampleType", measures=alpha_meas)

p + geom_boxplot(data=p$data, aes(x=human, y=value, color=NULL), alpha=0.1)


#----------------------------------------------------------------------------#
                    # Curvas de Rango-Abundancia #
# Crear gráficas de rango-abundancia para cada tipo de muestra 

soil_subset <- subset_samples(GP, SampleType == "Soil")
ab_soils <- taxa_sums(soil_subset)
ab_soils_ord  <-sort(ab_soils, decreasing = TRUE)

Feces_subset <- subset_samples(GP, SampleType == "Feces")
ab_feces <- taxa_sums(Feces_subset)
ab_feces_ord <- sort(ab_feces, decreasing = TRUE)

Skin_subset <- subset_samples(GP, SampleType == "Skin")
ab_skin <- taxa_sums(Skin_subset)
ab_skin_ord <- sort(ab_skin, decreasing = TRUE)

Tongue_subset <- subset_samples(GP, SampleType == "Tongue")
ab_tongue <- taxa_sums(Tongue_subset)
ab_tongue_ord <- sort(ab_tongue, decreasing = TRUE)

Freshwater_subset <- subset_samples(GP, SampleType == "Freshwater")
ab_fw <- taxa_sums(Freshwater_subset)
ab_fw_ord <- sort(ab_fw, decreasing = TRUE)

Freshwater_creek_subset <- subset_samples(GP, SampleType == "Freshwater (creek)")
ab_fw_cr <- taxa_sums(Freshwater_creek_subset)
ab_fw_cr_ord <- sort(ab_fw_cr, decreasing = TRUE)

Ocean_subset <- subset_samples(GP, SampleType == "Ocean")
ab_ocean <- taxa_sums(Ocean_subset)
ab_ocean_ord <- sort(ab_ocean, decreasing = TRUE)

Mock_subset <- subset_samples(GP, SampleType == "Mock")
ab_mock <- taxa_sums(Mock_subset)
ab_mock_ord <- sort(ab_mock, decreasing = TRUE)

par(mgp = c(0, 1, -0.5))
barplot(ab_soils_ord[1:200], 
        main = "Rank-Abundance", 
        xlab = "Taxones (ordenados de mayor a menor)", 
        ylab = "Abundancia total", 
        col = "lightseagreen", 
        las = 2,  
        cex.names = 0.5,
        cex.axis = 0.7,
        width = 1)


        





#----------------------------------------------------------------------------#
                         # Perfil taxonómico #
# Crear gráfico apilado de abundancia a nivel de Phylum
plot_bar(cinco, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

# Mostrar solo los 5 phyla más abundantes
GP_phylum <- tax_glom(GP, taxrank = "Phylum")
phylum_ab <- taxa_sums(GP_phylum)
cinco <- names(sort(phylum_ab, decreasing = TRUE))[1:5]
seleccion <- prune_taxa(taxa_names(GP_phylum) %in% cinco, GP_phylum)
plot_bar(seleccion, fill = "Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")


# comentar resultados biológicos





#----------------------------------------------------------------------------#
                        # Diversidad Beta #

# Calcular distancia Bray-Curtis
# Realizar PCoA
# Visualizar con:
#   Colores por tipo de muestra
#   Elipses de confianza del 95%
#   Incluir stress plot
#   Realizar PERMANOVA para diferencias entre grupos
#   Interpretar resultados en contexto ecológico
























