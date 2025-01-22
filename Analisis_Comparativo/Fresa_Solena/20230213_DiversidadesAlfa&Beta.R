### AHORA QUEREMOS VER DIFERENCIAS ENTRE DIFERENTES NIVELES TAXONOMICOS

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data),1,6)
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1/metadata.csv",header =  FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
fresa_kraken@sam_data$Sample<-row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data)<-c('Treatment','Samples')
samples_to_remove <- c("MP2079","MP2080","MP2088","MP2109","MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
percentages_df <- psmelt(percentages_fil)
## De aqui en adelante solo trabajaremos con los datos filtrados, 

## Grafico de barras de abundancia 
## Trazando las muestras en el eje x y las abundancias en el eje y
## Los valores de abundancia para cada OTUen cada muestra se apilan en el orden de mayor a menor, separados por una fina línea horizontal.
plot_bar(fresa_kraken_fil,fill="Treatment")
ggsave("FresaKraken_StackBar.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
plot_bar(percentages_fil,fill="Treatment")

ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='Phylum' ,color='Treatment'))  +
  scale_colour_manual(values=c('white','black')) +
  geom_bar(aes(), stat="identity", position="stack") +
  #scale_x_discrete(limits = rev(levels(percentages_df$Category))) +
  labs(title = "Abundance", x='Sample', y='Abundance', color = 'Treatment') +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=8, face = "bold"),
        legend.text=element_text(size=6),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
ggsave("FresaKraken_fil_StackBar.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")


## Queremos explorar nuestras muestras a diferentes niveles taxonómicos específicos  
## A nivel de Kingdom
percentages_glom_kingdom <- tax_glom(percentages_fil, taxrank = 'Kingdom')
#View(percentages_glom_kingdom@tax_table@.Data)
## A nivel de Phylum
percentages_glom_phylum <- tax_glom(percentages_fil, taxrank = 'Phylum')
#View(percentages_glom_phylum@tax_table@.Data)
percentages_df_phylum <- psmelt(percentages_glom_phylum)
str(percentages_df_phylum)

absolute_glom_phylum <- tax_glom(physeq = fresa_kraken_fil, taxrank = "Phylum")
absolute_df_phylum <- psmelt(absolute_glom_phylum)
str(absolute_df_phylum)

absolute_df_phylum$Phylum <- as.factor(absolute_df_phylum$Phylum)
phylum_colors_abs <- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(absolute_df_phylum$Phylum)))

absolute_plot <- ggplot(data= absolute_df_phylum, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors_abs) +
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

percentages_df_phylum$Phylum <- as.factor(percentages_df_phylum$Phylum)
phylum_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(percentages_df_phylum$Phylum)))
relative_plot <- ggplot(data=percentages_df_phylum, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = phylum_colors_rel) +
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

absolute_plot
ggsave("absolute_plot.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
relative_plot
ggsave("relative_plot.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")


## Usaremos un comando llamado unique() para explorar cuántos filos y reinos tenemos. 
unique(fresa_kraken_fil@tax_table@.Data[,"Kingdom"])
unique(fresa_kraken_fil@tax_table@.Data[,"Phylum"])
## con esto podemos ver cuantos "Eukaryota" tenemos en "Kingdom"
sum(fresa_kraken_fil@tax_table@.Data[,"Kingdom"] == "Eukaryota")
## y cuantos "Bacteria"
sum(fresa_kraken_fil@tax_table@.Data[,"Kingdom"] == "Bacteria")

## Diversidad Beta

## veremos aqui solo el reino Eucariota
merge_Eukaryota<-subset_taxa(fresa_kraken_fil,Kingdom=="Eukaryota")
## sacamos  las abundancias relativas 
percentages_Eukaryota <- transform_sample_counts(merge_Eukaryota, function(x) x*100 / sum(x) )
percentages_Eukaryota_df <- psmelt(percentages_Eukaryota)
## beta diversidad de Eukaryota
meta_ord_Eukaryota <- ordinate(physeq = percentages_Eukaryota, method = "NMDS", distance = "bray")  
plot_ordination(physeq = percentages_Eukaryota, ordination = meta_ord_Eukaryota, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)

## Eukaryota por Phylum
ggplot(data= percentages_Eukaryota_df, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack") +
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

merge_Eukaryota_Phylum<-tax_glom(merge_Eukaryota,taxrank="Phylum")
## sacamos  las abundancias relativas 
percentages_Eukaryota_Phylum <- transform_sample_counts(merge_Eukaryota_Phylum, function(x) x*100 / sum(x) )
percentages_Eukaryota_Phylum_df <- psmelt(percentages_Eukaryota_Phylum)
meta_ord_Eukaryota_Phylum <- ordinate(physeq = percentages_Eukaryota_Phylum, method = "NMDS", distance = "bray")  
plot_ordination(physeq = percentages_Eukaryota_Phylum, ordination = meta_ord_Eukaryota_Phylum, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)

## Eukaryota por Family
ggplot(data= percentages_Eukaryota_df, aes(x=Sample, y=Abundance, fill=Family))+ 
  geom_bar(aes(), stat="identity", position="stack") +
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

merge_Eukaryota_Family<-tax_glom(merge_Eukaryota,taxrank="Family")
## sacamos  las abundancias relativas 
percentages_Eukaryota_Family <- transform_sample_counts(merge_Eukaryota_Family, function(x) x*100 / sum(x) )
percentages_Eukaryota_Family_df <- psmelt(percentages_Eukaryota_Family)
meta_ord_Eukaryota_Family <- ordinate(physeq = percentages_Eukaryota_Family, method = "NMDS", distance = "bray")  
plot_ordination(physeq = percentages_Eukaryota_Family, ordination = meta_ord_Eukaryota_Family, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)

## Eukaryota por Genus
ggplot(data= percentages_Eukaryota_df, aes(x=Sample, y=Abundance, fill=Genus))+ 
         geom_bar(aes(), stat="identity", position="stack") +
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


merge_Eukaryota_Genus<-tax_glom(merge_Eukaryota,taxrank="Genus")
## sacamos  las abundancias relativas 
percentages_Eukaryota_Genus <- transform_sample_counts(merge_Eukaryota_Genus, function(x) x*100 / sum(x) )
percentages_Eukaryota_Genus_df <- psmelt(percentages_Eukaryota_Genus)
meta_ord_Eukaryota_Genus <- ordinate(physeq = percentages_Eukaryota_Genus, method = "NMDS", distance = "bray")  
plot_ordination(physeq = percentages_Eukaryota_Genus, ordination = meta_ord_Eukaryota_Genus, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)

## tomando diferentes porcentajes de abundancia
percentages_Eukaryota_Genus_df$Genus[percentages_Eukaryota_Genus_df$Abundance < 10.0] <- "Genus < 0.5% abund." #no deja poner ningun numero diferente a 0.5%
percentages_Eukaryota_Genus_df$Genus <- as.factor(percentages_Eukaryota_Genus_df$Genus)

genus_colors_rel <- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(percentages_Eukaryota_Genus_df$Genus)))
relative_plot <- ggplot(data=percentages_Eukaryota_Genus_df, aes(x=Sample, y=Abundance, fill=Genus ,color = Treatment))+ 
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = phylum_colors_rel) +
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
relative_plot

relative_plot <- ggplot(data=percentages_Eukaryota_Genus_df, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = phylum_colors_rel) +
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
relative_plot

## sacando la beta diversidad con aglomerado de 10%
meta_ord_Eukaryota_Genus <- ordinate(physeq = percentages_Eukaryota_Genus, method = "NMDS", distance = "bray")  
plot_ordination(physeq = percentages_Eukaryota_Genus, ordination = meta_ord_Eukaryota_Genus, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)






## Veremos aqui solo el reino Bacteriano
merge_Bacteria<-subset_taxa(fresa_kraken_fil,Kingdom=="Bacteria")
## sacamos  las abundancias relativas 
percentages_Bacteria <- transform_sample_counts(merge_Bacteria, function(x) x*100 / sum(x) )
percentages_Bacteria_df <- psmelt(percentages_Bacteria)
## beta diversidad de Bacteria
meta_ord_Bacteria <- ordinate(physeq = percentages_Bacteria, method = "NMDS", distance = "bray")   
plot_ordination(physeq = percentages_Bacteria, ordination = meta_ord_Bacteria, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)

## Bacterias por Phylum
ggplot(data= percentages_Bacteria_df, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack") +
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

merge_Bacteria_Phylum<-tax_glom(merge_Bacteria,taxrank="Phylum")
## sacamos  las abundancias relativas 
percentages_Bacteria_Phylum <- transform_sample_counts(merge_Bacteria_Phylum, function(x) x*100 / sum(x) )
percentages_Bacteria_Phylum_df <- psmelt(percentages_Bacteria_Phylum)
meta_ord_Bacteria_Phylum <- ordinate(physeq = percentages_Bacteria_Phylum, method = "NMDS", distance = "bray")  
plot_ordination(physeq = percentages_Bacteria_Phylum, ordination = meta_ord_Bacteria_Phylum, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)

## Bacterias por Family
ggplot(data= percentages_Bacteria_df, aes(x=Sample, y=Abundance, fill=Family))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

merge_Bacteria_Family<-tax_glom(merge_Bacteria,taxrank="Family")
## sacamos  las abundancias relativas 
percentages_Bacteria_Family <- transform_sample_counts(merge_Bacteria_Family, function(x) x*100 / sum(x) )
percentages_Bacteria_Family_df <- psmelt(percentages_Bacteria_Family)
meta_ord_Bacteria_Family <- ordinate(physeq = percentages_Bacteria_Family, method = "NMDS", distance = "bray")  
plot_ordination(physeq = percentages_Bacteria_Family, ordination = meta_ord_Bacteria_Family, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)

## Bacterias por Genero
ggplot(data= percentages_Bacteria_df, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack") +
  theme(text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

merge_Bacteria_Genus<-tax_glom(merge_Bacteria,taxrank="Genus")

## sacamos  las abundancias relativas 
percentages_Bacteria_Genus <- transform_sample_counts(merge_Bacteria_Genus, function(x) x*100 / sum(x) )
percentages_Bacteria_Genus_df <- psmelt(percentages_Bacteria_Genus)
meta_ord_Bacteria_Genus <- ordinate(physeq = percentages_Bacteria_Genus, method = "NMDS", distance = "bray")  
plot_ordination(physeq = percentages_Bacteria_Genus, ordination = meta_ord_Bacteria_Genus, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)


















