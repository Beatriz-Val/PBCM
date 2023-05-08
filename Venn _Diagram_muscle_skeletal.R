# Venn Diagram - Muscle Skeletal results from GSEA - Transcription Factor Target

# load ggvenn package
library("ggvenn")

#Open Sets

foxo1_02 <- read.csv("C:/Users/beatr/OneDrive - FCT NOVA/Ambiente de Trabalho/3º ano - 2022,2023/Projeto/Resultados/Final/Expression/FOXO3/V$FOXO1_02.csv", sep = ",", header = T, row.names = 1)
foxo3_01 <- read.csv("C:/Users/beatr/OneDrive - FCT NOVA/Ambiente de Trabalho/3º ano - 2022,2023/Projeto/Resultados/Final/Expression/FOXO3/V$FOXO3_01.csv", sep = ",", header = T, row.names = 1)
foxo4_02 <- read.csv("C:/Users/beatr/OneDrive - FCT NOVA/Ambiente de Trabalho/3º ano - 2022,2023/Projeto/Resultados/Final/Expression/FOXO3/V$FOXO4_02.csv", sep = ",", header = T, row.names = 1)
foxo1_01 <- read.csv("C:/Users/beatr/OneDrive - FCT NOVA/Ambiente de Trabalho/3º ano - 2022,2023/Projeto/Resultados/Final/Expression/FOXO3/V$FOXO1_01.csv", sep = ",", header = T, row.names = 1)

set1 <- foxo1_02["Gene.Symbol"]
set2 <- foxo3_01["Gene.Symbol"]
set3 <- foxo4_02["Gene.Symbol"]
set4 <- foxo1_01["Gene.Symbol"]

set1= set1[,1]
set2= set2[,1]
set3= set3[,1]
set4= set4[,1]

#Create list as input
H <-list("V$FOXO1_02" = set1, "V$FOXO3_01" = set2, "V$FOXO4_02" = set3, "V$FOXO1_01" = set4)

fill_colors <- c("red", "blue", "green", "orange")

ggvenn(H,c("V$FOXO1_02", "V$FOXO3_01", "V$FOXO4_02", "V$FOXO1_01"), show_percentage=F,
       fill_color=c("red","orange", "blue", "green"))



