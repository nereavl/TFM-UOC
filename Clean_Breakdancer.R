
# Rscript para limpiar el output de breakdancer

# Arguments:

args <- commandArgs(trailingOnly=TRUE)
strain <- args[1]  
chromosome <- args[2]
min_length <- as.numeric(args[3])
min_quality <- as.numeric(args[4])

suppressMessages(library(dplyr))

theoretical_coords <- read.table("inversions_theory.bed") 

raw_data <- read.table(paste("./seqs/4_BreakDancer/OF74_breakdancer_", chromosome, ".txt", sep=""), sep="\t", comment.char="#", col.names=c("chr1", "pos1", "strand1", "chr2", "pos2", "strand2", "Tipo_SV", "Longitud", "score", "n_Reads", "num_Reads_orig", "Valor"), as.is=TRUE)

# Print some interesting basic data
paste("Se detectan",nrow(raw_data), "variantes cromosómicas")
paste("Se detectan",sum(raw_data$Tipo_SV == "INV"), "inversiones")

filter_BD <- function(BD, coords){
  my_coords <- coords[coords$V1 %in% BD$chr1,]
  
  matches1 <- data.frame(matrix(NA, nrow=nrow(BD), ncol=nrow(my_coords)))
  matches2 <- data.frame(matrix(NA, nrow=nrow(BD), ncol=nrow(my_coords)))

  for (i in 1:nrow(my_coords)){
    matches1[i] <- BD$pos1 %in% my_coords[i,2]:my_coords[i,3]
	matches2[i] <- BD$pos2 %in% my_coords[i,2]:my_coords[i,3]
  }
  
  BD_filtered <- BD[rowSums(matches1) > 0 & rowSums(matches2) > 0 &
					abs(BD$pos1-BD$pos2) > min_length,]  
  BD_filtered <- BD_filtered %>% select(-c("chr2","num_Reads_orig","Valor"))
  BD_filtered <- BD_filtered %>% filter(score > min_quality)

  return(BD_filtered)
} 

INV_data <- filter_BD(raw_data, theoretical_coords)

paste("Se seleccionan", nrow(INV_data), "variantes en las coordenadas esperadas")

# Print the data

write.table(INV_data, paste("./seqs/5_Clean_BreakDancer/", strain, "_SVclean_", chromosome, ".txt", sep=""))
write.table(raw_data, paste("./seqs/5_Clean_BreakDancer/", strain, "_SVall_", chromosome, ".txt", sep=""))