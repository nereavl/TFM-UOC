
# Rscript para limpiar el output de breakdancer

# Arguments:

args <- commandArgs(trailingOnly=TRUE)
strain <- args[1]  
chromosome <- args[2]
min_length <- as.numeric(args[3])
min_quality <- as.numeric(args[4])

suppressMessages(library(StructuralVariantAnnotation))
suppressMessages(library(dplyr))

theoretical_coords <- read.table("inversions_theory.bed") 

vcf <- readVcf(paste("./seqs/4_Manta/", chromosome, "/results/variants/diploidSV.vcf", sep=""))
bpgr <- breakpointRanges(vcf, nominalPosition=FALSE)
				   
bedpe <- breakpointgr2bedpe(bpgr)

# Print some interesting basic data
paste("Se detectan",nrow(bedpe), "variantes cromosómicas")
paste("Se detectan",sum(simpleEventType(bpgr)=="INV")/2, "inversiones")

# Filter the inversions
filter_BD <- function(BD, coords){
  my_coords <- coords[coords$V1 %in% BD$chrom1,]
  
  matches1 <- data.frame(matrix(NA, nrow=nrow(BD), ncol=nrow(my_coords)))
  matches2 <- data.frame(matrix(NA, nrow=nrow(BD), ncol=nrow(my_coords)))

  for (i in 1:nrow(my_coords)){
  
    matches1[i] <- (BD$start1+BD$end1)%/%2 %in% my_coords[i,2]:my_coords[i,3]
	matches2[i] <- (BD$start2+BD$end2)%/%2 %in% my_coords[i,2]:my_coords[i,3]
  
  }
  
  BD_filtered <- BD[rowSums(matches1) > 0 & rowSums(matches2) > 0 &
					abs((BD$start2+BD$end2)%/%2 - (BD$start1+BD$end1)%/%2) > abs(min_length),]  
  BD_filtered <- BD_filtered %>% select(-c("chrom2"))  

  return(BD_filtered)
} 

INV_data <- filter_BD(bedpe, theoretical_coords)

paste("Se seleccionan", nrow(INV_data), "variantes en las coordenadas esperadas")

# Print the data

write.table(INV_data, paste("./seqs/5_Clean_Manta/", strain, "_SVclean_", chromosome, ".txt", sep=""))
write.table(bedpe, paste("./seqs/5_Clean_Manta/", strain, "_SVall_", chromosome, ".txt", sep=""))