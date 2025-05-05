library(bio3d)
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/01_structure_contacts/data/pdbs/raw/20230923/")

# ---
prep_file <- read.pdb("8ic0_prep1.pdb")
clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = FALSE, fix.chain = FALSE, rm.h = TRUE)
write_file <- write.pdb(clean_file, file = "8ic0_clean.pdb")
rm(prep_file, clean_file, write_file)



