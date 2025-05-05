library(bio3d)
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/01_structure_contacts/data/pdbs/raw/20220718/")

# Some PDB files did not need "cleaned" as they were already free from hydrogen
# and chemokine numbered as chain A and receptor as chain B
# # ---
# prep_file <- read.pdb("7f1r_prep.pdb")
# clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = FALSE, fix.chain = FALSE, rm.h = TRUE)
# write_file <- write.pdb(clean_file, file = "7f1r_clean.pdb")
# rm(prep_file, clean_file, write_file)
# 
# # ---
# prep_file <- read.pdb("7f1q_prep.pdb")
# clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = FALSE, fix.chain = FALSE, rm.h = TRUE)
# write_file <- write.pdb(clean_file, file = "7f1q_clean.pdb")
# rm(prep_file, clean_file, write_file)
# 
# # ---
# prep_file <- read.pdb("7xbx_prep2.pdb")
# clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = FALSE, fix.chain = FALSE, rm.h = TRUE)
# write_file <- write.pdb(clean_file, file = "7xbx_clean.pdb")
# rm(prep_file, clean_file, write_file)
# 
# # ---
# prep_file <- read.pdb("7sk3_prep2.pdb")
# clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = FALSE, fix.chain = FALSE, rm.h = TRUE)
# write_file <- write.pdb(clean_file, file = "7sk3_clean.pdb")
# rm(prep_file, clean_file, write_file)

# ---
prep_file <- read.pdb("6meo_prep2.pdb")
clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = FALSE, fix.chain = FALSE, rm.h = TRUE)
write_file <- write.pdb(clean_file, file = "6meo_clean.pdb")
rm(prep_file, clean_file, write_file)



