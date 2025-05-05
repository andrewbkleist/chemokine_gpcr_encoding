library(bio3d)
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/01_structure_contacts/data/pdbs/raw/20211002/")

clean.pdb.custom <- function(PREPPDB, OUTNAME){
  prep_file <- read.pdb(PREPPDB)
  clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = FALSE, fix.chain = FALSE, rm.h = TRUE)
  write_file <- write.pdb(clean_file, file = OUTNAME)
  rm(prep_file, clean_file, write_file)
}


# ---
clean.pdb.custom("7o7f_prep2.pdb", "7o7f_prep.pdb")

# ---
prep_file <- read.pdb("7f1t_prep1.pdb")
clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = TRUE, fix.chain = TRUE, rm.h = TRUE)
write_file <- write.pdb(clean_file, file = "7f1t_prep.pdb")
rm(prep_file, clean_file, write_file)

# ---
prep_file <- read.pdb("af_ccl20_ccr6_prep1.pdb")
clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = TRUE, fix.chain = TRUE, rm.h = TRUE)
write_file <- write.pdb(clean_file, file = "af_ccl20_ccr6_prep.pdb")
rm(prep_file, clean_file, write_file)

# ---
prep_file <- read.pdb("af_cxcl8_cxcr2_prep1.pdb")
clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = TRUE, fix.chain = TRUE, rm.h = TRUE)
write_file <- write.pdb(clean_file, file = "af_cxcl8_cxcr2_prep.pdb")
rm(prep_file, clean_file, write_file)

# ---
prep_file <- read.pdb("af_ccl3_ccr5_prep1.pdb")
clean_file <- clean.pdb(prep_file, consecutive = FALSE, force.renumber  = TRUE, fix.chain = TRUE, rm.h = TRUE)
write_file <- write.pdb(clean_file, file = "af_ccl3_ccr5_prep.pdb")
rm(prep_file, clean_file, write_file)





