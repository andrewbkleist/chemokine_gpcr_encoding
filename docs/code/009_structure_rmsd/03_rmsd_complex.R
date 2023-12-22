# Calculates pairwise RMSD of all chemokines-GPCR complexes by first aligning
# receptor residues then calculating chemokine RMSD 

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# define PDBs and lookup column names
pdbfiles <- c("data/structure/raw/pdb_files/5uiw_ck_clean.pdb",
              "data/structure/raw/pdb_files/7o7f_clean.pdb",
              "data/structure/raw/pdb_files/7f1r_clean.pdb",
              "data/structure/raw/pdb_files/zheng_model_clean.pdb",
              "data/structure/raw/pdb_files/7vl9_clean.pdb",
              "data/structure/raw/pdb_files/7xa3_clean.pdb",
              "data/structure/raw/pdb_files/7f1t_clean.pdb",
              "data/structure/raw/pdb_files/6wwz_clean.pdb",
              "data/structure/raw/pdb_files/6lfo_clean.pdb",
              "data/structure/raw/pdb_files/ngo_model_clean.pdb",
              "data/structure/raw/pdb_files/7sk3_clean.pdb",
              "data/structure/raw/pdb_files/7xbx_clean.pdb",
              "data/structure/raw/pdb_files/4rws_ck_clean.pdb",
              "data/structure/raw/pdb_files/4xt1_ck_clean.pdb",
              "data/structure/raw/pdb_files/5wb2_clean.pdb",
              "data/structure/raw/pdb_files/8ic0_clean.pdb")
lookup.ck.gnccn <- c("ccn_5uiw_ck",
                     "ccn_7o7f_ck",
                     "ccn_7f1r_ck",
                     "ccn_zheng_ck",
                     "ccn_7vl9_ck",
                     "ccn_7xa3_ck",
                     "ccn_7f1t_ck",
                     "ccn_6wwz_ck",
                     "ccn_6lfo_ck",
                     "ccn_ngo_ck",
                     "ccn_7sk3_ck",
                     "ccn_7xbx_ck",
                     "ccn_4rws_ck",
                     "ccn_4xt1_ck",
                     "ccn_5wb2_ck",
                     "ccn_8ic0_ck")
lookup.ck.resno <- c("clean_5uiw_ck",
                     "clean_7o7f_ck",
                     "clean_7f1r_ck",
                     "clean_zheng_ck",
                     "clean_7vl9_ck",
                     "clean_7xa3_ck",
                     "clean_7f1t_ck",
                     "clean_6wwz_ck",
                     "clean_6lfo_ck",
                     "clean_ngo_ck",
                     "clean_7sk3_ck",
                     "clean_7xbx_ck",
                     "clean_4rws_ck",
                     "clean_4xt1_ck",
                     "clean_5wb2_ck",
                     "clean_8ic0_ck")
lookup.ckr.gnccn <- c("bw_5uiw_ckr",
                      "bw_7o7f_ckr",
                      "bw_7f1r_ckr",
                      "bw_zheng_ckr",
                      "bw_7vl9_ckr",
                      "bw_7xa3_ckr",
                      "bw_7f1t_ckr",
                      "bw_6wwz_ckr",
                      "bw_6lfo_ckr",
                      "bw_ngo_ckr",
                      "bw_7sk3_ckr",
                      "bw_7xbx_ckr",
                      "bw_4rws_ckr",
                      "bw_4xt1_ckr",
                      "bw_5wb2_ckr",
                      "bw_8ic0_ckr")
lookup.ckr.resno <- c("clean_5uiw_ckr",
                      "clean_7o7f_ckr",
                      "clean_7f1r_ckr",
                      "clean_zheng_ckr",
                      "clean_7vl9_ckr",
                      "clean_7xa3_ckr",
                      "clean_7f1t_ckr",
                      "clean_6wwz_ckr",
                      "clean_6lfo_ckr",
                      "clean_ngo_ckr",
                      "clean_7sk3_ckr",
                      "clean_7xbx_ckr",
                      "clean_4rws_ckr",
                      "clean_4xt1_ckr",
                      "clean_5wb2_ckr",
                      "clean_8ic0_ckr")


# loop all-by-all comparisons of CA RMSD
data <- as.data.frame(NULL)
for (i in 1:length(pdbfiles)){
  for(j in 1:length(pdbfiles)){
    temp <- GetPairwiseComplexCalphaRMSD(pdbfiles[i], pdbfiles[j],
                                         lookup.ck.gnccn[i], lookup.ck.resno[i],
                                         lookup.ckr.gnccn[i], lookup.ckr.resno[i],
                                         lookup.ck.gnccn[j], lookup.ck.resno[j],
                                         lookup.ckr.gnccn[j], lookup.ckr.resno[j])
    data <- rbind(data, temp)
  }
}
rm(temp)

# remove self-by-self comparisons
data <- data %>% filter(file1 != file2)

# remove duplicates
# see https://stackoverflow.com/questions/25297812/pair-wise-duplicate-removal-from-dataframe
# you may be tempted to run the below two lines on the "data" df alone...
# DO NOT DO THIS - it only works for 2-column dfs, and "scrambles" columns 
# with > 2 columns, hence the "left_join" solution below
temp <- data %>% dplyr::select(file1, file2)
temp <- data.frame(t(apply(temp,1,sort)))
temp <- temp[!duplicated(temp),]
colnames(temp) <- c("file1","file2")
temp$temp <- c("temp")
data <- left_join(data, temp)
data <- data %>% dplyr::filter(temp == "temp") %>% dplyr::select(-temp)
rm(temp)


# write csv
# write_csv(data, "data/structure/processed/pairwise_complex_rmsd.csv")
# WRITTEN 20231031
