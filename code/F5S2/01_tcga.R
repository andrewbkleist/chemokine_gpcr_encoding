source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) TCGA INTERFACE -----------------------------------------------------------
# CHEMOKINE
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("source_gnccn") %>% unique() %>%
  filter(!is.na(source_gnccn))
interface <- interface$source_gnccn
data <- read_csv("data/variant/tcga/processed/CK_TCGA_COUNTS.csv") %>%
  filter(ccn %in% c(interface))
data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(protein, "_", ccn, "_", uniqueid)) 

data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[order(data$cancer_mut_count)])

data %>%
  top_n(10, cancer_mut_count) %>%
  # arrange(cancer_mut_count) %>%
  # arrange(match(gene_gnccn, unique(gene_gnccn))) %>%
  ggplot(aes(cancer_mut_count, gene_gnccn)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# RECEPTOR
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("target_gnccn") %>% unique() %>%
  filter(!is.na(target_gnccn))
interface <- interface$target_gnccn
data <- read_csv("data/variant/tcga/processed/CKR_TCGA_COUNTS.csv") %>%
  separate(col = gn, into = c("a", "gn"), sep = "gn", remove = FALSE) %>% dplyr::select(-a) %>%
  filter(gn %in% c(interface))
data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(protein, "_", gn, "_", uniqueid)) 

data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[order(data$cancer_mut_count)])

data %>%
  top_n(10, cancer_mut_count) %>%
  # arrange(cancer_mut_count) %>%
  # arrange(match(gene_gnccn, unique(gene_gnccn))) %>%
  ggplot(aes(cancer_mut_count, gene_gnccn)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# (2) TCGA GENE -----------------------------------------------------------
# CHEMOKINE
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("source_gnccn") %>% unique() %>%
  filter(!is.na(source_gnccn))
interface <- interface$source_gnccn
data <- read_csv("data/variant/tcga/processed/CK_TCGA_COUNTS.csv")
data <- data %>%
  mutate(cancer_mut_count = case_when(
    is.na(cancer_mut_count) ~ 0,
    !(is.na(cancer_mut_count)) ~ cancer_mut_count
  ))
data <- data %>%
  filter(ccn %in% c(interface)) %>%
  dplyr::group_by(protein) %>%
  dplyr::summarise(n = sum(cancer_mut_count)) %>%
  ungroup()

data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(protein, "_", uniqueid)) 

data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[order(data$n)])

data %>%
  # top_n(10, n) %>%
  # arrange(cancer_mut_count) %>%
  # arrange(match(gene_gnccn, unique(gene_gnccn))) %>%
  ggplot(aes(n, gene_gnccn)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# ggsave(filename = "tcga_ck_gene.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 3,
#        height = 6)

# RECEPTOR
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("target_gnccn") %>% unique() %>%
  filter(!is.na(target_gnccn))
interface <- interface$target_gnccn
data <- read_csv("data/variant/tcga/processed/CKR_TCGA_COUNTS.csv")
data <- data %>%
  mutate(cancer_mut_count = case_when(
    is.na(cancer_mut_count) ~ 0,
    !(is.na(cancer_mut_count)) ~ cancer_mut_count
  ))
data <- data %>%
  separate(col = gn, into = c("a", "gn"), sep = "gn", remove = FALSE) %>% dplyr::select(-a) %>%
  filter(gn %in% c(interface)) %>%
  dplyr::group_by(protein) %>%
  dplyr::summarise(n = sum(cancer_mut_count)) %>%
  ungroup()

data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(protein, "_", uniqueid)) 

data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[order(data$n)])

data %>%
  # top_n(10, n) %>%
  # arrange(cancer_mut_count) %>%
  # arrange(match(gene_gnccn, unique(gene_gnccn))) %>%
  ggplot(aes(n, gene_gnccn)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# ggsave(filename = "tcga_ckr_gene.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 3,
#        height = 6)
