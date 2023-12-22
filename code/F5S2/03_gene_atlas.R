source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) GENEATLAS INTERFACE ------------------------------------------------------
# CHEMOKINE
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("source_gnccn") %>% unique() %>%
  filter(!is.na(source_gnccn))
interface <- interface$source_gnccn
data <- read_csv("data/variant/geneatlas/processed/CK_GENEATLAS.csv") %>%
  filter(P_value < 1e-8) %>%
  filter(ccn %in% c(interface)) 

# uniqueid
data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(gene_symbol, "_", ccn, "_", AA_Consequences, "_", Trait_Description, "_", uniqueid)) 
data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[rev(order(data$P_value))])

data %>%
  ggplot(aes(-log(P_value), gene_gnccn)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# RECEPTOR
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("target_gnccn") %>% unique() %>%
  filter(!is.na(target_gnccn))
interface <- interface$target_gnccn
data <- read_csv("data/variant/geneatlas/processed/CKR_GENEATLAS.csv") %>%
  filter(P_value < 1e-8) %>%
  filter(gn %in% c(interface)) 

# uniqueid
data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(gene_symbol, "_", gn, "_", AA_Consequences, "_", Trait_Description, "_", uniqueid)) 
data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[rev(order(data$P_value))])

data %>%
  ggplot(aes(-log(P_value), gene_gnccn)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# (2) GENEATLAS GENE -----------------------------------------------------------
# CHEMOKINE
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("source_gnccn") %>% unique() %>%
  filter(!is.na(source_gnccn))
interface <- interface$source_gnccn
data <- read_csv("data/variant/geneatlas/processed/CK_GENEATLAS.csv") %>%
  filter(P_value < 1e-8) %>%
  filter(ccn %in% c(interface)) %>%
  dplyr::count(Gene)

# uniqueid
data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(Gene, "_", uniqueid)) 
data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[(order(data$n))])

data %>%
  ggplot(aes(n, gene_gnccn)) +
  geom_bar(stat = "identity") +
  xlim(0,8) +
  theme_minimal()

# ggsave(filename = "geneatlas_ck_gene.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 3,
#        height = 3)

# RECEPTOR
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("target_gnccn") %>% unique() %>%
  filter(!is.na(target_gnccn))
interface <- interface$target_gnccn
data <- read_csv("data/variant/geneatlas/processed/CKR_GENEATLAS.csv") %>%
  filter(P_value < 1e-8) %>%
  filter(gn %in% c(interface)) %>%
  dplyr::count(Gene)

# uniqueid
data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(Gene, "_", uniqueid)) 
data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[(order(data$n))])

data %>%
  ggplot(aes(n, gene_gnccn)) +
  geom_bar(stat = "identity") +
  xlim(0,8) +
  theme_minimal()

# ggsave(filename = "geneatlas_ckr_gene.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 3,
#        height = 3)

