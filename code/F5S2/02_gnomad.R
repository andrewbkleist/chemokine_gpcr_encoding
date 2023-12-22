source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) GNOMAD INTERFACE ---------------------------------------------------------
# CHEMOKINE
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("source_gnccn") %>% unique() %>%
  filter(!is.na(source_gnccn))
interface <- interface$source_gnccn
data <- read_csv("data/variant/gnomad/processed/CK_GNOMAD_TABLE.csv") %>%
  dplyr::select(gene_symbol, ccn, aln_pos, Consequence, Allele.Number, Allele.Count, Allele.Frequency) %>%
  filter(ccn %in% c(interface))
data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(gene_symbol, "_", ccn, "_", Consequence, "_", uniqueid)) 

data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[order(data$Allele.Frequency)])

data %>%
  top_n(10, Allele.Frequency) %>%
  arrange(Allele.Frequency) %>%
  arrange(match(gene_gnccn, unique(gene_gnccn))) %>%
  ggplot(aes(Allele.Frequency, gene_gnccn)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# ggsave(filename = "gnomad_ck_position.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 4.5,
#        height = 5)

# RECEPTOR
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("target_gnccn") %>% unique() %>%
  filter(!is.na(target_gnccn))
interface <- interface$target_gnccn
data <- read_csv("data/variant/gnomad/processed/CKR_GNOMAD_TABLE.csv") %>%
  dplyr::select(gene_symbol, gn, aln_pos, Consequence, Allele.Number, Allele.Count, Allele.Frequency) %>%
  separate(col = gn, into = c("a", "gn"), sep = "gn", remove = FALSE) %>%
  dplyr::select(-a) %>%
  filter(gn %in% c(interface))
data$uniqueid <- 1:nrow(data)
data <- data %>% 
  mutate(gene_gnccn = paste0(gene_symbol, "_", gn, "_", Consequence, "_", uniqueid)) 

data$gene_gnccn <- factor(data$gene_gnccn, levels = data$gene_gnccn[order(data$Allele.Frequency)])

data %>%
  top_n(10, Allele.Frequency) %>%
  arrange(Allele.Frequency) %>%
  arrange(match(gene_gnccn, unique(gene_gnccn))) %>%
  ggplot(aes(Allele.Frequency, gene_gnccn)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# ggsave(filename = "gnomad_ckr_position.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 4.5,
#        height = 5)

# (2) GNOMAD GENE -----------------------------------------------------------
# CHEMOKINE
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("source_gnccn") %>% unique() %>%
  filter(!is.na(source_gnccn))
interface <- interface$source_gnccn
data <- read_csv("data/variant/gnomad/processed/CK_GNOMAD_TABLE.csv") %>%
  filter(ccn %in% c(interface)) %>%
  dplyr::select(gene_symbol, ccn, aln_pos, Consequence, Allele.Number, Allele.Count, Allele.Frequency) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(n = sum(Allele.Count)) %>%
  ungroup()

data$gene_symbol <- factor(data$gene_symbol, levels = data$gene_symbol[order(data$n)])

data %>%
  ggplot(aes(n, gene_symbol)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# ggsave(filename = "gnomad_ck_gene.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 3,
#        height = 6)

# RECEPTOR
interface <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("target_gnccn") %>% unique() %>%
  filter(!is.na(target_gnccn))
interface <- interface$target_gnccn
data <- read_csv("data/variant/gnomad/processed/CKR_GNOMAD_TABLE.csv") %>%
  separate(col = gn, into = c("a", "gn"), sep = "gn", remove = FALSE) %>% dplyr::select(-a) %>%
  filter(gn %in% c(interface)) %>%
  dplyr::select(gene_symbol, gn, aln_pos, Consequence, Allele.Number, Allele.Count, Allele.Frequency) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(n = sum(Allele.Count)) %>%
  ungroup()

data$gene_symbol <- factor(data$gene_symbol, levels = data$gene_symbol[order(data$n)])

data %>%
  ggplot(aes(n, gene_symbol)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# ggsave(filename = "gnomad_ckr_gene.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 3,
#        height = 6)
