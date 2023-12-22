source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

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
  dplyr::summarise(allele_cnt_gene = sum(Allele.Count), allele_no_gene = sum(Allele.Number)) %>%
  ungroup()
data <- data %>% mutate(allele_freq = allele_cnt_gene/allele_no_gene)
data$gene_symbol <- factor(data$gene_symbol, levels = data$gene_symbol[order(data$allele_freq)])

data %>%
  ggplot(aes(allele_freq, gene_symbol)) +
  geom_bar(stat = "identity") +
  theme_minimal()


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
  dplyr::summarise(allele_cnt_gene = sum(Allele.Count), allele_no_gene = sum(Allele.Number)) %>%
  ungroup()
data <- data %>% mutate(allele_freq = allele_cnt_gene/allele_no_gene)

data$gene_symbol <- factor(data$gene_symbol, levels = data$gene_symbol[order(data$allele_freq)])

data %>%
  ggplot(aes(allele_freq, gene_symbol)) +
  geom_bar(stat = "identity") +
  theme_minimal()
