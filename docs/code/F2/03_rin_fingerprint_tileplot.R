source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

  # import contacts (remove 6meo)
  rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% filter(file != "6meo")
  rin <- rin %>% dplyr::select(file, source_gnccn, target_gnccn) 
  rin <- rin %>% unite(rin, c(source_gnccn, target_gnccn), sep = "_")
  rin$n <- c(1)
  
  # spread
  rin <- rin %>% pivot_wider(names_from = rin, values_from = n)
  
  # gather, then replace NULL with "0"
  rin <- rin %>% pivot_longer(cols = c("NTc.Cm10_2x60":"b2b3.12_NTr.Cm8"), names_to = "rin", values_to =  "n")
  rin$n <- as.character((rin$n))
  rin <- rin %>% mutate(n = case_when(
    n == "NULL" ~ "0",
    n == "1" ~ n
  ))
  rin$n <- as.numeric((rin$n))
  
  # matrix plot
  order.ck  <- c("7vl9",  # CCL15:CCR1
                 "7xa3",  # CCL2:CCR2
                 "7f1t",  # CCL3:CCR5
                 "7f1r",  # CCL5:CCR5
                 "zheng", # CCL5:CCR5 (model)
                 "5uiw",  # 5P7:CCR5
                 "7o7f",  # 6P4:CCR5
                 "6wwz",  # CCL20:CCR6
                 "8ic0",  # CXCL8:CXCR1
                 "6lfo",  # CXCL8:CXCR2
                 "ngo",   # CXCL12:CXCR4
                 "4rws",  # vMIPII:CXCR4
                 "7sk3",  # CXCL12:ACKR3
                 "7xbx",  # CX3CL1:CX3CR1
                 "4xt1",  # CX3CL1:US28
                 "5wb2")  # CX3CL1.35:US28
  rin$file <- factor(rin$file, levels = order.ck)
  
  # order based on appearance
  rin <- rin %>% filter(n > 0) %>% add_count(rin) %>% dplyr::select(-n)
  rin <- rin %>% arrange(nn, file)
  
  rin$rin  <- as.factor(rin$rin)
  order.res_par <- unique(rin$rin)
  levels(rin$rin)
  rin$rin <- factor(rin$rin, levels = rev(order.res_par))
  
  # plot
  rin %>%
    ggplot() + 
    geom_tile(aes(file, rin), fill = "gray25")+
    theme_minimal() +
    theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  
  # ggsave(filename = "rin_tileplot_fingerprint.pdf",
  #        plot = last_plot(), path = "output/F2/",
  #        width = 4,
  #        height = 5)
  