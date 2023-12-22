source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import data
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(file != "6meo")
  
# plot
rin %>%
  ggplot(aes(all_non_ackr_para_ckr, all_para_ck))  +
  geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
  xlim(0,1) +
  ylim(0,1) +
  theme_minimal()

# note that "missing values" correspond to contacts made with GPCR positions
# in regions after 4x65 and prior to 45x50 which are not structurally 
# equivalent and thus do not have associated conservation scores; 
# this is because GPCRdb alignment "left adjusts" these residues 
# without gaps; see Methods for details; there is one additional GPCR position 
# (2x67 in CCR1 that makes a contact in the CCR1 complex) but is not defined in 
# our alignment; this is due to D97 being relabeled as 2x67 in updated GPCRdb 
# alignment but unlabeled in the original GPCRdb alignment. Note that this 
# corresponds to ECL1.1 in CGN scheme here.
# 2x67 and not ECL1.1 was assigned since CCL15-CCR1 complex was added recently 
# and after initial CGN assignments were made using the "old" GPCRdb numbering 
# scheme; regardless of name this is a structurally-non-equivalent position
# like the 4x65-45x50 (not inclusive) positions

# ggsave(filename = "paracons_by_paracons.pdf",
#        plot = last_plot(), path = "output/F2/",
#        width = 4.5,
#        height = 4)

  

  
  
  