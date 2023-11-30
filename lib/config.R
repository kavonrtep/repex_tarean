#/usr/bin/env Rscript
TANDEM_RANKS = c(
      "Putative satellites (high confidence)" =  1,
      "Putative satellites (low confidence)" = 2,
      "Putative LTR elements" = 3,
      "rDNA" = 4,
      "Other" = 0
)
# inverted - key value
RANKS_TANDEM = names(TANDEM_RANKS)
names(RANKS_TANDEM) = TANDEM_RANKS
