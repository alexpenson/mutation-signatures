library(data.table)
library(ggplot2)
old <- theme_set(theme_minimal(20))

# Annotate maf with Stratton Plot bin
add_mut_tri <- function(maf) {

  if (!"TriNuc" %in% names(maf)) {
    if ("Ref_Tri" %in% names(maf)) {
      maf[, TriNuc := Ref_Tri]
    } else {
      stop("must have either Ref_Tri or TriNuc column")
    }
  }

  if (!'t_var_freq' %in% names(maf)) maf[!t_ref_count %in% c(NA,'.') & !t_alt_count %in% c(NA,'.'),
                                         t_var_freq := as.numeric(t_alt_count)/(as.numeric(t_alt_count)+as.numeric(t_ref_count))]

  maf[, c('TriNuc_CT', 'Tumor_Seq_Allele2_CT', 'Mut_Tri') := 'X']
  maf[Variant_Type == "SNP" & !is.na(Reference_Allele) & !is.na(t_var_freq) & !is.na(TriNuc),
      TriNuc_CT := ifelse(Reference_Allele %in% c("G", "A"),
                          as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(TriNuc))),
                          TriNuc)]

  maf[Variant_Type == "SNP" & !is.na(Reference_Allele) & !is.na(t_var_freq) & !is.na(Tumor_Seq_Allele2),
      Tumor_Seq_Allele2_CT := ifelse(Reference_Allele %in% c("G", "A"),
                                     as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(Tumor_Seq_Allele2))),
                                     Tumor_Seq_Allele2)]

  maf[Variant_Type == "SNP" & !is.na(TriNuc_CT) & !is.na(Tumor_Seq_Allele2_CT),
      Mut_Tri := paste0(substr(TriNuc_CT, 1, 2),
                        Tumor_Seq_Allele2_CT,
                        substr(TriNuc_CT, 3, 3))]

  maf
}

stratton_plot <- function(maf){
  if (!"Mut_Tri" %in% names(maf)) maf <- add_mut_tri(maf)

  adam_trinuc_colnames <-
    c("ACAA", "ACAC", "ACAG", "ACAT", "CCAA", "CCAC", "CCAG", "CCAT",
      "GCAA", "GCAC", "GCAG", "GCAT", "TCAA", "TCAC", "TCAG", "TCAT",
      "ACGA", "ACGC", "ACGG", "ACGT", "CCGA", "CCGC", "CCGG", "CCGT",
      "GCGA", "GCGC", "GCGG", "GCGT", "TCGA", "TCGC", "TCGG", "TCGT",
      "ACTA", "ACTC", "ACTG", "ACTT", "CCTA", "CCTC", "CCTG", "CCTT",
      "GCTA", "GCTC", "GCTG", "GCTT", "TCTA", "TCTC", "TCTG", "TCTT",
      "ATAA", "ATAC", "ATAG", "ATAT", "CTAA", "CTAC", "CTAG", "CTAT",
      "GTAA", "GTAC", "GTAG", "GTAT", "TTAA", "TTAC", "TTAG", "TTAT",
      "ATCA", "ATCC", "ATCG", "ATCT", "CTCA", "CTCC", "CTCG", "CTCT",
      "GTCA", "GTCC", "GTCG", "GTCT", "TTCA", "TTCC", "TTCG", "TTCT",
      "ATGA", "ATGC", "ATGG", "ATGT", "CTGA", "CTGC", "CTGG", "CTGT",
      "GTGA", "GTGC", "GTGG", "GTGT", "TTGA", "TTGC", "TTGG", "TTGT"
    )
  maf <- maf[Mut_Tri %in% adam_trinuc_colnames]
  maf[, Mut_Tri := factor(Mut_Tri, levels = adam_trinuc_colnames)]
  ggplot(maf[Mut_Tri %in% adam_trinuc_colnames],
         aes(Mut_Tri
             ,
             fill = factor(
               paste0(
                 substr(Mut_Tri, 2, 2),
                 " > ",
                 substr(Mut_Tri, 3, 3)
               ))
             )
         ) +
    geom_bar() +
    facet_grid(Tumor_Sample_Barcode~.) +
    scale_fill_manual("", values = c("#1EBFF0", "#050708", "#E62725", "#CBCACB", "#A1CF64", "#EDC8C5")) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          strip.text.y = element_text(angle=0),
          axis.text.x = element_text(size = rel(0.2), angle = 90),
#          axis.title.y=element_text(vjust=1),
          panel.margin = unit(0, "cm")) +
    xlab("") +
    ylab("Fraction of Mutations") +
      theme(legend.text=element_text(size=16, family="Courier"))
  # ggsave(file=paste0(gsub(".txt$", "_stratton_plot.pdf", summ_file)),
  #                    width = 15, height = 3)
}
