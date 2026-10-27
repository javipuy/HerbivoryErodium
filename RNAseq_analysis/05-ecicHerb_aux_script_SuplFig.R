##############################################
### E. CICUTARIUM HERBIVORY: FIGURE SUPPL. ###
##############################################

#### Import topGO results from azacitidine comparisons ####

library(ggplot2)
library(scales)
library(forcats)
library(gridExtra)

# load enriched GO terms per comparison
setwd("/path/to/project/directory/topGO/") ## SET CUSTOM PATH TO WORK DIRECTORY

goterms <- read.csv("GOtermPlotSummary.csv", sep = ';', header = T)

# generate subsets
goterms.UNvsDEdown <- goterms[which(goterms$Comparison == "herbAA_down"), ]
goterms.deUNvsdeMHup <- goterms[which(goterms$Comparison == "deUNdeMH_up"), ]
goterms.deMHvsdeHSMH <- goterms[which(goterms$Comparison == "deMHdeHS_MH"), ]

# generate a funtion to print GO term plots for specific GO term subsets from modified topGO output
goTermPlot <- function(goTermSubset) {
  mainPlot <- 
    ggplot(goTermSubset, 
           aes(x = fct_reorder(Term, weightFDR, .desc = TRUE), 
               y = -log10(weightFDR), 
               size = log2(Significant/Expected), 
               fill = -log10(weightFDR))) +
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5)) +
    scale_fill_continuous(low = 'steelblue', high = 'tomato') +
    xlab('') + ylab('Enrichment score') +
    labs(title = paste('Enriched GO terms')) +
    # Theme and general format
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      # Axis format
      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
      axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
      axis.title = element_text(size = 12, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.line = element_line(colour = 'black'),
      #Legend format
      legend.key = element_blank(), # removes the border
      legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
      legend.text = element_text(size = 14, face = "bold"), # Text size
      title = element_text(size = 14, face = "bold")) +
    coord_flip()
  return(mainPlot)
}

# Figure S4
goplot.UNvsDEdown <- goTermPlot(goterms.UNvsDEdown)

goplot.UNvsDEdown2 <- ggplotGrob(goplot.UNvsDEdown)
goplot.UNvsDEdown2$widths <- goplot.UNvsDEdown2$widths
library(gridExtra)
grid.arrange(goplot.UNvsDEdown2)

# Figure demethylated
goplot.deUNvsdeMHup <- goTermPlot(goterms.deUNvsdeMHup)
goplot.deMHvsdeHSMH <- goTermPlot(goterms.deMHvsdeHSMH)

goplot.deUNvsdeMHup2 <- ggplotGrob(goplot.deUNvsdeMHup)
goplot.deUNvsdeMHup2$widths <- goplot.deUNvsdeMHup2$widths
goplot.deMHvsdeHSMH2 <- ggplotGrob(goplot.deMHvsdeHSMH)
goplot.deMHvsdeHSMH2$widths <- goplot.deMHvsdeHSMH2$widths

library(gridExtra)
grid.arrange(goplot.deUNvsdeMHup2)
grid.arrange(goplot.deMHvsdeHSMH2)

