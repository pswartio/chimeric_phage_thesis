library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)

readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

mumgp = readDelta("/Users/pimswart/References/ecolibl21de3versuscp0536021.delta")

mumgp %>% head %>% kable

filterMum <- function(df, minl=1000, flanks=1e1){
  coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
    summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
    ungroup %>% arrange(desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
  merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
    mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

mumgp.filt = filterMum(mumgp, minl=1e2)
mumgp.filt %>% head %>% kable






# Define prophage regions for Genome 1
prophage_genome1 <- data.frame(
  start = c(1608803, 1909486, 2543444,2765522,2959247,3804165),  # Add your prophage start coordinates for Genome 1
  end = c(1626129, 1930783, 2561466,2777007,2998366,3821325)     # Add your prophage end coordinates for Genome 1
)

# Define prophage regions for Genome 2
prophage_genome2 <- data.frame(
  start = c(521230,748334,873642,1395977,1570151),   # Add your prophage start coordinates for Genome 2
  end = c(564358,800959,910549,1411816,1603239)     # Add your prophage end coordinates for Genome 2
)

# Modify the plot code
output_plot <- ggplot(mumgp.filt, aes(x = rs, xend = re, y = qs, yend = qe, colour = strand)) +
  geom_segment() +
  geom_point(alpha = 0.1) +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 180, size = 5),
    legend.position = c(0.99, 0.01),
    legend.justification = c(1, 0),
    strip.background = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.ticks.y = element_blank()
  ) +
  xlab('E_coli_bl21_DE3') +
  ylab('CP053602.1') +
  scale_colour_brewer(palette = 'Set1') +
  # Add a layer for highlighting prophage regions in Genome 1
  annotate(
    "rect",
    xmin = prophage_genome1$start, xmax = prophage_genome1$end,
    ymin = -Inf, ymax = Inf,
    fill = "seagreen", alpha = 0.2
  ) +
  # Add a layer for highlighting prophage regions in Genome 2
  annotate(
    "rect",
    ymin = prophage_genome2$start, ymax = prophage_genome2$end,
    xmin = Inf, xmax = -Inf,
    fill = "purple", alpha = 0.2
  )

ggsave("outputt_pplot.png", plot, width = 18.74, height = 16.01, units = "cm", dpi = 300)
