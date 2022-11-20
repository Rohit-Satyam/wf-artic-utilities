if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos='http://cran.us.r-project.org',dependencies = TRUE)
packages <- c('dplyr',"BiocManager",'stringr','pacman','tidyr','ggplot2','R3port','ggpubr','plotly','data.table','argparse','ggsci','htmlwidgets','egg')
if (!suppressMessages(require("Biostrings", quietly = TRUE,warn.conflicts = FALSE))) BiocManager::install('Biostrings')

install.packages(setdiff(packages, rownames(installed.packages())), repos='http://cran.us.r-project.org',dependencies = TRUE)
#pacman::p_load(c(packages,'Biostrings'), character.only = T,verbose=F,quietly=TRUE)
t <- c(packages,'S4Vectors','BiocGenerics','Biostrings')
invisible(lapply(t, function(xxx) suppressMessages(library(xxx,character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE))))

parser <- ArgumentParser(description='Produce Coverage Plot',add_help=T)
parser$add_argument('-b', "--bed",help="Provide BED file of Masked regions", default=NULL,nargs=1)
parser$add_argument('-fs', "--fastcat_summary",help="Provide Fastcat Stats file", default=NULL,nargs=1)
parser$add_argument('-gl', "--genome_len",help="Provide genome length", default=29903,nargs=1)
parser$add_argument('-o',"--output_prefix",help="Output File Name Prefix", default="Summary", nargs=1)


args <- parser$parse_args()

## Reading and plotting the BED file
bed <- read.csv(args$bed, sep = "\t", header = F)
bed$V4 <- bed$V3 - bed$V2
tmp <- data.table::setorder(aggregate(bed$V4, by=list(group=bed$V1), FUN=sum), -x)
tmp$coverage <- (args$genome_len - tmp$x ) * 100/args$genome_len 

tmp <- tmp %>% mutate(qc_status = case_when(
  coverage > 90 ~ "Good",
  coverage <90 & coverage > 80 ~ "Mediocre",
  coverage < 80 ~ "GISAID Bad",
  coverage < 50 ~ "Unsubmittable"
)
)
bed$V1 <- factor(bed$V1, levels = unique(tmp$group))
bed$V5 <- tmp$qc_status[match(bed$V1, tmp$group)]

colnames(bed) <- c("barcode","start","end","length","Status")
bed$Status <- factor(bed$Status, levels = c("GISAID Bad","Mediocre","Good"))

## Ploting range plot
g <- ggplot(bed) + 
  geom_linerange(aes(ymin = start, ymax = end, x = barcode, color=Status),
                                size = 2, alpha = 1)+ 
  coord_flip()+
  theme_article()+
  scale_color_uchicago()+
  scale_fill_uchicago()

p1 <- ggplotly(g)

n_cases <- length(unique(bed$Status))
for (i in 1:n_cases) {
  p1$x$data[[i]]$name <- levels(bed$Status)[i]
  p1$x$data[[i]]$legendgroup <- levels(bed$Status)[i]
  p1$x$data[[i]]$showlegend <- FALSE
}

# Reading the fastcat stats file
file <- read.csv(args$fastcat_summary, sep = '\t')
file$filename <- stringr::str_split(file$filename,pattern = "/", simplify = T) %>% .[,ncol(.)-1]
df <- file %>% group_by(.,filename) %>% summarise(count = n())
df$Status <- tmp$qc_status[match(df$filename, tmp$group)]
df2 <- df[order(match(df$filename,tmp$group)),]
df2$Status <- factor(df2$Status, levels = c("GISAID Bad","Mediocre","Good"))


g2 <- ggbarplot(df2, "filename", "count",
                fill = "Status", color = "Status",label = F,
                xlab = FALSE,  ylab = F, ggtheme = theme_void())+
  coord_flip()+
  theme(axis.ticks.y = element_blank(),axis.text.y=element_blank(), 
        axis.line.y = element_blank())+
  scale_color_uchicago()+
  scale_fill_uchicago()
p2 <- ggplotly(g2)

g3 <- plotly::subplot(p1,p2)
saveWidget(g3, paste0(args$output_prefix,".html"), selfcontained = T, libdir = "lib")
write.csv(tmp,"temp.cov.csv", row.names = F, quote = F)
