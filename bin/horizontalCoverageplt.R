packages <- c('dplyr','tidyr','ggplot2','scales','R3port','ggpubr',"patchwork",'data.table','argparse','ggsci','R3port')
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
invisible(suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE)))

parser <- ArgumentParser(description='Produce Coverage Plot',add_help=T)
parser$add_argument('-l', "--list",help="Provide coverage bed files", default=NULL,nargs='+')
parser$add_argument('-t', "--temp",help="Provide Barcode vs Coverage file", default=NULL,nargs=1)
parser$add_argument('-o',"--output_prefix",help="Output File Name Prefix", default="Summary", nargs=1)
parser$add_argument('-w',"--width",help="Plot Width", default=2000, nargs=1)
parser$add_argument('-hi',"--height",help="Plot Height", default=1000, nargs=1)
parser$add_argument('-s',"--subset",help="Subset only mediocre and GISAID Bad quality samples", default=TRUE,nargs=1)

args <- parser$parse_args()

## Reading mosdepth files
beds <- args$list
#t<- list.files(".","per-base.bed.gz$", full.names = T)
df <- lapply(beds, function(x) {read.csv(x,header = F, sep = '\t')})
names(df) <- basename(sub('\\.per-base.bed.gz$', '', beds) )

cov <- read.csv(args$temp)
cov<-subset(cov, cov$group %in% names(df))

if(args$subset){
cov <- subset(cov ,cov$qc_status != "Good")
}
l <- lapply(cov$group, function(x){
 ggplot(df[[x]]) +
  geom_rect(aes(xmin = V2, xmax = V3, ymin = 0, ymax = V4), col="blue", fill="blue")+geom_hline(aes(yintercept = 20))+
    theme_minimal()+theme(axis.text=element_text(size=17),
                          axis.title=element_text(size=17,face="bold"))+xlab(paste(x, "Cov:",percent(round(cov[cov$group==x,]$coverage)/100)))
})

names(l) <- cov$group
t <- wrap_plots(l, ncol = 6)
html_plot(t,paste0(args$output_prefix,".html"), pwidth = as.numeric(args$width), pheight = as.numeric(args$height), title = "Coverage plot for Mediocre and Bad Quality Samples",show = F)
