
# clear workspace
rm(list = ls(all.names = TRUE))

# list packages to load
packages <- c("devtools","conflicted","here","magrittr", "dplyr", "tidyr", "ggplot2","ggpubr","RColorBrewer", "ciftiTools","tableone", "data.table", "reshape2","neuroCombat")

# install packages if not yet installed
# note: ciftiTools install fails if R is started without enough memory on cluster (try 16G)
all_packages <- rownames(installed.packages())
installed_packages <- packages %in% all_packages
if (any(installed_packages == FALSE)){install.packages(packages[!installed_packages])}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# use the filter function from dplyr, not stats
conflict_prefer("filter", "dplyr")

# get path to project repo directory
project <- here()
print(paste("Project directory:", project))

# read rsfa
rsfa <- read.csv(file.path(project,"cabnp_key_rsfa.csv"))
rsfa_cortex <- rsfa[1:360,]
rsfa_subcort <- rsfa[361:nrow(rsfa),]

# read netho
netho <- read.csv(file.path(project,"cabnp_key_netho.csv"))
netho_cortex <- netho[1:360,]
netho_subcort <- netho[361:nrow(netho),]

# correlations
cor_cortex <- cor.test(x=rsfa_cortex$rsfa_diff_beta, y=netho_cortex$netho_diff_beta, na.action="na.omit")
cor_cortex
cor_subcort <- cor.test(x=rsfa_subcort$rsfa_diff_beta, y=netho_subcort$netho_diff_beta, na.action="na.omit")
cor_subcort


# merge rsfa and netho
merged <- merge(x=rsfa, y=netho, by="INDEX")
merged$both_sig <- !is.na(merged$netho_diff_beta_parasig) & !is.na(merged$rsfa_diff_beta_parasig)

for (f in 1:nrow(merged)){
  if(sum(!is.na(merged[f,c("netho_diff_beta_parasig","rsfa_diff_beta_parasig")]))==2){
    merged[f,"Significance"] <- "Both"
  }else if(!is.na(merged[f,"netho_diff_beta_parasig"])){
    merged[f,"Significance"] <- "LC"
  }else if(!is.na(merged[f,"rsfa_diff_beta_parasig"])){
    merged[f,"Significance"] <- "RSFA"
  }else{
    merged[f,"Significance"] <- NA
  }
} 
        

# plot lm
# cortex
text <- paste0("r = ",signif(cor_cortex$estimate, digits=3), ", p = ",signif(cor_cortex$p.value, digits=3))
cort_plot <- ggplot(data=merged[1:360,], aes(x=rsfa_diff_beta, y=netho_diff_beta))+
  geom_point(aes(fill=Significance), alpha=0.7, shape=21, color="gray20")+
  scale_fill_manual(limits = c(NA,"RSFA","LC","Both"),values=c("gray90","blue","yellow","green"))+
  geom_smooth(method = "lm", color="red")+
  annotate("text",y=-0.6, x=0.3, label=text)+
  xlab(expression(paste("RSFA group difference (",beta,")")))+
  ylab(expression(paste("NetHo group difference (",beta,")")))+
  ggtitle("Comparing RSFA vs NetHo group differences: cortex")+
  theme_classic()+
  #theme(axis.text=element_text(color="gray60"), legend.position="none")
  theme(axis.text=element_text(color="gray60"))
cort_plot
ggsave("~/Dropbox/PhD/conferences/2023_SOBP/poster/rsfa_netho_cortex.pdf", plot=cort_plot, device="pdf", width=5.9, height=4)

# sc
text <- paste0("r = ",signif(cor_subcort$estimate, digits=3), ", p = ",signif(cor_subcort$p.value, digits=3))
subcort_plot <- ggplot(data=merged[361:nrow(merged),], aes(x=rsfa_diff_beta, y=netho_diff_beta))+
  geom_point(aes(fill=Significance), alpha=0.7, shape=21, color="gray20")+
  scale_fill_manual(limits = c(NA,"RSFA","LC","Both"),values=c("gray90","blue","yellow","green"))+
  geom_smooth(method = "lm", color="red")+
  annotate("text",y=-0.47, x=0.64, label=text)+
  xlab(expression(paste("RSFA group difference (",beta,")")))+
  ylab(expression(paste("NetHo group difference (",beta,")")))+
  ggtitle("Comparing RSFA vs NetHo group differences: subcortex")+
  theme_classic()+
  theme(axis.text=element_text(color="gray60"))
subcort_plot
ggsave("~/Dropbox/PhD/conferences/2023_SOBP/poster/rsfa_netho_subcortex.pdf", plot=subcort_plot, device="pdf", width=5.9, height=4)



