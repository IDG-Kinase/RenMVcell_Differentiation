#source('http://depot.sagebase.org/CRAN.R')
#pkgInstall("synapseClient")
#install.packages("rstudioapi")

library(synapseClient)
library(plyr);library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
`%nin%`<-Negate(`%in%`)

##################################################################################################################T
# Get user input ------------
##################################################################################################################T
#User input
Synapse_UID<-""
Synapse_password<-""
dir_project<-""
synapse_ids<-c("syn14619049")#c("synID_1","synID_2",...,"synID_n") ## requires user input
variable_names<-c("Phospho_baseline_g_day") ## optional user input

##################################################################################################################T
# make connection to synapse, set directories ------------
##################################################################################################################T
synapseLogin(Synapse_UID, Synapse_password) ## requires user input
c.time<-Sys.time()%>%as.character()%>%gsub("-","",.)%>%gsub(" ","_",.)%>%gsub(":","",.)
#dir_project<-"/Users/nienke/Dropbox (HMS)/CBDM-SORGER/Collaborations/Dark_kinome/" ##requires user input
paste0(dir_project,paste0(rstudioapi::getActiveDocumentContext()$path%>%dirname()%>%basename()))
dir_c_script<-paste0(dir_project,paste0(rstudioapi::getActiveDocumentContext()$path%>%dirname()%>%basename()))
dir_log<-paste0(dir_project,paste0(rstudioapi::getActiveDocumentContext()$path%>%dirname()%>%basename()),"_SynapseLog")
tempdir<-paste0(dir_project,"temp_local_synapse_files",Sys.time()%>%as.character()%>%gsub("-","",.)%>%gsub(" ","_",.)%>%gsub(":","",.))
dir.create(tempdir, showWarnings = TRUE, recursive = TRUE)
#synapse_ids<-c("syn12074621","syn12074611")#c("synID_1","synID_2",...,"synID_n") ## requires user input
#variable_names<-NA#optional: change to c("myname_1","myname_2",...,"myname_n") ## optional user input
##################################################################################################################T
# get info of files copied in for user to check ------------
##################################################################################################################T
file_info.l<-list()
i=0
for(syn_id in synapse_ids){
  i=i+1
  syn_file<-synGet(syn_id, downloadFile =F, downloadLocation = tempdir, ifcollision = "overwrite.local")
  c.info<-list()
  c.info$file_name<-syn_file[[1]]$name
  c.info$R_variable_name<-ifelse(is.na(variable_names)==T,tools::file_path_sans_ext(c.info$file_name),variable_names[i])
  c.info$id<-syn_file[[1]]$id%>%as.character(.)
  c.info$descreption<-syn_file[[1]]$description
  c.info$url<-syn_file@synapseWebUrl
  c.info$version_number<-syn_file[[1]]$versionNumber
  c.info$versions_label<-syn_file[[1]]$versionLabel
  c.info$DateTime_modified<-syn_file[[1]]$modifiedOn
  file_info.l[[i]]<-c.info%>%as.data.frame(.,stringsAsFactors=F)
}
file_info<-file_info.l%>%bind_rows(.)

writeLines(
  paste0("to be loaded:\n ",
         file_info$file_name%>%toString()%>%gsub(",","\n",.),
         "\n===========================================================
         \nplease check 'file_info' before continuing loading process"))

##################################################################################################################T
# download files from Synapse ------------
##################################################################################################################T

for(syn_id in synapse_ids){
  syn_file<-synGet(syn_id, downloadFile =T, downloadLocation = tempdir, ifcollision = "overwrite.local")
  c.info<-list()
  c.info$file_name<-syn_file[[1]]$name
  c.info$R_variable_name<-ifelse(is.na(variable_names)==T,tools::file_path_sans_ext(c.info$file_name),variable_names[i])
  c.info$id<-syn_file[[1]]$id
  c.info$descreption<-syn_file[[1]]$description
  c.info$url<-syn_file@synapseWebUrl
  c.info$version_number<-syn_file[[1]]$versionNumber
  c.info$versions_label<-syn_file[[1]]$versionLabel
  c.info$DateTime_modified<-syn_file[[1]]$modifiedOn
  file_info.l[[i]]<-c.info%>%as.data.frame(.,stringsAsFactors=F)
}
file_info<-file_info.l%>%bind_rows(.)
synapseLogout()

##################################################################################################################T
# write logfile ------------
##################################################################################################################T
setwd(dir_c_script)
paste0(dir_log)%nin%list.files()
if(dir_log%nin%list.files()){
  dir.create(dir_log, showWarnings = TRUE, recursive = TRUE)}

currentscript_insert<-paste0("_",
                             rstudioapi::getActiveDocumentContext()$path%>%basename()%>%gsub(".R","",.),
                             "_")
logfile_name<-Sys.time()%>%as.character()%>%gsub("-","",.)%>%
  gsub(" ",currentscript_insert,.)%>%gsub(":","",.)%>%paste0(.,"GMT")%>%
  paste0("FilesUsed_",.,".csv")

#Sys.time()%>%as.character()%>%gsub("-","",.)%>%
#  gsub(" ",currentscript_insert,.)%>%gsub(":","",.)%>%paste0(.,"GMT")%>%
#  paste0("FilesUsed_",.,".csv")

setwd(dir_log)
write.csv(file_info,file = logfile_name,row.names = F)

##################################################################################################################T
# load data into R ------------
##################################################################################################################T
##!!modify to allow more datatypes
setwd(tempdir)
file_info<-file_info%>%arrange(file_name)
files<-ifelse(list.files()%>%sort()==file_info$file_name%>%as.character(),
              list.files(),"something went wrong")
i=0
for(c.file in files){
  i=i+1
  c.varname<-file_info$R_variable_name[i]%>%as.character()
  c.file<-read.csv(c.file, stringsAsFactors = F) #!!modify to allow more datatypes
  assign(c.varname,c.file)
}
##################################################################################################################T
# remove excess variables ------------
##################################################################################################################T
vars_to_keep<-c(c("%nin%","file_info",file_info$R_variable_name),
                "Synapse_UID","Synapse_password", "dir_project","tempdir",
                "dir_log","currentscript_insert","dir_c_script")
rm(list = ls()[ls()%nin%vars_to_keep])

##################################################################################################################T

##################################################################################################################T
# determine AUC phosphorylation over time ------------
##################################################################################################################T
Phospho_baseline_g_gene<-Phospho_baseline_g_day%>%
  group_by(Gene_Symbol, gene_id, uniprot_id)%>%
  summarize(variance=var(score_mean,na.rm=T),
            mean=mean(score_mean,na.rm=T),
            max=max(score_mean,na.rm=T),
            min=min(score_mean,na.rm=T),
            is_kinase=unique(is_kinase),
            is_dark_kinase=unique(is_dark_kinase))%>%
  ungroup()%>%as.data.frame()

# calculate AUC
AUC_split<-dlply(Phospho_baseline_g_day,.(Gene_Symbol,gene_id,uniprot_id),c)

#genes<-Phospho_baseline%>%select(Gene_Symbol,gene_id,uniprot_id)%>%unique()
AUC_table<-list()
for(gene_index in 1:length(AUC_split)){
  c.gene<-AUC_split[[gene_index]]%>%as.data.frame()
  n_timepoints<-dim(c.gene)[1]
  AUC_accum=0
  for(time_index in 2:n_timepoints){
    A<-c.gene$day[time_index]-c.gene$day[time_index-1]
    B<-min(c.gene$score_mean[time_index],c.gene$score_mean[time_index-1])
    C<-0.5*abs(c.gene$score_mean[time_index]-c.gene$score_mean[time_index-1])
    AUC<-A*(B+C)
    AUC_accum<-AUC+AUC_accum
  }
  c.result<-list()
  c.result$Gene_Symbol<-unique(c.gene$Gene_Symbol)
  c.result$gene_id<-unique(c.gene$gene_id)
  c.result$uniprot_id<-unique(c.gene$uniprot_id)
  c.result$is_kinase<-unique(c.gene$is_kinase)
  c.result$is_dark_kinase<-unique(c.gene$is_dark_kinase)
  c.result$AUC<-AUC_accum
  c.result$score_Day_0<-c.gene[c.gene$day==0,]$score_mean
  c.result$score_Day_15<-c.gene[c.gene$day==15,]$score_mean
  c.result$n_timepoints<-n_timepoints
  c.result<-as.data.frame(c.result)
  AUC_table[[gene_index]]<-c.result
  print(paste0(gene_index,"-",length(AUC_split)))
}

AUCtable_bound<-AUC_table%>%bind_rows()
AUCtable_bound<-AUCtable_bound%>%
  mutate(target_type=
           ifelse(is.na(is_dark_kinase)==F,"dark_kinase",
                  ifelse(is.na(is_kinase)==F,"kinase",
                         ifelse(is.na(is_kinase)==T & is.na(is_dark_kinase)==T,"other","error"))))
head(AUCtable_bound)

#setwd(dir_c_script)
#write.csv(AUCtable_bound,file="RenMV_baseline_proteomics.csv",row.names = F)
#AUCtable_bound<-read.csv("RenMV_baseline_proteomics.csv")

AUC_plot<-
  ggplot(AUCtable_bound, aes(x=AUC))+
  geom_histogram()+
  geom_rug(data=AUCtable_bound%>%filter(is_kinase==TRUE),aes(colour=target_type))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 ),
        axis.text.y = element_text(angle = 0, hjust = 0.5 ),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black") )+
  labs(title="distribution AUC protein abundance", x="AUC",y="Count")
setwd(dir_c_script)
ggsave(AUC_plot,filename = "ReNMV_AUC_protein_abundance.pdf",
       device = "pdf",height = 6, width=8)

AUC_plot_2D_FC<-
  ggplot(AUCtable_bound, aes(x=AUC))+
  geom_point(aes(y=log2(score_Day_0/10)), colour="red", alpha=0.3)+
  geom_point(aes(y=log2(score_Day_15/10)),colour="blue",alpha=0.3)+
  geom_point(data=AUCtable_bound%>%filter(target_type=="dark_kinase"),
             aes(y=log2(score_Day_0/10)),colour="black",size=1.5, alpha=1)+
  geom_point(data=AUCtable_bound%>%filter(target_type=="dark_kinase"),
             aes(y=log2(score_Day_0/10)),colour="orange",size=1, alpha=1)+
  geom_point(data=AUCtable_bound%>%filter(target_type=="dark_kinase"),
             aes(y=log2(score_Day_15/10)),colour="black",size=1.5,alpha=1)+
  geom_point(data=AUCtable_bound%>%filter(target_type=="dark_kinase"),
             aes(y=log2(score_Day_15/10)),colour="lightblue",size=1,alpha=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 ),
        axis.text.y = element_text(angle = 0, hjust = 0.5 ),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black") )+
  labs(title="AUC versus Relative Expression", 
       x="AUC",
       y="FoldChange\nover mean expression")

AUC_plot_2D_FC
setwd(dir_c_script)
ggsave(AUC_plot_2D_FC,filename = "ReNMV_AUC_vs_RelativeExpression.pdf",
       device = "pdf",height = 6, width=8)

# define max upreg and downreg gene for three types --
# label those w/ different background colors

##################################################################################################################T
# define reference kinase & reference gene ------------
##################################################################################################################T
#reference kinase: the most extreme kinase example
#reference gene: the most extreme gene example
#references: variance_min, variance_max, max_mean (over time), min_mean(over time) --> min_min, min_max, max_min, max_max
#data seems normalized on mean --> use genes with max and min variance
## shouldn't take min and max --> gives outliers mostly --> take area under curve
## max AUC is only genes going up over time --> take max AUC for Day0>Day15 and Day15<Day0 respectively --> most Day0>Day15 are noise
## reference light kinases gives best reference overview

AUCtable_bound%>%arrange(desc(AUC))%>%filter(target_type=="other" & score_Day_0>score_Day_15)%>%head
AUCtable_bound%>%arrange(desc(AUC))%>%filter(target_type=="other" & score_Day_0<score_Day_15)%>%head
AUCtable_bound%>%arrange(desc(AUC))%>%filter(target_type=="kinase" & score_Day_0>score_Day_15)%>%head
AUCtable_bound%>%arrange(desc(AUC))%>%filter(target_type=="kinase" & score_Day_0<score_Day_15)%>%head
AUCtable_bound%>%arrange(desc(AUC))%>%filter(target_type=="dark_kinase" & score_Day_0>score_Day_15)%>%head
AUCtable_bound%>%arrange(desc(AUC))%>%filter(target_type=="dark_kinase" & score_Day_0<score_Day_15)%>%head

ref_gene<-c("VEPH1","MTSS1L")
ref_kinase<-c("EGFR","HIPK2")
ref_dark_kinase<-c("MAPK4","NEK1")

plot_ref_other<-
  ggplot(Phospho_baseline_g_day%>%filter(Gene_Symbol%in% ref_gene),
         aes(x=day, y=score_mean, colour=interaction(Gene_Symbol,uniprot_id)))+
  geom_errorbar(aes(
    ymin=score_mean-score_sd, 
    ymax=score_mean+score_sd))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits=c(0,30),expand=c(0,0))+
  scale_color_manual("Symbol.UniprotID",values = c("#a6cee3","#b2df8a","#33a02c"))+ # use "#1f78b4" for dark kinase of interest
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 ),
        axis.text.y = element_text(angle = 0, hjust = 0.5 ),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black") )+
  labs(title="Reference Genes",x="Day",y="Relative phosphorylation protein")
plot_ref_other

plot_ref_kinase<-
  ggplot(Phospho_baseline_g_day%>%filter(Gene_Symbol%in% ref_kinase),
         aes(x=day, y=score_mean, colour=interaction(Gene_Symbol,uniprot_id)))+
  geom_errorbar(aes(
    ymin=score_mean-score_sd, 
    ymax=score_mean+score_sd))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits=c(0,30),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0.15))+
  scale_color_manual("Symbol.UniprotID",values = c("#a6cee3","#b2df8a"))+ # use "#1f78b4" for dark kinase of interest
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 ),
        axis.text.y = element_text(angle = 0, hjust = 0.5 ),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black") )+
  labs(title="Reference Kinases",x="Day",y="Relative phosphorylation protein")
plot_ref_kinase

plot_ref_darkkinase<-
  ggplot(Phospho_baseline_g_day%>%filter(Gene_Symbol%in% ref_dark_kinase),
         aes(x=day, y=score_mean, colour=interaction(Gene_Symbol,uniprot_id)))+
  geom_errorbar(aes(
    ymin=score_mean-score_sd, 
    ymax=score_mean+score_sd))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits=c(0,30),expand=c(0,0))+
  scale_color_manual("Symbol.UniprotID",values = c("#a6cee3","#b2df8a"))+ # use "#1f78b4" for dark kinase of interest
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 ),
        axis.text.y = element_text(angle = 0, hjust = 0.5 ),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black") )+
  labs(title="Reference Dark Kinases",x="Day",y="Relative phosphorylation protein")
plot_ref_darkkinase

#setwd(dir_c_script)
#ggsave(plot_ref_kinase,filename = "ReNMV_refplot_kinase.pdf",
#       device = "pdf",height = 6, width=8)
#ggsave(plot_ref_darkkinase,filename = "ReNMV_refplot_darkkinase.pdf",
#       device = "pdf",height = 6, width=8)
#ggsave(plot_ref_other,filename = "ReNMV_refplot_other.pdf",
#       device = "pdf",height = 6, width=8)

##################################################################################################################T
# plot timecourses ------------
##################################################################################################################T
DarkKinase_phospho<-Phospho_baseline_g_day%>%
  filter(is_dark_kinase==TRUE)%>%
  dlply(.,.(Gene_Symbol,gene_id,uniprot_id),c)

for(n in 1){#:length(DarkKinase_phospho)){
  c.df<-DarkKinase_phospho[[n]]%>%as.data.frame()
  c.symbol<-c.df$Gene_Symbol%>%unique()
  c.uniprotID<-c.df$uniprot_id%>%unique() 
  legend_df<-data_frame(symbol_uniprot=c("EGFR.P00533","HIPK2.Q9H2X6",paste0(c.symbol,".",c.uniprotID)),colors=c("#b2df8a","#a6cee3","black"))%>%
    arrange(symbol_uniprot)
  c.values<-legend_df$colors
  c.filename<-paste0("Example_fig_ReNcell_",c.symbol,"_",c.uniprotID,".pdf")
  
  plot_DK_phospho<-
    ggplot(Phospho_baseline_g_day%>%filter(Gene_Symbol%in% ref_kinase),
           aes(x=day, y=score_mean, colour=interaction(Gene_Symbol,uniprot_id)))+
    geom_errorbar(aes(
      ymin=score_mean-score_sd, 
      ymax=score_mean+score_sd))+
    geom_point()+
    geom_line()+
    geom_errorbar(data=c.df,
                  aes(ymin=score_mean-score_sd, ymax=score_mean+score_sd))+
    geom_point(data=c.df)+
    geom_line(data=c.df,size=1)+
    scale_color_manual("Symbol.UniprotID", values = c.values)+ 
    scale_y_continuous(limits=c(0,30),expand=c(0,0))+
    scale_x_continuous(expand=c(0,0.15))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5 ),
          axis.text.y = element_text(angle = 0, hjust = 0.5 ),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black") )+
    labs(title=paste0(c.symbol," phosphorylation\ncompared to EGFR & HIPK2 phosphorylation"),
         x="Day",y="Relative phosphorylation protein")
  
  setwd(paste0(dir_c_script,"/Example_fig_ReNcell/"))
  #ggsave(plot_DK_phospho,file=c.filename,device = "pdf", height = 6,width = 8,
  #       useDingbats=FALSE)
  print(paste0(n,"-",length(DarkKinase_phospho)))
}
plot_DK_phospho

##################################################################################################################T
# remove tempdir  ------------
##################################################################################################################T
setwd(dir_project)
unlink(tempdir,recursive =T)




