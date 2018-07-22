#source('http://depot.sagebase.org/CRAN.R')
#pkgInstall("synapseClient")
#install.packages("rstudioapi")

library(synapseClient)
library(plyr);library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
`%nin%`<-Negate(`%in%`)


##! to do: divide this script into a) create database table, b) create figures

##################################################################################################################T
# Get user input ------------
##################################################################################################################T
#User input
Synapse_UID<-""
Synapse_password<-""
dir_project<-""
synapse_ids<-c("syn11701291", "syn12617467","syn13549817")#c("synID_1","synID_2",...,"synID_n") ## requires user input
variable_names<-NA#optional: change to c("myname_1","myname_2",...,"myname_n") ## optional user input

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
# Reformat data ------------
##################################################################################################################T
map_unitprot2geneID_RenCells_full<-
  map_unitprot2geneID_RenCells_full%>%mutate(uniprot_id=gsub(" ","",uniprot_id))
Dark_kinases<-Table_001_all_kinases[Table_001_all_kinases$in_IDG_darkkinases==TRUE,]

Phospho_baseline<-Table3_phosphoMS_baseline%>%
  filter(Gene_Symbol!="")%>%
  mutate(uniprot_id=gsub("-.*","",Uniprot_Id))%>%
  merge(.,map_unitprot2geneID_RenCells_full,by="uniprot_id")%>%
  select( - X,-HGNC_Gene_Name,- Site_Position, 
          -Motif,-kmeans_cluster_name, 
          -set1_max_score,-set1_max_score_1,-set1_max_score_2,-set1_max_score_3,
          -set2_max_score, -set2_max_score_1, -set2_max_score_2, -set2_max_score_3)%>%
  gather(parameter,score,Set1_Day_0:Set2_Day_7)%>%arrange(Gene_Symbol)%>%
  mutate(set=gsub("_..._..*","",parameter),
         day=gsub("Set._Day_*","",parameter))%>%
  arrange(Gene_Symbol,day)%>%
  mutate(is_kinase=ifelse(gene_id %in% Table_001_all_kinases$gene_id,TRUE,NA))%>%
  mutate(is_dark_kinase=
           ifelse(gene_id %in% Dark_kinases$gene_id,TRUE,NA))

map_unitprot2geneID_RenCells_full%>%filter(gene_id %nin% Phospho_baseline$gene_id)%>%dim()==c(0,2)

head(Phospho_baseline)

unique(Phospho_baseline$day)
Phospho_baseline$day<-as.numeric(Phospho_baseline$day)
symbols<-unique(Phospho_baseline$Gene_Symbol)

Phospho_baseline_g_day<-Phospho_baseline%>%
  group_by(day, Gene_Symbol, gene_id,uniprot_id)%>%
  summarise(score_mean=mean(score, na.rm=T),
            score_sd=sd(score,na.rm=T),
            is_kinase=unique(is_kinase),
            is_dark_kinase=unique(is_dark_kinase))%>%
  arrange(Gene_Symbol)%>%ungroup()%>%as.data.frame()%>%
  select(gene_id,Gene_Symbol, uniprot_id, day, score_mean,score_sd,is_kinase,is_dark_kinase)

setwd(dir_c_script)
write.csv(Phospho_baseline_g_day,file ="phosphoproteomics_RenMVcellDifferentiation_databasetable.csv",row.names = F )

##################################################################################################################T
# user input for upload  ------------
##################################################################################################################T
file_info_upload<-
  data.frame(Rvariable_name=c("Phospho_baseline_g_day"),
             filenames=c("Table_006_RenMVcell_phosphoproteomics.csv"), 
             file_description=c("Summary table phosphoproteomics of RenMVcell differentiation 15 days"),
             synapse_folder_id=c("syn12176087"),
             synapse_folder_url=c("https://www.synapse.org/Portal.html#!Synapse:syn12176087"),
             date_time=c(Sys.time()%>%as.character()%>%gsub("-","",.)%>%paste0(.,"GMT")),
             Rscript=c(currentscript_insert),
             local_copy_name=c("phosphoproteomics_RenMVcellDifferentiation_databasetable.csv"),
             stringsAsFactors = F)
#synapse_folder<-"syn12176087"

##################################################################################################################T
# create logfile input  ------------
##################################################################################################################T
logfile_name_upload<-Sys.time()%>%as.character()%>%gsub("-","",.)%>%
  gsub(" ",currentscript_insert,.)%>%gsub(":","",.)%>%paste0(.,"GMT")%>%
  paste0("FilesUploaded_",.,".csv")
setwd(dir_log)
write.csv(file_info_upload,logfile_name_upload,row.names = F)

##################################################################################################################T
# upload to synapse  ------------
##################################################################################################################T
synapseLogin(Synapse_UID, Synapse_password)
setwd(tempdir)
for(i in 1:dim(file_info_upload)[1]){
  write.csv(get(file_info_upload$Rvariable_name[i]),file=file_info_upload$filenames[i],row.names = F)
  c.file<-File(paste0(tempdir,"/",file_info_upload$filenames[i]),
               parentId=file_info_upload$synapse_folder_id[i], synapseStore = T)
  synSetAnnotations(c.file)<-list(description = file_info_upload$file_description[i])
  file<-synStore(c.file)
  print(paste0(i,"-",dim(file_info_upload)[1]))
}
synapseLogout()
##################################################################################################################T
# remove tempdir  ------------
##################################################################################################################T
setwd(dir_project)
unlink(tempdir,recursive =T)



