####Ionomic Known Gene List####
##Identifying Orthologs##
library(readr)
library(plyr)
library(dplyr)
library(data.table)
library(doParallel)

#########functions##########
orthoPerOrg<-function(AdataFrame,organismName){
  lilList<-paste(as.character(AdataFrame$GeneID[AdataFrame$OrthologSpecies==organismName]), collapse = ",")
}
noDups<-function(AdataFrame){
  UnDupped<-data.frame(GeneID=unique(AdataFrame$GeneID))
  Inferred<-foreach(id=UnDupped$GeneID, .combine = rbind) %do% {
    somthing<-paste(as.character(AdataFrame$Inferred[AdataFrame$GeneID==id]),collapse = ",")
  }
  Elements<-foreach(id=UnDupped$GeneID, .combine = rbind) %do% {
    something<-data.frame(elems=as.character(AdataFrame$Elements[AdataFrame$GeneID==id]))
    nocommas<-as.data.table(unlist(strsplit(as.character(something$elems), split = ",")))
    returning<-paste(as.character(unique(nocommas$V1)), collapse = ",")
  }
  UnDupped<-cbind(UnDupped,Inferred,Elements)
  return(UnDupped)
}
surprisinglySimple<-function(...){mapply(rbind, ..., SIMPLIFY = F)}
geneFilePrep<-function(filePath){
  comparison<-fread(filePath, sep = ",", stringsAsFactors = FALSE)
  if(!length(grep("Homolog",colnames(comparison)))==0){
    setnames(comparison,old = c("Homolog > Gene . Primary Identifier","Homolog > Gene > Organism Name",
                                "Homolog > Ortholog _ Gene . Primary Identifier","Homolog > Ortholog _ Gene > Organism Name"),
             new = c("Gene1id","Gene1org","Gene2id","Gene2org"))
    comparison[,Gene1org := gsub(" ","",Gene1org),]
    comparison[,Gene2org := gsub(" ","",Gene2org),]
  }else{
    setnames(comparison, old = c("Gene name","Organism name","Ortholog gene name","Ortholog organism name"),
             new = c("Gene1id","Gene1org","Gene2id","Gene2org"))
    comparison[,Gene1org := gsub("^([A-Z])","\\1.",Gene1org),]
    comparison[,Gene2org := gsub("^([A-Z])","\\1.",Gene2org),]
  }
  return(comparison)
}
############################
####beginning of pipeline###

#dir.create("./knownIonomicsGenesWOrthologs")
listBase <- read_csv("IonomicsKnownGenes/ionomics_known_genes_input.csv")
listBase<-listBase[order(listBase$Species, listBase$GeneID),]
###checking for duplicate entries
listBase[which(duplicated(listBase$GeneID)),]


###removes whitespaces in lists so they collapse correctly
listBase$Elements<-gsub(" ","",listBase$Elements, fixed = TRUE)

##this order may be important, depending on how ionomicskonwngenes gets sorted
inferredSpecies<-sort(unique(c(gsub("^([^_]+)_.*","\\1",list.files("./data/phytozome/current/", pattern="^[^_]+_[^_]+.csv$", recursive = F)),
                          gsub("^[^_]+_([^_]+).csv","\\1",list.files("./data/phytozome/current/", pattern="^[^_]+_[^_]+.csv$", recursive = F)))))
possibleOrthologs<-c(as.character(unique(listBase$Species)),inferredSpecies[!inferredSpecies %in% unique(listBase$Species)])
###with a short list of species, listing them is the best way to do this. But with a larger list you could search and grep the file names of your phytozome dir
###my phytozome dir includes species not on the list so this way I can make sure it's just the species I want for KIG

###finding orthologs for the primary genes
####don't worry about the warnings generated in this loop, it has to do w the removal of mitochondira or chloroplast genes in the phytozome lists, they will be removed
OrthologLists<-foreach(i=unique(listBase$Species), .packages = c('plyr','dplyr','reader')) %do% {
  OrgSub<-listBase[listBase$Species==i,]
  List_Holder<-foreach(f=possibleOrthologs[!possibleOrthologs==i], .packages = c('plyr','dplyr','reader'), .combine = surprisinglySimple) %do% {
    ###I name my files depending on wich species is the homolog species vs ortholog species: Soy_corn.csv vs corn_soy.csv
    ###pattern string can be edited to search for your own naming convention
    file<-list.files("./data/phytozome/current/", pattern = paste0(f,"_",i,"|",i,"_",f), full.names = T, recursive = F)[1]
    if(is.na(file)){   ##checks to make sure the file exists before proceeding w pipeline
      print(paste0("No file found for ",i," and ",f))
    }else{
      geneFile<-geneFilePrep(file)
      Frame_Maker<-foreach(iD=OrgSub$GeneID, .packages = c('plyr','dplyr','reader')) %do% {
        if(!length(grep(i,geneFile$Gene1org[1]))==0){    #checking to see if the species id's match the first or second species
          base<-data.frame(OrthologSpecies=rep(f,length(which(geneFile$Gene1id==iD))), GeneID=geneFile$Gene2id[geneFile$Gene1id==iD])
          if(nrow(base)!=0){base$Inferred<-iD; base$Elements<-OrgSub$Elements[OrgSub$GeneID==iD]}
          return(base)
        }else{
          base<-data.frame(OrthologSpecies=rep(f,length(which(geneFile$Gene2id==iD))), GeneID=geneFile$Gene1id[geneFile$Gene2id==iD])
          if(nrow(base)!=0){base$Inferred<-iD; base$Elements<-OrgSub$Elements[OrgSub$GeneID==iD]}
          return(base)
        }
      }
    }
    return(Frame_Maker)
  }
  names(List_Holder)<-OrgSub$GeneID
  assign(paste0(i,"List"), List_Holder)
}
ConciseOrthologList<-unlist(OrthologLists, recursive = FALSE)
names(OrthologLists)<-unique(listBase$Species)
inferredOrthologs<-possibleOrthologs[which(!possibleOrthologs %in% listBase$Species)]

InferredOrthologsIDs<-foreach(i=possibleOrthologs, .packages = c('plyr','dplyr','reader')) %do% {
  backside<-foreach(q=ConciseOrthologList, .combine = rbind) %do% {
    line<-q[q$OrthologSpecies==i,]
    #return(line)
  }
}
InferredOrthologsIDs<-lapply(InferredOrthologsIDs,noDups)
names(InferredOrthologsIDs)<-possibleOrthologs

n=1
Inferredframes<-foreach(i=InferredOrthologsIDs, .packages = c('plyr','dplyr','reader')) %do% {
  tableHolder<-data.frame(Species=possibleOrthologs[n], GeneID=i$GeneID, GeneName=NA, Elements=i$Elements, `Citation(s) - DOI only`=NA, 
                          `Closest ortholog species`=NA, Tissue=NA,Comments=NA, `Primary/Inferred`="Inferred", 
                          `Inferred from?`= i$Inferred)
  n=n+1
  return(tableHolder)
}
names(Inferredframes)<-possibleOrthologs

####finding orthologs for the inferred genes
InferredLists<-foreach(i=possibleOrthologs, .packages = c('plyr','dplyr','reader')) %do% {
  OrgSub<-Inferredframes[[i]]
  #files<-list.files("./data/phytozome/current/", recursive = FALSE, full.names = TRUE, pattern = paste0(i,"_"))
  List_Holder<-foreach(f=possibleOrthologs[!possibleOrthologs==i], .packages = c('plyr','dplyr','reader'), .combine = surprisinglySimple) %do% {
    file<-list.files("./data/phytozome/current/", pattern = paste0(f,"_",i,"|",i,"_",f), full.names = T, recursive = F)[1]
    if(is.na(file)){   ##checks to make sure the file exists before proceeding w pipeline
      print(paste0("No file found for ",i," and ",f))
    }else{
      geneFile<-geneFilePrep(file)
      Frame_Maker<-foreach(iD=OrgSub$GeneID, .packages = c('plyr','dplyr','reader')) %do% {
        if(!length(grep(i,geneFile$Gene1org[1]))==0){
          base<-data.frame(OrthologSpecies=rep(f, length(which(geneFile$Gene1id==iD))), GeneID=geneFile$Gene2id[geneFile$Gene1id==iD])
        }else{
          base<-data.frame(OrthologSpecies=rep(f, length(which(geneFile$Gene2id==iD))), GeneID=geneFile$Gene1id[geneFile$Gene2id==iD])
        }
      }
    }
    return(Frame_Maker)
  }
  names(List_Holder)<-OrgSub$GeneID
  assign(paste0(i,"List"), List_Holder)
}
names(InferredLists)<-possibleOrthologs

###combining the primary and inferred lists for the organisms with no primary entries
###have to keep the primary and inferred entries for an organism seperate or they won't merge back to the table correctly
Overall<-c(OrthologLists,InferredLists[inferredOrthologs])
PrimaryInferred<-InferredLists[!possibleOrthologs %in% inferredOrthologs]
###making a list of orthologs in each species for all the primary and non-primary, inferred organisms
n=1
listings<-foreach(i=Overall, .packages = c('plyr','dplyr','reader'), .combine = 'rbind') %do% {
  HomologIDs<-names(i)
  SomeTable<-foreach(p=possibleOrthologs, .combine = 'cbind') %do% {
    as.data.frame(sapply(i, orthoPerOrg, organismName=p))
  }
  SomeTable[,n]<-HomologIDs
  colnames(SomeTable)<-possibleOrthologs
  n=n+1
  return(SomeTable)
}
###squishing down to a list for all the primary-inferred organisms
n=1
listings2<-foreach(i=PrimaryInferred, .packages = c('plyr','dplyr','reader'), .combine = 'rbind') %do% {
  HomologIDs<-names(i)
  SomeTable<-foreach(p=possibleOrthologs, .combine = 'cbind') %do% {
    as.data.frame(sapply(i, orthoPerOrg, organismName=p))
  }
  SomeTable[,n]<-HomologIDs
  colnames(SomeTable)<-possibleOrthologs
  n=n+1
  return(SomeTable)
}

parsedlistings<-foreach(i=possibleOrthologs, .packages = c('plyr','dplyr','reader')) %do% {
  OrgSubHolder<-listings[c(1:length(Overall[[i]])),]
  listings<-listings[-c(1:length(Overall[[i]])),]
  return(OrgSubHolder)
}
parsedlistings2<-foreach(i=unique(listBase$Species), .packages = c('plyr','dplyr','reader')) %do% {
  OrgSubHolder<-listings2[c(1:length(PrimaryInferred[[i]])),]
  listings2<-listings2[-c(1:length(PrimaryInferred[[i]])),]
  return(OrgSubHolder)
}
names(parsedlistings)<-possibleOrthologs
names(parsedlistings2)<-unique(listBase$Species)

#####splitting the table into seperate files for each organism - inferred stuff will need to be added onto this
OrthologTables<-foreach(i=unique(listBase$Species), .packages = c('plyr','dplyr','reader')) %do% {
  TableHolder<-listBase[listBase$Species==i,]
  TableHolder$`Primary/Inferred`<-"Primary"
  TableHolder$`Inferred from?`<-"NA"
  TableHolder<-cbind(TableHolder,parsedlistings[[i]])
  return(TableHolder)
}
names(OrthologTables)<-unique(listBase$Species)
okInferrs<-c(parsedlistings2,parsedlistings[-c(1:length(unique(listBase$Species)))])
InferredTables<-foreach(i=possibleOrthologs, .packages = c('plyr','dplyr','reader')) %do% {
  TableHolder<-Inferredframes[[i]]
  TableHolder<-cbind(TableHolder,distinct(okInferrs[[i]]))
}
names(InferredTables)<-possibleOrthologs
#####making tables for the inferred orthologs#####
summaryTable<-data.frame(Org=character(),GeneNum=integer(),NumPrimary=integer(),`NumPrimary/Inferred`=integer(),NumInferred=integer(),
                         NoOrthologs=integer())
Combining<-foreach(i=unique(listBase$Species), .packages = c('plyr','dplyr','reader')) %do% {
  Ortho<-OrthologTables[[i]]
  Infer<-InferredTables[[i]]
  colnames(Infer)<-colnames(Ortho)
  for(iD in Infer$GeneID){
    if(iD %in% Ortho$GeneID){
      Ortho$`Primary/Inferred`[Ortho$GeneID==iD]<-"Primary/Inferred"
      Ortho$`Inferred from?`[Ortho$GeneID==iD]<-as.character(Infer$`Inferred from?`[Infer$GeneID==iD])
    }else{
      Ortho<-rbind(Ortho,Infer[Infer$GeneID==iD,])
    }
  }
  setnames(Ortho, old = possibleOrthologs, new = gsub("(.+)","\\1 orthologs", possibleOrthologs))
  Ortho<-Ortho[ , -which(names(Ortho)==paste0(i," orthologs"))]
  summaryTableEntry<-data.frame(Org=i,GeneNum=nrow(Ortho),NumPrimary=length(grep("^Primary$",Ortho$`Primary/Inferred`)),
                           `NumPrimary/Inferred`=length(grep("Primary/Inferred",Ortho$`Primary/Inferred`)),
                           NumInferred=length(grep("^Inferred$",Ortho$`Primary/Inferred`)),
                           NoOrthologs=length(which(apply(Ortho[,grep("orthologs",colnames(Ortho))],1,function(x) if(sum(is.na(x)|x=="")==length(possibleOrthologs)){return(TRUE)} else{FALSE}))))
  summaryTable<-rbind(summaryTable,summaryTableEntry)
  write.table(Ortho, file = paste0("./IonomicsKnownGenes/knownIonomicsGenesWOrthologs/",i,"_knownIonomicsGenesWOrthologs.csv"), row.names = FALSE,col.names = TRUE, sep = ",")
  return(Ortho)
}
for(i in inferredOrthologs){
  setnames(InferredTables[[i]], old = possibleOrthologs, new = gsub("(.+)","\\1 orthologs", possibleOrthologs))
  setnames(InferredTables[[i]], old=c("Primary.Inferred"), new=c("Primary/Inferred"))
  InferredTables[[i]]<-InferredTables[[i]][ , -which(names(InferredTables[[i]])==paste0(i," orthologs"))]
  summaryTableEntry<-data.frame(Org=i,GeneNum=nrow(InferredTables[[i]]),NumPrimary=length(grep("^Primary$",InferredTables[[i]]$`Primary/Inferred`)),
                                `NumPrimary/Inferred`=length(grep("Primary/Inferred",InferredTables[[i]]$`Primary/Inferred`)),
                                NumInferred=length(grep("^Inferred$",InferredTables[[i]]$`Primary/Inferred`)),
                                NoOrthologs=length(which(apply(InferredTables[[i]][,grep("orthologs",colnames(InferredTables))],1,function(x) if(sum(is.na(x)|x=="")==length(possibleOrthologs)){return(TRUE)} else{FALSE}))))
  summaryTable<-rbind(summaryTable,summaryTableEntry)
  write.table(InferredTables[[i]], file = paste0("./IonomicsKnownGenes/knownIonomicsGenesWOrthologs/",i,"_knownIonomicsGenesWOrthologs.csv"), row.names = FALSE,col.names = TRUE, sep = ",")
}
write.table(summaryTable, file = "./IonomicsKnownGenes/CurrentIonomicsSummaryTable.csv", row.names = FALSE, col.names = TRUE, sep = ",")
