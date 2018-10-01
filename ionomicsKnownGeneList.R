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
  setnames(comparison,c("Gene1id","Gene1org","Gene2id","Gene2org","Relationship","Gene2end","Gene2start","Gene1end","Gene1start","Gene1chr",
                        "Gene2chr","Gene1Defline","Gene2Defline"))
  comparison[,Gene1chr := gsub("Chr|Chr_","",Gene1chr),]
  comparison[,Gene2chr := gsub("Chr|Chr_","",Gene2chr),]
  comparison[,Gene1org := gsub(" ","",Gene1org),]
  comparison[,Gene2org := gsub(" ","",Gene2org),]
  ##make chromosome and base pair columns numeric (this will convert some non-numeric chromsome IDs to NA (with warnings), 
  #they aren't the major chromosomes anyways)
  for(col in c("Gene1chr","Gene2chr","Gene1start","Gene2start","Gene1end","Gene2end")) set(comparison, j=col, 
                                                                                           value=as.numeric(comparison[[col]]))
  return(comparison)
}
############################
dir.create("./knownIonomicsGenesWOrthologs")
ionomicsKnownGenes <- read_csv("IonomicsKnownGenes/ionomics_known_genes_input.csv")
ionomicsKnownGenes<-ionomicsKnownGenes[order(ionomicsKnownGenes$Species, ionomicsKnownGenes$GeneID),]
ionomicsKnownGenesbase<-ionomicsKnownGenes[-88,]
tagalongs<-ionomicsKnownGenes[88,]

###removes whitespaces in lists so they collapse correctly
ionomicsKnownGenesbase$Elements<-gsub(" ","",ionomicsKnownGenesbase$Elements, fixed = TRUE)

##this order may be important, depending on how ionomicskonwngenes gets sorted
possibleOrthologs<-c("A.thaliana","M.truncatula","O.sativa","Z.mays","G.max","S.bicolor","S.viridisearly-release","S.italica")
#sort(possibleOrthologs)
###finding orthologs for the primary genes
OrthologLists<-foreach(i=unique(ionomicsKnownGenesbase$Species), .packages = c('plyr','dplyr','reader')) %do% {
  OrgSub<-ionomicsKnownGenesbase[ionomicsKnownGenesbase$Species==i,]
  files<-list.files("./data/phytozome/current/", recursive = FALSE, full.names = TRUE, pattern = paste0(i,"_"))
  List_Holder<-foreach(f=files, .packages = c('plyr','dplyr','reader'), .combine = surprisinglySimple) %do% {
    geneFile<-geneFilePrep(f)
    Frame_Maker<-foreach(iD=OrgSub$GeneID, .packages = c('plyr','dplyr','reader')) %do% {
      base<-data.frame(OrthologSpecies=geneFile$Gene2org[geneFile$Gene1id==iD], GeneID=geneFile$Gene2id[geneFile$Gene1id==iD])
      if(nrow(base)!=0){base$Inferred<-iD; base$Elements<-OrgSub$Elements[OrgSub$GeneID==iD]}
      return(base)
    }
    return(Frame_Maker)
  }
  names(List_Holder)<-OrgSub$GeneID
  assign(paste0(i,"List"), List_Holder)
}
ConciseOrthologList<-unlist(OrthologLists, recursive = FALSE)
names(OrthologLists)<-unique(ionomicsKnownGenesbase$Species)
inferredOrthologs<-possibleOrthologs[which(!possibleOrthologs %in% ionomicsKnownGenesbase$Species)]

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
  files<-list.files("./data/phytozome/current/", recursive = FALSE, full.names = TRUE, pattern = paste0(i,"_"))
  List_Holder<-foreach(f=files, .packages = c('plyr','dplyr','reader'), .combine = surprisinglySimple) %do% {
    geneFile<-geneFilePrep(f)
    Frame_Maker<-foreach(iD=OrgSub$GeneID, .packages = c('plyr','dplyr','reader')) %do% {
      base<-data.frame(OrthologSpecies=geneFile$Gene2org[geneFile$Gene1id==iD], GeneID=geneFile$Gene2id[geneFile$Gene1id==iD])
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
parsedlistings2<-foreach(i=unique(ionomicsKnownGenesbase$Species), .packages = c('plyr','dplyr','reader')) %do% {
  OrgSubHolder<-listings2[c(1:length(PrimaryInferred[[i]])),]
  listings2<-listings2[-c(1:length(PrimaryInferred[[i]])),]
  return(OrgSubHolder)
}
names(parsedlistings)<-possibleOrthologs
names(parsedlistings2)<-unique(ionomicsKnownGenesbase$Species)

#####splitting the table into seperate files for each organism - inferred stuff will need to be added onto this
OrthologTables<-foreach(i=unique(ionomicsKnownGenesbase$Species), .packages = c('plyr','dplyr','reader')) %do% {
  TableHolder<-ionomicsKnownGenesbase[ionomicsKnownGenesbase$Species==i,]
  TableHolder$`Primary/Inferred`<-"Primary"
  TableHolder$`Inferred from?`<-"NA"
  TableHolder<-cbind(TableHolder,parsedlistings[[i]])
  return(TableHolder)
}
names(OrthologTables)<-unique(ionomicsKnownGenesbase$Species)
okInferrs<-c(parsedlistings2,parsedlistings[-c(1:length(unique(ionomicsKnownGenesbase$Species)))])
InferredTables<-foreach(i=possibleOrthologs, .packages = c('plyr','dplyr','reader')) %do% {
  TableHolder<-Inferredframes[[i]]
  TableHolder<-cbind(TableHolder,distinct(okInferrs[[i]]))
}
names(InferredTables)<-possibleOrthologs
#####making tables for the inferred orthologs#####
Combining<-foreach(i=unique(ionomicsKnownGenesbase$Species), .packages = c('plyr','dplyr','reader')) %do% {
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
  #write.table(Ortho, file = paste0("./IonomicsKnownGenes/",i,"_knownIonomicsGenesWOrthologs.csv"), row.names = FALSE,col.names = TRUE, sep = ",")
  return(Ortho)
}
for(i in inferredOrthologs){
  setnames(InferredTables[[i]], old = possibleOrthologs, new = gsub("(.+)","\\1 orthologs", possibleOrthologs))
  #write.table(InferredTables[[i]], file = paste0("./IonomicsKnownGenes/",i,"_knownIonomicsGenesWOrthologs.csv"), row.names = FALSE,col.names = TRUE, sep = ",")
}