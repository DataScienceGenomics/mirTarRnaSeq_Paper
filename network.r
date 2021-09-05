##MakenetworkPlot
library(EnsDb.Hsapiens.v86)
library(PANTHER.db)
library(biomaRt)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(tidyr)
library(data.table)
library(ggalluvial)
library(stringr)
library(readxl)
library(visNetwork)


rC<-as.data.frame(results[["corrs"]])
unique_Gene<-unique(as.data.frame(rC$V2))
write.table(unique_Gene,"~/Desktop/unique_genes_correlation_part2_Sig.txt",sep="\t",
            quote=FALSE,row.names = FALSE,header=FALSE)
names(unique_Gene)<-c("symbol")

ReactCharc<-read.table("~/Desktop/annot_react.txt",sep="\t",header=TRUE)
Joined_anntSymb<-full_join(ReactCharc,unique_Gene,by="symbol")
write.table(Joined_anntSymb,"~/Desktop/Joined_anntSymb.txt",sep="\t",
            quote=FALSE,row.names = FALSE)
#Tera<-Joined_anntSymb
#Tera$Pathway[which(Tera$Pathway == "")] <- "Not Identified"

Joined_reRead<-read.table("~/Desktop/Joined_anntSymb.txt", sep="\t", header=TRUE)
Tera2<-Joined_reRead%>%select(Pathway,symbol)
Tera2$Pathway[which(is.na(Tera2$Pathway))] <- "Not Identified"
#&36 characterized and 752 genes uncharacterized which is normal considering 50% of 
#genes in the genome are uncharacterized.
Tera2
mRNAres<-as.data.frame(results[["mrna"]])
mRNAres$symbol<-rownames(mRNAres)
names(mRNAres)<-c("TimePoint1vsTimePoint2","TimePoint2vsTimePoint3","TimePoint1vsTimePoint3")
data_long <- gather(mRNAres, TimePoint, FoldChange, TimePoint1vsTimePoint2:TimePoint1vsTimePoint3, factor_key=TRUE)
data_long
names(data_long)<-c("symbol","TimePoint","FoldChange")
joinAlluvi2Conn4<-dplyr::inner_join(Tera2,data_long,by="symbol")
alluFin<-unique(joinAlluvi2Conn4)
is_alluvia_form(as.data.frame(alluFin), axes = 1:3, silent = TRUE)
alluFinNet<-alluFin
uPW<-as.data.frame(unique(alluFinNet$symbol))
uPW1<-uPW[-c(9,12),]
uPW1<-as.data.frame(uPW1)
names(uPW1)<-c("symbol")

#UniquePathwaysInOldAnalysis

uPK<-as.data.frame(unique(alluFinNet$Pathway))
uPK1<-uPK[-c(9,12),]
uPK1<-as.data.frame(uPK1)
names(uPK1)<-c("Pathway")
#ParseReactomeFile:
fKidCol <- rep("text", 15)##All  text
Rf<-read_excel("~/Desktop/ReactomePart2.xlsx",col_types=fKidCol)
RfP<-Rf%>%select(`Pathway name`,`Submitted entities found`)
names(RfP)<-c("Pathway","genes")
awesome <- lapply(1:nrow(RfP), function(i) data.frame(Pathway=RfP[i, "Pathway", drop=TRUE], 
                                                       symbol=unlist(strsplit(RfP[i, "genes", drop=TRUE], ';')))) %>% bind_rows()

#SubsetForGenesCharacterized
#Make Nodes
ljsymb<-inner_join(uPW1,awesome,by="symbol")
symb<-as.data.frame(awesome$symbol)
names(symb)<-c("Nodes")
Path<-as.data.frame(awesome$Pathway)
names(Path)<-c("Nodes")
ljsymb2<-as.data.frame(rbind(symb,Path))
names(ljsymb2)<-c("id")
#Make Edges
EdgeInit<-awesome%>%select("symbol","Pathway")
names(EdgeInit)<-c("from","to")
#Find unique edges
uED<-as.data.frame(unique(EdgeInit$to))
uED1<-as.data.frame(uED)
names(uED1)<-c("uniqueP")
#ImportPathSubPath
PathSub<-read.table("~/Desktop/PathSubPath.txt",sep="\t",header=TRUE)
names(PathSub)<-c("from","to")
FinEdge<-rbind(EdgeInit,PathSub)
#DrawTheNetwork
ljsymb3<-as.data.frame(unique(ljsymb2))
nodes<-as.data.frame(ljsymb3)
edges<-as.data.frame(FinEdge)
#ForAllGens
visNetwork(nodes, edges)

#Now Do this for only sign Reactome Pathways
PKidCol <- rep("text", 3)##All  text
PSigRE<-read_excel("~/Desktop/PSignificantReactomeEnrichmentCorrelation.xlsx",col_types=PKidCol)

names(PSigRE)<-c("SubPathwayPathway","Pathway","genes")
awesomePSig <- lapply(1:nrow(PSigRE), function(i) data.frame(Pathway=PSigRE[i, "Pathway", drop=TRUE], 
                                                      SubPathway=PSigRE[i, "SubPathwayPathway", drop=TRUE],
                                                      symbol=unlist(strsplit(PSigRE[i, "genes", drop=TRUE], ';')))) %>% bind_rows()
#MakeEdges
edge1<-awesomePSig%>%select(SubPathway,Pathway)
names(edge1)<-c("from","to")
edge2<-awesomePSig%>%select(symbol,SubPathway)
names(edge2)<-c("from","to")
bindEdge<-rbind(edge1,edge2)
names(bindEdge)<-c("from","to")
#MakeNodes
nodeA<-awesomePSig%>%select(SubPathway,Pathway)
names(nodeA)<-c("id","group")
nodeB<-awesomePSig%>%select(Pathway)
nodeB$pathway2<-nodeB$Pathway
names(nodeB)<-c("id","group")
nodeC<-awesomePSig%>%select(symbol,Pathway)
names(nodeC)<-c("id","group")

bindedNode<-as.data.frame(rbind(nodeA,nodeB,nodeC))
names(bindedNode)<-c("id","group")

bindedNodeFin<-as.data.frame(unique(bindedNode))
bindedNodeFin$id<-as.character(bindedNodeFin$id)
bindedNodeFin$group<-as.character(bindedNodeFin$group)
bindedNodeFin$label<-as.character(bindedNodeFin$id)

bindedNodeFin <- bindedNodeFin[order(bindedNodeFin$id, bindedNodeFin$group, decreasing=TRUE),]
bindedNodeFinClean <- bindedNodeFin[!duplicated(bindedNodeFin$id),]
bindedNodeFinClean

visNetwork(bindedNodeFinClean, bindEdge,height = "1200px", width = "250%") %>% 
  visOptions(selectedBy = "label", highlightNearest = TRUE, 
             nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = TRUE)

#For Immune Genes Only

ISigRE<-read.table("~/Desktop/immune_import_visGenes.txt",sep="\t",header=TRUE)
ISigRE$SubPathwayPathway<-as.character(ISigRE$SubPathwayPathway)
ISigRE$Pathway<-as.character(ISigRE$Pathway)
ISigRE$genes<-as.character(ISigRE$genes)

names(ISigRE)<-c("SubPathwayPathway","Pathway","genes")
awesomeISig <- lapply(1:nrow(ISigRE), function(i) data.frame(Pathway=ISigRE[i, "Pathway", drop=TRUE], 
                                                             SubPathway=ISigRE[i, "SubPathwayPathway", drop=TRUE],
                                                             symbol=unlist(strsplit(ISigRE[i, "genes", drop=TRUE], ';')))) %>% bind_rows()
#MakeEdges
edge1<-awesomeISig%>%select(SubPathway,Pathway)
names(edge1)<-c("from","to")
edge2<-awesomeISig%>%select(symbol,SubPathway)
names(edge2)<-c("from","to")
bindEdge<-rbind(edge1,edge2)
names(bindEdge)<-c("from","to")
#MakeNodes
nodeA<-awesomeISig%>%select(SubPathway,Pathway)
names(nodeA)<-c("id","group")
nodeB<-awesomeISig%>%select(Pathway)
nodeB$pathway2<-nodeB$Pathway
names(nodeB)<-c("id","group")
nodeC<-awesomeISig%>%select(symbol,Pathway)
names(nodeC)<-c("id","group")

bindedNode<-as.data.frame(rbind(nodeA,nodeB,nodeC))
names(bindedNode)<-c("id","group")

bindedNodeFin<-as.data.frame(unique(bindedNode))
bindedNodeFin$id<-as.character(bindedNodeFin$id)
bindedNodeFin$group<-as.character(bindedNodeFin$group)
bindedNodeFin$label<-as.character(bindedNodeFin$id)

bindedNodeFin <- bindedNodeFin[order(bindedNodeFin$id, bindedNodeFin$group, decreasing=TRUE),]
bindedNodeFinClean <- bindedNodeFin[!duplicated(bindedNodeFin$id),]
bindedNodeFinClean

visNetwork(bindedNodeFinClean, bindEdge,height = "1200px", width = "250%") %>% 
  visOptions(selectedBy = "label", highlightNearest = TRUE, 
             nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = TRUE)

#Add Legened
visNetwork(bindedNodeFinClean, bindEdge,height = "700px", width = "100%") %>% visLegend() %>%
  visOptions(selectedBy = "label", highlightNearest = TRUE, 
             nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = TRUE)


##ExampleCode
nb <- 10
nodes <- data.frame(id = 1:nb, label = paste("Label", 1:nb),
                    group = sample(LETTERS[1:3], nb, replace = TRUE), value = 1:nb,
                    title = paste0("<p>", 1:nb,"<br>Tooltip !</p>"), stringsAsFactors = FALSE)

edges <- data.frame(from = trunc(runif(nb)*(nb-1))+1,
                    to = trunc(runif(nb)*(nb-1))+1,
                    value = rnorm(nb, 10), label = paste("Edge", 1:nb),
                    title = paste0("<p>", 1:nb,"<br>Edge Tooltip !</p>"))
visNetwork(nodes, edges,height = "700px", width = "100%") %>% visLegend() %>%
  visOptions(selectedBy = "label", highlightNearest = TRUE, 
              nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = FALSE)



