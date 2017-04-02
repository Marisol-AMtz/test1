###############
## LIBRARIES ##
###############
library(shiny)
#library(ggplot2)
library(Gviz)
library(rtracklayer)
library(GenomicFeatures)
#library(Rsamtools)
library(biomaRt)
options(ucscChromosomeNames=FALSE)
options(scipen=999)
#####################
### Load datasets ###
#####################
load("./resources/Mphg/PeaksAnnotations_A_object.Rdata")
load("./resources/Mphg/PeaksAnnotations_B_object.Rdata")
load("./resources/Mphg/PeakClassificationData.Rdata")
load("./resources/Mphg/ReducedMHC_AltPeaks.Rdata")
load("./resources/TxChr6_genesymbol.Rdata")

HaplotypeId2Symbol<-data.frame(Symbol=c("PGF","COX","APD","DBB","MANN","SSTO","QBL","MCF","Alt"),haplotype=c("chr6","chr6_GL000251v2_alt","chr6_GL000250v2_alt","chr6_GL000252v2_alt","chr6_GL000253v2_alt","chr6_GL000256v2_alt","chr6_GL000255v2_alt","chr6_GL000254v2_alt","Alt"))
##################
## Prepare Data ##
##################

## ANNOTATED PEAK TABLE
MHC_Annot<-rbind(MHC_Annot_A,MHC_Annot_B)
DeepExam<-MHC_Annot[,c("seqnames","width","Name","distancetoFeature","gene_symbol")]
DeepExam<-DeepExam[with(DeepExam, order(gene_symbol, Name)), ]


## REDUCED ALT PEAKS Annotated
Reduced_MHC_Annot_B$Condition<-"B";Reduced_MHC_Annot_A$Condition<-"A"
Reduced_MHC_Annot<-rbind(Reduced_MHC_Annot_A,Reduced_MHC_Annot_B)
Reduced_MHC_Annot<-GRanges(seqnames="6",ranges=IRanges(start=Reduced_MHC_Annot$start,end=Reduced_MHC_Annot$end), peak= paste(as.vector(Reduced_MHC_Annot$sample),as.vector(Reduced_MHC_Annot$Condition),sep = "_"), Haplotype="Alt")

## LIFTED PEAK TABLE 
QryByPosTable<-PeakTable
b1<-subset(QryByPosTable, LStart==-2)
b1$start<-b1$Start; b1$end<-b1$End
QryByPosTable<-subset(QryByPosTable, LStart>0)
QryByPosTable$start<-QryByPosTable$LStart; QryByPosTable$end<-QryByPosTable$LEnd
QryByPosTable<-rbind(b1,QryByPosTable)
QryByPosTable<-GRanges(seqnames="6",ranges=IRanges(start=QryByPosTable$start,end=QryByPosTable$end), peak=QryByPosTable$PeakName,Haplotype=QryByPosTable$Haplotype)


## LOAD BIOMART
bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
fm <- Gviz:::.getBMFeatureMap()

## LOAD COORD AXIS
Gaxis <- GenomeAxisTrack()

#################
### Functions ###
#################
## Collapse ALT
reducealt<-function(qsample,sampledf){
  sampledf<-subset(sampledf, grepl(qsample,peak))
  altdf<-subset(sampledf,Haplotype !="chr6" )
  names(altdf)<-NULL
  altdf<-reduce(altdf)
  altdf<-as.data.frame(altdf)
  altdf$peak<-qsample
  altdf$Haplotype<-"Alt"
  return(altdf)
}
###   ###   ###
## Collapse ALT table
reducealt_info<-function(qsample,sampledf){
  sampledf<-subset(sampledf, grepl(qsample,PeakName) & Haplotype !="chr6")
  if(dim(sampledf)[1]==0){return()}
  sampledf<-GRanges(seqnames="Alt",ranges=IRanges(start=sampledf$LStart,end=sampledf$LEnd),gene=sampledf$gene_symbol,peakname=sampledf$PeakName,distanceF=sampledf$distancetoFeature)
  names(sampledf)<-NULL
  samplered<-reduce(sampledf)
  names(samplered)<-NULL
  listrange<-lapply(samplered, SubsetSubjectByOverlaps,sbj=sampledf)
  listrange<-do.call("rbind", listrange)
  samplered$gene<-as.vector(listrange$Feature_Annotated)
  samplered$Peaks_Contained<-as.vector(listrange$Peaks_Contained)
  samplered$sample<-qsample
  samplered$DisFeat<-as.vector(listrange$Dist_Feat)
  samplered<-as.data.frame(samplered)
  return(samplered)
}


SubsetSubjectByOverlaps<-function(qry,sbj){
  subjecthitrange<-sbj[subjectHits(findOverlaps(qry,sbj))]
  gene<-as.vector(subjecthitrange$gene)
  gene<-paste(unique(sort(gene)),collapse=";")
  peaksnames<-as.vector(subjecthitrange$peakname)
  peaksnames<-paste(unique(sort(peaksnames)),collapse=";")
  distanceF<-as.vector(subjecthitrange$distanceF)
  distanceF<-paste(unique(sort(distanceF)),collapse=";")
  info<-data.frame(Feature_Annotated=gene, Peaks_Contained=peaksnames,Dist_Feat=distanceF)
  return(info)
}

###   ###   ###
## Query by Coordinates
GetATACtrackPOS<-function(MySample, Myrange,Mylett,MyDF){
  caselifted<-subset(MyDF, grepl(paste0(MySample,"_",Mylett),peak) )
  caselifted<-subsetByOverlaps(caselifted,Myrange)
  if(length(caselifted)==0){return()}
  mybg<-"white"
  if(Mylett=="B"){mybg<-"gray96"}
  temp<-as.vector(HaplotypeId2Symbol$Symbol)
  names(temp)<-as.vector(HaplotypeId2Symbol$haplotype)
  tempcol<-c("cadetblue3","chartreuse3","hotpink3","salmon","steelblue3","darkorchid3","darkgoldenrod2","khaki1","orange")
  names(tempcol)<-as.vector(HaplotypeId2Symbol$haplotype)
  caselifted$group<-as.vector(temp[as.vector(caselifted$Haplotype)])
  caselifted$colour<-as.vector(tempcol[as.vector(caselifted$Haplotype)])
  
  dtrack<- AnnotationTrack(range=caselifted,name =MySample, chromosome = "6",cex.title=.56,cex.axis=.28,background.panel =mybg ,group=as.vector(caselifted$group), groupAnnotation = "group",cex.group=1.2, fill=as.vector(caselifted$colour))
  return(dtrack)
}


###################
### START SHINY ###
###################
shinyServer(function(input, output){
      PlotHeight=reactive(return(as.numeric(input$mywindowsize)))
      output$Peaks<-renderPlot({
        req(input$plotme)
        isolate({
          chrplot<-"6"
          SAMPLES<-as.list(input$Sample)
          if(length(input$Haps)==1){
            if(as.vector(input$Haps)=="PGF"){
              #PGF selected
              QryByPosTable<-PeakTable
              QryByPosTable<-subset(QryByPosTable,Haplotype=="chr6")
              row.names(QryByPosTable)<-NULL
              QryByPosTable<-GRanges(seqnames = chrplot,ranges=IRanges(start=QryByPosTable$Start, end= QryByPosTable$End),peak=QryByPosTable$PeakName,Haplotype=QryByPosTable$Haplotype)
            }else{
              #One Alt selected, do not lift
              chrplot<-as.vector(subset(HaplotypeId2Symbol,Symbol==as.vector(input$Haps))$haplotype)
              QryByPosTable<-PeakTable
              QryByPosTable<-subset(QryByPosTable,Haplotype==chrplot)
              row.names(QryByPosTable)<-NULL
              QryByPosTable<-GRanges(seqnames = chrplot,ranges=IRanges(start=QryByPosTable$Start, end= QryByPosTable$End),peak=QryByPosTable$PeakName,Haplotype=QryByPosTable$Haplotype)
            }
          }else{
            QryByPosTable<-PeakTable
            b1<-subset(QryByPosTable, LStart==-2)
            b1$start<-b1$Start; b1$end<-b1$End
            QryByPosTable<-subset(QryByPosTable, LStart>0)
            QryByPosTable$start<-QryByPosTable$LStart; QryByPosTable$end<-QryByPosTable$LEnd
            QryByPosTable<-rbind(b1,QryByPosTable)
            row.names(QryByPosTable)<-NULL
            
            #Multiple Haplotypes selected, lifting coords
            
            temp<-as.vector(HaplotypeId2Symbol[HaplotypeId2Symbol$Symbol%in%as.vector(input$Haps),]$haplotype)
            QryByPosTable<-subset(QryByPosTable,Haplotype%in%temp)
            QryByPosTable<-GRanges(seqnames=chrplot,ranges=IRanges(start=QryByPosTable$start,end=QryByPosTable$end), peak=QryByPosTable$PeakName,Haplotype=QryByPosTable$Haplotype)
            if(input$redcd=="yes"){
              smpls<-as.list(paste0(as.vector(unlist(SAMPLES)),"_A"))
              altrgA<-lapply(smpls,FUN=reducealt,sampledf=QryByPosTable)
              altrgA<- do.call("rbind", altrgA)
              output$redcd<-renderTable(tail(as.data.frame(QryByPosTable),5))
              
              smpls<-as.list(paste0(as.vector(unlist(SAMPLES)),"_B"))
              altrgB<-lapply(smpls,FUN=reducealt,sampledf=QryByPosTable)
              altrgB<- do.call("rbind", altrgB)
              altrg<-rbind(as.data.frame(subset(QryByPosTable,Haplotype=="chr6")),altrgA,altrgB)
              QryByPosTable<-altrg
              row.names(QryByPosTable)<-NULL
              QryByPosTable<-GRanges(seqnames=chrplot,ranges=IRanges(start=QryByPosTable$start,end=QryByPosTable$end), peak=QryByPosTable$peak,Haplotype=QryByPosTable$Haplotype)

            }
          }
          
          qryrange<-GRanges(seqnames=chrplot,ranges=IRanges(start=as.numeric(input$Start),end=as.numeric(input$End)))
          SampleTracks_A<-lapply(SAMPLES,FUN = GetATACtrackPOS, Myrange=qryrange,Mylett="A",QryByPosTable)
          SampleTracks_B<-lapply(SAMPLES,FUN = GetATACtrackPOS, Myrange=qryrange,Mylett="B",QryByPosTable)
          TracksBySample<-list()
          n<-1
          for( i in seq(1,length(SAMPLES))){
            TracksBySample[n]<-SampleTracks_A[i]; TracksBySample[n+1]<-SampleTracks_B[i]
            n=n+2
          }
          TracksBySample<-TracksBySample[!sapply(TracksBySample, is.null)]
          
          if(chrplot=="6"){
            biomartTrack <- BiomartGeneRegionTrack(genome="hg38", chromosome="6", start = start(qryrange), end = end(qryrange),name = "ENSEMBL", biomart=bm, featureMap=fm)
            plotTracks(c(Gaxis,TracksBySample),from =start(qryrange),to =end(qryrange),chromosome =chrplot,fontcolor.feature = "darkblue",transcriptAnnotation="symbol", stacking="squish",background.title = "white",fontcolor.title="black",cex.title=.8,cex.group=.7)
            
          }else{
            btrack_FLAG=0
            Txdb_hg38_genes<-subset(TxChr6,seqnames==chrplot)
            Txdb_hg38_genes<-GRanges(seqnames = Txdb_hg38_genes$seqnames, ranges=IRanges(start=Txdb_hg38_genes$start,end=Txdb_hg38_genes$end,names=Txdb_hg38_genes$GeneID),strand=Txdb_hg38_genes$strand,gene=Txdb_hg38_genes$gene_symbol)
            grtrack<-GeneRegionTrack(Txdb_hg38_genes,chromosome = chrplot, name = "Annotation", transcriptAnnotation = "gene", shape='arrow',labelPos="below",fill="lightsteelblue1",cex.label=3.3,fontsize.group=15.3)
          
            plotTracks(c(Gaxis,TracksBySample,grtrack),from =start(qryrange),to =end(qryrange),chromosome =chrplot,fontcolor.feature = "darkblue",transcriptAnnotation="gene", stacking="squish",background.title = "white",fontcolor.title="black",cex.title=.8,cex.group=.7)
          }
          })
        })
      output$Peaks.ui <- renderUI({plotOutput("Peaks", height = PlotHeight())})
      ###########################################################
       output$Genes<-renderPlot({ 
         req(input$plotme)
         isolate({
           qryrange<-GRanges(seqnames="6",ranges=IRanges(start=as.numeric(input$Start),end=as.numeric(input$End)))
           biomartTrack <- BiomartGeneRegionTrack(genome="hg38", chromosome="6", start = start(qryrange), end = end(qryrange),name = "ENSEMBL", biomart=bm, featureMap=fm)
           
           if((length(input$Haps)!=1 ) || (as.vector(input$Haps)=="PGF")){
           plotTracks(c(Gaxis,biomartTrack),from=start(qryrange),to=end(qryrange),chromosome = "6",fontcolor.feature = "darkblue",genome="hg38",transcriptAnnotation="symbol", stacking="squish",genome="hg38",background.title = "white",fontcolor.title="black",cex.title=.8,cex.group=.7)
           }
         })
       })
       ###########################################################
       ###########################################################
  
    output$TablePeaks<-renderTable({
        req(input$querytable)
        isolate({
          AltPeaksTab<-AllAlternativePeaks
          row.names(AltPeaksTab)<-NULL
          if(!("All"%in%input$HapsT)){
            qchr<-as.vector(subset(HaplotypeId2Symbol,Symbol%in%as.vector(input$HapsT))$haplotype)
            AltPeaksTab<-subset(AltPeaksTab,Haplotype%in%qchr)
            row.names(AltPeaksTab)<-NULL
          }
          
          if(grepl("Macrophages",input$State)){
            tmp<-as.vector(AltPeaksTab$PeakName)
            tmp<-unlist(strsplit(tmp,"_peak"))
            tmp<-tmp[seq(1,length(tmp),2)]
            tmp2<-rep("Naive",length(tmp))
            tmp2[grepl("_B",tmp)]<-"IFNg"
            names(tmp2)<-tmp
            AltPeaksTab$Condition<-as.vector(tmp2)
            if(input$State!="Macrophages-All"){
              # separate in active or inactive
              if(input$State=="Macroghages-Naive"){stt="Naive"
              }else{stt="IFNg"}
              AltPeaksTab<-subset(AltPeaksTab,Condition==stt)
              row.names(AltPeaksTab)<-NULL
            }else{ p("DO NOT SUPPOR LCL DATA YET")}
          }
          
          if(!("All"%in%input$Cats)){
            tmp<-input$Cats
            tmp<-tmp[tmp!="All"]
            AltPeaksTab<-subset(AltPeaksTab,Category%in%tmp)
            row.names(AltPeaksTab)<-NULL
          }
          if(!(input$StartT%in% c("-",""))){
            coords<-as.numeric(input$StartT)
            if(input$coordsys=="yes"){AltPeaksTab<-subset(AltPeaksTab,LStart>=coords)
            }else{AltPeaksTab<-subset(AltPeaksTab,Start>=coords)}
            row.names(AltPeaksTab)<-NULL
          }
          if(input$EndT!="-"){
            coords<-as.numeric(input$EndT)
            if(input$coordsys=="yes"){ AltPeaksTab<-subset(AltPeaksTab,LEnd<=coords)
            }else{AltPeaksTab<-subset(AltPeaksTab,End<=coords)}
            row.names(AltPeaksTab)<-NULL
          }
            #Cross with Peak Annotation
            tmp<-MHC_Annot[,c("Name","GeneID","gene_symbol","distancetoFeature")]
            colnames(tmp)<-c("PeakName","GeneID","gene_symbol","distancetoFeature")
            AltPeaksTab<-merge(AltPeaksTab,tmp,by="PeakName",all.x=T)
            row.names(AltPeaksTab)<-NULL
            
          if(input$qrygene!="All"){
            if(grepl("ENSG0",input$qrygene,ignore.case = T)){
              AltPeaksTab<-subset(AltPeaksTab,grepl(input$qrygene,GeneID,,ignore.case = T))
            }else{ AltPeaksTab<-subset(AltPeaksTab,grepl(input$qrygene,gene_symbol,ignore.case = T))}
            row.names(AltPeaksTab)<-NULL
          }
            #AltPeaksTab<-subset(AltPeaksTab,Category!="Failed")
            #row.names(AltPeaksTab)<-NULL

          if(input$collapsetable== "yes" ){
             if(grepl("Macrophages",input$State)){
                 samples2<-as.vector(AltPeaksTab$PeakName)
                 samples2<-unlist(strsplit(samples2,"_peak_"))
                 samples2<-as.vector(samples2[seq(1,length(samples2),2)])
                 samples2<-as.list(as.vector(unique(samples2)))
                 AltPeaksTab<-lapply(samples2,FUN=reducealt_info,sampledf=AltPeaksTab)
                 AltPeaksTab<-do.call("rbind", AltPeaksTab)
                 output$dimtable<-renderText("Number of peaks:")
                 AltPeaksTab<-AltPeaksTab[,c("seqnames","start","end","gene","DisFeat","Peaks_Contained")]
                
                 colnames(AltPeaksTab)<-c("Haplotype","Lifted_Start","Lifted_End","Annotated_Feature","Distance_To_Feature","Peaks_Contained")
                 row.names(AltPeaksTab)<-NULL
                 
                 return(AltPeaksTab)
             }else{tmp<-("LCL data not supported yet")}
          }else{
            AltPeaksTab<-AltPeaksTab[,c("PeakName","Haplotype","Start","End","LStart","LEnd","HeterozygousPGFTo","Category","GeneID","gene_symbol","distancetoFeature")]
            colnames(AltPeaksTab)<-c("Peak_Name","Haplotype","Start","End","Lifted_Start","Lifted_End","Matching_PGF_Peak","Category","Annotated_Feature(ENSG)","Annotated_Feature","Distance_to_Feature")
            row.names(AltPeaksTab)<-NULL
            return(AltPeaksTab)
          }
          
          output$dimtable<-renderText("Number of peaks: ",paste0(dim(AltPeaksTab)[1]))
          
        })
    })
})
