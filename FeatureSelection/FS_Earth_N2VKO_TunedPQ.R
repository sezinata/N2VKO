library(RANKS)
library("e1071")
library(caret)
library(ggplot2)

library(clusterSim)
library(plyr)
library(earth)

datasets <- c("IntAct", "NCBI", "STRING")

diseases <- c( "prostateCancer","breastCancer","diabetesMellitus","obesity","lungCancer","alzheimer", "colorectalCancer")
#p q tuned values
p_q_max<-matrix(c(1,4,4,2, 0.25,1,0.5,4, 0.25,1,4,2, 0.25,0.5,1,2, 2,4,1,0.25, 2,2,4,4, 4,2,4,4
),nrow=7,ncol=4,byrow = TRUE)
memory.limit(size=50000)
sized<-length(datasets)
k=5
dis=1
z=1
j=1
path1<-"D:\\"
pathWrite<-"D:\\"

for(dis in 1:7){
  disease=diseases[dis]
  p_q_max_disease<-p_q_max[dis,]
  labelfile<-paste0(path1,"2018OMIMDiseaseLabels/")
  
  for(z in 1:2){
    dataset=datasets[z]
    print(paste0("Currently handling: ",disease," ",dataset))
    indicep=z*2-1
    indiceq=z*2  #for intact first two columns gives p-q values, for ncbi the last 2 cols gives p-q values
    pvalue<- p_q_max_disease[indicep]
    qvalue<- p_q_max_disease[indiceq]
    
    dataset.name <-paste0(path1, dataset,"MatrixforSVMvector.txt")
    keyw = read.table(dataset.name,header =FALSE,sep="\t")
    mapread <-paste0(path1, dataset ,"KeywordsandIndices.txt")
    map= read.table(mapread,header =FALSE,sep="\t")
    keyw <- as.matrix(keyw)
    columnNames=as.vector(map$V1)
    names<-append("Length",columnNames)
    colnames(keyw)<-names
    var_names <- colnames(keyw)
    keywords<-list()
    for(i in var_names) {
      keywords[[i]]= keyw[,i]  #each keyword individually a list element with its values for all proteins
    }
    
    dataset.name <- paste0(path1,"node2vec-master\\emb\\TunePQ\\node2vec5_5_10\\", dataset,"-",pvalue,"-",qvalue,".emb")
    data.name = read.table(dataset.name,header =FALSE,sep=" ",skip=1) #for original embedding results skips first line
    Kemb <- as.matrix(data.name)
    Kemb <-Kemb[order(Kemb[,1]),]  #sort is AUTOMATICALLY DONE HERE
    keywordsEmbeddingRows<-keyw[Kemb[,1]+1,1:length(keywords)]#+1 EMBEDDINGS START WITH 0 INDEX,VALUES WITH SAME INDICES ARE CHOSEN  
    K<-cbind(Kemb,keywordsEmbeddingRows)

    n <- nrow(K)
    p <- numeric(n)
    firstcolumn<-K[,1]  
    for (j in 1:k){
      importance<-NULL
      RESULTSTEMP<-NULL
      testIndices <- paste0(path1,"TestIndices/", dataset,disease, j,"indicesofTestSet.txt")  
      testIndices= read.table(testIndices,header =FALSE,sep="\n")
      testIndices<-as.array(testIndices$V1-1)# node2vec ids/indices starts with 0
      
      xtrainWithindices<-K[ !(firstcolumn%in% testIndices),]  # use !(%in)
      xtrain<-xtrainWithindices[,-1]#removes first column
      
      xtestWithindices<-K[ firstcolumn%in% testIndices,]
      xtest<-xtestWithindices[,-1]
      
      
      allLabels <-paste0(labelfile , dataset , "SVMclassfor",disease,".txt")
      allLabel= read.table(allLabels,header =FALSE,sep="\n")
      allIndices=c()
      for (i in 1:length(allLabel$V1)) { 
        allIndices[i]=i-1
      }
      allLabel<-cbind(allLabel,allIndices)
      testlabelswithindices<-allLabel[ allLabel$allIndices%in%xtestWithindices[,1] ,]
      testlabel<-testlabelswithindices[,1]
      testlabel<-as.matrix(testlabel)
      
      trainlabelswithindices<-allLabel[ allLabel$allIndices%in%xtrainWithindices[,1] ,]
      trainlabel<-trainlabelswithindices[,1]
      

      ytrain <- as.matrix(trainlabel)
      nclasses <- ncol(ytrain)
      trainIndices<-xtrainWithindices[,1]
      trainIndices_vector<-as.vector(trainIndices)
      
      trainMean <- apply(xtrain,2,mean)
      trainSd <- apply(xtrain,2,sd)
      xtrain<- sweep(sweep(xtrain, 2L, trainMean), 2, trainSd, "/") # using the default, "-" to subtract mean column-wise   
      # ## centered AND scaled
      xtest<-sweep(sweep(xtest, 2L, trainMean), 2, trainSd, "/")
      trainsize<-nrow(xtrain)
      merged<-rbind(xtrain,xtest)
      merged<-merged[,!(colSums(is.na(merged)))]
      xtrain<-merged[1:trainsize,]
      xtest<-merged[(trainsize+1):nrow(merged),]

      #>>>>>>>>>>REMOVE ZERO VARIANCE VARIABLES
      zv <- apply(xtrain, 2, function(x) length(unique(x)) == 1)
      dfr <- xtrain[, !zv]
      #>>>>>>>>>>>>>>>>>
      correlationMatrix <- cor(dfr[,1:ncol(dfr)],use="complete.obs")
      # summarize the correlation matrix
      #print(correlationMatrix)
      # find attributes that are highly corrected (ideally >0.75)
      highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
      # print indexes of highly correlated attributes
      print(highlyCorrelated)
      
       xtrain<-xtrain[,-highlyCorrelated ]

      if (nrow(xtrain) != nrow(ytrain)) 
        stop("train data and label matrices do not agree")
      colnames(ytrain) <- "Labels"
      DATA=cbind(xtrain, ytrain)
      DATA<-as.data.frame(DATA)
      names(DATA)
      xtest<-as.data.frame(xtest)
      xtrain<-as.data.frame(xtrain)
      
      var_names<-colnames(xtrain)
      marsModel <- earth(Labels~. , data = DATA)  # fit lm() model
      importance <- evimp (marsModel) # calculate relative importance scaled to 100
      importance<-as.matrix.data.frame(importance)
      loop=1
      for (loop in 1:length(importance[,1])){  #found keywords length. 3 columns are added as a row to the resulttemp data frame
        resultkeywordname=var_names[importance[loop,1]]
        resultgcv=importance[loop,4] #Generalized cross validation
        resultrss=importance[loop,6]  #residual sum of squares 
        RESULTSTEMP<-rbind(RESULTSTEMP,data.frame(resultkeywordname,resultgcv,resultrss),stringsAsFactors=FALSE)#data.frame usage is better than cbind. cbind merges double values as char arrays
      }
      RESULTSorted<-RESULTSTEMP[order(RESULTSTEMP$resultgcv,decreasing = TRUE),]
      
      path<- paste0(pathWrite,"earthVarimpImportance/")
      dir.create(path,recursive = TRUE,showWarnings = FALSE)
      
      write.table(RESULTSorted, file = paste0(path,disease,dataset,j,"KeywordsEmbeddingsUsefulSubset.txt"),append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE, sep = "\t")
      
      
    }
  }
}




