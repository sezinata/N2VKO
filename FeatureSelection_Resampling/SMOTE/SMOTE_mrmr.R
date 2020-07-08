library(RANKS) 
library("PRROC")
library("AUC")
library("e1071")
library(caret)
library(ggplot2)
library(ggrepel)
library(doParallel) 
library('stringr')
library(clusterSim)
library('DMwR')

networks<-c("PPI","PPIK")

datasets <- c("IntAct", "NCBI", "STRING")

diseases <- c( "prostateCancer","breastCancer","diabetesMellitus","obesity","lungCancer","alzheimer", "colorectalCancer")
#p q max values 
p_q_max<-matrix(c(1,4,4,2, 0.25,1,0.5,4, 0.25,1,4,2, 0.25,0.5,1,2, 2,4,1,0.25, 2,2,4,4, 4,2,4,4
),nrow=7,ncol=4,byrow = TRUE)
memory.limit(size=50000)
sized<-length(datasets)
features<-array()
featuresFolders<-array()
ifeat<-1
for(to in seq(10, 500, by = 10) ){
  features[ifeat]<- paste0("mRMRTOP",to)
  featuresFolders[ifeat]<-paste0("mRMRTOP",to)
  ifeat<-ifeat+1
}
OversamplingMethods<-c("ROSE","SMOTE")
OversamplingMethod<-OversamplingMethods[2]
parameters<-as.numeric(c('200','300','400','500','600','700','800','900','1000'))
parametersFolder<-c('200','300','400','500','600','700','800','900','1000')
sized<-length(datasets)
k=5
pv=1
q=1
net=1
feat=1
dis=1
z=1
j=1
net=1
empty=0
para=1


path1<-"D:\\temp\\"
pathRead<-"D:\\temp\\FS\\"
pathWrite<-"D:\\temp\\OversamplingRuns\\"
for(dis in 1:7){
  disease=diseases[dis]
  p_q_max_disease<-p_q_max[dis,]
  labelfile<-paste0(path1,"2018OMIMDiseaseLabels/")

  for(feat in 1:length(features)){
    feature<-features[feat]
    featuresFolder<-featuresFolders[feat]
    gc()
    
    for(z in 1:2){
      dataset=datasets[z]
      print(paste0("Currently handling: ",disease," ",dataset," ",featuresFolder))
      
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

      dataset.name <- paste0(path1,"node2vec-master\\emb\\TunePQ\\node2vec5_5_10\\", dataset,"-",pvalue,"-",qvalue,".emb")

        #ppi
      data.name = read.table(dataset.name,header =FALSE,sep=" ",skip=1) #for original embedding results skips first line
      
      #.........................................
      
      result<-array()  
      for(para in 1:9){
        parameter<-parameters[para]
        Oversampling<-paste0(OversamplingMethod,parametersFolder[para])
        sum=0
        for (j in 1:k){
          KeywordsEmbeddingsconcatfs <-paste0(pathRead,feature,"/",disease,dataset,j,'KeywordsEmbeddingsUsefulSubset.txt')
          selectedKeywords = read.table(KeywordsEmbeddingsconcatfs,header =FALSE,sep="\t")
          selectedkeywordsnoNAs<-str_replace_all(selectedKeywords[,1],"`","")#just cleaning
          var_names <- as.vector(selectedkeywordsnoNAs)
          selectedkeywordsnoNAs<-as.array(na.omit(var_names))
          
          if(length(selectedkeywordsnoNAs)>0){
            ####CONCAT for whole matrix#############
            Kemb <- as.matrix(data.name)
            Kemb <-Kemb[order(Kemb[,1]),]  #sort is AUTOMATICALLY DONE HERE, 
            keywordsEmbeddingRows<-keyw[Kemb[,1]+1,1:ncol(keyw)]#+1 EMBEDDINGS START WITH 0 INDEX,VALUES WITH SAME INDICES ARE CHOSEN  
            K<-cbind(Kemb,keywordsEmbeddingRows)
            filteresKeywords=as.data.frame( K[,selectedkeywordsnoNAs])  #each keyword individually a list element with its values for all proteins
            colnames(filteresKeywords)<-selectedkeywordsnoNAs
            n <- nrow(K)
            firstcolumn<-K[,1]  # node2vec ids/indices starts with 0
            filteresKeywords<-cbind(firstcolumn,filteresKeywords)         
            testIndices <- paste0(path1,"TestIndices/", dataset,disease, j,"indicesofTestSet.txt")  

            testIndices= read.table(testIndices,header =FALSE,sep="\n")
            testIndices<-as.array(testIndices$V1-1)
            
            xtrainWithindices<-filteresKeywords[ !(firstcolumn%in% testIndices),]  # use !(%in)
            xtrain<-xtrainWithindices[,-1]#removes first column
            
            xtestWithindices<-filteresKeywords[ firstcolumn%in% testIndices,]
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
            xtrain<-as.matrix(xtrain)
            xtest<-as.matrix(xtest)
            trainIndices<-xtrainWithindices[,1]
            trainIndices_vector<-as.vector(trainIndices)
            trainMean <- apply(xtrain,2,mean) #2 for columns, 1 for rows
            trainSd <- apply(xtrain,2,sd)
            xtrain<- sweep(sweep(xtrain, 2L, trainMean), 2, trainSd, "/") # using the default "-" to subtract mean column-wise   
            # ## centered AND scaled
            xtest<-sweep(sweep(xtest, 2L, trainMean), 2, trainSd, "/")

            trainsize<-nrow(xtrain)
            merged<-as.matrix(rbind(xtrain,xtest))  # added as.matrix for just one column results, it turns it to vector
            merged<-as.matrix(merged[,!(colSums(is.na(merged)))]) # added as.matrix for just one column results, it turns it to vector
            xtrain<-as.matrix(merged[1:trainsize,])
            xtest<-merged[(trainsize+1):nrow(merged),]
            
            if (nrow(xtrain) != nrow(ytrain)) 
              stop("train data and label matrices do not agree")
            
            colnames(ytrain) <- "Labels"
            path<- paste0(pathWrite,feature,"/LogisticRegression/",Oversampling,"/")
            dir.create(path,recursive = TRUE,showWarnings = FALSE)
            if(ncol(xtrain)==1){
              colnames(xtrain)<-'xtrain'
              DATA=cbind(xtrain, ytrain)
              DATA<-as.data.frame(DATA)
              names(DATA)
              DATA$Labels = as.factor(DATA$Labels)
              data.smote<-SMOTE(Labels~xtrain, data = DATA, perc.over = parameter, k = 5,  learner = NULL) #K NEAREST NEIGHBOR
              data.smote$Labels = as.matrix(as.numeric(data.smote$Labels))
              
              model <-lm(Labels~xtrain, data =  data.smote)
              new <- data.frame(xtrain = xtest) 
              glmpredict<- predict(model, newdata = new)
              ntrain <- data.frame(xtrain = xtrain) 
              colnames(ntrain)<-'xtrain'
              glmpredictrain<-predict(model, newdata=ntrain)
              label_test <- factor(testlabel)
              sum<-sum+auc(roc(glmpredict,label_test))#BestOversampling 
              
            }
            else{
              DATA=cbind(xtrain, ytrain)
              DATA<-as.data.frame(DATA)
              names(DATA)
              DATA$Labels = as.factor(DATA$Labels)
              data.smote<-SMOTE(Labels~., data = DATA, perc.over = parameter, k = 5,  learner = NULL) #K NEAREST NEIGHBOR
              data.smote$Labels = as.matrix(as.numeric(data.smote$Labels))
              
              xtest<-as.data.frame(xtest)
              xtrain<-as.data.frame(xtrain)
              
              model <- glm(Labels~., data = data.smote, family = "gaussian")
              glmpredict<- predict(model, newdata = xtest, type = "response")
              glmpredictrain<-predict(model, newdata=xtrain,  type="response")
              label_test <- factor(testlabel)
              sum<-sum+auc(roc(glmpredict,label_test))#BestOversampling  
              
            }
            write.table(glmpredict, file = paste0(path,disease,dataset, j,"Probabilities.txt"),append = FALSE, row.names = FALSE, sep = " ")
            
            write.table(glmpredictrain, file = paste0(path,disease,dataset, j,"trainProbabilities.txt"),append = FALSE, row.names = FALSE, sep = " ")
            write(testIndices, file = paste0(path,disease,dataset, j,"indicesofTestSet.txt"),append = FALSE, sep = "\n")
            write(testlabel, file = paste0(path,disease,dataset, j,"Labels.txt"),append = FALSE, sep = "\n")
            write.table(ytrain, file = paste0(path,disease,dataset, j,"TrainLabels.txt"),append = FALSE,quote=FALSE, row.names = FALSE, sep = " ")
          }
          else{
            empty=empty+1
            empties<-as.data.frame(rbind(empties,cbind(disease,dataset,j,feature)))
            print(empties)
            
          }
        }  
        result[para]<-sum/k
        
        DATA <- matrix(c(dataset,result[para],disease,parameter,feature),ncol=5,byrow=TRUE)
        
        write.table(DATA, file = paste0(pathWrite,"MRMRSMOTEParameterScores.txt"),quote=FALSE,append = TRUE,col.names = FALSE, row.names = FALSE, sep = " ")
        
      }
    }
  }
}