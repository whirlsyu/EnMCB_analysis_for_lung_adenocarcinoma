# tested environment:
# those functions are compatible with 3.6.0 or higher

require('survivalROC')
require('prognosticROC')

get_the_datasets<-function(type){
  if (!file.exists('annotation_met.txt')) download.file('ftp://server.genelibs.com/en_mcb/annotation_met.txt',destfile = 'annotation_met.txt')      #download the annotation data
  if (type == 'training/testing'){
    if (!file.exists('training_testing_data.RData')){
      #Considering the data size up to 1.7Gb, the data downloading could be a time consuming process.
      #download the data sets
      download.file('ftp://server.genelibs.com/en_mcb/training_testing_precedure/training_testing_data.RData',
                    destfile = 'training_testing_data.RData' )
    }
    load('training_testing_data.RData',envir = globalenv())
  }else if (type == 'training/validation/testing'){
    if (!file.exists('training_validation_testing_data.RData')){
      #Considering the data size up to 1.7Gb, the data downloading could be a time consuming process.
      #download the data sets
      download.file('ftp://server.genelibs.com/en_mcb/training_testing_precedure/training_validation_testing_data.RData',
                    destfile = 'training_validation_testing_data.RData' )
    }
    load('training_validation_testing_data.RData',envir = globalenv())
  }else{
    cat("unable to recognize the type code, please check it.")
  }
}

auc_cal_cv<- function(single_one,train_set,y_surv_train,nfold=10,seed=NULL) {
    library(survival)
    
    #remove all the NA and zero survival data
    na_or_zero_data<-(is.na(y_surv_train[,1])|y_surv_train[,1]==0)
    train_set = train_set[,!na_or_zero_data]
    y_surv_train = y_surv_train[!na_or_zero_data]
    
    cox<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'cox',silent = T,seed = seed)
    svr<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'svm',silent = T,seed = seed)
    enet<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'enet',silent = T,seed = seed)
    em<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'ensemble',silent = T,seed = seed)
    
    auc_COX_cv<-auc_roc(y_surv = y_surv_train,cox$MCB_matrix)
    auc_SVR_cv<-auc_roc(y_surv = y_surv_train,svr$MCB_matrix)
    auc_eNet_cv<-auc_roc(y_surv = y_surv_train,enet$MCB_matrix)
    auc_em_cv<-auc_roc(y_surv = y_surv_train,em$MCB_matrix)
  
  return(c(single_one,
           auc_COX_cv,auc_SVR_cv,auc_eNet_cv,auc_em_cv))
}

auc_cal_cv_tvt <- function(single_one,
                           train_set,validation_set,test_set,
                           y_surv_train,y_surv_validation,y_surv_test,
                           seed=9){
  require(survival)
  
  #remove all the NA and zero survival data
  na_or_zero_data<-(is.na(y_surv_train[,1])|y_surv_train[,1]==0)
  train_set = train_set[,!na_or_zero_data]
  y_surv_train = y_surv_train[!na_or_zero_data]
  
  na_or_zero_data<-(is.na(y_surv_validation[,1])|y_surv_validation[,1]==0)
  validation_set = validation_set[,!na_or_zero_data]
  y_surv_validation = y_surv_validation[!na_or_zero_data]
  
  na_or_zero_data<-(is.na(y_surv_test[,1])|y_surv_test[,1]==0)
  test_set = test_set[,!na_or_zero_data]
  y_surv_test = y_surv_test[!na_or_zero_data]
  
  
  cox<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'cox',silent = T,seed = seed)
  svr<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'svm',silent = T,seed = seed)
  enet<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'enet',silent = T,seed = seed)
  em<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'ensemble',silent = T,seed = seed)
  
  auc_COX_train_cv<-auc_roc(y_surv = y_surv_train,cox$MCB_matrix)
  auc_SVR_train_cv<-auc_roc(y_surv = y_surv_train,svr$MCB_matrix)
  auc_eNet_train_cv<-auc_roc(y_surv = y_surv_train,enet$MCB_matrix)
  auc_em_train_cv<-auc_roc(y_surv = y_surv_train,em$MCB_matrix)
  
  
  vt_data<-cbind(validation_set,test_set)
  vt_survival<-c(y_surv_validation,y_surv_test)
  vt_flag<-c(rep('validation',ncol(validation_set)),
             rep('test',ncol(test_set)))
  
  cox = metricMCB(t(single_one),
                      training_set = train_set,
                      Surv = y_surv_train,
                      testing_set = vt_data,
                      Surv.new = vt_survival,
                      Method = 'cox',silent = T)
  
  svr = metricMCB(t(single_one),
                      training_set = train_set,
                      Surv = y_surv_train,
                      testing_set = vt_data,
                      Surv.new = vt_survival,
                      Method = 'svm',silent = T) 
 
  enet = metricMCB(t(single_one),
                  training_set = train_set,
                  Surv = y_surv_train,
                  testing_set = vt_data,
                  Surv.new = vt_survival,
                  Method = 'enet',silent = T)  
  
  em = EnMCB::ensemble_model(t(single_one),
                   training_set = train_set,
                   Surv_training = y_surv_train,
                   testing_set = vt_data,
                   Surv_testing = vt_survival)
  prid_em<- EnMCB::ensemble_prediction(em,train_set)
  test_and_validation_em<- EnMCB::ensemble_prediction(em,vt_data)
  
  auc_COX_train<-auc_roc(y_surv = y_surv_train,cox$MCB_cox_matrix_training)
  auc_SVR_train<-auc_roc(y_surv = y_surv_train,svr$MCB_svm_matrix_training)
  auc_eNet_train<-auc_roc(y_surv = y_surv_train,enet$MCB_enet_matrix_training)
  auc_em_train<-auc_roc(y_surv = y_surv_train,prid_em)
  
  
  auc_COX_validation<-auc_roc(y_surv = y_surv_validation,cox$MCB_cox_matrix_test_set[vt_flag == 'validation'])
  auc_SVR_validation<-auc_roc(y_surv = y_surv_validation,svr$MCB_svm_matrix_test_set[vt_flag == 'validation'])
  auc_eNet_validation<-auc_roc(y_surv = y_surv_validation,enet$MCB_enet_matrix_test_set[vt_flag == 'validation'])
  auc_em_validation<-auc_roc(y_surv = y_surv_validation,test_and_validation_em[vt_flag == 'validation'])

  auc_COX_test<-auc_roc(y_surv = y_surv_test,cox$MCB_cox_matrix_test_set[vt_flag == 'test'])
  auc_SVR_test<-auc_roc(y_surv = y_surv_test,svr$MCB_svm_matrix_test_set[vt_flag == 'test'])
  auc_eNet_test<-auc_roc(y_surv = y_surv_test,enet$MCB_enet_matrix_test_set[vt_flag == 'test'])
  auc_em_test<-auc_roc(y_surv = y_surv_test,test_and_validation_em[vt_flag == 'test'])
  
  return(c(single_one,
           auc_COX_train_cv,auc_SVR_train_cv,auc_eNet_train_cv,auc_em_train_cv,
           auc_COX_train,auc_SVR_train,auc_eNet_train,auc_em_train,
           auc_COX_validation,auc_SVR_validation,auc_eNet_validation,auc_em_validation,
           auc_COX_test,auc_SVR_test,auc_eNet_test,auc_em_test))
}


auc_cal_cv_tt <- function(single_one,all_data_train,data_test,y_surv_train,y_surv_test,nfold=10,seed=NULL) {
  #remove all the NA and zero survival data
  na_or_zero_data<-(is.na(y_surv_train[,1])|y_surv_train[,1]==0)
  all_data_train = all_data_train[,!na_or_zero_data]
  y_surv_train = y_surv_train[!na_or_zero_data]
  
  na_or_zero_data<-(is.na(y_surv_test[,1])|y_surv_test[,1]==0)
  data_test = data_test[,!na_or_zero_data]
  y_surv_test = y_surv_test[!na_or_zero_data]
  
  set.seed(seed)
  sp<-sample(1:length(y_surv_train),replace = F)
  order_sp<-order(sp)
  y_surv_train = y_surv_train[sp]
  all_data_train = all_data_train[,sp]
  
  train<-1:length(y_surv_train)
  folds <- cut(seq(1,length(train)),breaks=nfold,labels=FALSE)
  survival_predictions<-NULL
  for (i in seq(unique(folds))) {
    rz<- which(folds==i,arr.ind=TRUE)
    em<-ensemble_model(t(single_one),all_data_train[,train[-rz]],y_surv_train[-rz])
    survival_prediction<-ensemble_prediction(em,all_data_train[,train[rz]],mutiple_results = T)
    survival_predictions<-cbind(survival_predictions,survival_prediction)
  }
  
  em_all<-ensemble_model(t(single_one),all_data_train,y_surv_train)
  pre_test<-ensemble_prediction(em_all,data_test,mutiple_results = T)
  
  auc_COX_test<-auc_roc(y_surv = y_surv_test,pre_test['cox',])
  auc_SVR_test<-auc_roc(y_surv = y_surv_test,pre_test['svm',])
  auc_eNet_test<-auc_roc(y_surv = y_surv_test,pre_test['enet',])
  auc_em_test<-auc_roc(y_surv = y_surv_test,pre_test['ensemble',])
  
  auc_COX_cv<-auc_roc(y_surv = y_surv_train,survival_predictions['cox',])
  auc_SVR_cv<-auc_roc(y_surv = y_surv_train,survival_predictions['svm',])
  auc_eNet_cv<-auc_roc(y_surv = y_surv_train,survival_predictions['enet',])
  auc_em_cv<-auc_roc(y_surv = y_surv_train,survival_predictions['ensemble',])
  
  return(c(single_one,
           auc_COX_cv,auc_SVR_cv,auc_eNet_cv,auc_em_cv,
           auc_COX_test,auc_SVR_test,auc_eNet_test,auc_em_test))
}

auc_roc <- function(y_surv,predictions) {
  survivalROC::survivalROC.C(Stime = y_surv[,1],
                           status = y_surv[,2],
                           marker = predictions,
                           predict.time = 5,
                           span =0.25*length(y_surv[,1])^(-0.20))$AUC
}

#'@title ROC calculation for multiple clinical data
#'@author Xin Yu
#'@description Demo matrix for methylation matrix.
#'@param test_frame clinical data as data.frame
#'@param y_surv Surv object created by survival package
#'@param file_name the file name used for saving the results
#'@param ntime time point for time dependent ROC (year).
#'@return This function will only used for AUC calculation.
#'@examples 

ROC_multiple_clinical<-function(test_frame,y_surv,file_name="title",ntime=5){
  require("prognosticROC") #version 0.7
  require("ggplot2")
  require("plotROC")
  sroclong_all<-NULL
  for (n in 1:ncol(test_frame)) {
    ROC_res= survivalROC::survivalROC(Stime=y_surv[,1],
                                      status=y_surv[,2],
                                      marker =as.numeric(test_frame[,n]),
                                      lambda = NULL,
                                      predict.time = ntime,method = "NNE",span =0.25*length(y_surv)^(-0.20))#
    sroclong_all<-ROCdata_save(sroclong_all,ROC_res,mark = paste(ntime,"year AUC at",colnames(test_frame)[n],round(ROC_res$AUC,2),collapse = " "))
  }
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # The palette with black:
  # cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # To use for fills, add
  #scale_fill_manual(values=cbPalette)
  # To use for line and point colors, add
  gg<-ggplot2::ggplot(sroclong_all, aes(x = FPF, y = TPF, label = c, color = group))+
    coord_cartesian(ylim=c(0, 1.05),xlim=c(0, 1.05),expand=FALSE)+
    geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_light())+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    theme(legend.position=c(0.8,0.3), legend.justification=c(1,1),legend.title = element_blank())+
    geom_abline(slope=1, colour="black")+scale_colour_manual(values=cbPalette)
  ggplot2::ggsave(filename = paste("Time ROC of ",file_name,".jpeg",sep=""),plot = gg,device ="jpeg" ,
         path = getwd(),dpi = 300,units = "in",width = 5, height = 4.5,
         limitsize=F)
}

ROCdata_save<-function(origin=NULL,perf,mark="none"){
  sroclong<-data.frame(TPF = perf$TP, FPF =perf$FP, 
                       c = perf$cut.values, 
                       group = rep(mark, length(perf$TP)) )
  sroclong<-rbind(origin,sroclong)
  sroclong
}
