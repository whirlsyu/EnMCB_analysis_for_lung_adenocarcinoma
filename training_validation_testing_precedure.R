#analyze data using training, validation, and testing procedure

#include the private functions
source('functions.R')

#For GSE39279 in GEO database, one can download the data using GEOquery package in bioconductor:
library(GEOquery)
GSE39279 <- getGEO('GSE39279',GSEMatrix=TRUE)[[1]]
#the function will return a expressionset in R
#You can get the clinical data using:
pData(GSE39279)

#For the lung adenocarcinoma data in TCGA database,
#One can get the lung adenocarcinoma methylation 450k data using RTCGA package
library(RTCGA)
downloadTCGA(cancerTypes = 'LUAD', 'LUAD-FFPE.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3', destDir ='your_dir', date = "2016-01-28")
#Alternatively, you can also use our ftp mirror server for RTCGA files at:
#ftp://server.genelibs.com/raw_data_tcga/
#This ftp server may fail to load with an IP outside mainland China due to some regulation of our network operators, you may need a proxy.

#For clinical data in TCGA, you can download it using:
downloadTCGA(cancerTypes = 'LUAD', 'LUAD.Merge_Clinical.Level_1', destDir ='your_dir', date = "2016-01-28")
#We also recommend you read this reference:
#Liu J, Lichtenberg T, et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell 173:400-416 e411. 10.1016/j.cell.2018.02.052

#You can get the subsets (training, validation, and testing sets) from the original datasets by barcode and GSM code saved in the file of tvt_names.txt.

#Also you may want to get the data automatically using the following built in function:
get_the_datasets('training/validation/testing')
#after downloading, the data will be loaded automatically, which contained the training, validation and testing sets.
#Those data were obtained from the TCGA and GEO data sets.
#data_combat_train : training set, obtained from the TCGA (LUAD) databases.
#data_combat_validation : validation set  51 samples from the GEO (GSE39279) databases.
#data_combat_test : testing set 104 samples from the GEO (GSE39279) databases.
#y_surv_train : Surv object for the training set, DFS
#y_surv_validation : Surv object for the validation set, DFS
#y_surv_test : Surv object for the testing set, DFS

#If you want do the further resampling, you need to combine 
#the survival data as well as the methylation matrix by following (for example):
#y_surv_for_resampling <- c(y_surv_validation,y_surv_test)
#methylation_data_for_resampling <- cbind(data_combat_validation, data_combat_test)
#And then sample them by following:
#resampling <- sample(colnames(methylation_data_for_resampling),floor(ncol(methylation_data_for_resampling)*0.33))
#y_surv_resampled_validation_set <- y_surv_for_resampling[resampling]
#methylation_data_resampled_validation_set <- methylation_data_for_resampling[,resampling]


total_res<-EnMCB::IdentifyMCB(data_combat_train)
# It should be printed like this:
# Statistics ( 31728  MCBs in total):
#  chr1 : total MCBs: 3080  Mean Length: 195.6571  (Range:  2 2348 )
# chr10 : total MCBs: 1713  Mean Length: 179.2691  (Range:  2 2075 )
# chr11 : total MCBs: 2000  Mean Length: 176.279  (Range:  2 2135 )
# chr12 : total MCBs: 1479  Mean Length: 165.0892  (Range:  2 2056 )
# chr13 : total MCBs: 806  Mean Length: 189.9082  (Range:  2 2528 )
# chr14 : total MCBs: 888  Mean Length: 206.9167  (Range:  2 2333 )
# chr15 : total MCBs: 949  Mean Length: 168.9494  (Range:  2 2606 )
# chr16 : total MCBs: 1359  Mean Length: 192.9801  (Range:  2 2143 )
# chr17 : total MCBs: 1964  Mean Length: 173.9308  (Range:  2 2499 )
# chr18 : total MCBs: 450  Mean Length: 205.9578  (Range:  2 1284 )
# chr19 : total MCBs: 1652  Mean Length: 170.1713  (Range:  2 2303 )
# chr2 : total MCBs: 2153  Mean Length: 189.8834  (Range:  2 2730 )
# chr20 : total MCBs: 743  Mean Length: 177.8964  (Range:  2 2036 )
# chr21 : total MCBs: 284  Mean Length: 197.3415  (Range:  2 1678 )
# chr22 : total MCBs: 480  Mean Length: 195.3646  (Range:  2 1543 )
# chr3 : total MCBs: 1316  Mean Length: 185.0479  (Range:  2 2384 )
# chr4 : total MCBs: 1194  Mean Length: 180.1985  (Range:  2 1625 )
# chr5 : total MCBs: 1536  Mean Length: 196.7161  (Range:  2 2341 )
# chr6 : total MCBs: 2594  Mean Length: 161.3204  (Range:  2 3176 )
# chr7 : total MCBs: 2219  Mean Length: 207.0717  (Range:  2 2924 )
# chr8 : total MCBs: 1435  Mean Length: 189.2767  (Range:  2 1696 )
# chr9 : total MCBs: 386  Mean Length: 204.1554  (Range:  2 1566 )
# chrX : total MCBs: 1046  Mean Length: 436.7247  (Range:  2 4282 )
# chrY : total MCBs: 2  Mean Length: 158  (Range:  29 287 )



#if you want do global comparison, you just use total_res$MCBinformation as total_res_select_filtered
#and skip this selection section, however, skip this selection will make the calculation very time consuming (4-15 days).
#selection procedure
#MCBs which have more than 5 CpG sites were retained.

total_res_select<-total_res$MCBinformation[as.numeric(total_res$MCBinformation[,'CpGs_num']) >=5,]

#feature selection based on the mean value of CpG sites in MCB using L1 penalty 
res_mean_all<-apply(total_res_select,1, function(x,data_set=data_combat_train){
  CpGs<-strsplit(x['CpGs']," ")[[1]]
  colMeans(data_set[CpGs,])
})
rz=!(is.na(y_surv_train)|y_surv_train[,1]==0)
library(glmnet)
set.seed(6)
cvg<-cv.glmnet(res_mean_all[rz,],y_surv_train[rz],family='cox',type.measure = 'deviance')
plot(cvg)
total_res_select_filtered<- total_res_select[which(coef(cvg, s=0.01)!=0),]
#end the selection procedure

#get the AUC results for models in training, validation and testing sets
res_cv_3<-NULL
for (i in seq(nrow(total_res_select_filtered))) {
  cat(i," ")
  single_in_total_res_select_filtered<-total_res_select_filtered[i,]
  #Several models can't be built due to the errors in survivalsvm package, this can be improved further.
  try(res_cv_3<-rbind(res_cv_3,auc_cal_cv_tvt(single_in_total_res_select_filtered,
                                              data_combat_train,
                                              data_combat_validation,
                                              data_combat_test,
                                              y_surv_train,y_surv_validation,y_surv_test)),silent=TRUE)
}
colnames(res_cv_3)<-c(colnames(total_res_select_filtered),
                      'COX_train_cv','SVR_train_cv','eNet_train_cv','em_train_cv',
                      'COX_train','SVR_train','eNet_train','em_train',
                      'COX_validation','SVR_validation','eNet_validation','em_validation',
                      'COX_test','SVR_test','eNet_test','em_test')


#calculation of the Rank Product value

RP_data_cal<-res_cv_3[,9:16]
RP_data_cal<-matrix(data=apply(RP_data_cal, 2, function(x){(length(x)-rank(x,ties.method='first')+1)}),nrow = nrow(RP_data_cal),ncol = ncol(RP_data_cal))
RP_data_cal<-data.frame(RP_data_cal,RP=apply(RP_data_cal, 1, function(x)sum(log(x))))
rownames(RP_data_cal)<-res_cv_3[,'MCB_no']
res_cv_3_RP<-data.frame(res_cv_3,RP=RP_data_cal$RP)
write.csv(res_cv_3_RP,file = 'res_cv_3_RP.csv')

# select the best mcb
# results may showed that the best one is MCB 29147
# The MCB_no may different in the future due to the annoation file for 450K may updated in minfi package, however, general results should be the same.
# MCB 29147 is identical to that of MCB 29016 in our paper,
# 29147 and 29016 are the same MCB with the CpG sites of "cg01957585 cg11323985 cg26118821 cg12082271 cg22276811"
# we used parallel computation using clusters in our lab, this give different MCB_no.
MCBblocks_selected='29147'

single_res<-t(as.matrix(total_res_select_filtered[total_res_select_filtered[,'MCB_no'] %in% MCBblocks_selected,]))

#build the model based on the five CpGs (cg01957585 cg11323985 cg26118821 cg12082271 cg22276811)
em<-ensemble_model(single_res = single_res,training_set = data_combat_train,Surv_training = y_surv_train)

#use the model to predict the responses in training set
ep<-ensemble_prediction(ensemble_model = em,predition_data = data_combat_train,multiple_results = T)

#calculation of AUC, plot and save the results in pdf files
ROC_multiple_clinical(test_frame = data.frame(cox=ep['cox',],
                                             svm=ep['svm',],
                                             eNet=ep['enet',],
                                             ensemble=ep['ensemble',]
),
y = y_surv_train,file_name = "multiple model training"
)

#use the model to predict the response in validation set
ep_v<-ensemble_prediction.m(ensemble_model = em,predition_data = data_combat_validation)

#calculation of AUC, plot and save the results in pdf files
ROC_multiple_clinical(test_frame = data.frame(cox=ep_v['cox',],
                                             svm=ep_v['svm',],
                                             eNet=ep_v['enet',],
                                             ensemble=ep_v['ensemble',]
),
y = y_surv_validation,file_name = "multiple model validation"
)

#use the model to predict the response in testing set
ep_t<-ensemble_prediction(ensemble_model = em,predition_data = data_combat_test, multiple_results = T)

#calculation of AUC, plot and save the results in pdf files
ROC_multiple_clinical(test_frame = data.frame(cox=ep_t['cox',],
                                             svm=ep_t['svm',],
                                             eNet=ep_t['enet',],
                                             ensemble=ep_t['ensemble',]
),
y = y_surv_test,file_name = "multiple model test"
)

#draw survival curve for training set
EnMCB::draw_survival_curve(ep['ensemble',],
                    living_days = y_surv_train[,1],
                    living_events = y_surv_train[,2],
                    write_name = "ensemble model training")

#draw survival curve for testing set
EnMCB::draw_survival_curve(ep_t['ensemble',],
                    living_days = y_surv_test[,1],
                    living_events = y_surv_test[,2],
                    write_name = "ensemble model testing",threshold = median(ep['ensemble',]))
