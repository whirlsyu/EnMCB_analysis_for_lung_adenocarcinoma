#analyze data using training and testing procedure

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

#For clinical data in TCGA, you can download it using:
downloadTCGA(cancerTypes = 'LUAD', 'LUAD.Merge_Clinical.Level_1', destDir ='your_dir', date = "2016-01-28")
#We also recommend you read this reference:
#Liu J, Lichtenberg T, et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell 173:400-416 e411. 10.1016/j.cell.2018.02.052

#You can get the subsets (training and testing sets) from the original datasets by barcode and GSM code saved in the file of tt_names.txt. 

#Also you may want to get the data automatically using the following built in function:
get_the_datasets('training/testing')
#after downloading, the data will be loaded automatically, which contained the training and testing sets.
#Those data were obtained from the TCGA and GEO data sets.
#data_combat_train : training set sampled from the TCGA (LUAD) and GEO (GSE39279) databases.
#data_combat_test : testing set sampled from the TCGA (LUAD) and GEO (GSE39279) databases.
#y_surv_train : Surv object for the training set, DFS
#y_surv_test : Surv object for the testing set, DFS
#Note that the data sets were pre-sampled from the database for reproducing the results. 
#One may use sample() function in R if you want do the sampling.


mcb_new_res<-EnMCB::IdentifyMCB(data_combat_train)
# Statistics ( 31251  MCBs in total):
# chr1 : total MCBs: 3059  Mean Length: 190.5567  (Range:  2 2230 )
# chr10 : total MCBs: 1684  Mean Length: 168.633  (Range:  2 2075 )
# chr11 : total MCBs: 1955  Mean Length: 171.268  (Range:  2 1822 )
# chr12 : total MCBs: 1489  Mean Length: 159.7414  (Range:  2 2056 )
# chr13 : total MCBs: 791  Mean Length: 185.7573  (Range:  2 2528 )
# chr14 : total MCBs: 867  Mean Length: 205.3645  (Range:  2 2333 )
# chr15 : total MCBs: 913  Mean Length: 162.7536  (Range:  2 2905 )
# chr16 : total MCBs: 1325  Mean Length: 194.6196  (Range:  2 2143 )
# chr17 : total MCBs: 1937  Mean Length: 164.0754  (Range:  2 2499 )
# chr18 : total MCBs: 433  Mean Length: 206.5866  (Range:  2 1284 )
# chr19 : total MCBs: 1640  Mean Length: 168.8378  (Range:  2 2303 )
# chr2 : total MCBs: 2103  Mean Length: 184.6562  (Range:  2 3253 )
# chr20 : total MCBs: 728  Mean Length: 178.9354  (Range:  2 2036 )
# chr21 : total MCBs: 280  Mean Length: 182.0786  (Range:  2 1678 )
# chr22 : total MCBs: 475  Mean Length: 182.4042  (Range:  2 1543 )
# chr3 : total MCBs: 1307  Mean Length: 176.9725  (Range:  2 1973 )
# chr4 : total MCBs: 1207  Mean Length: 169.193  (Range:  2 1625 )
# chr5 : total MCBs: 1521  Mean Length: 191.0039  (Range:  2 2341 )
# chr6 : total MCBs: 2526  Mean Length: 156.2637  (Range:  2 2596 )
# chr7 : total MCBs: 2193  Mean Length: 203.7852  (Range:  2 2924 )
# chr8 : total MCBs: 1422  Mean Length: 188.5682  (Range:  2 1876 )
# chr9 : total MCBs: 363  Mean Length: 205.2645  (Range:  2 1566 )
# chrX : total MCBs: 1031  Mean Length: 448.4656  (Range:  2 4282 )
# chrY : total MCBs: 2  Mean Length: 158  (Range:  29 287 )


#MCBs which have more than 5 CpG sites were retained.
mcb_new_res_select<-mcb_new_res$MCBinformation[as.numeric(mcb_new_res$MCBinformation[,'CpGs_num']) >=5,]

#feature selection based on the mean value of CpG sites in MCB using L1 penalty 
mcb_new_res_select_mean_all<-apply(mcb_new_res_select,1, function(x,data_set=data_combat_train){
  CpGs<-strsplit(x['CpGs']," ")[[1]]
  colMeans(data_set[CpGs,])
})
rz=!(is.na(y_surv_train)|y_surv_train[,1]==0)

library(glmnet)
set.seed(6)
cvg<-glmnet::cv.glmnet(mcb_new_res_select_mean_all[rz,],y_surv_train[rz],family='cox' ,type.measure = 'deviance')
plot(cvg)
mcb_new_res_select_filtered<- mcb_new_res_select[which(coef(cvg, s=0.01)!=0),]

new_res_cv<-NULL
for (i in seq(nrow(mcb_new_res_select_filtered))) {
  single_in_mcb_new_res_select_filtered<-mcb_new_res_select_filtered[i,]
  try(new_res_cv<-rbind(new_res_cv,auc_cal_cv_tt(single_in_mcb_new_res_select_filtered,
                                                 data_combat_train,
                                                 data_combat_test,
                                                 y_surv_train,
                                                 y_surv_test,
                                                 nfold = 10,seed = 6)))
}
colnames(new_res_cv)<-c(colnames(mcb_new_res_select_filtered),'COX_cv','SVR_cv','eNet_cv','em_cv',
                        'COX_test','SVR_test','eNet_test','em_test')

#calculation of the Rank Product value
RP_data_cal<-new_res_cv[,9:11]
RP_data_cal<-matrix(data=apply(RP_data_cal, 2, function(x){(length(x)-rank(x,ties.method='first')+1)}),nrow = nrow(RP_data_cal),ncol = ncol(RP_data_cal))
RP_data_cal<-data.frame(RP_data_cal,RP=apply(RP_data_cal, 1, function(x)sum(log(x))))
rownames(RP_data_cal)<-new_res_cv[,'MCB_no']
res_cv_RP<-data.frame(mcb_new_res_select_filtered[mcb_new_res_select_filtered[,'MCB_no'] %in% rownames(RP_data_cal),],new_res_cv[,9:16],RP=RP_data_cal$RP)
#save the results
write.csv(res_cv_RP,file = 'res_cv_RP.csv')
