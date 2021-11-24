# No evidence for association between polygenic risk of Multiple Sclerosis and MRI phenotypes in ~30,000 healthy UK Biobank participants

A quick replication effort to see if the results reported [here](https://journals.sagepub.com/doi/full/10.1177/13524585211034826) can be reproduced in the UKB imaging cohort.

Steps:

0. Generate polygenic risk scores
1. Divide the UKB cohort into training (not included in MRI) and testing (MRI data available) sets
2. Determine the optimal MS PRS in the training cohort
3. Check that the optimal PRS performs ok in the MRI cohort
4. Correlate PRS with MRI imaging-derived phenotypes in the MRI cohort


## 0. Generate polygenic risk scores
We generated polygenic risk scores in UKB using the IMSGC meta-analysis discovery stage summary statistics. We used the clumping-and-thresholding method to create a variety of scores. For a full description of our methods please see the [code](https://github.com/benjacobs123456/MS_UKB_PRS) and read the [manuscript](https://nn.neurology.org/content/8/4/e1007).

## 1. Divide the UKB cohort into two sets - the non-imaging cohort (to select the best PRS, i.e. training set), and the imaging cohort (validation or testing set)
````R
#############################################
#               Load packages
#############################################
library(gridExtra)
library(Hmisc)
library(reshape2)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(rcompanion)
library(ROCR)
library(RNOmni)

setwd("/data/Wolfson-UKBB-Dobson/UKB_PRS_MRI")

# read in data
source_of_report_data = read_tsv("../ukb_pheno_17_03_21/ukb_pheno_nocognitive_17032021.tsv.gz",col_types=cols_only(
  `Source of report of G35 (multiple sclerosis).0.0` = col_character(),
  EID = col_double()
  ))

# define MS status
source_of_report_data = source_of_report_data %>% mutate(MS_status = ifelse(!is.na(`Source of report of G35 (multiple sclerosis).0.0`),1,0))

# read in remainder of phenotype data
ukb_pheno = read_tsv("../ukb_pheno_911/ukb_pheno_final_MS_1301")

# get drb1*15 dose
drb15_doses = list()
for (i in 1:nrow(ukb_pheno)){
  drb = strsplit(ukb_pheno$`HLA imputation values.0.0`[i][1],split=",")[[1]][282]
  drb15_doses[[i]] = c(ukb_pheno$EID[i],drb)
}
drb15_doses_df = data.frame(do.call("rbind",drb15_doses))
colnames(drb15_doses_df)=c("EID","DRB1_15")

# impose posterior threshold
drb15_doses_df = drb15_doses_df %>%
mutate(DRB1_15 = as.numeric(as.character(DRB1_15))) %>%
mutate(DRB1_15 = ifelse(DRB1_15<0.7,0,DRB1_15)) %>%
mutate(DRB1_15 = ifelse(DRB1_15>=0.7 & DRB1_15 <1.4,1,DRB1_15)) %>%
mutate(DRB1_15 = ifelse(DRB1_15>=1.4,2,DRB1_15))

# calculate 'risk score' for DRB1*15 using OR 3.92 from Moutsianas 2015
drb_doses_df = drb15_doses_df %>% mutate(DRB_risk_score = DRB1_15*log(3.92))
drb_doses_df = drb_doses_df %>% mutate(EID=as.numeric(as.character(EID)))

# filter cols
ukb_pheno = ukb_pheno %>% select(EID,`Age at recruitment.0.0`,Sex.0.0,`Ethnic background.0.0`,contains("rincipal components"),contains("ethnic grouping"))

# combine
ukb_pheno = ukb_pheno %>% left_join(source_of_report_data,by="EID")
ukb_pheno = ukb_pheno %>% left_join(drb_doses_df,by="EID")

# read in mri data
mri_data = read_tsv("mri.tsv")

# exclude participants who have withdrawn
withdrawn = read_tsv("../helper_progs_and_key/excluded_indivs",col_names=FALSE)
ukb_pheno = ukb_pheno %>% filter(!EID %in% withdrawn$X1)

table(ukb_pheno$MS_status)
````

Total number of cases and controls:

Controls     |     Cases
------------ | --------------
485809       |   2405


````R

# Remove non-European participants
ukb_pheno = ukb_pheno %>% filter(`Genetic ethnic grouping.0.0` =="Caucasian")
table(ukb_pheno$MS_status)
````

After exclusion of non-European participants:

Controls     |     Cases
-----------  | --------------
407478       |   2102


````R
# Filter highly related individuals
# filter relatedness
# nb this is taking care to exclude the non-ms control from each pair to boost case numbers
kin = read_table2("/data/Wolfson-UKBB-Dobson/helper_progs_and_key/ukb43101_rel_s488282.dat")
highly_related = kin %>% filter(Kinship>0.0884) %>% filter(ID1 %in% ukb_pheno$EID) %>% filter(ID2 %in% ukb_pheno$EID)
highly_related = highly_related %>% left_join((ukb_pheno %>% filter(EID %in% highly_related$ID1) %>% select(EID,MS_status)  %>% rename(ID1 = EID,MS_status_ID1 = MS_status)),by="ID1") %>%
left_join((ukb_pheno %>% filter(EID %in% highly_related$ID2) %>% select(EID,MS_status) %>% rename(ID2 = EID, MS_status_ID2 = MS_status)),by="ID2")

exclusion = bind_rows(highly_related %>% filter(MS_status_ID1==1 & MS_status_ID2==0) %>% select(ID2) %>% rename(EID = ID2),
highly_related %>% filter(MS_status_ID1==0 & MS_status_ID2==1) %>% select(ID1) %>% rename(EID = ID1),
highly_related %>% filter(MS_status_ID1==0 & MS_status_ID2==0) %>% select(ID1) %>% rename(EID = ID1))
ukb_pheno = ukb_pheno %>% filter(!EID %in% exclusion$EID)
table(ukb_pheno$MS_status)
````   

After exclusion of non-European participants:

Controls     |     Cases
------------ | --------------
376587       |   2102

````R
# define non-MRI dataset
mri_data = mri_data %>% filter(`Brain MRI measurement completed.2.0`==1)
ukb_pheno_nomri = ukb_pheno %>% filter(!EID %in% mri_data$EID)
table(ukb_pheno_nomri$MS_status)
````

In the non-MRI dataset:

Controls     |     Cases
------------ | --------------
346547       |  1978

````R
ukb_pheno_mri = ukb_pheno %>% filter(EID %in% mri_data$EID)
table(ukb_pheno_mri$MS_status)
````

And in the MRI dataset:

Controls     |     Cases
------------ | --------------
30040        |  124

## 2. Select the best PRS in the UKB non-imaging cohort
Now we will test the PRS in the non-MRI cohort and select the score with the best explanatory power for MS susceptibility in this cohort.
````R
# collate PRS scores
pvals = c("1","0.8","0.6","0.4","0.2","0.1","0.05","0.005","0.0005","0.00005","0.000005","0.0000005","0.00000005")
r2_vector = c(0.1,0.2,0.4,0.6,0.8)

for(r2 in r2_vector){
  for(pval in pvals){
    scores = list()
    for(chr in 1:22){
      message("P value: ",pval)
      message("Clumping R2: ",r2)
      message("Chr: ",chr)
      file = paste0("/data/scratch/hmy117/chr",chr,"_pval",pval,"_r2_",r2,".sscore")
      if(file.exists(file)){
        df = read_table2(file)
        df = df %>% select(2,6)
        scores[[chr]] = df
      }
    }
    overall_score = do.call("rbind",scores)

    overall_score = overall_score %>% group_by(IID) %>% summarise(PRS = sum(SCORE1_SUM))
    colnames(overall_score)=c("EID","PRS")
    write_tsv(overall_score,paste0("/data/Wolfson-UKBB-Dobson/ms_prs/summarised_PRS_results_pval",pval,"r2",r2))
    message("Done")
  }
}


# Define scoring function
score_prs = function(pval,r2,hla="no",mri="no"){
  # read in prs
  filename = paste0("/data/Wolfson-UKBB-Dobson/ms_prs/summarised_PRS_results_pval",pval,"r2",r2)
  prs = read_table2(filename)
  colnames(prs) = c("EID","PRS")

  dataset = if(mri=="no"){
      ukb_pheno_nomri
    } else if(mri=="yes"){
      ukb_pheno_mri
    }

  prs_tuning_dataset = dataset %>% select(`Age at recruitment.0.0`,`Sex.0.0`,`Genetic principal components.0.1`,`Genetic principal components.0.2`,`Genetic principal components.0.3`,`Genetic principal components.0.4`,MS_status,EID,DRB_risk_score)
  prs_tuning_dataset$EID = as.numeric(as.character(prs_tuning_dataset$EID))
  prs_tuning_dataset = prs_tuning_dataset %>% left_join(prs,by="EID")
  prs_tuning_dataset = prs_tuning_dataset %>% filter(!(is.na(PRS)))
  prs_tuning_dataset = prs_tuning_dataset %>% mutate(combined_score = PRS+DRB_risk_score)

  prs_tuning_dataset$PRS=if(hla=="no"){
    message("Excluding HLA")
    rankNorm(prs_tuning_dataset$PRS)
  } else if(hla=="yes"){
    message("Including HLA")
    rankNorm(prs_tuning_dataset$combined_score)
  }

  null_model = glm(data=prs_tuning_dataset,
                     MS_status~`Age at recruitment.0.0`+
                       `Sex.0.0`+
                       `Genetic principal components.0.1`+
                       `Genetic principal components.0.2`+
                       `Genetic principal components.0.3`+
                       `Genetic principal components.0.4`,
                     family=binomial(link="logit"))

  prs_model = glm(data=prs_tuning_dataset,
                    MS_status~`Age at recruitment.0.0`+
                      `Sex.0.0`+
                      `Genetic principal components.0.1`+
                      `Genetic principal components.0.2`+
                      `Genetic principal components.0.3`+
                      `Genetic principal components.0.4`+
                      PRS,
                    family=binomial(link="logit"))
  nag = nagelkerke(prs_model,null_model)$Pseudo.R.squared.for.model.vs.null[3]
  message("Full output:")
  print(nagelkerke(prs_model,null_model))
  message("P value: ",pval)
  message("Clumping R2: ",r2)
  message("Nagelkerke's Pseudo-R2: ",nag)
  return(nag)
}

# Loop scoring function over different PRS
pvals = c("1","0.8","0.6","0.4","0.2","0.1","0.05","0.005","0.0005","0.00005","0.000005","0.0000005","0.00000005")
r2 = c(0.1,0.2,0.4,0.6,0.8)
params = expand.grid(pvals,r2)

nohla_nags = mapply(score_prs,pval=params[,1],r2=params[,2],hla="no",mri="no")
hla_nags = mapply(score_prs,pval=params[,1],r2=params[,2],hla="yes",mri="no")

params = data.frame(params)
params$Nagelkerke_Pseudo_R2 = hla_nags
colnames(params)=c("Pval","R2","Nagelkerke_Pseudo_R2")
hla_scores = params
params = data.frame(params)
params$Nagelkerke_Pseudo_R2 = nohla_nags
colnames(params)=c("Pval","R2","Nagelkerke_Pseudo_R2")
nohla_scores = params

nohla_scores$MHC = "MHC excluded"
hla_scores$MHC = "MHC included"
params = bind_rows(nohla_scores,hla_scores)
training_params = params

nagel_plot = ggplot(params,aes(factor(as.numeric(as.character(Pval))),Nagelkerke_Pseudo_R2,fill=factor(R2)))+
  facet_wrap(~MHC)+
  geom_col(position=position_dodge(),col="black")+
  theme_classic()+
  scale_fill_brewer(palette = "Set2")+
  labs(x="P value threshold",y=expression(Nagelkerke~Pseudo-R^{"2"}),fill=expression(Clumping~R^{"2"}~parameter))+
  theme(text=element_text(size=12))

png("nagel_plot.png",height=8,width=16,res=300,units="in")
nagel_plot
dev.off()

# save results
write_tsv(training_params,"training_results.tsv")
write_tsv(ukb_pheno_nomri,"ukb_pheno_nomri.tsv")
write_tsv(ukb_pheno_mri,"ukb_pheno_mri.tsv")
````
Here's a plot of the prediction pseudo-R2 values for each PRS in the training set (roughly equivalent to variance explained)
![Nagel plot](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/nagel_plot.png)

The best PRS explain ~ 1.3% (with MHC) and ~1.2% (without MHC) of MS liability in the training (non-MRI) cohort:

## 3. Validate the PRS in the MRI cohort
This step checks that these optimal PRS actually have some predictive/discriminative power to distinguish MS from controls in the MRI dataset. These results are a rough check as there are only 164 MS cases in the MRI cohort.


````R
# Even though there are only 124 MS cases in the MRI cohort, let's just check the PRS is working reasonably well.

# read in results
ukb_pheno_nomri = read_tsv("ukb_pheno_nomri.tsv")
ukb_pheno_mri = read_tsv("ukb_pheno_mri.tsv")

print("printing nagelkerke for best non-HLA PRS in testing set")
score_prs(pval="0.05",r2=0.8,hla="no",mri="yes")
print("printing nagelkerke for best HLA PRS in testing set")
score_prs(pval="0.0005",r2=0.8,hla="yes",mri="yes")


````
Clearly the PRS works well in discriminating MS cases from controls in the MRI dataset (or as well as can be expected for PRS):

Parameter  |  PRS (no HLA)   | PRS (with HLA)
--- | --- | ---
Nagelkerke's Pseudo-R2 | 0.0134 |  0.015
Likelihood ratio P value | 5.70e-06 | 1.43e-06

Now we'll check that with increasing PRS deciles, MS risk increases as we would expect.
````R
# read in best PRS

hla_prs = read_table2("/data/Wolfson-UKBB-Dobson/ms_prs/summarised_PRS_results_pval0.0005r20.8")
prs = read_table2("/data/Wolfson-UKBB-Dobson/ms_prs/summarised_PRS_results_pval0.05r20.8")

# HLA PRS
ukb_pheno_mri$EID = as.numeric(as.character(ukb_pheno_mri$EID))
ukb_pheno_mri = ukb_pheno_mri %>% left_join(hla_prs,by="EID")
ukb_pheno_mri = ukb_pheno_mri %>% filter(!(is.na(PRS)))
ukb_pheno_mri = ukb_pheno_mri %>% mutate(combined_score = PRS+DRB_risk_score)
ukb_pheno_mri = ukb_pheno_mri %>% mutate(HLA_PRS = rankNorm(ukb_pheno_mri$combined_score)) %>% select(-PRS)

# add in non-HLA
ukb_pheno_mri = ukb_pheno_mri %>% left_join(prs,by="EID")
ukb_pheno_mri = ukb_pheno_mri %>% filter(!(is.na(PRS)))
ukb_pheno_mri = ukb_pheno_mri %>% mutate(NOHLA_PRS = rankNorm(ukb_pheno_mri$PRS))


# decile plot
ukb_pheno_mri$prs_decile = cut2(ukb_pheno_mri$HLA_PRS,g=10)
prs_decile_model = glm(data=ukb_pheno_mri,
            MS_status~`Age at recruitment.0.0`+
              `Sex.0.0`+
              `Genetic principal components.0.1`+
              `Genetic principal components.0.2`+
              `Genetic principal components.0.3`+
              `Genetic principal components.0.4`+
              prs_decile,
            family=binomial(link="logit"))

tbl = data.frame(summary(prs_decile_model)$coefficients[-c(1:7),])
tbl$decile=c(2:10)
tbl$or=exp(tbl$Estimate)
tbl$lower_ci=exp(tbl$Estimate-1.96*tbl$Std..Error)
tbl$upper_ci=exp(tbl$Estimate+1.96*tbl$Std..Error)
tbl

prs_decile_plot=ggplot(tbl,aes(decile,or))+
geom_errorbar(aes(x=decile,ymin=lower_ci,ymax=upper_ci,width=0.2))+
geom_point(size=3,shape=22,fill="black")+
theme_classic()+
theme(legend.position="none",text=element_text(size=12))+
labs(x="PRS Decile",y="OR for MS (vs lowest decile)")+
annotate("text",label="MHC Included",x = -Inf, y = Inf, hjust = -0.5, vjust = 1,size=8)+
ylim(c(0,10))


ukb_pheno_mri$prs_decile = cut2(ukb_pheno_mri$NOHLA_PRS,g=10)
prs_decile_model = glm(data=ukb_pheno_mri,
            MS_status~`Age at recruitment.0.0`+
              `Sex.0.0`+
              `Genetic principal components.0.1`+
              `Genetic principal components.0.2`+
              `Genetic principal components.0.3`+
              `Genetic principal components.0.4`+
              prs_decile,
            family=binomial(link="logit"))

tbl = data.frame(summary(prs_decile_model)$coefficients[-c(1:7),])
tbl$decile=c(2:10)
tbl$or=exp(tbl$Estimate)
tbl$lower_ci=exp(tbl$Estimate-1.96*tbl$Std..Error)
tbl$upper_ci=exp(tbl$Estimate+1.96*tbl$Std..Error)

no_hla_prs_decile_plot=ggplot(tbl,aes(decile,or))+
geom_errorbar(aes(x=decile,ymin=lower_ci,ymax=upper_ci,width=0.2))+
geom_point(size=3,shape=22,fill="black")+
theme_classic()+
theme(legend.position="none",text=element_text(size=12))+
labs(x="PRS Decile",y="OR for MS (vs lowest decile)")+
annotate("text",label="MHC Excluded",x = -Inf, y = Inf, hjust = -0.5, vjust = 1,size=8)+
ylim(c(0,10))



png("decile_plot.png",height=8,width=8,res=300,units="in")
grid.arrange(no_hla_prs_decile_plot,prs_decile_plot,nrow=1)
dev.off()
````
Here is the decile plot:
![decile plot](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/decile_plot.png)

Now let's see the distributions of PRS:
````R
hist_nohla = ggplot(ukb_pheno_mri,aes(NOHLA_PRS,fill=factor(MS_status)))+
geom_density(alpha=0.5)+
scale_fill_brewer(palette ="Set2",labels=c("Controls","MS"))+
theme_classic()+
theme(text=element_text(size=12))+
labs(x="PRS",y="Density",fill="MS status")+
annotate("text",label="MHC Excluded",x = -Inf, y = Inf, hjust = -0.5, vjust = 1,size=5)

hist_hla = ggplot(ukb_pheno_mri,aes(HLA_PRS,fill=factor(MS_status)))+
geom_density(alpha=0.5)+
scale_fill_brewer(palette ="Set2",labels=c("Controls","MS"))+
theme_classic()+
theme(text=element_text(size=12))+
labs(x="PRS",y="Density",fill="MS status")+
annotate("text",label="MHC Included",x = -Inf, y = Inf, hjust = -0.5, vjust = 1,size=5)


png("prs_hist.png",height=8,width=8,res=300,units="in")
grid.arrange(hist_nohla,hist_hla,nrow=1)
dev.off()
````
Here is the PRS density plot:
![density_plot](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/prs_hist.png)

Now the calibration plots:

````R
# roc analysis
library(ROCR)
hlaprs_model = glm(data=ukb_pheno_mri,
            MS_status~`Age at recruitment.0.0`+
              `Sex.0.0`+
              `Genetic principal components.0.1`+
              `Genetic principal components.0.2`+
              `Genetic principal components.0.3`+
              `Genetic principal components.0.4`+
              HLA_PRS,
            family=binomial(link="logit"))

nohla_prs_model = glm(data=ukb_pheno_mri,
            MS_status~`Age at recruitment.0.0`+
              `Sex.0.0`+
              `Genetic principal components.0.1`+
              `Genetic principal components.0.2`+
              `Genetic principal components.0.3`+
              `Genetic principal components.0.4`+
              NOHLA_PRS,
            family=binomial(link="logit"))

null_model = glm(data=ukb_pheno_mri,
            MS_status~`Age at recruitment.0.0`+
              `Sex.0.0`+
              `Genetic principal components.0.1`+
              `Genetic principal components.0.2`+
              `Genetic principal components.0.3`+
              `Genetic principal components.0.4`,
            family=binomial(link="logit"))

make_prediction = function(x){
predict(x,newdata=ukb_pheno_mri,type="response")
}

preds = data.frame("hlaprs_model"=make_prediction(hlaprs_model),
"nohla_prs_model"=make_prediction(nohla_prs_model),
"null_model"=make_prediction(null_model),
"prs_decile"=ukb_pheno_mri$prs_decile,
"MS_status"=ukb_pheno_mri$MS_status)

preds_summary = preds %>% group_by(prs_decile) %>% summarise(mean(hlaprs_model,na.rm=TRUE),mean(nohla_prs_model,na.rm=TRUE),mean(null_model,na.rm=TRUE))
tbl = table(ukb_pheno_mri$prs_decile,ukb_pheno_mri$MS_status)
obs_risk = tbl[,2]/rowSums(tbl)

pred_df = data.frame(preds_summary,obs_risk)
pred_df$prs_decile = c(1:10)
library(reshape2)
pred_df = melt(pred_df,id="prs_decile")
pred_df$prs_decile = factor(pred_df$prs_decile)

calib_plot = ggplot(pred_df,aes(prs_decile,value,col=variable,group=variable))+geom_point(col="black",shape=22,alpha=0.3)+geom_line()+labs(col="Source of risk estimate",x="PRS decile",y="MS risk (probability scale)")+scale_color_brewer(palette="Set2",labels=c("HLA PRS model","Non-HLA PRS model","Null model (Age, Sex, PCs)","Observed risk"))+theme_classic()

png("calibration_plot.png",height=8,width=8,res=300,units="in")
calib_plot
dev.off()
````
Here is the calibration plot:
![calib](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/calibration_plot.png)

And finally the ROC curves:
````R
# auc

preds = preds %>% na.omit()
preds$MS=factor(preds$MS_status)
predictions = prediction(list(preds$hlaprs_model,preds$nohla_prs_model,preds$null_model),
list(preds$MS,preds$MS,preds$MS))
roc.perf = ROCR::performance(predictions, measure = "tpr", x.measure = "fpr")

hla_prs_model_auc = data.frame(x=roc.perf@x.values[[1]],y=roc.perf@y.values[[1]],model="HLA PRS model")
nohla_prs_model_auc = data.frame(x=roc.perf@x.values[[2]],y=roc.perf@y.values[[2]],model="Non-HLA PRS model")
null_model_auc = data.frame(x=roc.perf@x.values[[3]],y=roc.perf@y.values[[3]],model="Null model")

df = bind_rows(hla_prs_model_auc,nohla_prs_model_auc,null_model_auc)
auc.perf = performance(predictions, measure = "auc")
auc.perf@y.values

aucplot=ggplot(df,aes(x,y,col=model))+
geom_line()+
scale_color_brewer(palette="Set2")+
labs(x="False Positive Rate",y="True Positive Rate",col="Model")+
theme_classic()+
geom_abline()+
annotate("text",x=0.8,y=0.3,label="AUC",hjust=0)+
annotate("text",x=1,y=0.25,label=paste0("HLA PRS model: ",round(auc.perf@y.values[[1]],3)),hjust=1)+
annotate("text",x=1,y=0.2,label=paste0("Non-HLA PRS model: ",round(auc.perf@y.values[[2]],3)),hjust=1)+
annotate("text",x=1,y=0.15,label=paste0("Null model: ",round(auc.perf@y.values[[3]],3)),hjust=1)

png("discrimination_plot.png",height=8,width=8,res=300,units="in")
aucplot
dev.off()
````
Here are the ROC curves:
![roc](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/discrimination_plot.png)


## 4. Correlation between PRS and MRI metrics

Now let's look at the MRI data itself. First we need to make sure that we exclude people with Alzheimer's, PD, and cerebral small vessel disease as we are interested in as healthy a cohort as possible.  

````R
# First exclude strokes and other neurodegenerative diseases
# read in data
source_of_report_data = read_tsv("../ukb_pheno_17_03_21/ukb_pheno_nocognitive_17032021.tsv.gz",col_types=cols_only(
  `Source of report of G20 (parkinson's disease).0.0` = col_character(),
  `Source of report of G30 (alzheimer's disease).0.0` = col_character(),
  `Source of report of G45 (transient cerebral ischaemic attacks and related syndromes).0.0` = col_character(),
  EID = col_double()))

source_of_report_data = source_of_report_data %>%
mutate(pd_status = ifelse(!is.na(`Source of report of G20 (parkinson's disease).0.0`),1,0)) %>%
mutate(ad_status = ifelse(!is.na(`Source of report of G30 (alzheimer's disease).0.0`),1,0)) %>%
mutate(stroke_status = ifelse(!is.na(`Source of report of G45 (transient cerebral ischaemic attacks and related syndromes).0.0`),1,0))

table(source_of_report_data$pd_status)
table(source_of_report_data$ad_status)
table(source_of_report_data$stroke_status)

any_neuro_disease = source_of_report_data %>% filter(pd_status==1 | ad_status ==1 | stroke_status ==1)

# filter out neurodegenerative diseases
ukb_pheno_mri = ukb_pheno_mri %>% filter(!EID %in% any_neuro_disease)
table(ukb_pheno_mri$MS_status)
````

The MRI dataset after excluding people with AD, PD, and cerebral SVD:
Controls     |     Cases
------------ | --------------
29988        |  124


Now we'll do some basic sense checks.
````R
# combine MRI and rest of phenotype data
mri = ukb_pheno_mri %>% left_join(mri_data,by="EID")

# sense checks
age_vs_brainvol = ggplot(mri,aes(`Age at recruitment.0.0`,`Volume of brain, grey+white matter (normalised for head size).2.0`))+geom_smooth()+geom_point()
png("age_brainvol.png",height=8,width=8,res=300,units="in")
age_vs_brainvol
dev.off()

age_vs_brainvol_bysex = ggplot(mri,aes(`Age at recruitment.0.0`,`Volume of brain, grey+white matter (normalised for head size).2.0`))+geom_smooth()+geom_point()+facet_wrap(~Sex.0.0)
png("age_brainvol.png",height=8,width=8,res=300,units="in")
age_vs_brainvol
dev.off()

````
1. Brain volume is lower in older people:
![brainvol](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/age_brainvol.png)

````R
# Do pwMS have higher WM lesion volume
png("wm_lesion_plot_ms_v_control.png",height=8,width=8,res=300,units="in")
ggplot(mri,aes(factor(MS_status),`Total volume of white matter hyperintensities (from T1 and T2_FLAIR images).2.0`,fill=factor(MS_status)))+
stat_summary(fun=median,geom="crossbar")+
geom_violin(alpha=0.5)+
theme_bw()+
labs(x="MS status",y="Normalised White Matter Lesion Volume",fill="MS status")+
scale_fill_brewer(palette="Set2",labels=c("Controls","MS"))+
theme(axis.text.x=element_blank())
dev.off()
````
2. People with MS have higher T2 hyperintensity volume (presumably this is picking up old MS lesions).
Interestingly there seem to be some controls (a small number) with very high lesion load, but the average T2 lesion load is clearly higher in MS, as we would expect.
![ms_lesions](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/wm_lesion_plot_ms_v_control.png)


````R
fa_vars = mri %>% select(EID,`Age at recruitment.0.0`,Sex.0.0,MS_status,contains("Mean FA")) %>%
 select(EID,`Age at recruitment.0.0`,Sex.0.0,MS_status,contains(".2.0")) %>%
 select(-contains("Weighted"))


overall_coefs = data.frame()
make_model = function(x){
  fa_vars_no_missing = fa_vars %>% filter(!is.na(fa_vars[[x]]))
  model = glm(data=fa_vars_no_missing, rankNorm(fa_vars_no_missing[[x]]) ~ `Age at recruitment.0.0`+ `Sex.0.0`+ MS_status)
  coefs = summary(model)$coefficients[4,]
  overall_coefs <<- bind_rows(overall_coefs,coefs)
}
sapply(colnames(fa_vars)[-c(1:4)],make_model)

overall_coefs = as.tbl(overall_coefs)
overall_coefs$variable = colnames(fa_vars)[-c(1:4)]

png("ms_vs_control_fa.png",height=8,width=8,res=300,units="in")
ggplot(overall_coefs,aes(Estimate,variable))+
geom_point(aes(size=-log10(`Pr(>|t|)`)))+
geom_errorbarh(aes(y=variable,xmin=Estimate-1.96*`Std. Error`,xmax=Estimate+1.96*`Std. Error`),height=0.1)+
geom_vline(xintercept=0,alpha=0.5)+
theme_bw()+
labs(x="Effect of MS status on normalised FA variable")
dev.off()

````
3. MS is strongly associated with widespread FA changes, consistent with [previous findings in normal-appearing cortex vs non-neurologic controls](https://academic.oup.com/brain/article/142/7/1921/5511700)
![ms_fa_lesion](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/ms_vs_control_fa.png)

So the data look good. Now we can try to replicate the main findings we are interested in. First we'll look at white matter hyperintensities

````R
# WM hyperintensities
# first let's remove people with MS
mri_noms = mri %>% filter(MS_status==0)

# normalise WM lesions
mri_noms = mri_noms %>% filter(!is.na(`Total volume of white matter hyperintensities (from T1 and T2_FLAIR images).2.0`))
mri_noms$norm_wm_lesions = rankNorm(mri_noms$`Total volume of white matter hyperintensities (from T1 and T2_FLAIR images).2.0`)

hla_prs = ggplot(mri_noms,aes(HLA_PRS,norm_wm_lesions))+
geom_point()+labs(x="HLA PRS",y="Normalised WM lesion volume")+theme_bw()

nohla_prs = ggplot(mri_noms,aes(NOHLA_PRS,norm_wm_lesions))+
geom_point()+labs(x="Non-HLA PRS",y="Normalised WM lesion volume")+theme_bw()

png("wmlesions_prs_plot.png",res=300,units="in",height=8,width=8)
grid.arrange(hla_prs,nohla_prs,nrow=1)
dev.off()
````

There is no obvious correlation between PRS and normalised WM lesion volume.
![WM lesion plot](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/wmlesions_prs_plot.png)

There is no statistical evidence of such an association either:
````R
# Does MS PRS predict WM lesion volume in healthy people
hla_model = glm(data=mri_noms,norm_wm_lesions~`Age at recruitment.0.0`+Sex.0.0+`Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+`Volume of brain, grey+white matter (normalised for head size).2.0`+`Volume of ventricular cerebrospinal fluid (normalised for head size).2.0` +HLA_PRS)
summary(hla_model)$coefficients["HLA_PRS",]
Estimate   Std. Error      t value     Pr(>|t|)
-0.004833947  0.005107555 -0.946430771  0.343937218

nohla_model = glm(data=mri_noms,norm_wm_lesions~`Age at recruitment.0.0`+Sex.0.0+`Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+`Volume of brain, grey+white matter (normalised for head size).2.0`+`Volume of ventricular cerebrospinal fluid (normalised for head size).2.0` +NOHLA_PRS)
summary(nohla_model)$coefficients["NOHLA_PRS",]

Estimate  Std. Error     t value    Pr(>|t|)
0.003524292 0.005110308 0.689643836 0.490424002
````

As a sensitivity analysis, let's repeat this analysis over a range of PRS parameters:
````R
pvals = c("1","0.8","0.6","0.4","0.2","0.1","0.05","0.005","0.0005","0.00005","0.000005","0.0000005","0.00000005")
r2_vector = c(0.1,0.2,0.4,0.6,0.8)
pval_summary = list()
for(hla in c("yes","no")){
  for(pval in pvals){
    for(r2 in r2_vector){
      message("calculating WM association with PRS")
      message("R2 = ",r2)
      message("P = ",pval)
      prs = read_table2(paste0("/data/Wolfson-UKBB-Dobson/ms_prs/summarised_PRS_results_pval",pval,"r2",r2))
      colnames(prs)=c("EID","PRS")
      mri_noms = mri_noms %>% select(-contains("PRS"))
      mri_noms = mri_noms %>% left_join(prs,by="EID")
      mri_noms = mri_noms %>% filter(!(is.na(PRS)))
      mri_noms = mri_noms %>% mutate(combined_score = PRS+DRB_risk_score)

      mri_noms$PRS=if(hla=="no"){
        message("Excluding HLA")
        rankNorm(mri_noms$PRS)
      } else if(hla=="yes"){
        message("Including HLA")
        rankNorm(mri_noms$combined_score)
      }


      # make model
      prs_model = glm(data=mri_noms,norm_wm_lesions~`Age at recruitment.0.0`+Sex.0.0+`Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+`Volume of brain, grey+white matter (normalised for head size).2.0`+`Volume of ventricular cerebrospinal fluid (normalised for head size).2.0` + PRS)
      pval_res = summary(prs_model)$coefficients["PRS","Pr(>|t|)"]
      beta_res = summary(prs_model)$coefficients["PRS","Estimate"]
      se_res = summary(prs_model)$coefficients["PRS","Std. Error"]
      pval_summary[[length(pval_summary)+1]] = c(as.numeric(pval),as.numeric(r2),as.numeric(pval_res),as.numeric(beta_res),as.numeric(se_res))
      message("P value for this model is ",pval_res)
    }
  }
}

res_df = data.frame(do.call("rbind",pval_summary))
colnames(res_df) = c("PRS_P","PRS_R2","P_val","Beta","SE")
res_df$hla = c(rep("HLA",nrow(res_df)/2),rep("NO HLA",nrow(res_df)/2))
res_df$fdr = p.adjust(res_df$P_val,method="fdr")
res_df$sig = ifelse(res_df$fdr<0.1,"*","NS")

png("sensitivity_wm_prs_plot.png",res=300,units="in",height=8,width=16)
ggplot(res_df,aes(factor(PRS_P),factor(PRS_R2),label=round(P_val,2),fill=-log10(P_val)))+
geom_tile()+
scale_fill_viridis_c()+
geom_text()+
theme_bw()+
facet_wrap(~hla)+
labs(x="PRS P value threshold",y="PRS clumping R2 threshold")
dev.off()
````
None of these associations surpass an FDR threshold of 10%.
![PRS WM lesion plot](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/sensitivity_wm_prs_plot.png)
So there is no association between the MS-PRS and WM hyperintensities in healthy (non-MS) UKB controls.

What about regional mean Fractional Anisotropy?
````R

fa_vars = mri_noms %>% select(EID,`Age at recruitment.0.0`,Sex.0.0,MS_status,HLA_PRS,NOHLA_PRS,DRB_risk_score,`Genetic principal components.0.1`,`Genetic principal components.0.2`,`Genetic principal components.0.3`,`Genetic principal components.0.4`,contains("Mean FA")) %>%
 select(EID,`Age at recruitment.0.0`,Sex.0.0,MS_status,,HLA_PRS,NOHLA_PRS,DRB_risk_score,`Genetic principal components.0.1`,`Genetic principal components.0.2`,`Genetic principal components.0.3`,`Genetic principal components.0.4`,contains(".2.0")) %>%
 select(-contains("Weighted"))


# HLA PRS
overall_coefs = data.frame()
make_model = function(x){
  fa_vars_no_missing = fa_vars %>% filter(!is.na(fa_vars[[x]]))
  model = glm(data=fa_vars_no_missing, rankNorm(fa_vars_no_missing[[x]]) ~ `Age at recruitment.0.0`+ `Sex.0.0`+ `Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+ HLA_PRS)
  coefs = summary(model)$coefficients["HLA_PRS",]
  overall_coefs <<- bind_rows(overall_coefs,coefs)
}
sapply(colnames(fa_vars)[-c(1:11)],make_model)

overall_coefs = as.tbl(overall_coefs)
overall_coefs$variable = colnames(fa_vars)[-c(1:11)]

hla_fa_plot = ggplot(overall_coefs,aes(Estimate,variable))+
geom_point(aes(size=-log10(`Pr(>|t|)`)))+
geom_errorbarh(aes(y=variable,xmin=Estimate-1.96*`Std. Error`,xmax=Estimate+1.96*`Std. Error`),height=0.1)+
geom_vline(xintercept=0,alpha=0.5)+
theme_bw()+
labs(x="Effect of HLA PRS on normalised FA variable")
hla_coefs = overall_coefs

# NOHLA PRS
overall_coefs = data.frame()
make_model = function(x){
  fa_vars_no_missing = fa_vars %>% filter(!is.na(fa_vars[[x]]))
  model = glm(data=fa_vars_no_missing, rankNorm(fa_vars_no_missing[[x]]) ~ `Age at recruitment.0.0`+ `Sex.0.0`+ `Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+ NOHLA_PRS)
  coefs = summary(model)$coefficients["NOHLA_PRS",]
  overall_coefs <<- bind_rows(overall_coefs,coefs)
}
sapply(colnames(fa_vars)[-c(1:11)],make_model)

overall_coefs = as.tbl(overall_coefs)
overall_coefs$variable = colnames(fa_vars)[-c(1:11)]

nonhla_fa_plot = ggplot(overall_coefs,aes(Estimate,variable))+
geom_point(aes(size=-log10(`Pr(>|t|)`)))+
geom_errorbarh(aes(y=variable,xmin=Estimate-1.96*`Std. Error`,xmax=Estimate+1.96*`Std. Error`),height=0.1)+
geom_vline(xintercept=0,alpha=0.5)+
theme_bw()+
labs(x="Effect of Non-HLA PRS on normalised FA variable")
nohla_coefs = overall_coefs

png("fa_prs_plot.png",res=300,units="in",height=8,width=16)
grid.arrange(hla_fa_plot,nonhla_fa_plot,nrow=1)
dev.off()

overall_coefs = bind_rows(nohla_coefs,hla_coefs)
overall_coefs$fdr = p.adjust(overall_coefs$`Pr(>|t|)`)

````
The answer is no. None of these associations pass the multiple testing threshold (Alpha = 0.05, Bonferroni correction).
![fa_plot](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/fa_prs_plot.png)
