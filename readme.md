# No evidence for association between polygenic risk of Multiple Sclerosis and MRI phenotypes in ~30,000 healthy UK Biobank participants

A quick replication effort to see if the results reported [here](https://journals.sagepub.com/doi/full/10.1177/13524585211034826) can be reproduced in the UKB imaging cohort.

Steps:
0. Generate polygenic risk scores
1. Divide the UKB cohort into training (not included in MRI) and testing (MRI data available) sets
2. Determine the optimal MS PRS in the training cohort



# 0. Generate polygenic risk scores
We generated polygenic risk scores in UKB using the IMSGC meta-analysis discovery stage summary statistics. We used the clumping-and-thresholding method to create a variety of scores. For a full description of our methods please see the [code](https://github.com/benjacobs123456/MS_UKB_PRS) and read the [manuscript](https://nn.neurology.org/content/8/4/e1007).

# 1. Divide the UKB cohort into two sets - the non-imaging cohort (to select the best PRS, i.e. training set), and the imaging cohort (validation or testing set)
````R
#############################################
#               Load packages
#############################################

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

  Source of report of G20 (parkinson's disease)
  	Source of report of G30 (alzheimer's disease)
    Source of report of G45 (transient cerebral ischaemic attacks and related syndromes)

# define MS status
source_of_report_data = source_of_report_data %>% mutate(MS_status = ifelse(!is.na(`Source of report of G35 (multiple sclerosis).0.0`),1,0))

# read in remainder of phenotype data
ukb_pheno = read_tsv("../ukb_pheno_911/ukb_pheno_final_MS_1301")
ukb_pheno = ukb_pheno %>% select(EID,`Age at recruitment.0.0`,Sex.0.0,`Ethnic background.0.0`,contains("rincipal components"),contains("ethnic grouping"))

# combine
ukb_pheno = ukb_pheno %>% left_join(source_of_report_data,by="EID")

# read in mri data
mri_data = read_tsv("mri_data.tsv.gz")

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
------------ | --------------
407478       |   2102

````R
# Filter highly related individuals
# filter relatedness
# nb this is taking care to exclude the non-ms control from each pair to boost case number
kin = read_table2("/data/Wolfson-UKBB-Dobson/helper_progs_and_key/ukb43101_rel_s488282.dat")
highly_related = kin %>% filter(Kinship>0.0442) %>% filter(ID1 %in% ukb_pheno$EID) %>% filter(ID2 %in% ukb_pheno$EID)
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
336766       |   2102
````R
# define non-MRI dataset
ukb_pheno_nomri = ukb_pheno %>% filter(!EID %in% mri_data$EID)
table(ukb_pheno_nomri$MS_status)
````
In the non-MRI dataset:
Controls     |     Cases
------------ | --------------
308905       |  1962
````R
ukb_pheno_mri = ukb_pheno %>% filter(EID %in% mri_data$EID)
table(ukb_pheno_mri$MS_status)
````
And in the non-MRI dataset:
In the non-MRI dataset:
Controls     |     Cases
------------ | --------------
27861        |  140

# 2. Select the best PRS in the UKB non-imaging cohort
Now we will test the PRS in the non-MRI cohort and select the score with the best explanatory power for MS susceptibility in this cohort.
````R
# Define scoring function
score_prs = function(pval,r2){
# read in prs
filename = paste0("/data/Wolfson-UKBB-Dobson/ms_prs/summarised_PRS_results_pval",pval,"r2",r2)
prs = read_table2(filename)
colnames(prs) = c("EID","PRS")
prs_tuning_dataset = ukb_pheno_nomri %>% select(`Age at recruitment.0.0`,`Sex.0.0`,`Genetic principal components.0.1`,`Genetic principal components.0.2`,`Genetic principal components.0.3`,`Genetic principal components.0.4`,MS_status,EID)
prs_tuning_dataset$EID = as.numeric(as.character(prs_tuning_dataset$EID))
prs_tuning_dataset = prs_tuning_dataset %>% left_join(prs,by="EID")
prs_tuning_dataset = prs_tuning_dataset %>% filter(!(is.na(PRS)))
prs_tuning_dataset$PRS=rankNorm(prs_tuning_dataset$PRS)

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
return(nag)
}

# Loop scoring function over different PRS
pvals = c("1","0.8","0.6","0.4","0.2","0.1","0.05","0.00000005")
r2 = c(0.2,0.4,0.6,0.8)
params = expand.grid(pvals,r2)

nags = mapply(score_prs,pval=params[,1],r2=params[,2])

params = data.frame(params)
params$Nagelkerke_Pseudo_R2 = nags
colnames(params)=c("Pval","R2","Nagelkerke_Pseudo_R2")
hla_scores = params

# repeat with hla excluded
score_prs_nohla = function(pval,r2){
# read in prs
filename = paste0("/data/Wolfson-UKBB-Dobson/ms_prs/nomhc_summarised_PRS_results_pval",pval,"r2",r2)
prs = read_table2(filename)
colnames(prs) = c("EID","PRS")
prs_tuning_dataset = ukb_pheno_nomri %>% select(`Age at recruitment.0.0`,`Sex.0.0`,`Genetic principal components.0.1`,`Genetic principal components.0.2`,`Genetic principal components.0.3`,`Genetic principal components.0.4`,MS_status,EID)
prs_tuning_dataset$EID = as.numeric(as.character(prs_tuning_dataset$EID))
prs_tuning_dataset = prs_tuning_dataset %>% left_join(prs,by="EID")
prs_tuning_dataset = prs_tuning_dataset %>% filter(!(is.na(PRS)))
prs_tuning_dataset$PRS=rankNorm(prs_tuning_dataset$PRS)

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
return(nag)
}

pvals = c("1","0.8","0.6","0.4","0.2","0.1","0.05","0.00000005")
r2 = c(0.2,0.4,0.6,0.8)
params = expand.grid(pvals,r2)

nags = mapply(score_prs_nohla,pval=params[,1],r2=params[,2])

params = data.frame(params)
params$Nagelkerke_Pseudo_R2 = nags
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

png("nagel_plot.png",height=8,width=12,res=300,units="in")
nagel_plot
dev.off()
````
The best PRS explain ~ 3.4% (with MHC) and ~1.5% (without MHC) of MS liability in the training (non-MRI) cohort.

````R
params %>% arrange(desc(Nagelkerke_Pseudo_R2))
         Pval  R2 Nagelkerke_Pseudo_R2          MHC
1         0.4 0.8           0.03444300 MHC included
2         0.6 0.8           0.03409030 MHC included
3         0.8 0.8           0.03384930 MHC included
4           1 0.8           0.03383730 MHC included
5         0.2 0.8           0.03358050 MHC included
6         0.1 0.8           0.03274790 MHC included
7        0.05 0.8           0.03205990 MHC included
8        0.05 0.6           0.03117320 MHC included
9         0.4 0.6           0.03072340 MHC included
10        0.1 0.6           0.03053970 MHC included
11        0.2 0.6           0.03027150 MHC included
12        0.6 0.6           0.02984240 MHC included
13        0.8 0.6           0.02938830 MHC included
14          1 0.6           0.02936110 MHC included
15       0.05 0.4           0.02727050 MHC included
16        0.1 0.4           0.02474370 MHC included
17        0.4 0.4           0.02372710 MHC included
18        0.2 0.4           0.02349400 MHC included
19        0.6 0.4           0.02277380 MHC included
20        0.8 0.4           0.02222430 MHC included
21          1 0.4           0.02221840 MHC included
22 0.00000005 0.8           0.02217370 MHC included
23 0.00000005 0.6           0.02183440 MHC included
24 0.00000005 0.4           0.02118930 MHC included
25 0.00000005 0.2           0.01906000 MHC included
26       0.05 0.2           0.01664680 MHC included
27        0.1 0.2           0.01385020 MHC included
28        0.4 0.2           0.01337980 MHC included
29        0.2 0.2           0.01309840 MHC included
30        0.6 0.2           0.01293740 MHC included
31        0.8 0.2           0.01265030 MHC included
32          1 0.2           0.01263120 MHC included
33       0.05 0.8           0.01254230 MHC excluded
34        0.4 0.8           0.01243330 MHC excluded
35        0.6 0.8           0.01194590 MHC excluded
36        0.2 0.8           0.01190250 MHC excluded
37        0.1 0.8           0.01189930 MHC excluded
38        0.8 0.8           0.01168780 MHC excluded
39          1 0.8           0.01168190 MHC excluded
40        0.4 0.6           0.01037480 MHC excluded
41        0.6 0.6           0.00991984 MHC excluded
42       0.05 0.6           0.00973293 MHC excluded
43        0.8 0.6           0.00966366 MHC excluded
44          1 0.6           0.00965404 MHC excluded
45        0.2 0.6           0.00929053 MHC excluded
46        0.1 0.6           0.00907070 MHC excluded
47        0.4 0.4           0.00811307 MHC excluded
48        0.6 0.4           0.00779066 MHC excluded
49          1 0.4           0.00752575 MHC excluded
50        0.8 0.4           0.00752307 MHC excluded
51       0.05 0.4           0.00734313 MHC excluded
52        0.2 0.4           0.00697575 MHC excluded
53        0.1 0.4           0.00667549 MHC excluded
54 0.00000005 0.8           0.00564719 MHC excluded
55        0.4 0.2           0.00554609 MHC excluded
56        0.6 0.2           0.00544342 MHC excluded
57        0.8 0.2           0.00529467 MHC excluded
58          1 0.2           0.00528087 MHC excluded
59       0.05 0.2           0.00482235 MHC excluded
60        0.2 0.2           0.00469673 MHC excluded
61 0.00000005 0.6           0.00455533 MHC excluded
62        0.1 0.2           0.00425953 MHC excluded
63 0.00000005 0.4           0.00351158 MHC excluded
64 0.00000005 0.2           0.00248418 MHC excluded
````

3. Validate the PRS in the MRI cohort

````R
# Even though there are only 140 MS cases in the MRI cohort, let's just check the PRS is working reasonably well.

# choose best prs based on r2 and read it in
prs = read_table2("/data/Wolfson-UKBB-Dobson/ms_prs/summarised_PRS_results_pval0.4r20.8")
colnames(prs)=c("EID","PRS")
nohla_prs = read_table2("/data/Wolfson-UKBB-Dobson/ms_prs/nomhc_summarised_PRS_results_pval0.05r20.8")
colnames(nohla_prs)=c("EID","PRS")

ukb_pheno_mri = ukb_pheno_mri %>% select(-contains("PRS"))
ukb_pheno_mri = ukb_pheno_mri %>% left_join(prs,by="EID")
ukb_pheno_mri = ukb_pheno_mri %>% filter(!(is.na(PRS)))
ukb_pheno_mri$HLA_PRS=rankNorm(ukb_pheno_mri$PRS)
ukb_pheno_mri = ukb_pheno_mri %>% select(-PRS)
ukb_pheno_mri = ukb_pheno_mri %>% left_join(nohla_prs,by="EID")
ukb_pheno_mri = ukb_pheno_mri %>% filter(!(is.na(PRS)))
ukb_pheno_mri$NOHLA_PRS=rankNorm(ukb_pheno_mri$PRS)

# nagelkerke
null_model = glm(data=ukb_pheno_mri,
               MS_status~`Age at recruitment.0.0`+
                 `Sex.0.0`+
                 `Genetic principal components.0.1`+
                 `Genetic principal components.0.2`+
                 `Genetic principal components.0.3`+
                 `Genetic principal components.0.4`,
               family=binomial(link="logit"))

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

print("printing nagelkerke for best PRS in testing set")
print(nagelkerke(hlaprs_model,null_model))
````
Clearly the PRS works well in discriminating MS cases from controls in the MRI dataset (or as well as can be expected for PRS):
Parameter  |  PRS (no HLA)   | PRS (with HLA)
--- | --- | ---
Nagelkerke's Pseudo-R2 | 0.0161361 |  0.03975300
Likelihood ratio P value | 2.3576e-07 | 4.8362e-16

````R
library(Hmisc)

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
ylim(c(0,25))


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
tbl
write_csv(tbl,"nohla_decile_ORs.csv")

no_hla_prs_decile_plot=ggplot(tbl,aes(decile,or))+
geom_errorbar(aes(x=decile,ymin=lower_ci,ymax=upper_ci,width=0.2))+
geom_point(size=3,shape=22,fill="black")+
theme_classic()+
theme(legend.position="none",text=element_text(size=12))+
labs(x="PRS Decile",y="OR for MS (vs lowest decile)")+
annotate("text",label="MHC Excluded",x = -Inf, y = Inf, hjust = -0.5, vjust = 1,size=8)+
ylim(c(0,25))


library(gridExtra)

png("decile_plot.png",height=8,width=8,res=300,units="in")
grid.arrange(no_hla_prs_decile_plot,prs_decile_plot,nrow=1)
dev.off()


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

4. Correlation between PRS and MRI metrics
Now let's look at the MRI data itself.

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

#
mri = ukb_pheno_mri %>% left_join(mri_data,by="EID") %>% filter(!is.na(wm_lesion_vol))

# First let's do a sense check and see if pwMS have higher WM lesion volume
png("wm_lesion_plot_ms_v_control.png",height=8,width=8,res=300,units="in")
ggplot(mri,aes(factor(MS_status),wm_lesion_vol,fill=factor(MS_status)))+
stat_summary(fun=median,geom="crossbar")+
geom_violin(alpha=0.5)+
theme_bw()+
labs(x="MS status",y="Normalised White Matter Lesion Volume",fill="MS status")+
scale_fill_brewer(palette="Set2",labels=c("Controls","MS"))+
theme(axis.text.x=element_blank())
dev.off()

# Let's normalise WM lesion volume before we make any regression models
mri$norm_wm_lesions = rankNorm(mri$wm_lesion_vol)

# does MS status predict WM lesion volume?
model = glm(data=mri,norm_wm_lesions~
  `Age at recruitment.0.0`+Sex.0.0+`Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+csf_vol_normalised_headvol+grey_and_white_vol_normalised_headvol+MS_status)
anova(model,test="Chisq")

# The answer is yes:
# Likelihood ration pval 2.2e-16
summary(model)

or = exp(summary(model)$coefficients[10,1])
lowerci = exp(summary(model)$coefficients[10,1]-1.96*summary(model)$coefficients[10,2])
upperci = exp(summary(model)$coefficients[10,1]+1.96*summary(model)$coefficients[10,2])
paste(round(or,3),round(lowerci,3),round(upperci,3))

# OR 3.288, 95% CI 2.85 - 3.794

# Does MS PRS predict WM lesion volume in healthy people
model = glm(data=mri,norm_wm_lesions~`Age at recruitment.0.0`+Sex.0.0+`Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+csf_vol_normalised_headvol+grey_and_white_vol_normalised_headvol+HLA_PRS)
summary(model)

mhc = ggplot(ms_mri,aes(HLA_PRS,norm_wm_lesions))+geom_point()+theme_classic()+
labs(x="MHC PRS",y="Normalised WM T2/FLAIR lesion load")
nomhc = ggplot(ms_mri,aes(NOHLA_PRS,norm_wm_lesions))+geom_point()+theme_classic()+
labs(x="Non-MHC PRS",y="Normalised WM T2/FLAIR lesion load")

png("mri_prs_plot.png",res=300,units="in",height=8,width=8)
grid.arrange(mhc,nomhc,nrow=1)
dev.off()

# repeat with controls
mri_data = read_table2("mri_data.tsv")
overall_mri = selected_vars %>% filter(EID %in% mri_data$EID)
overall_mri = overall_mri %>% left_join(mri_data,by="EID") %>% filter(!is.na(wm_lesion_vol))
overall_mri$norm_wm_lesions = rankNorm(overall_mri$wm_lesion_vol)


mhc = ggplot(overall_mri,aes(HLA_PRS,norm_wm_lesions,col=factor(MS_status)))+geom_point(alpha=0.5)+theme_classic()+
labs(x="MHC PRS",y="Normalised WM T2/FLAIR lesion load",col="MS status")+scale_color_discrete(labels=c("Controls","MS cases"))
nomhc = ggplot(overall_mri,aes(NOHLA_PRS,norm_wm_lesions,col=factor(MS_status)))+geom_point(alpha=0.5)+theme_classic()+
labs(x="Non-MHC PRS",y="Normalised WM T2/FLAIR lesion load",col="MS status")+scale_color_discrete(labels=c("Controls","MS cases"))

png("healthy_indivs_mri_prs_plot.png",res=300,units="in",height=8,width=8)
grid.arrange(mhc,nomhc,nrow=1)
dev.off()

model = glm(data=overall_mri,norm_wm_lesions~`Age at recruitment.0.0`+Sex.0.0+`Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+csf_vol_normalised_headvol+grey_and_white_vol_normalised_headvol+HLA_PRS)
summary(model)

model = glm(data=overall_mri,norm_wm_lesions~`Age at recruitment.0.0`+ Sex.0.0+`Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+csf_vol_normalised_headvol+grey_and_white_vol_normalised_headvol+NOHLA_PRS)
summary(model)
