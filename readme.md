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
ukb_pheno_nomri = ukb_pheno %>% filter(!EID %in% mri_data$EID)
table(ukb_pheno_nomri$MS_status)
````

In the non-MRI dataset:

Controls     |     Cases
------------ | --------------
340537       |  1938

````R
ukb_pheno_mri = ukb_pheno %>% filter(EID %in% mri_data$EID)
table(ukb_pheno_mri$MS_status)
````

And in the MRI dataset:

Controls     |     Cases
------------ | --------------
36050        |  164

## 2. Select the best PRS in the UKB non-imaging cohort
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
Here's a plot of the prediction pseudo-R2 values for each PRS in the training set (roughly equivalent to variance explained)
![Nagel plot](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/nagel_plot.png)

The best PRS explain ~ 3.3% (with MHC) and ~1.2% (without MHC) of MS liability in the training (non-MRI) cohort:

````R
params %>% arrange(desc(Nagelkerke_Pseudo_R2))
        Pval  R2 Nagelkerke_Pseudo_R2          MHC
1         0.4 0.8           0.03366990 MHC included
2         0.6 0.8           0.03334440 MHC included
3         0.8 0.8           0.03311840 MHC included
4           1 0.8           0.03309640 MHC included
5         0.2 0.8           0.03281090 MHC included
6         0.1 0.8           0.03205260 MHC included
7        0.05 0.8           0.03139980 MHC included
8        0.05 0.6           0.03059800 MHC included
9         0.4 0.6           0.03008720 MHC included
10        0.1 0.6           0.02993970 MHC included
11        0.2 0.6           0.02959850 MHC included
12        0.6 0.6           0.02925860 MHC included
13        0.8 0.6           0.02882640 MHC included
14          1 0.6           0.02878650 MHC included
15       0.05 0.4           0.02680780 MHC included
16        0.1 0.4           0.02425730 MHC included
17        0.4 0.4           0.02326720 MHC included
18        0.2 0.4           0.02295200 MHC included
19        0.6 0.4           0.02235550 MHC included
20        0.8 0.4           0.02183290 MHC included
21          1 0.4           0.02181220 MHC included
22 0.00000005 0.8           0.02165780 MHC included
23 0.00000005 0.6           0.02134920 MHC included
24 0.00000005 0.4           0.02068890 MHC included
25 0.00000005 0.2           0.01867150 MHC included
26       0.05 0.2           0.01638110 MHC included
27        0.1 0.2           0.01358170 MHC included
28        0.4 0.2           0.01312180 MHC included
29        0.2 0.2           0.01280330 MHC included
30        0.6 0.2           0.01268200 MHC included
31       0.05 0.8           0.01245110 MHC excluded
32        0.8 0.2           0.01240420 MHC included
33          1 0.2           0.01237700 MHC included
34        0.4 0.8           0.01215340 MHC excluded
35        0.1 0.8           0.01169500 MHC excluded
36        0.6 0.8           0.01168870 MHC excluded
37        0.2 0.8           0.01161880 MHC excluded
38        0.8 0.8           0.01144460 MHC excluded
39          1 0.8           0.01142960 MHC excluded
40        0.4 0.6           0.01015730 MHC excluded
41        0.6 0.6           0.00973436 MHC excluded
42       0.05 0.6           0.00964942 MHC excluded
43        0.8 0.6           0.00949456 MHC excluded
44          1 0.6           0.00947578 MHC excluded
45        0.2 0.6           0.00905918 MHC excluded
46        0.1 0.6           0.00890722 MHC excluded
47        0.4 0.4           0.00797541 MHC excluded
48        0.6 0.4           0.00767666 MHC excluded
49        0.8 0.4           0.00742282 MHC excluded
50          1 0.4           0.00741685 MHC excluded
51       0.05 0.4           0.00730127 MHC excluded
52        0.2 0.4           0.00681032 MHC excluded
53        0.1 0.4           0.00656468 MHC excluded
54 0.00000005 0.8           0.00566601 MHC excluded
55        0.4 0.2           0.00543246 MHC excluded
56        0.6 0.2           0.00532709 MHC excluded
57        0.8 0.2           0.00518394 MHC excluded
58          1 0.2           0.00516490 MHC excluded
59       0.05 0.2           0.00475391 MHC excluded
60        0.2 0.2           0.00457437 MHC excluded
61 0.00000005 0.6           0.00454180 MHC excluded
62        0.1 0.2           0.00417019 MHC excluded
63 0.00000005 0.4           0.00352574 MHC excluded
64 0.00000005 0.2           0.00248654 MHC excluded
````

## 3. Validate the PRS in the MRI cohort
This step checks that these optimal PRS actually have some predictive/discriminative power to distinguish MS from controls in the MRI dataset. These results are a rough check as there are only 164 MS cases in the MRI cohort.


````R
# Even though there are only 164 MS cases in the MRI cohort, let's just check the PRS is working reasonably well.

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
print(nagelkerke(nohla_prs_model,null_model))

````
Clearly the PRS works well in discriminating MS cases from controls in the MRI dataset (or as well as can be expected for PRS):

Parameter  |  PRS (no HLA)   | PRS (with HLA)
--- | --- | ---
Nagelkerke's Pseudo-R2 | 0.0134 |  0.0394
Likelihood ratio P value | 2.54e-07 | 9.98e-19

Now we'll check that with increasing PRS deciles, MS risk increases as we would expect.
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
35991        |  163


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
summary(hla_model)$coefficients[10,]
summary(hla_model)$coefficients[10,]
    Estimate   Std. Error      t value     Pr(>|t|)
-0.005647618  0.004715386 -1.197700098  0.231042586

nohla_model = glm(data=mri_noms,norm_wm_lesions~`Age at recruitment.0.0`+Sex.0.0+`Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+`Volume of brain, grey+white matter (normalised for head size).2.0`+`Volume of ventricular cerebrospinal fluid (normalised for head size).2.0` +NOHLA_PRS)
summary(nohla_model)$coefficients[10,]

   Estimate  Std. Error     t value    Pr(>|t|)
0.004583730 0.004707508 0.973706143 0.330209763
````

So there is no association between either PRS and WM hyperintensities.

What about regional mean Fractional Anisotropy?
````R

fa_vars = mri_noms %>% select(EID,`Age at recruitment.0.0`,Sex.0.0,MS_status,contains("Mean FA"),HLA_PRS,NOHLA_PRS,`Genetic principal components.0.1`,`Genetic principal components.0.2`,`Genetic principal components.0.3`,`Genetic principal components.0.4`) %>%
 select(EID,`Age at recruitment.0.0`,Sex.0.0,MS_status,HLA_PRS,NOHLA_PRS,`Genetic principal components.0.1`,`Genetic principal components.0.2`,`Genetic principal components.0.3`,`Genetic principal components.0.4`,contains(".2.0")) %>%
 select(-contains("Weighted"))

# HLA PRS
overall_coefs = data.frame()
make_model = function(x){
  fa_vars_no_missing = fa_vars %>% filter(!is.na(fa_vars[[x]]))
  model = glm(data=fa_vars_no_missing, rankNorm(fa_vars_no_missing[[x]]) ~ `Age at recruitment.0.0`+ `Sex.0.0`+ `Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+ HLA_PRS)
  coefs = summary(model)$coefficients[8,]
  overall_coefs <<- bind_rows(overall_coefs,coefs)
}
sapply(colnames(fa_vars)[-c(1:10)],make_model)

overall_coefs = as.tbl(overall_coefs)
overall_coefs$variable = colnames(fa_vars)[-c(1:10)]

hla_fa_plot = ggplot(overall_coefs,aes(Estimate,variable))+
geom_point(aes(size=-log10(`Pr(>|t|)`)))+
geom_errorbarh(aes(y=variable,xmin=Estimate-1.96*`Std. Error`,xmax=Estimate+1.96*`Std. Error`),height=0.1)+
geom_vline(xintercept=0,alpha=0.5)+
theme_bw()+
labs(x="Effect of HLA PRS on normalised FA variable")

# NOHLA PRS
overall_coefs = data.frame()
make_model = function(x){
  fa_vars_no_missing = fa_vars %>% filter(!is.na(fa_vars[[x]]))
  model = glm(data=fa_vars_no_missing, rankNorm(fa_vars_no_missing[[x]]) ~ `Age at recruitment.0.0`+ `Sex.0.0`+ `Genetic principal components.0.1`+`Genetic principal components.0.2`+`Genetic principal components.0.3`+`Genetic principal components.0.4`+ NOHLA_PRS)
  coefs = summary(model)$coefficients[8,]
  overall_coefs <<- bind_rows(overall_coefs,coefs)
}
sapply(colnames(fa_vars)[-c(1:10)],make_model)

overall_coefs = as.tbl(overall_coefs)
overall_coefs$variable = colnames(fa_vars)[-c(1:10)]

nonhla_fa_plot = ggplot(overall_coefs,aes(Estimate,variable))+
geom_point(aes(size=-log10(`Pr(>|t|)`)))+
geom_errorbarh(aes(y=variable,xmin=Estimate-1.96*`Std. Error`,xmax=Estimate+1.96*`Std. Error`),height=0.1)+
geom_vline(xintercept=0,alpha=0.5)+
theme_bw()+
labs(x="Effect of Non-HLA PRS on normalised FA variable")

png("fa_prs_plot.png",res=300,units="in",height=8,width=16)
grid.arrange(hla_fa_plot,nonhla_fa_plot,nrow=1)
dev.off()
````
The answer is no. None of these associations pass the multiple testing threshold (Alpha = 0.05, Bonferroni correction).
![fa_plot](https://github.com/benjacobs123456/PRS_UKB_MRI/blob/main/fa_prs_plot.png)
