####Multicolinearity####
#initialize and load the libraries
library(gplots)
library(RColorBrewer)
library(pheatmap)

#reading the csv file with the data
df<- read.csv("/Users/..." ,header=TRUE)

# remove certain unwanted columns
dff1 <- df[colMeans(is.na(df)) <= 0.1 & colMeans((df == 0), na.rm = T) <= 0.1]
dff1<-df

#replacing NA values with '0'
dff1[is.na(dff1)] <- 0

# correlation analysis
library(ggplot2)
library(reshape2)
library(corrplot)
library(Hmisc)

matriz_cor <-cor(data.matrix(dff1[,1:ncol(dff1)]))
write.csv(matriz_cor, file="/Users/...")

for (i in 1:nrow(matriz_cor)){
  correlations <-  which((abs(matriz_cor[i,]) > 0.5) & (matriz_cor[i,] != 1))
  
  if(length(correlations)> 0){
    print(colnames(dff1)[i+1])
    print(correlations)
  }
}

corrplot(matriz_cor, method="circle", type = "lower", tl.cex = .8, tl.col = "black", order="hclust")





####Pec Survival####
##Load Dataset##
MESA <- data.frame(read.csv(file = '/Users/...'))

##Impute Median Values##
detach("package:randomForestSRC", unload = TRUE)
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
dat <- as.data.frame(imputed$data)

#Load Packages
library(survival)
library(pec)
library(prodlim)

#COXPH Regression
coxR <- coxph(Surv(duration,event)~Race,data=dat,x=TRUE,y=TRUE)

coxA <- coxph(Surv(duration,event)~ASCVD,data=dat,x=TRUE,y=TRUE)

coxB <- coxph(Surv(duration,event)~Age+Race+Gender+Diabetes+Hypertension+Hyperlipidemia+
                Statin+Smoking+Metabolic_Syndrome+Mean_SBP+Mean_DBP+LDL+HDL+Chol+Trig+
                Alb_Creat_Ratio+Pack_Years+BMI+Walking_Min_Wk+Education+Income+FH_Heart
              ,data=dat,x=TRUE,y=TRUE)

coxBI <- coxph(Surv(duration,event)~Age+Race+Gender+Diabetes+Hypertension+Hyperlipidemia+
                 Statin+Smoking+Metabolic_Syndrome+Mean_SBP+Mean_DBP+LDL+HDL+Chol+Trig+
                 Alb_Creat_Ratio+Pack_Years+BMI+Walking_Min_Wk+Education+Income+FH_Heart+
                 Homocysteine+IL6+Plasmin_Antiplasmin+Fibrinogen_Antigen+C_Reactive_Protein+
                 D_Dimer+Factor_VIII,data=dat,x=TRUE,y=TRUE)

coxBCT <- coxph(Surv(duration,event)~Age+Race+Gender+Diabetes+Hypertension+Hyperlipidemia+
                  Statin+Smoking+Metabolic_Syndrome+Mean_SBP+Mean_DBP+LDL+HDL+Chol+Trig+
                  Alb_Creat_Ratio+Pack_Years+BMI+Walking_Min_Wk+Education+Income+FH_Heart+
                  LV_EF+Left_Ventricular_Area+Pericardial_Fat+Coronary_Calcium
                ,data=dat,x=TRUE,y=TRUE)

coxALL <- coxph(Surv(duration,event)~Age+Race+Gender+Diabetes+Hypertension+Hyperlipidemia+
                  Statin+Smoking+Metabolic_Syndrome+Mean_SBP+Mean_DBP+LDL+HDL+Chol+Trig+
                  Alb_Creat_Ratio+Pack_Years+BMI+Walking_Min_Wk+Education+Income+FH_Heart+
                  Homocysteine+IL6+Plasmin_Antiplasmin+Fibrinogen_Antigen+
                  C_Reactive_Protein+D_Dimer+Factor_VIII+Left_Ventricular_Area+
                  LV_EF+Pericardial_Fat+Coronary_Calcium,data=dat,x=TRUE,y=TRUE)

#Concordance Index#
ApparrentCindex <- pec::cindex(list("COXPH Race"=coxR,
                                     "COXPH ASCVD"=coxA,
                                     "COXPH Baseline"=coxB,
                                     "COXPH Baseline + Inflammatory"=coxBI,
                                     "COXPH Baseline + CT"=coxBCT,
                                     "COXPH All Features"=coxALL),
                                formula=Surv(duration,event)~Age+Race+Gender+Diabetes+Hypertension+Hyperlipidemia+
                                  Statin+Smoking+Metabolic_Syndrome+Mean_SBP+Mean_DBP+LDL+HDL+Chol+Trig+
                                  Alb_Creat_Ratio+Pack_Years+BMI+Walking_Min_Wk+Education+Income+FH_Heart+
                                  Homocysteine+IL6+Plasmin_Antiplasmin+Fibrinogen_Antigen+
                                  C_Reactive_Protein+D_Dimer+Factor_VIII+Left_Ventricular_Area+
                                  LV_EF+Pericardial_Fat+Coronary_Calcium,data=dat,
                                eval.times=seq(0,6000,1), pred.times=seq(0,6000,1))

plot(ApparrentCindex, legend = FALSE)
write.csv(ApparrentCindex$AppCindex, file = "/Users/...")

#Concordance Index with bootstrapping#
set.seed(100)
bcvCindex <- pec::cindex(list("COXPH Race"=coxR,
                              "COXPH ASCVD"=coxA,
                              "COXPH Baseline"=coxB,
                              "COXPH Baseline + Inflammatory"=coxBI,
                              "COXPH Baseline + CT"=coxBCT,
                              "COXPH All Features"=coxALL),
                         formula=Surv(duration,event)~Age+Race+Gender+Diabetes+Hypertension+Hyperlipidemia+
                           Statin+Smoking+Metabolic_Syndrome+Mean_SBP+Mean_DBP+LDL+HDL+Chol+Trig+
                           Alb_Creat_Ratio+Pack_Years+BMI+Walking_Min_Wk+Education+Income+FH_Heart+
                           Homocysteine+IL6+Plasmin_Antiplasmin+Fibrinogen_Antigen+
                           C_Reactive_Protein+D_Dimer+Factor_VIII+Left_Ventricular_Area+
                           LV_EF+Pericardial_Fat+Coronary_Calcium,data=dat,splitMethod="bootcv",B=5,
                         eval.times=seq(0,6000,1), pred.times=seq(0,6000,1))

plot(bcvCindex, legend = FALSE)
write.csv(bcvCindex$AppCindex, file = "/Users/...")





####COXPH Regression####
##Load Packages##
library(tidyverse)
library(reshape2)
library(xgboost)
library(randomForest)
library(rfUtilities)
library(caret)
library(survival)
library(survminer)

##Load Dataset##
MESA <- data.frame(read.csv(file = '/Users/...'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

#Multivariate Analysis#
res.cox <- coxph(Surv(duration, event) ~ Age+Race+Gender+Diabetes+Hypertension+Hyperlipidemia+
                   Statin+Smoking+Metabolic_Syndrome+Mean_SBP+Mean_DBP+LDL+HDL+Chol+Trig+
                   Alb_Creat_Ratio+Pack_Years+BMI+Walking_Min_Wk+Education+Income+FH_Heart+
                   Homocysteine+IL6+Plasmin_Antiplasmin+Fibrinogen_Antigen+
                   C_Reactive_Protein+D_Dimer+Factor_VIII+LV_EF+Left_Ventricular_Area+
                   Pericardial_Fat+Coronary_Calcium,
                 data = final)
summary(res.cox)

#Hazard Ratios#
surv_object <- Surv(time = final$duration, event = final$event)
fit.coxph <- coxph(surv_object ~ Age+Race+Gender+Diabetes+Hypertension+Hyperlipidemia+
                     Statin+Smoking+Metabolic_Syndrome+Mean_SBP+Mean_DBP+LDL+HDL+Chol+Trig+
                     Alb_Creat_Ratio+Pack_Years+BMI+Walking_Min_Wk+Education+Income+FH_Heart+
                     Homocysteine+IL6+Plasmin_Antiplasmin+Fibrinogen_Antigen+
                     C_Reactive_Protein+D_Dimer+Factor_VIII+LV_EF+Left_Ventricular_Area+
                     Pericardial_Fat+Coronary_Calcium, 
                   data = final)
ggforest(fit.coxph, data = final)





####Non-Categorical NRI####
library(survival)
library(nricens)

##Load Dataset##
MESA <- data.frame(read.csv(
  file = '/Users/...'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)
final$event=as.integer(final$event)
final$duration=as.integer(final$duration)

# detach mlr package because of 'plotCalibration' conflict
detach("package:mlr", unload = TRUE)

## predciting the event of 'death'
time  = final$duration
event = final$event

z.std <- as.matrix(final[1])
z.new <- as.matrix(final[2])

#5475 days = 15 year cutoff
nricens(time = time, event = event, z.std = z.std, z.new = z.new, t0 = 5475,  updown = "diff", cut = 0, point.method = "ipw",
        niter = 50, msg = TRUE)
