# RBD_Clinical_Analyses
Clinical analysis for the RBD GWAS - June 2020

## Preparing files
* Age at onset (AAO) and age at diagnoses (AADx) provided by clinicians from each centre. If missing, sample is not included in any age analyses.    
* Time to converstion (TTC) is calculated for both AAO and AADx. If not converted, time since onset or diagnosis is used.  
* If conversion has not been indicates YES or NO, sample is considered "unknown" and not used in conversion analyses. 
* Date of last follow up used as indicated unless YEAR OF FOLLOW UP - YEAR OF BIRTH - YEAR OF RBD ONSET/DIAGNOSIS was negative. In this case, 2018 was used (last collection of clinical data).  
* Genotypes were collected from imputed GWAS data - only common variants were included (+ known variant N370S).  

## Effect of variants on AAO/AADx 
**Load libraries & data.**  
Data includes two files:  
* covar = complete covariate csv created as described above
* variants = csv containing only the GWAS significant variants *in the same sample order as the covariate csv.*  
```R
library(devtools)
install_github("dgrtwo/broom")
library(broom)

covar = read.csv("cases_analysesInR.csv", header = T)
variants = read.csv("VARIANTSfor_cases_analysesInR.csv")

attach(covar)
````
**AAO Summary:**
```R
min(AAO, na.rm = T) # 18
max(AAO, na.rm = T) # 85
mean(AAO, na.rm = T) # 60.3
sd(AAO, na.rm = T) # 10.11

has_aao = subset(covar, AAO != "NA")
nrow(has_aao) # 699
````
**AADx Summary:**
```R
min(AADx, na.rm = T) # 28
max(AADx, na.rm = T) # 89
mean(AADx, na.rm = T) # 65.38
sd(AADx, na.rm = T) # 8.49

has_aadx = subset(covar, AADx != "NA")
nrow(has_aadx) # 660
```
**Test GWAS significant variants on AGE AT ONSET.**  
Using covariates sex and PC1-5:
```R
# create a function to run regression
aao_fun = function(x) {
  fit = glm(covar$AAO ~ x + covar$sex + covar$PC1 + covar$PC2 + covar$PC3 + covar$PC4 + covar$PC5)
  return(fit) 
}

# create an empty matrix to put results
z = matrix(, nrow = ncol(variants), ncol = 5, dimnames = list(colnames(variants), 
                                                              c("variant", "estimate","se","statistic","p")))

# run regression on all variants
for (i in 1:ncol(variants))
{ f = aao_fun(variants[,i])
  z[i,1] = coef(summary(f))[2,1]
  z[i,2] = coef(summary(f))[2,2]
  z[i,3] = coef(summary(f))[2,3]
  z[i,4] = coef(summary(f))[2,4]
  z[i,5] = ""}

# save results
write.table(z, file="AAO_regression.txt", col.names=T, row.names=T, sep="\t", quote=F)
````

## Effect of variants on rate and type of conversion  
* Time to conversion is calculated as the difference between AADx of RBD and AADx of overt neurodegeneration.  
* All who converted to *any* overt neurodegeneration are included (summary of types below).  
* For those who have not converted, time to conversion is the difference between their AADx RBD and age at last follow-up. *If last follow-up is ambiguous, their age in 2018 was used (year of last collection of clinical data).  




