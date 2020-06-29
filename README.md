# RBD_Clinical_Analyses
Clinical analysis for the RBD GWAS - June 2020

## Preparing files
* Age at onset (AAO) and age at diagnoses (AADx) provided by clinicians from each centre. If missing, sample is not included in any age analyses.    
* Time to converstion (TTC) is calculated for both AAO and AADx. If not converted, time since onset or diagnosis is used.  
* If conversion has not been indicates YES or NO, sample is considered "unknown" and not used in conversion analyses. 
* Date of last follow up used as indicated unless YEAR OF FOLLOW UP - YEAR OF BIRTH - YEAR OF RBD ONSET/DIAGNOSIS was negative. In this case, 2018 was used (last collection of clinical data).  
* Genotypes were collected from imputed GWAS data - only common variants were included (+ known variant N370S).  

## Clinical analysis (on cases only). 
**Load libraries & data.**  
Data includes two files:  
* data = complete covariate csv created as described above
* variants = csv containing only the GWAS significant variants *in the same sample order as the covariate csv.*  
```R
library(devtools)
install_github("dgrtwo/broom")
library(broom)

data = read.csv("cases_analysesInR.csv", header = T)
variants = read.csv("VARIANTSfor_cases_analysesInR.csv")

attach(data)
````
**AAO Summary:**
```R
min(AAO, na.rm = T) # 18
max(AAO, na.rm = T) # 85
mean(AAO, na.rm = T) # 60.3
sd(AAO, na.rm = T) # 10.11
```
**Test GWAS significant variants on AGE AT ONSET.**  
Using covariates sex and PC1-5:
```R
aao_fun = function(x) {
  fit = glm(data$AAO ~ x + data$sex + data$PC1 + data$PC2 + data$PC3 + data$PC4 + data$PC5)
  return(fit) 
}

z = matrix(, nrow = ncol(variants), ncol = 5, dimnames = list(colnames(variants), 
                                                              c("variant", "estimate","se","statistic","p")))

for (i in 1:ncol(variants))
{ f = aao_fun(variants[,i])
  z[i,1] = coef(summary(f))[2,1]
  z[i,2] = coef(summary(f))[2,2]
  z[i,3] = coef(summary(f))[2,3]
  z[i,4] = coef(summary(f))[2,4]
  z[i,5] = "" }

write.table(z, file="AAO_regression.txt", col.names=T, row.names=T, sep="\t", quote=F)
````



