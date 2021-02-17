#!/usr/bin/env Rscript

# Goal: association analysis of SNVs with Brugada syndrome using logistic regressions under a recessive, dominant and multiplicative model of inheritance.
# NOTE: Input file is a ".txt" file containing all cases and control individuals (rows) with their genotipes (columns), differently coded depending on the inheritance model to be tested.
    ## For the recessive model: genotypes must be coded as 0 or 1 depending on the number of alleles carried for the test SNV allele (ex: no test allele or heterozygous=0, homozygous for test allele=1).
    ## For the dominant model: genotypes must be coded as 0 or 1 depending on the number of alleles carried for the test SNV allele (ex: no test allele=0, heterozygous or homozygous for test allele=1).
    ## For the multiplicative model: genotypes must be coded as 0, 1, or 2 depending on the number of alleles carried for the test SNV allele (ex: no test allele=0, heterozygous=1, homozygous for test allele=2).

# NOTE: Remember to change output name for every model tested.

args = commandArgs(trailingOnly=TRUE)

# Test if there is at least one argument: if not, return an error with usage information
if (length(args)==0) {
  stop("At least one argument must be supplied (input file)\nUsage: Rscript logistic_regression_SNVs.R <INPUT_FILE>\n", call.=FALSE)
}

input_file = args[1]

#Logistic regression for Haplotypes
library(epiDisplay)

mydata <- read.table(input_file, header=TRUE, sep="\t")
rs6801957 <-as.numeric(mydata$rs6801957)
rs6799257 <-as.numeric(mydata$rs6799257)
rs9836859 <-as.numeric(mydata$rs9836859)
rs6790396 <-as.numeric(mydata$rs6790396)
rs9874633 <-as.numeric(mydata$rs9874633)
rs10428132 <-as.numeric(mydata$rs10428132)
rs10428168 <-as.numeric(mydata$rs10428168)

gender <- as.factor(mydata$Gender)

df <- data.frame("SNV", "OR", "Pvalue", "BIC")
if ("rs6801957" %in% names(mydata)) {
    mrs6801957 <- glm(mydata$Population ~ rs6801957 + gender, family="binomial")
    results_mrs6801957 <- try(logistic.display(mrs6801957), silent=TRUE)
    if ("try-error" %in% class(results_mrs6801957)) {
        OR_mrs6801957 = "NA"
        pvalue_mrs6801957 = "NA"
        BIC_rs6801957 <- BIC(mrs6801957)
         
    } else {
        OR_mrs6801957 <- results_mrs6801957$table[,2][1]
        summary_mrs6801957 <-summary(mrs6801957)
        pvalue_mrs6801957 <- summary_mrs6801957$coefficients[,4][2]
        BIC_rs6801957 <- BIC(mrs6801957) 
    }
    linea <- c("6801957", OR_mrs6801957, pvalue_mrs6801957, BIC_rs6801957)
    df[nrow(df) +1,] = linea
}


if ("rs6799257" %in% names(mydata)) {
    mrs6799257 <- glm(mydata$Population ~ rs6799257 + gender, family="binomial")
    results_mrs6799257 <- try(logistic.display(mrs6799257), silent=TRUE)
    if ("try-error" %in% class(results_mrs6799257)) {
        OR_mrs6799257 = "NA"
        pvalue_mrs6799257 = "NA"
        BIC_rs6799257 <- BIC(mrs6799257) 
         
    } else {
        OR_mrs6799257 <- results_mrs6799257$table[,2][1]
        summary_mrs6799257 <-summary(mrs6799257)
        pvalue_mrs6799257 <- summary_mrs6799257$coefficients[,4][2]
        BIC_rs6799257 <- BIC(mrs6799257)        
    }
    linea <- c("rs6799257", OR_mrs6799257, pvalue_mrs6799257, BIC_rs6799257)
    df[nrow(df) +1,] = linea
}

if ("rs9836859" %in% names(mydata)) {
    mrs9836859 <- glm(mydata$Population ~ rs9836859 + gender, family="binomial")
    results_mrs9836859 <- try(logistic.display(mrs9836859), silent=TRUE)
    if ("try-error" %in% class(results_mrs9836859)) {
        OR_mrs9836859 = "NA"
        pvalue_mrs9836859 = "NA"
        BIC_rs9836859 <- BIC(mrs9836859) 
         
    } else {
        OR_mrs9836859 <- results_mrs9836859$table[,2][1]
        summary_mrs9836859 <-summary(mrs9836859)
        pvalue_mrs9836859 <- summary_mrs9836859$coefficients[,4][2]
        BIC_rs9836859 <- BIC(mrs9836859)   
    }
    linea <- c("rs9836859", OR_mrs9836859, pvalue_mrs9836859, BIC_rs9836859)
    df[nrow(df) +1,] = linea    
}

if ("rs6790396" %in% names(mydata)) {
    mrs6790396 <- glm(mydata$Population ~ rs6790396 + gender, family="binomial")
    results_mrs6790396 <- try(logistic.display(mrs6790396), silent=TRUE)
    if ("try-error" %in% class(results_mrs6790396)) {
        OR_mrs6790396 = "NA"
        pvalue_mrs6790396 = "NA"
        BIC_rs6790396 <- BIC(mrs6790396) 
         
    } else {
        OR_mrs6790396 <- results_mrs6790396$table[,2][1]
        summary_mrs6790396 <-summary(mrs6790396)
        pvalue_mrs6790396 <- summary_mrs6790396$coefficients[,4][2]
        BIC_rs6790396 <- BIC(mrs6790396) 
           
    }
    linea <- c("rs6790396", OR_mrs6790396, pvalue_mrs6790396, BIC_rs6790396)
    df[nrow(df) +1,] = linea  
}

if ("rs9874633" %in% names(mydata)) {
    mrs9874633 <- glm(mydata$Population ~ rs9874633 + gender, family="binomial")
    results_mrs9874633 <- try(logistic.display(mrs9874633), silent=TRUE)
    if ("try-error" %in% class(results_mrs9874633)) {
        OR_mrs9874633 = "NA"
        pvalue_mrs9874633 = "NA"
        BIC_rs9874633 <- BIC(mrs9874633) 
          
    } else {
        OR_mrs9874633 <- results_mrs9874633$table[,2][1]
        summary_mrs9874633 <-summary(mrs9874633)
        pvalue_mrs9874633 <- summary_mrs9874633$coefficients[,4][2]
        BIC_rs9874633 <- BIC(mrs9874633) 
           
    }
    linea <- c("rs9874633", OR_mrs9874633, pvalue_mrs9874633, BIC_rs9874633)
    df[nrow(df) +1,] = linea      
}

if ("rs10428132" %in% names(mydata)) {
    mrs10428132 <- glm(mydata$Population ~ rs10428132 + gender, family="binomial")
    results_mrs10428132 <- try(logistic.display(mrs10428132), silent=TRUE)
    if ("try-error" %in% class(results_mrs10428132)) {
        OR_mrs10428132 = "NA"
        pvalue_mrs10428132 = "NA"
        BIC_rs10428132 <- BIC(mrs10428132) 
         
    } else {
        OR_mrs10428132 <- results_mrs10428132$table[,2][1]
        summary_mrs10428132 <-summary(mrs10428132)
        pvalue_mrs10428132 <- summary_mrs104281326$coefficients[,4][2]
        BIC_rs10428132 <- BIC(mrs10428132)     
    }
    linea <- c("rs10428132", OR_mrs10428132, pvalue_mrs10428132, BIC_rs10428132)
    df[nrow(df) +1,] = linea      
}

if ("rs10428168" %in% names(mydata)) {
    mrs10428168 <- glm(mydata$Population ~ rs10428168 + gender, family="binomial")
    results_mrs10428168 <- try(logistic.display(mrs10428168), silent=TRUE)
    if ("try-error" %in% class(results_mrs10428168)) {
        OR_mrs10428168 = "NA"
        pvalue_mrs10428168 = "NA"
        BIC_rs10428168 <- BIC(mrs10428168) 
         
    } else {
        OR_mrs10428168 <- results_mrs10428168$table[,2][1]
        summary_mrs10428168 <-summary(mrs10428168)
        pvalue_mrs10428168 <- summary_mrs10428168$coefficients[,4][2]
        BIC_rs10428168 <- BIC(mrs10428168) 
    }
    linea <- c("rs10428168", OR_mrs10428168, pvalue_mrs10428168, BIC_rs10428168)
    df[nrow(df) +1,] = linea    
}

write.table(df, "logistic_regression_SNVs.txt", sep="\t") #Remember to change output name for every model tested.