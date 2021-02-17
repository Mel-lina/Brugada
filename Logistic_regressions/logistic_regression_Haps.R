#!/usr/bin/env Rscript

# Goal: association analysis of haplotypes with Brugada syndrome using logistic regressions under a recessive, dominant and multiplicative model of inheritance.
# NOTE: Input file is a ".txt" file containing all cases and control individuals (rows) with their genotipes (columns), differently coded depending on the inheritance model to be tested.
    ## For the recessive model: genotypes must be coded as 0 or 1 depending on the number of alleles carried for the test haplotype (ex: no test haplotype or heterozygous=0, homozygous for test haplotype=1).
    ## For the dominant model: genotypes must be coded as 0 or 1 depending on the number of alleles carried for the test haplotype (ex: no test haplotype=0, heterozygous or homozygous for test haplotype=1).
    ## For the multiplicative model: genotypes must be coded as 0, 1, or 2 depending on the number of alleles carried for the test haplotype (ex: no test haplotype= 0, heterozygous=1, homozygous for test haplotype=2).

# NOTE: Remember to change output name for every model tested.

args = commandArgs(trailingOnly=TRUE)

# Test if there is at least one argument: if not, return an error with usage information
if (length(args)==0) {
  stop("At least one argument must be supplied (input file)\nUsage: Rscript logistic_regression_Haps.R <INPUT_FILE>\n", call.=FALSE)
}

input_file = args[1]

#Logistic regression for Haplotypes
library(epiDisplay)

mydata <- read.table(input_file, header=TRUE, sep="\t")
Hap1 <-as.numeric(mydata$Hap1)
Hap2 <-as.numeric(mydata$Hap2)
Hap3 <-as.numeric(mydata$Hap3)
Hap4 <-as.numeric(mydata$Hap4)
Hap5 <-as.numeric(mydata$Hap5)
Hap6 <-as.numeric(mydata$Hap6)
Hap7 <-as.numeric(mydata$Hap7)
Hap8 <-as.numeric(mydata$Hap8)
Hap9 <-as.numeric(mydata$Hap9)
Hap10 <-as.numeric(mydata$Hap10)
Hap11 <-as.numeric(mydata$Hap11)
Hap12 <-as.numeric(mydata$Hap12)
Hap13 <-as.numeric(mydata$Hap13)
Hap14 <-as.numeric(mydata$Hap14)
Hap15 <-as.numeric(mydata$Hap15)
gender <- as.factor(mydata$Gender)

df <- data.frame("Haplotype", "OR", "Pvalue", "BIC")
if ("Hap1" %in% names(mydata)) {
    mHap1 <- glm(mydata$Population ~ Hap1 + gender, family="binomial")
    results_mHap1 <- try(logistic.display(mHap1), silent=TRUE)
    if ("try-error" %in% class(results_mHap1)) {
        OR_mHap1 = "NA"
        pvalue_mHap1 = "NA"
        BIC_Hap1 <- BIC(mHap1)
         
    } else {
        OR_mHap1 <- results_mHap1$table[,2][1]
        summary_mHap1 <-summary(mHap1)
        pvalue_mHap1 <- summary_mHap1$coefficients[,4][2]
        BIC_Hap1 <- BIC(mHap1) 
    }
    linea <- c("Hap1", OR_mHap1, pvalue_mHap1, BIC_Hap1)
    df[nrow(df) +1,] = linea
}


if ("Hap2" %in% names(mydata)) {
    mHap2 <- glm(mydata$Population ~ Hap2 + gender, family="binomial")
    results_mHap2 <- try(logistic.display(mHap2), silent=TRUE)
    if ("try-error" %in% class(results_mHap2)) {
        OR_mHap2 = "NA"
        pvalue_mHap2 = "NA"
        BIC_Hap2 <- BIC(mHap2) 
         
    } else {
        OR_mHap2 <- results_mHap2$table[,2][1]
        summary_mHap2 <-summary(mHap2)
        pvalue_mHap2 <- summary_mHap2$coefficients[,4][2]
        BIC_Hap2 <- BIC(mHap2)        
    }
    linea <- c("Hap2", OR_mHap2, pvalue_mHap2, BIC_Hap2)
    df[nrow(df) +1,] = linea
}

if ("Hap3" %in% names(mydata)) {
    mHap3 <- glm(mydata$Population ~ Hap3 + gender, family="binomial")
    results_mHap3 <- try(logistic.display(mHap3), silent=TRUE)
    if ("try-error" %in% class(results_mHap3)) {
        OR_mHap3 = "NA"
        pvalue_mHap3 = "NA"
        BIC_Hap3 <- BIC(mHap3) 
         
    } else {
        OR_mHap3 <- results_mHap3$table[,2][1]
        summary_mHap3 <-summary(mHap3)
        pvalue_mHap3 <- summary_mHap3$coefficients[,4][2]
        BIC_Hap3 <- BIC(mHap3)   
    }
    linea <- c("Hap3", OR_mHap3, pvalue_mHap3, BIC_Hap3)
    df[nrow(df) +1,] = linea    
}

if ("Hap4" %in% names(mydata)) {
    mHap4 <- glm(mydata$Population ~ Hap4 + gender, family="binomial")
    results_mHap4 <- try(logistic.display(mHap4), silent=TRUE)
    if ("try-error" %in% class(results_mHap4)) {
        OR_mHap4 = "NA"
        pvalue_mHap4 = "NA"
        BIC_Hap4 <- BIC(mHap4) 
         
    } else {
        OR_mHap4 <- results_mHap4$table[,2][1]
        summary_mHap4 <-summary(mHap4)
        pvalue_mHap4 <- summary_mHap4$coefficients[,4][2]
        BIC_Hap4 <- BIC(mHap4) 
           
    }
    linea <- c("Hap4", OR_mHap4, pvalue_mHap4, BIC_Hap4)
    df[nrow(df) +1,] = linea  
}

if ("Hap5" %in% names(mydata)) {
    mHap5 <- glm(mydata$Population ~ Hap5 + gender, family="binomial")
    results_mHap5 <- try(logistic.display(mHap5), silent=TRUE)
    if ("try-error" %in% class(results_mHap5)) {
        OR_mHap5 = "NA"
        pvalue_mHap5 = "NA"
        BIC_Hap5 <- BIC(mHap5) 
          
    } else {
        OR_mHap5 <- results_mHap5$table[,2][1]
        summary_mHap5 <-summary(mHap5)
        pvalue_mHap5 <- summary_mHap5$coefficients[,4][2]
        BIC_Hap5 <- BIC(mHap5) 
           
    }
    linea <- c("Hap5", OR_mHap5, pvalue_mHap5, BIC_Hap5)
    df[nrow(df) +1,] = linea      
}

if ("Hap6" %in% names(mydata)) {
    mHap6 <- glm(mydata$Population ~ Hap6 + gender, family="binomial")
    results_mHap6 <- try(logistic.display(mHap6), silent=TRUE)
    if ("try-error" %in% class(results_mHap6)) {
        OR_mHap6 = "NA"
        pvalue_mHap6 = "NA"
        BIC_Hap6 <- BIC(mHap6) 
         
    } else {
        OR_mHap6 <- results_mHap6$table[,2][1]
        summary_mHap6 <-summary(mHap6)
        pvalue_mHap6 <- summary_mHap6$coefficients[,4][2]
        BIC_Hap6 <- BIC(mHap6)     
    }
    linea <- c("Hap6", OR_mHap6, pvalue_mHap6, BIC_Hap6)
    df[nrow(df) +1,] = linea      
}

if ("Hap7" %in% names(mydata)) {
    mHap7 <- glm(mydata$Population ~ Hap7 + gender, family="binomial")
    results_mHap7 <- try(logistic.display(mHap7), silent=TRUE)
    if ("try-error" %in% class(results_mHap7)) {
        OR_mHap7 = "NA"
        pvalue_mHap7 = "NA"
        BIC_Hap7 <- BIC(mHap7) 
         
    } else {
        OR_mHap7 <- results_mHap7$table[,2][1]
        summary_mHap7 <-summary(mHap7)
        pvalue_mHap7 <- summary_mHap7$coefficients[,4][2]
        BIC_Hap7 <- BIC(mHap7) 
    }
    linea <- c("Hap7", OR_mHap7, pvalue_mHap7, BIC_Hap7)
    df[nrow(df) +1,] = linea    
}

if ("Hap8" %in% names(mydata)) {
    mHap8 <- glm(mydata$Population ~ Hap8 + gender, family="binomial")
    results_mHap8 <- try(logistic.display(mHap8), silent=TRUE)
    if ("try-error" %in% class(results_mHap8)) {
        OR_mHap8 = "NA"
        pvalue_mHap8 = "NA"
        BIC_Hap8 <- BIC(mHap8) 
         
    } else {
        OR_mHap8 <- results_mHap8$table[,2][1]
        summary_mHap8 <-summary(mHap8)
        pvalue_mHap8 <- summary_mHap8$coefficients[,4][2]
        BIC_Hap8 <- BIC(mHap8)  
    }
    linea <- c("Hap8", OR_mHap8, pvalue_mHap8, BIC_Hap8)
    df[nrow(df) +1,] = linea  
}

if ("Hap9" %in% names(mydata)) {
mHap9 <- glm(mydata$Population ~ Hap9 + gender, family="binomial")
results_mHap9 <- try(logistic.display(mHap9), silent=TRUE)
    if ("try-error" %in% class(results_mHap9)) {
        OR_mHap9 = "NA"
        pvalue_mHap9 = "NA"
        BIC_Hap9 <- BIC(mHap9) 
         
    } else {
        OR_mHap9 <- results_mHap9$table[,2][1]
        summary_mHap9 <-summary(mHap9)
        pvalue_mHap9 <- summary_mHap9$coefficients[,4][2]
        BIC_Hap9 <- BIC(mHap9)            
    }
    linea <- c("Hap9", OR_mHap9, pvalue_mHap9, BIC_Hap9)
    df[nrow(df) +1,] = linea  
}

if ("Hap10" %in% names(mydata)) {
    mHap10 <- glm(mydata$Population ~ Hap10 + gender, family="binomial")
    results_mHap10 <- try(logistic.display(mHap10), silent=TRUE)
    if ("try-error" %in% class(results_mHap10)) {
        OR_mHap10 = "NA"
        pvalue_mHap10 = "NA"
        BIC_Hap10 <- BIC(mHap10) 
         
    } else {
        OR_mHap10 <- results_mHap10$table[,2][1]
        summary_mHap10 <-summary(mHap10)
        pvalue_mHap10 <- summary_mHap10$coefficients[,4][2]
        BIC_Hap10 <- BIC(mHap10) 
    }
    linea <- c("Hap10", OR_mHap10, pvalue_mHap10, BIC_Hap10)
    df[nrow(df) +1,] = linea  
}

if ("Hap11" %in% names(mydata)) {
    mHap11 <- glm(mydata$Population ~ Hap11 + gender, family="binomial")
    results_mHap11 <- try(logistic.display(mHap11), silent=TRUE)
    if ("try-error" %in% class(results_mHap11)) {
        OR_mHap11 = "NA"
        pvalue_mHap11 = "NA"
        BIC_Hap11 <- BIC(mHap11) 
         
    } else {
        OR_mHap11 <- results_mHap11$table[,2][1]
        summary_mHap11 <-summary(mHap11)
        pvalue_mHap11 <- summary_mHap11$coefficients[,4][2]
        BIC_Hap11 <- BIC(mHap11) 
    }
    linea <- c("Hap11", OR_mHap11, pvalue_mHap11, BIC_Hap11)
    df[nrow(df) +1,] = linea  
}

if ("Hap12" %in% names(mydata)) {
    mHap12 <- glm(mydata$Population ~ Hap12 + gender, family="binomial")
    results_mHap12 <- try(logistic.display(mHap12), silent=TRUE)
    if ("try-error" %in% class(results_mHap12)) {
        OR_mHap12 = "NA"
        pvalue_mHap12 = "NA"
        BIC_Hap12 <- BIC(mHap12) 
            
    } else {
        OR_mHap12 <- results_mHap12$table[,2][1]
        summary_mHap12-summary(mHap12)
        pvalue_mHap12 <- summary_mHap12$coefficients[,4][2]
        BIC_Hap12 <- BIC(mHap12)         
    }
    linea <- c("Hap12", OR_mHap12, pvalue_mHap12, BIC_Hap12)
    df[nrow(df) +1,] = linea  
}

if ("Hap13" %in% names(mydata)) {
    mHap13 <- glm(mydata$Population ~ Hap13 + gender, family="binomial")
    results_mHap13 <- try(logistic.display(mHap13), silent=TRUE)
    if ("try-error" %in% class(results_mHap13)) {
        OR_mHap13 = "NA"
        pvalue_mHap13 = "NA"
        BIC_Hap13 <- BIC(mHap13) 
         
    } else {
        OR_mHap13 <- results_mHap13$table[,2][1]
        summary_mHap13 <-summary(mHap13)
        pvalue_mHap13 <- summary_mHap13$coefficients[,4][2]
        BIC_Hap13 <- BIC(mHap13)       
    }
    linea <- c("Hap13", OR_mHap13, pvalue_mHap13, BIC_Hap13)
    df[nrow(df) +1,] = linea  
}

if ("Hap14" %in% names(mydata)) {
    mHap14 <- glm(mydata$Population ~ Hap14 + gender, family="binomial")
    results_mHap14 <- try(logistic.display(mHap14), silent=TRUE)
    if ("try-error" %in% class(results_mHap14)) {
        OR_mHap14 = "NA"
        pvalue_mHap14 = "NA"
        BIC_Hap14 <- BIC(mHap14) 
         
    } else {
        OR_mHap14 <- results_mHap14$table[,2][1]
        summary_mHap14 <-summary(mHap14)
        pvalue_mHap14 <- summary_mHap14$coefficients[,4][2]
        BIC_Hap14 <- BIC(mHap14)  
    }
    linea <- c("Hap14", OR_mHap14, pvalue_mHap14, BIC_Hap14)
    df[nrow(df) +1,] = linea  
}

if ("Hap15" %in% names(mydata)) {
    mHap15 <- glm(mydata$Population ~ Hap15 + gender, family="binomial")
    results_mHap15 <- try(logistic.display(mHap15), silent=TRUE)
    if ("try-error" %in% class(results_mHap15)) {
        OR_mHap15 = "NA"
        pvalue_mHap15 = "NA"
        BIC_Hap15 <- BIC(mHap15) 
         
    } else {
        OR_mHap15 <- results_mHap15$table[,2][1]
        summary_mHap15 <-summary(mHap15)
        pvalue_mHap15 <- summary_mHap15$coefficients[,4][2]
        BIC_Hap15 <- BIC(mHap15)    
    }
    linea <- c("Hap15", OR_mHap15, pvalue_mHap15, BIC_Hap15)
    df[nrow(df) +1,] = linea  
}  
write.table(df, "logistic_regression_Haps∫.txt", sep="\t") #Remember to change output name for every model tested.
