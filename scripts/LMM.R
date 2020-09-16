## Regression Problem: Data has LACTO GRP (clearance) as reference and we are looking at cytokines of those who remained CST3/4 (persistent) and 
## CST3/4  for again at 3rd visit ( recurrent grp). Control for age, sti, psa and genital exam. If you can also do False descovery rate (fdr) on those 
## significant


## Load required packages
library(readxl) 
library("nnet")

# Define cytokine column names
allCytokineColnames <- c("IL.1a","IL.1b","IL.6","IL.12p40","IL.12p70","IL.18", "MIF","TNF.a","TNF.b","TRAIL","IL.8","IL.16","CTACK","Eotaxin","IP.10","IFN.a2","GROa","MCP.1","MCP.3","MIG","MIP.1a", "MIP.1b","RANTES","IL.3","IL.7",  "IL.9","b.NGF","FGF.basic","G.CSF","GM.CSF",  "HGF", "LIF","M.CSF", "PDGF.bb", "SCF","SCGF.b", "SDF.1a", "VEGF","IL.2Ra", "IL.2",  "IL.4", "IL.5","IL.13","IL.15","IL.17","IFN.g","IL.1ra", "IL.10")

# Import data
LinearModelData <- read_excel("metadata/LinearModelData.xlsx")

# Explore dataset
View(LinearModelData)
names(LinearModelData)
str(LinearModelData)

# Redefine group var
LinearModelData$Group <- ifelse(LinearModelData$Group == "Lacto group", "LG", ifelse(LinearModelData$Group == "BV group", "BG", "BA"))

# Check age variable
boxplot(LinearModelData$Age, main="Age", sub=paste("Outlier rows: ", boxplot.stats(LinearModelData$Age)$out)) 

## Correlation to establish relationship
# A low correlation (-0.2 < x < 0.2) probably suggests that much of variation of the response variable (Y) is unexplained by the predictor (X). In 
# that case, you should probably look for better explanatory variables.

# Define factors
LinearModelData$Group <- factor(LinearModelData$Group)
levels(LinearModelData$Group)

## Specify reference group of MLR
LinearModelData$Group <- relevel(LinearModelData$Group, ref = "LG") 
levels(LinearModelData$Group)

head(LinearModelData$PSA)
head(LinearModelData$Any_STI)
head(LinearModelData$`0_Genital_exam`)
head(LinearModelData$Age)

LinearModelData$PSA <- factor(LinearModelData$PSA)
LinearModelData$Any_STI <- factor(LinearModelData$Any_STI)
LinearModelData$`0_Genital_exam` <- factor(LinearModelData$`0_Genital_exam`)

str(LinearModelData)


# Using multinomial logistic regression as response variable is categorical
MLR <- multinom(Group ~ b.NGF, data = LinearModelData)
summary(MLR)

## Calculate Z score and p-Value for the variables in the model.
z <- summary(MLR)$coefficients/summary(MLR)$standard.errors
z

p <- (1 - pnorm(abs(z), 0, 1))*2
p

#https://stats.stackexchange.com/questions/157875/logistic-regression-coefficients-and-expcoefficients-meaning-and-relationship
exp(coef(MLR))

regCytokines <- function(ModelData, ListOfCytokineNames){
  for(cname in ListOfCytokineNames){
    #f <- paste0("Group ~ ", noquote(deparse(substitute(cname))))
    #MLR <- multinom(as.formula(f), data = ModelData)
    MLR_code <- paste0("Group ~ ", cname)
    #MLR_code <- paste0("MLR <- multinom(Group ~ ",cname, ", data=ModelData)")
    #eval(parse(text=MLR_code))
    MLR <- multinom(as.formula(MLR_code), data = ModelData)
    print(summary(MLR))
  
    ## Calculate Z score and p-Value for the variables in the model.
    cat("Z - scores\n\n")
    z <- summary(MLR)$coefficients/summary(MLR)$standard.errors
    print(z)
    cat("\n\n")
    
    cat("P-values")
    p <- (1 - pnorm(abs(z), 0, 1))*2
    print(p)
    cat("\n\n")
    
    cat("Coefficient log ratios")
    print(exp(coef(MLR)))
    cat("\n\n") 
    
    ## With interactions
    cat("With interactions")
    MLRI_code <- paste0("Group ~ ", cname, " + Age + PSA + `0_Genital_exam` + Any_STI")

    MLRI <- multinom(as.formula(MLRI_code), data = ModelData)
    print(summary(MLRI))
    cat("\n\n") 
  
    # Calculate Z score and p-Value for the variables in the model.
    cat("Z - scores (interaction)")
    z <- summary(MLRI)$coefficients/summary(MLRI)$standard.errors
    print(z)
    cat("\n\n") 
  
    cat("P-values (interaction)")
    p <- (1 - pnorm(abs(z), 0, 1))*2
    print(p)
    cat("\n\n") 
  
    cat("Coefficient log ratios (interaction)")
    print(exp(coef(MLRI)))
    cat("\n\n") 
  }
}

regCytokines(LinearModelData, allCytokineColnames)

MLR_all <- multinom(Group ~ ., data = LinearModelData[, c(4,5:7,11:59)])

summary(MLR_all)

## Calculate Z score and p-Value for the variables in the model.
z_all <- summary(MLR_all)$coefficients/summary(MLR_all)$standard.errors
z_all

p_all <- (1 - pnorm(abs(z_all), 0, 1))*2
p_all

exp(coef(MLR))


## Check fitted values
# head(fitted(MLR_all))

MLR_b.NGF <- multinom(Group ~ b.NGF + Age + PSA + `0_Genital_exam` + Any_STI, data = LinearModelData[, c(4,5:7,11:59)])

summary(MLR_b.NGF)

## Calculate Z score and p-Value for the variables in the model.
z_b.NGF <- summary(MLR_b.NGF)$coefficients/summary(MLR_b.NGF)$standard.errors
z_b.NGF

p_b.NGF <- (1 - pnorm(abs(z_b.NGF), 0, 1))*2
p_b.NGF

exp(coef(MLR))

