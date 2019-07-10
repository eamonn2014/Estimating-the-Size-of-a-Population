---
title: "t_4.2-1.1.Rmd"
author: "Eamonn O'Brien"
date: '`r format(Sys.time(), "%d %B, %Y")`'
knit: (function(inputFile, encoding) { 
      out_dir <- '/view/obrieea1_view/vob/CAIN457/pool/pool_002/report/pgm_a/';
          rmarkdown::render(inputFile,
                        encoding=encoding, 
      output_file=file.path(out_dir,'t_4.2-1.1.html'))})
output:  
  html_document: 
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    theme: united
    highlight: tango
    fig_caption: yes
    fig_height: 6
    fig_width: 8
    number_sections: no
  word_document: 
header-includes:
- \usepackage{eso-pic,graphicx,transparent}
- \usepackage{graphicx}
- \usepackage{fancyhdr}
- \pagestyle{fancy}
- \setlength\headheight{22pt}
- \fancyfoot[RO]{Novartis IGE025E3401 AWARE}
- \usepackage{lastpage}
- \cfoot{Page \thepage\ of \pageref{LastPage}}
---


\AddToShipoutPictureFG{
  \AtPageCenter{% or \AtTextCenter
    \makebox[0pt]{\rotatebox[origin=c]{45}{%
      \scalebox{5}{\texttransparent{0.3}{ }}%
    }}
  }
}


```{r set-options, echo=FALSE, cache=FALSE, warning = FALSE}

#/***************************************************************************#
#Filename           :      t_4.2-1.1.Rmd
#Author             :      obrieea1
#Date               :      11-Mar-2019
#R                  :      3.2.3 (2015-12-10) 
#Platform           :      x86_64-pc-linux-gnu (64-bit) 
#Project/Study      :      CAIN457\pub_2
#Description        :      Analysis 
#Assumptions:   
#Input              :      adhoc3 dataset  
#Output             :      A Rmd file is outputed to GPS and checked in to share with team
#Macros used        :      none
#---------------------------------------------------------------------------
#MODIFICATION HISTORY: 
#    <DD-MON-YYYY>,
#    <Description> 
#    
#***************************************************************************/

        rm(list=ls())

        set.seed(123)
        startTime<-proc.time()
        library(knitr)
        options(width=120)
        opts_chunk$set(comment = "", warning = FALSE, message = FALSE,
                       echo = FALSE, tidy = FALSE, size="tiny",  cache=FALSE,
                       progress=TRUE,
                         fig.width=7, fig.height=3.5,
                       cache.path = 'program_Cache/',
                       fig.path='figure/')
         
        knitr::knit_hooks$set(inline = function(x) {
          knitr:::format_sci(x, 'md')
        })
         
        
        options(scipen=999)  # remove scientific notation
        
        
        # create an R file of the code!
        # https://stackoverflow.com/questions/26048722/knitr-and-tangle-code-without-execution
        
         knit_hooks$set(purl = function(before, options) {
           if (before) return()
           input  = current_input()  # filename of input document
           output = paste(tools::file_path_sans_ext(input), 'R', sep = '.')
           if (knitr:::isFALSE(knitr:::.knitEnv$tangle.start)) {
           assign('tangle.start', TRUE, knitr:::.knitEnv)
           unlink(output)
         }
           cat(options$code, file = output, sep = '\n', append = TRUE)
         })

          # function to calculate the mode     
          getmode <- function(v) {
               uniqv <- unique(v)
               uniqv[which.max(tabulate(match(v, uniqv)))]
          }
  
          
        # https://stackoverflow.com/questions/36868287/purl-within-knit-duplicate-label-error
          
          options(knitr.duplicate.label = 'allow') # to help with purl
 
          
```

```{r introduction}
 
      cat(" \nAnalysis\n\n
            Primary objective\n
            Main effects Cox model with #1 to #13 variables (including dose) and stratified by RCT\n
            Secondary exploratory objectives\n
            Main effects Cox model excluding placebo patients #1 to #9 variables (including dose) and stratified by RCT\n
            Interaction model variables #2 to #4 all interacting with dose and stratified by RCTs\n
            SAS Model\n

            Response is 1st record of AE Candida infection\n

            TRTDURE up until 1st record of AE Candida infection\n

            stratification variable:  RCT ('STUDYID')\n
     
            predictor variables:\n
            #1 Dose ('Dose')\n
            #2 weight ('WEIGHT')\n
            #3 Sex ('SEX')\n
            #4 smoking status ('CURSMH')\n
            #5 gastro intestinal disorders('GASTRDIS_STATUS')\n
            #6 infections and infestations ('INFECT_STATUS')\n
            #7 Glucocorticosteroids ('H02A')\n
            #8 Antibiotics ('J01')\n
            #9 proton pump inhibitors ('A02BC')\n
            #10 Antidiabetic agents ('A10')\n
            #11 immunosuppresive agents ('L04A')\n
            #12 medical history:gastrooesophageal reflux disease ('GRD_STATUS')\n
            #13 diabetes ('DIAMELL_STATUS')\n
          ")
  
```

#### Da vinci and GPS directories

```{r directories, echo=TRUE, cache=FALSE}

      # Da vinci directories
      wd <- Rhome <-'/home/obrieea1/Documents/CAIN457 RISK CALCULATOR/'
      wd.code <- paste0(Rhome,'CODE')                 
      wd.data <- paste0(Rhome,'DATA')
 
      # GPS directories
      gps <- '/view/obrieea1_view/vob/CAIN457/pool/pool_002/'
      path.script <- paste0(gps, '/report/pgm_a/')  
      path.import <- paste0(gps, '/report/data_a/')  # data is here
      path.export <- paste0(gps, '/report/export/pgm_saf/')    # outputs here?
 
      # rounding functions
      p2 <- function(x) {formatC(x, format="f", digits=2)}
      p3 <- function(x) {formatC(x, format="f", digits=3)}
      p4 <- function(x) {formatC(x, format="f", digits=4)}
      
      
```
  
#### some options 
  
```{r parameters, echo =TRUE, cache=FALSE}      
      
rez <- 800 # dots per inch for plots, more dots bigger file size, quality does not change that much

status<- "final" #"draft "
     
```
  
#### user library
  
```{r user libraries , echo =TRUE, cache=FALSE}


      
# see email for Install package to DaVinci R 3.4.3 08nov
  .libPaths('~/R-libraries-3.4.3')
  library('rtf')

 #preparing for bottome row of outputs
 path.short <- gsub(gps, "", path.script)
   

``` 

#### Call in the analysis R data adhoc3 and adsubgrp

```{r read in data, echo=TRUE, cache=FALSE}
      
 
    require(Hmisc)
    require(dplyr)

    adhoc3 <- NULL
    
    #load R version of analysis dataset
    load(file=paste(path.import,"adhoc3.rda",sep=""))
    r <- adhoc3

    
    
      r_adsubgrp <- NULL
      load(file=paste(path.import,"adsubgrp.rda",sep="")) # r_adsubgrp df is loaded for X01....Rmd
      dim(r_adsubgrp)
      foo <- r_adsubgrp
     # dplyr::glimpse(d)
      
 
    d <- f <- NULL  # ensure these are empty
    d <- r          # choose either the r or SAS datset here, they are the same.
   
```
#### Preparation for modelling

```{r examining data, echo=TRUE, cache=FALSE,results='asis'}


      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # function to tabulate by treatment 
     tabs <- function(var) {
       print(var)
       print(addmargins(with (d, table(eval(parse(text=var)), Treatment, useNA="ifany"))))
       cat("\nColumn percentages\n")
       pt <- prop.table(with (d, table(eval(parse(text=var)), Treatment)),2)*100  # column percentages  
       print(pt,digits=3)
       cat("\n")
     }
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
     # create response variable~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     d$e <- ifelse(d$SGRP07FL %in% "Y", 1, ifelse(d$SGRP08FL %in% "Y" , 0, NA)) # 4219 200
     # these numbers agree with outputs
     # higher AE of interest in treated but not by much ~5% vrs. ~4%
     # with 200 events we should only fit max 13 covariates..
     
     
     ##quick tabulation
     addmargins(table( d$STUDYID, d$e, d$Treatment ))
     
     
     tabs("e")


     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # Now look at predictor variables blinded to y
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
     #RCTstudies~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     d$STUDYID <- as.factor(d$STUDYID)
     unique(d$STUDYID) # 13 studies
     
     #Dose~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("Dose")

     #1 weight~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tapply(d$WEIGHT,d$Treatment,summary)

     #2 SEX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("SEX")
   
     #3 smoking status~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("CURSMH")
     
     #4 gastro intestinal disorders~~~~~~~~~~~~~~~~~~~~~~
     tabs("GASTRDIS_STATUS")
    
     #5 infections and infestations~~~~~~~~~~~~~~~~~~~~~~
     tabs("INFECT_STATUS")
      
     #6 Glucocorticosteroids~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("H02A")
     
     #7 Antibiotics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("J01")
      
     #8 proton pump inhibitors~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("A02BC")
     
     #9 Antidiabetic agents~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("A10")
     
     #10 immunosuppresive agents~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("L04A")
     
     #11 medical history:gastrooesophageal reflux disease
     tabs("GRD_STATUS")
     
     #12 diabetes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~?
     tabs("DIAMELL_STATUS")  # **check this, is this tabulation odd?**
     
     #13 neoplasm benign malignant unspecified~~~~~~~~~~~
     tabs("NEOPLASM_STATUS")
     
     #14 endocrine disorders~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("ENDOCDIS_STATUS")
     
     #15 antimycotic agents~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("J02")
     
     #16 prior biologic agents~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("THER_USE")
     
     #17 H2 receptor antagonists~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("A02BA")

     
     tabs("diab_status")
     #18 metabolic syndrome (not in dataset at moment due to number of missing)
     
     #19 is dose but this will have to be included~~~~~~~
     
     #20 Time since first diagnosis~~~~~~~~~~~~~~~~~~~~~~
     tapply(d$timediag,d$Treatment,summary)
     
     #21 comorbidities is this in the adhoc3 dataset?~~~~?
     
     #22 obesity is this in the adhoc3 dataset?~~~~~~~~~~?
     
     #23 BMI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tapply(d$BMI,d$Treatment,summary)
     
     #24 skin and subcutaneous tissue disorders~~~~~~~~~~
     tabs("SKINDIS_STATUS")
     
     #25 Metabolism and nutritional disorders~~~~~~~~~~~~
     tabs("METABDIS_STATUS")
     
     #26 Immune system disorders~~~~~~~~~~~~~~~~~~~~~~~~~
     tabs("IMMUNDIS_STATUS")

     #27 age~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tapply(d$AGE,d$Treatment,summary)

     #28 duration of treatment until AE of interest~~~~~
     tapply(d$TRTDURE,d$Treatment,summary)
     
     
  mymean <- function(x) {mean(x, na.rm=T)}
  mysd <- function(x) {sd(x, na.rm=T)}
  
  ss <- d %>% 
  group_by(Treatment, e) %>%  #variable rather than grp
  summarise_at(.vars = vars(TRTDURE), 
                
               funs(
                 
                 C.sd=mysd,
                 B.mean=mymean,
                 I.miss = sum(is.na(.)),
                 A.count = sum(!is.na(.)), 
                 E.Q1=quantile(.,probs=c(0.25) ,na.rm = TRUE, type=2),  # type 2 match SAS
                 D.median=quantile(.,probs=c(0.50) ,na.rm = TRUE, type=2),  # 
                 F.Q3=quantile(.,probs=c(0.75) ,na.rm = TRUE, type=2),
                 G.min=min(., na.rm=TRUE),
                 H.max=max(., na.rm=TRUE)
                 #H.perc = sum(is.na(.)) / length(.)*100
                 )
                 )
  
  
 
  cat("\nTrt duration\n")
  print(ss, digits =2)   
  
 
```
## Final model, the analysis is performed here

```{r final model,eval=TRUE,echo=TRUE}
require(rms)
require(Hmisc)

varList <- Cs(A02BC, J01,H02A,A10,J02,A02BA,GASTRDIS_STATUS,INFECT_STATUS, L04A)

for (i in varList) {
    d[[i]] <- relevel(as.factor(d[[i]]), ref = "NEVER")
}
 

 
 
  d$H02A <- plyr::mapvalues(d$H02A, from = c("CURRENT","BEFORE", "BEFORE/CURRENT"), to = c("CURRENT","PRIOR", "PRIOR/CURRENT"))
   
  d$H02A <- factor(d$H02A, levels = c("NEVER", "CURRENT", "PRIOR","PRIOR/CURRENT"))
  
     d$Dose <- as.factor(d$Dose)

     dd <<- datadist(d)        # needed to use global assigner due to error if otherwise
     options(datadist = "dd")

     # stratified Cox Model
     S <- Surv(d$TRTDURE, d$e)
     f <- NULL
  
     f2 <- cph(S ~ Dose + rcs(WEIGHT,4)+ SEX +CURSMH + GASTRDIS_STATUS + INFECT_STATUS + H02A +  strat(STUDYID), x=TRUE, y=TRUE,  data=d,surv=TRUE,time.inc=150)    #15 df
     
     print(f2)
          
     cox.zph(f2, transform = "log")  
     
    f2i <- cph(S ~ Dose *( rcs(WEIGHT,4)+ SEX +CURSMH +
                                 GASTRDIS_STATUS + INFECT_STATUS + H02A) +
                strat(STUDYID), x=TRUE, y=TRUE,  data=d)    #15 df

     
        lr<- lrtest(f2,f2i)
        lrp<-lr$stats["P"][[1]]
        
        
     #checking
        
        
     hazard.ratio.plot(d$WEIGHT,S, e=30, legendloc='ll',smooth = T,antilog=T)
     
     #~~~~~~~~~~~~~~~~
     
     mod <-f2
     
  # Assess PH by estimating log relative hazard over time
  
  

        
        pha <- cox.zph(f2, transform = "log")

        for(i in 1:15) {
        plot(pha[i])
        abline(h=f2$coefficients[i],col='red')
        par(mfrow=c(1,1))
        }

        
           
    
 
        gph.p <-pha$table[,3]["GLOBAL"][[1]]
        
        v <- pha$table[pha$table[,3]<0.05,]
        v<-rownames(v)
        v<- gsub("GLOBAL",NA,(v))
        v<-v[!is.na(v)]
        v

    
     
     ## sensitivity analyis nostrat?
 
          s1 <- cph(S ~ Dose + rcs(WEIGHT,4)+ SEX +CURSMH + GASTRDIS_STATUS + INFECT_STATUS + H02A  , x=TRUE, y=TRUE,  data=d)    #15 dff1
      
     
    ##sensitivity 2 statification by disease
          
          
          
            pop1 <- Cs(CAIN457A2302,CAIN457A2303,CAIN457A2304,CAIN457A2308,
             CAIN457A2309,CAIN457A2317,CAIN457A2313)
  
   pop2 <- Cs(CAIN457F2312,CAIN457F2318,CAIN457F2336,CAIN457F2342) 
  
   pop3 <- Cs(CAIN457F2310,CAIN457F2320)
  

d$pop <- ifelse(d$STUDYID %in% pop1, "Arthritis",
                ifelse(d$STUDYID %in% pop2,"Arthritis Psorasis",
                     ifelse(d$STUDYID %in% pop3,"Ankolysing Spondilitis",  NA  )))
          
          d$pop <-as.factor(d$pop)
          
          
      dd <<- datadist(d)        # needed to use global assigner due to error if otherwise
     options(datadist = "dd")

     # stratified Cox Model
     S <- Surv(d$TRTDURE, d$e)
     f <- NULL
 
     
     
   s2 <- cph(S ~ Dose + rcs(WEIGHT,4)+ SEX +CURSMH + GASTRDIS_STATUS + INFECT_STATUS + H02A + strat(pop) , x=TRUE, y=TRUE,  data=d)    #15 df
     
          
     

```

## t_4.2-1.1 regression table

```{r delete0,eval=TRUE,echo=TRUE}     
 



##basic approach
covariate.labels=c("Secukinumab dose (150mg)","Secukinumab dose (300mg)",
            "Patient weight (kg)","Patient weight (kg) (non linear 1)","Patient weight (kg) (non linear 2)",
                           "Patient sex (male)",
                           "Smoking status (yes)",
                           "Gastrointestinal disorders (current only)",
                           "Gastrointestinal disorders (prior only)",
                           "Gastrointestinal disorders (prior and current)",
                           "Infections and infestations (current only)",
                           "Infections and infestations (prior only)",
                           "Infections and infestations (prior and current)",
                           "Glucocorticosteroid use (current only)",
                           "Glucocorticosteroid use (prior only)",
                           "Glucocorticosteroid use (prior and current)"
                           )


# 
# 
#  f <- f2$sformula
#   
#       
#      
# 
#   # manage the model formula for incorporation into RTF output
#   ff <- as.character(f)
#   ff[3] <- gsub("\\<rcs\\>","", ff[3])
#   ff[3] <- gsub("[()]", "", ff[3])
#   ff[3] <- gsub('[0-9]+', '', ff[3])
#   ff[3] <- gsub(',','',ff[3])
#   ff[3] <- gsub('\\.','',ff[3])
#   ff[3] <- gsub('  ',' ',ff[3])
#   ff[3] <- gsub('  ',' ',ff[3])
#   ff[3] <- gsub('\\c\\>',' ',ff[3])
#   ff[3] <- gsub('  ',' ',ff[3])
#   ff[3] <- gsub('  ',' ',ff[3])
#   fff <- paste(ff[2],ff[1],ff[3], sep=" ")
#   fff <- toupper(fff)
#   fff <- gsub("AVISIT", "VISIT", fff)
#   fff<- as.character(fff)
#   fff<- tolower(fff)







#est <-  exp(f2$coefficients) #replace coefficents with odds ratios

HR.vector<-exp(f2$coef)

HRCI.vector<-exp(confint(f2))

beta = coef(f2)
B_SE = sqrt(diag(vcov(f2)))
pvalue =  pnorm(-abs(beta) / B_SE)  * 2

res <- cbind(covariate.labels, 
             exp( as.data.frame(beta)[,1]), 
              as.data.frame(HRCI.vector)[,1],
                    as.data.frame(HRCI.vector)[,2],
             as.data.frame(pvalue)[,1]
             )
                          

res<- as.data.frame(res)

res$V2<-as.numeric(levels(res$V2))[res$V2]
res$V3<-as.numeric(levels(res$V3))[res$V3]
res$V4<-as.numeric(levels(res$V4))[res$V4]
res$V5<-as.numeric(levels(res$V5))[res$V5]


star <-c("","","","","","","","","","","","","","*","","")



res$V2 <-p2(res$V2)
res$V3 <-p2(res$V3)
res$V4 <-p2(res$V4)
res$V2 <- paste0(res$V2,star)

res$V5 <- Hmisc::format.pval(res$V5,  eps = .001,digits=3)
res$CI <- paste0(" (",res$V3,", ",res$V4,")" )


res<- res[,c("covariate.labels","V2","CI","V5")]

names(res) <- c("Variable","Adjusted Hazard Ratio","95% Wald Confidence Interval", "    P-value")
print(res,center=T)
 

text1 <- "t_4.2-1.1.rtf"
figname3 <- "4.2-1.1"
namex3 <- "t_4.2-1.1.R                "
namex3 <- "t_4.2-1.1.R"

title1 <- "Cox proportional hazards regression coefficients of each predictor level (All patients)\n"
title2 <- "Global proportional hazard test P-value = "
title3 <- "; Likelihood ratio test P-value = "

foot1 <- "-Secukinumab dose = 'placebo', Patient sex = 'female', Smoking status = 'no' and for all other variables 'never' are the reference levels for the model predictors.\n"
foot1a <- "-An adjusted hazard ratio less than 1 indicates better outcomes compared to the reference level for the categorical predictors (i.e. all the predictors except patient weight).\n"
foot2 <- "-Patient weight is not assumed to be linear, the 'Patient weight (kg) (non linear x)' terms represent differences in cubes that restrict the function of patient weight to be linear beyond the outer knots, interpreting individual terms is not advisable.\n"
foot3 <- "-The presented global proportional hazard P-value is the result of the Grambsch-Therneau hypothesis test of the model assumption that hazards are proportional over the complete follow up time."
foot4 <- paste0(" -The fitted Cox proportional hazards (PH) model from which the coefficients are presented was determined a better fit to the data than one including interactions between dose and each of the predictors as evidenced by the likelihood ratio test P-value.\n" )
foot5 <-"The symbol '*' accompanying an adjusted hazard ratio denotes there is evidence of violation of the proportional hazards assumption for the predictor level. No corrective action has been taken, instead the 'principle of conservatism' is invoked, the P-value would likely be more impressive had the predictor exhibited proportionality.\n"
foot6 <-"-The fitted Cox PH model studies the follow up time and any adverse event within candida infection MedDRA high level term during secukinumab therapy (and placebo), in a pooled population including all psoriasis, psoriasis arthritis and ankylosing spondylitis patients using categorical predictors secukinumab dose and placebo, patient sex, smoking status, gastrointestinal disorders, infections and infestations other than candida infection (any event within candida infection MedDRA high level term) and glucocorticosteroid use. Patient weight in kg is a continuous predictor not assumed to be linear; modelled using restricted cubic splines. A stratification term comprising each study used in the analysis is included in the model.\n"
projcode='CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy\n'
  

  rtffile <- RTF(paste(path.export,text1,sep=""), font.size=11, omi=c(1,1,1,1.4), width=11, height=12)
 
  addText(rtffile, paste(projcode), bold=F, italic=F) 
  
  addText(rtffile, paste("Table",figname3, "(Page 1 of 1)\n"), bold=F,italic=F)
  
  addText(rtffile, paste(title1), bold=F, italic=F) 
  
  addParagraph(rtffile, paste0("              ")  )
   
  addParagraph(rtffile, paste0(title2,p3(gph.p),title3,p3(lrp))  )

  addTable(rtffile, res, col.justify=c('L','C','C','C'), col.widths=c(3.1,1.9,2.4,1.15 ) )
  
  addParagraph(rtffile, paste0("              ")  )
  
  addParagraph(rtffile, paste(foot4,foot6,foot1,foot2,foot1a,foot3,foot5))  
    
  # addText(rtffile, paste0("Source:",path.script,namex3,format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" ") )
    addText(rtffile, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),"                                                                        ",status,   sep=" ") )

  
   
  done(rtffile)  

  
  
```
## t_4.2-1.2 anova table

```{r delete1,eval=FALSE,echo=TRUE}
  
        ##basic approach
        covariate.labels=c("Secukinumab dose", 
                                   "Patient weight" ,
                                   "Patient sex",
                                   "Smoking status",
                                   "Gastrointestinal disorders",
                                   "Infections and infestations",
                                   "Glucocorticosteroid use"
                                   )
          res<-anova(f2)
          
          res<-res[c(-3,-9),]
          
          res<-as.data.frame(res)
          
          rownames(res) <- NULL
          
          res <- cbind(covariate.labels, res)
                        
          res<- as.data.frame(res)
          
          res$`Chi-Square`<-p2(res$`Chi-Square`)
          res$P <- Hmisc::format.pval(res$P,  eps = .001,digits=3)
          
          names(res) <- c("Variable","    Chi-square","      Degrees of freedom", "    P-value")

 


print(res,center=T)
 

text1 <- "t_4.2-1.2.rtf"
figname3 <- "4.2-1.2"
namex3 <- "t_4.2-1.2.R                "
namex3 <- "t_4.2-1.2.R"
projcode='CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy\n'

title1 <- "Cox proportional hazards regression for each predictor - ANOVA analysis (All patients)\n"
title2 <- "  "
title3 <- "  "

foot1 <- " -ANOVA denotes analysis of variance.\n"
foot2 <- "-The output presents a P-value for a hypothesis test for each predictor, testing if they are risk factors for the AE candid infection whilst adjusting for all other predictors.\n"
foot3 <- "-The expected value of a chi-square statistic is equal to the degrees of freedom under the null hypothesis of no association.\n"
foot4 <- "-The degrees of freedom are the number of regression parameters that uniquely define the hypothesis. For a factor with n levels there are n-1 degrees of freedom. The weight predictor requires k-1 degrees of freedom, here k is the total number of knots.\n"
 

  rtffile <- RTF(paste(path.export,text1,sep=""), font.size=11, omi=c(1,1,1,1.4), width=10.5, height=12)
 
  addText(rtffile, paste(projcode), bold=F, italic=F) 
  
  addText(rtffile, paste("Table",figname3, "(Page 1 of 1)\n"), bold=F,italic=F)
  
  addText(rtffile, paste(title1), bold=F, italic=F) 
  
  addParagraph(rtffile, paste0("              ")  )
   
  addParagraph(rtffile, paste0(title2,title3)  )

  addTable(rtffile, res, col.justify=c('L','C','C','C'), col.widths=c(2.9,1.5,2.4,1.15 ) )
  
  addParagraph(rtffile, paste0("              ")  )
  
  addParagraph(rtffile, paste(foot1,foot3,foot4,foot2, foot6))  
    
  #addText(rtffile, paste0("Source:",path.script,namex3,format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  
   addText(rtffile, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                               ",status,   sep=" ") )
  
   
  
  done(rtffile) 


```

## f_4.2-1.1  Forest plot

```{r delete2,eval=FALSE,echo=TRUE}


 myplot2 <- function() {
 
covariate.labels=c("Secukinumab dose (mg)", 
                                   "Patient weight (kg)" ,
                                   "Patient sex",
                                   "Smoking status",
                                   "Gastrointestinal disorders",
                                   "Infections and infestations",
                                   "Glucocorticosteroid use"
)
         
          label(d$Dose) <- covariate.labels[1]
          label(d$WEIGHT) <- covariate.labels[2]
          label(d$SEX) <- covariate.labels[3]
          label(d$CURSMH) <- covariate.labels[4]
          label(d$GASTRDIS_STATUS) <- covariate.labels[5]
          label(d$INFECT_STATUS) <- covariate.labels[6]
          label(d$H02A) <- covariate.labels[7]
         
       dd <<- datadist(d)        # needed to use global assigner due to error if otherwise
       options(datadist = "dd")
           ##secondary
     f2 <- cph(S ~ Dose + rcs(WEIGHT,4)+ SEX +CURSMH + GASTRDIS_STATUS + INFECT_STATUS + H02A +
                strat(STUDYID), x=TRUE, y=TRUE,  data=d,surv=TRUE,time.inc=150)    #15 df
     
    s<- (summary(f2,vnames='labels', Dose="0", WEIGHT=c(97,71),SEX="F"))
    
    plot( s,log=T,q=.95,at=c(0.1,.25,.5,1,2,4),
    col=rgb(red=.1,green=.1,blue=.8,alpha=c(.4,.2,.2)),
    col.points=rgb(red=.1,green=.1,blue=.8,alpha=.7), 
    pch=17, 
    lwd=if(length(q) == 1) 3 else 2 : (length(q) + 1), digits=4,
     main="Adjusted Hazard Ratio") 
    
  
  }

    
    
text1 <- "f_4.2-1.1.rtf"
figname3 <- "4.2-1.1"
namex3 <- "f_4.2-1.1.R                "
namex3 <- "f_4.2-1.1.R"

foot1 <- " - Forest plot of the adjusted hazard ratio with associated 95% confidence interval for each predictor in the model."
foot2 <- "The blue triangles denote hazard ratio point estimates, confidence intervals are denoted by blue bands.\n"
foot3 <- "- Between predictors, narrower bands indicate greater precision in the estimate.\n"
foot3 <-"- The adjusted hazard ratio axis is plotted on the log scale but labelled with anti-logs.\n"
foot3a <- "- Secukinumab dose = 'placebo', Patient sex = 'female', Smoking status = 'no' and for all other variables 'never' are the reference levels for the model predictors.\n"
foot4 <- "- The vertical dashed line represents the null value, a hazard ratio of 1. An adjusted hazard ratio less than 1 indicates better outcomes compared to the reference level for the categorical predictors. In the case of the continuous predictor weight, the inter-quartile range effect is presented. This is the average effect comparing two patients, one with a body weight equal to the lower quartile (25th percentile) and one with a body weight equal to the upper quartile (75th percentile) of the weight distribution and who are identical in all the other predictors.\n"

foot6a <-"-The fitted Cox PH model studies the follow up time and any adverse event within candida infection MedDRA high level term during secukinumab therapy (and placebo), in a pooled population including all psoriasis, psoriasis arthritis and ankylosing spondylitis patients using categorical predictors secukinumab dose and placebo, patient sex, smoking status, gastrointestinal disorders, infections and infestations other than candida infection (any event within candida infection MedDRA high level term) and glucocorticosteroid use. Patient weight in kg is a continuous predictor not assumed to be linear; modelled using restricted cubic splines. A stratification term comprising each study used in the analysis is included in the model."
 
  rtf3 <- RTF(paste(path.export,text1, sep=''), font.size=11, omi=c(1,1,1,1.4),
              width=11,  height=10)
 
  addText(rtf3, paste(projcode), bold=F, italic=F) 
  
  addText(rtf3, paste("Figure",figname3, "(Page 1 of 1)\n"), bold=F,italic=F)
  
  addText(rtf3, paste("Forest plot of effects (hazard ratio) of each predictor level (All patients)\n"), bold=F, italic=F) 
  
   
   

  # addParagraph(rtf3, 
  # paste0("Forest plot of effects (regression coefficients) of each predictor level (All patients)"))

  addPlot(rtf3, plot.fun=myplot2, width=8, height=4.5, res=rez)  
  
  addParagraph(rtf3, paste(foot1,foot2,foot3,foot3a,foot4, foot6a))
  
  # addText(rtf3, paste0("Source:",path.script,namex3,format(Sys.time(), "%b %d %Y, %H:%M:%S"), 
  #                      sep=" "))
  
  addText(rtf3, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),"                                                                       ",status,   sep=" ") )
  done(rtf3) 

  

```

## f_4.2-1.2 ANOVA plot

```{r delete3, eval=FALSE,echo=TRUE}
 
 
  myplot2 <- function() {
 
covariate.labels=c("Secukinumab dose", 
                                   "Patient weight" ,
                                   "Patient sex",
                                   "Smoking status",
                                   "Gastrointestinal disorders",
                                   "Infections and infestations",
                                   "Glucocorticosteroid use"
)
         
          label(d$Dose) <- covariate.labels[1]
          label(d$WEIGHT) <- covariate.labels[2]
          label(d$SEX) <- covariate.labels[3]
          label(d$CURSMH) <- covariate.labels[4]
          label(d$GASTRDIS_STATUS) <- covariate.labels[5]
          label(d$INFECT_STATUS) <- covariate.labels[6]
          label(d$H02A) <- covariate.labels[7]
         
       dd <<- datadist(d)        # needed to use global assigner due to error if otherwise
       options(datadist = "dd")
           ##secondary
     f2 <- cph(S ~ Dose + rcs(WEIGHT,4)+ SEX +CURSMH + GASTRDIS_STATUS + INFECT_STATUS + H02A +
                strat(STUDYID), x=TRUE, y=TRUE,  data=d,surv=TRUE,time.inc=150)    #15 df
     
     
       s<- anova(f2, vnames='labels' )
    plot(s, trans=NULL, margin=c('chisqminusdf','chisq','d.f.'), rm.totals = T ,ntrans=50 )
  
  }

    text1 <- "f_4.2-1.2.rtf"
figname3 <- "4.2-1.2"
namex3 <- "f_4.2-1.2.R                "
 namex3 <- "f_4.2-1.2.R"
 

  # add variable importance plot
 
  
   rtf3 <- RTF(paste(path.export,text1, sep=''), font.size=11, omi=c(1,1,1,1.4),
              width=11,  height=10)
 
   addText(rtf3, paste(projcode), bold=F, italic=F) 
   addText(rtf3, paste("Figure",figname3, "(Page 1 of 1)\n"), bold=F,italic=F)
   addText(rtf3, 
   paste("Ranking the importance of each predictor of candida infection in patients (All patients)\n"), bold=F, italic=F) 

  addPlot(rtf3, plot.fun=myplot2, width=8, height=4.5, res=rez)  
  
   addParagraph(rtf3, "- A graphical depiction of the analysis of variance table from Cox proportional hazards.\n- The importance of each predictor is presented as judged by the Wald chi-square statistic minus the predictor degrees of freedom (d.f.) for assessing the partial effect of each predictor in the model, adjusting for all other predictors.\n- The expected value of a chi-square statistic is equal to the d.f. (shown in rightmost column) under the null hypothesis of no association.\n- Predictors requiring a large number of parameters to achieve the chi-square value are penalized.\n- Larger values for a predictor denote greater relative importance.\n- The fitted Cox PH model studies the follow up time and any adverse event within candida infection MedDRA high level term during secukinumab therapy (and placebo), in a pooled population including all psoriasis, psoriasis arthritis and ankylosing spondylitis patients using categorical predictors secukinumab dose and placebo, patient sex, smoking status, gastrointestinal disorders, infections and infestations other than candida infection (any event within candida infection MedDRA high level term) and glucocorticosteroid use. Patient weight in kg is a continuous predictor not assumed to be linear; modelled using restricted cubic splines. A stratification term comprising each study used in the analysis is included in the model.\n" )
  # 
  # addText(rtf3, paste0("Source:",path.script,namex3,format(Sys.time(), "%b %d %Y, %H:%M:%S"), 
  #                      sep=" "))
  # 
  
    addText(rtf3, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),"                                                                        ",status,   sep=" ") )
  
  done(rtf3) 
  
  
```

## t_4.2-3.1 ANOVA table sensitivity analyis 1 no stratification term

```{r delete4,eval=FALSE,echo=TRUE}
  
        ##basic approach
        covariate.labels=c("Secukinumab dose", 
                                   "Patient weight" ,
                                   "Patient sex",
                                   "Smoking status",
                                   "Gastrointestinal disorders",
                                   "Infections and infestations",
                                   "Glucocorticosteroid use"
                                   )
          res<-anova(s1)
          
          res<-res[c(-3,-9),]
          
          res<-as.data.frame(res)
          
          rownames(res) <- NULL
          
          res <- cbind(covariate.labels, res)
                        
          res<- as.data.frame(res)
          
          res$`Chi-Square`<-p2(res$`Chi-Square`)
          res$P <- Hmisc::format.pval(res$P,  eps = .001,digits=3)
          
          names(res) <- c("Variable","    Chi-square","      Degrees of freedom", "    P-value")

 


print(res,center=T)
 

text1 <- "t_4.2-3.1.rtf"
figname3 <- "4.2-3.1"
namex3 <- "t_4.2-3.1.R                "
namex3 <- "t_4.2-3.1.R"
projcode='CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy\n'

title1 <- "Sensitivity Analysis (I) Cox proportional hazards regression for each predictor - ANOVA analysis (All patients)\n"
title2 <- "  "
title3 <- "  "

foot1 <- " -ANOVA denotes analysis of variance.\n"
foot2 <- "-The output presents a P-value for a hypothesis test for each predictor, testing if they are risk factors for the AE candid infection whilst adjusting for all other predictors.\n"
foot3 <- "-The expected value of a chi-square statistic is equal to the degrees of freedom under the null hypothesis of no association.\n"
foot4 <- "-The degrees of freedom are the number of regression parameters that uniquely define the hypothesis. For a factor with n levels there are n-1 degrees of freedom. The weight predictor requires k-1 degrees of freedom, here k is the total number of knots.\n"
foot6 <-"-The fitted Cox PH model studies the follow up time and any adverse event within candida infection MedDRA high level term during secukinumab therapy (and placebo), in a pooled population including all psoriasis, psoriasis arthritis and ankylosing spondylitis patients using categorical predictors secukinumab dose and placebo, patient sex, smoking status, gastrointestinal disorders, infections and infestations other than candida infection (any event within candida infection MedDRA high level term) and glucocorticosteroid use. Patient weight in kg is a continuous predictor not assumed to be linear; modelled using restricted cubic splines. A stratification term is not included in the model.\n"
 
  rtffile <- RTF(paste(path.export,text1,sep=""), font.size=11, omi=c(1,1,1,1.4), width=10.5, height=12)
 
  addText(rtffile, paste(projcode), bold=F, italic=F) 
  
  addText(rtffile, paste("Table",figname3, "(Page 1 of 1)\n"), bold=F,italic=F)
  
  addText(rtffile, paste(title1), bold=F, italic=F) 
  
  addParagraph(rtffile, paste0("              ")  )
   
  addParagraph(rtffile, paste0(title2,title3)  )

  addTable(rtffile, res, col.justify=c('L','C','C','C'), col.widths=c(2.9,1.5,2.4,1.15 ) )
  
  addParagraph(rtffile, paste0("              ")  )
  
  addParagraph(rtffile, paste(foot1,foot3,foot4,foot2, foot6))  
    
  #addText(rtffile, paste0("Source:",path.script,namex3,format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  
   addText(rtffile, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                               ",status,   sep=" ") )
  done(rtffile)    


```

## t_4.2-3.2 ANOVA table sensitivity analyis II  stratification term include 3 pops

```{r delete5,eval=FALSE,echo=TRUE}
  
        ##basic approach
        covariate.labels=c("Secukinumab dose", 
                                   "Patient weight" ,
                                   "Patient sex",
                                   "Smoking status",
                                   "Gastrointestinal disorders",
                                   "Infections and infestations",
                                   "Glucocorticosteroid use"
                                   )
          res<-anova(s2)
          
          res<-res[c(-3,-9),]
          
          res<-as.data.frame(res)
          
          rownames(res) <- NULL
          
          res <- cbind(covariate.labels, res)
                        
          res<- as.data.frame(res)
          
          res$`Chi-Square`<-p2(res$`Chi-Square`)
          res$P <- Hmisc::format.pval(res$P,  eps = .001,digits=3)
          
          names(res) <- c("Variable","    Chi-square","      Degrees of freedom", "    P-value")

 


print(res,center=T)
 

text1 <- "t_4.2-3.2.rtf"
figname3 <- "4.2-3.2"
namex3 <- "t_4.2-3.2.R                "
namex3 <- "t_4.2-3.2.R"
projcode='CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy\n'

title1 <- "Sensitivity Analysis (II) Cox proportional hazards regression for each predictor - ANOVA analysis (All patients)\n"
title2 <- "  "
title3 <- "  "

foot1 <- " -ANOVA denotes analysis of variance.\n"
foot2 <- "-The output presents a P-value for a hypothesis test for each predictor, testing if they are risk factors for the AE candid infection whilst adjusting for all other predictors.\n"
foot3 <- "-The expected value of a chi-square statistic is equal to the degrees of freedom under the null hypothesis of no association.\n"
foot4 <- "-The degrees of freedom are the number of regression parameters that uniquely define the hypothesis. For a factor with n levels there are n-1 degrees of freedom. The weight predictor requires k-1 degrees of freedom, here k is the total number of knots.\n"
foot6 <-"-The fitted Cox PH model studies the follow up time and any adverse event within candida infection MedDRA high level term during secukinumab therapy (and placebo), in a pooled population including all psoriasis, psoriasis arthritis and ankylosing spondylitis patients using categorical predictors secukinumab dose and placebo, patient sex, smoking status, gastrointestinal disorders, infections and infestations other than candida infection (any event within candida infection MedDRA high level term) and glucocorticosteroid use. Patient weight in kg is a continuous predictor not assumed to be linear; modelled using restricted cubic splines. A stratification term comprising patient disease (psoriasis, psoriasis arthritis and ankylosing spondylitis) is included in the model.\n"
 

  rtffile <- RTF(paste(path.export,text1,sep=""), font.size=11, omi=c(1,1,1,1.4), width=10.5, height=12)
 
  addText(rtffile, paste(projcode), bold=F, italic=F) 
  
  addText(rtffile, paste("Table",figname3, "(Page 1 of 1)\n"), bold=F,italic=F)
  
  addText(rtffile, paste(title1), bold=F, italic=F) 
  
  addParagraph(rtffile, paste0("              ")  )
   
  addParagraph(rtffile, paste0(title2,title3)  )

  addTable(rtffile, res, col.justify=c('L','C','C','C'), col.widths=c(2.9,1.5,2.4,1.15 ) )
  
  addParagraph(rtffile, paste0("              ")  )
  
  addParagraph(rtffile, paste(foot1,foot3,foot4,foot2, foot6))  
    
  #addText(rtffile, paste0("Source:",path.script,namex3,format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  
   addText(rtffile, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                               ",status,   sep=" ") )
  
   
  
  done(rtffile) 

```
  
## f_4.2-2.1 repeating KMcode treat v placebo to improve using log scale 
  
```{r delete6,eval=FALSE,echo=TRUE}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## plot function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

title2 <- "Kaplan Meier plot of dose (placebo 0mg, secukinumab 150mg and 300mg) candida infection free probability (All patients)"
title1 <-'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


  
surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x") {
  
  # Plot function  
  
  dt1 <- data
  
  myplot <- function() {
  
  d2 <- dt1
  
  
  S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  #  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
  fit <- do.call(survfit,
                 list(formula = S ~ Dose, data = d2))
  
  
      print(addmargins(table(d2$Dose, d2$Treatment)))
 
    print(addmargins(table(d2$Dose, d2$e)))
  
  # n<-fit$numevents
  # x<-fit$exposure
  # h<-fit$numevents/fit$exposure
  
     gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better readability 
    ylim = c(.85, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 8),
    palette = "jco",
    
        
    legend.title = "Secukinumab dose (mg)",
    legend.labs = c("0 (Placebo)", "150","300")
    )
     
    
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
  gg <- gg +  
  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000) ),
  labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
 cat("log times that n risk and n of events are calculated\n")
 print(unique(gg$data.survtable$time))
 cat("exponetiated times that n risk and n of events\n")
 print(exp(unique(gg$data.survtable$time)))
 print(gg$data.survtable)
  cat("\n")   
  
}

   
  
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   S <- Surv((dt1$TRTDURE), dt1$e)  
  # 
   ex<-survreg(S ~ Dose, dt1 , dist="exponential")
   exp(-coefficients(ex)[1])  # hazard treat
   exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz placebo
   exp(-coefficients(ex)[1] + -coefficients(ex)[3]) #haz placebo
   require(rms)
   f2 <- npsurv(S ~ Dose, dt1)
  #     
   n <- f2$numevents
   x <- f2$exposure/365.25 # change to years rather than days
   tapply(dt1$TRTDURE, dt1$Dose,sum,na.rm=TRUE)
   h <- n/x #f2$numevents/f2$exposure
   namez <- gsub(".*=","",names(h))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
  addText(rtf, " \n", bold=F, italic=F)
  
   r=7
  r2=3#5
  
    i=1
  addText(rtf, paste(text3,namez[i],"mg =",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),   text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
  i=2
 
  addText(rtf, paste(text3,namez[i],"mg =",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),   text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)    
  
  i=3
  addText(rtf, paste(text3,namez[i],"mg =",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),   text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  
 

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
       addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years);")
              addText(rtf, " The lower table depicts the number of patients at risk and the cumulative number of events in each arm.\n")
    addText(rtf, "- The time axis is plotted on the log scale but labelled with anti-logs.\n\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                         ",status,   sep=" ") )
  done(rtf) 

}





 


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
#id <- 28 # only change this and the file names and data used will all be altered accordingly
id<-1
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, 
              text1=namex2, 
              text2=title2 , 
              #text3='Number of events in treated=',text4='; hazard=', text5='Number of events placebo=', text6='')
               text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')
          

 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ Dose, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.85, 1)))
 
```
  
## f_4.2-2.2 Kaplan Meier plot of patient sex candida infection free probability (All patients)

```{r delete7,eval=FALSE,echo=TRUE}

title2 <- "Kaplan Meier plot of patient sex candida infection free probability (All patients)"
title1 <-'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x" ,var="x") {
  
  # Plot function  
  
  dt1 <- data
  
  myplot <- function() {
  
  d2 <- dt1
  
  #S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
   
  fit <- do.call(survfit,
                 list(formula = S ~ SEX, data = d2))  #eval(parse(text=var))
  
      print(addmargins(table(d2$SEX, d2$Treatment)))
 
    print(addmargins(table(d2$SEX, d2$e)))
  
  
    # n<-fit$numevents
    # x<-fit$exposure
    # h<-fit$numevents/fit$exposure
  
    gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better readability 
    ylim = c(.75, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 8),
    palette = "jco",
    
    legend.title = "Patient sex",
    legend.labs = c("Female", "Male")
    
    )
     
    
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
#  gg <- gg +  
#  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), #log(400),log(1000),log(2000) ),
 # labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
 cat("times that n risk and n of events are displayed\n")
 print(unique(gg$data.survtable$time))
 print(gg$data.survtable)
   cat("\n")    
     
}

 
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- Surv((dt1$TRTDURE), dt1$e)  
  
  ex<-survreg(S ~ SEX, dt1 , dist="exponential")
  exp(-coefficients(ex)[1])  # hazard treat
  exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz placebo
  
  require(rms)
  f2 <- npsurv(S ~ SEX, dt1)
      
  n <- f2$numevents
  x <- f2$exposure/365.25 # change to years rather than days
  tapply(dt1$TRTDURE, dt1$SEX,sum,na.rm=TRUE)
  h <- n/x #f2$numevents/f2$exposure
      
  namez <- gsub(".*=","",names(h))
  namez <- c("Female","Male")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
  addText(rtf, " \n", bold=F, italic=F)
  r=7
  r2=3
  
    i=1
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
  i=2
 
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
  

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
       addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years);")
              addText(rtf, " The lower table depicts the number of patients at risk and the cumulative number of events in each stratum.\n\n")
  #  addText(rtf, "- The time axis is plotted on the log scale but labelled with anti-logs.\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                            ",status,   sep=" ") )
  done(rtf) 

}


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
#id <- 28 # only change this and the file names and data used will all be altered accordingly
id<-2
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, var="SEX",
              text1=namex2, 
              text2=title2 , 
           text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')
          

 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ SEX, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.75, 1)))


```
  
## f_4.2-2.3 Kaplan Meier plot of smoking status candida infection free probability (All patients)

```{r delete8,eval=FALSE,echo=TRUE}

title2 <- "Kaplan Meier plot of smoking status candida infection free probability (All patients)"
title1 <-'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x" ,var="x") {
  
  # Plot function  
  
  dt1 <- data
  
  myplot <- function() {
  
  d2 <- dt1
  
  #S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
    

  
   
  fit <- do.call(survfit,
                 list(formula = S ~ CURSMH, data = d2))  #eval(parse(text=var))
  
  
        print(addmargins(table(d2$CURSMH, d2$Treatment)))
 
    print(addmargins(table(d2$CURSMH, d2$e)))
  
    # n<-fit$numevents
    # x<-fit$exposure
    # h<-fit$numevents/fit$exposure
  
    gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better readability 
    ylim = c(.75, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 8),
    palette = "jco",
        
    legend.title = "Smoking status",
    legend.labs = c("Non Smoker", "Smoker")
    )
     
    
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
#  gg <- gg +  
#  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), #log(400),log(1000),log(2000) ),
 # labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
  cat("times that n risk and n of events are displayed\n")
 print(unique(gg$data.survtable$time)) 
 print(gg$data.survtable)
     cat("\n")     
}

 
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- Surv((dt1$TRTDURE), dt1$e)  
  
  ex<-survreg(S ~ CURSMH, dt1 , dist="exponential")
  exp(-coefficients(ex)[1])  # hazard treat
  exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz placebo
  
  require(rms)
  f2 <- npsurv(S ~ CURSMH, dt1)
      
  n <- f2$numevents
  x <- f2$exposure/365.25 # change to years rather than days
  tapply(dt1$TRTDURE, dt1$CURSMH,sum,na.rm=TRUE)
  h <- n/x #f2$numevents/f2$exposure
  
  
  namez <- gsub(".*=","",names(h))
     namez = c("non smoker", "smoker")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
addText(rtf, " \n", bold=F, italic=F)
  r=7
  r2=3
  
    i=1
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),    text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
  i=2
 
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years)\n\n")
          #    addText(rtf, " The lower table depicts the number of patients at risk and the cumulative number of events in each arm.\n")
   # addText(rtf, "- The time axis is plotted on the log scale but labelled with anti-logs.\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                         ",status,   sep=" ") )
  done(rtf) 

}


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
id<-3
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, var="SEX",
              text1=namex2, 
              text2=title2 , 
                  text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')
          

 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ CURSMH, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.75, 1)))
  
```

## f_4.2-2.4 Kaplan Meier plot of patient gastrointestinal disorder candida infection free probability (All patients)

```{r delete9, eval=FALSE,echo=TRUE}

title2 <- "Kaplan Meier plot of patient gastrointestinal disorder candida infection free probability (All patients)"
title1 <-'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x" ,var="x") {
  
  # Plot function  
  
  dt1 <- data
  
  myplot <- function() {
  
  d2 <- dt1
  
  #S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
  
  
  
   
  fit <- do.call(survfit,
                 list(formula = S ~ GASTRDIS_STATUS, data = d2))  #eval(parse(text=var))
  
          print(addmargins(table(d2$GASTRDIS_STATUS, d2$Treatment)))
 
    print(addmargins(table(d2$GASTRDIS_STATUS, d2$e)))
  
  
    # n<-fit$numevents
    # x<-fit$exposure
    # h<-fit$numevents/fit$exposure
  
    gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better
    #readability 
    #legend.labs=c("NEVER","CURRENT","PRIOR","PRIOR/CURRENT"),
    legend.labs=c("Never","Current","Prior","Prior/current"),
    ylim = c(.75, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 7),
    palette = "jco",
        legend.title = "Gastrointestinal disorder"
    )
     

 
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
#  gg <- gg +  
#  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), #log(400),log(1000),log(2000) ),
 # labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
  cat("times that n risk and n of events are displayed\n")
 print(unique(gg$data.survtable$time)) 
 print(gg$data.survtable)
     cat("\n")     
     
}

 
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- Surv((dt1$TRTDURE), dt1$e)  
  
  ex<-survreg(S ~ GASTRDIS_STATUS, dt1 , dist="exponential")
  exp(-coefficients(ex)[1])  # hazard treat
  exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz placebo
  exp(-coefficients(ex)[1] + -coefficients(ex)[3]) #haz placebo
  exp(-coefficients(ex)[1] + -coefficients(ex)[4]) #haz placebo
  require(rms)
  f2 <- npsurv(S ~ GASTRDIS_STATUS, dt1)
      
  n <- f2$numevents
  x <- f2$exposure/365.25 # change to years rather than days
  tapply(dt1$TRTDURE, dt1$GASTRDIS_STATUS,sum,na.rm=TRUE)
  h <- n/x #f2$numevents/f2$exposure
  namez <- gsub(".*=","",names(h))
   namez=c("Never","Current","Prior","Prior/current")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)  #8.5
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
  addText(rtf, " \n", bold=F, italic=F)
  
   r=7
  r2=3
  
   i=2
 
  addText(rtf, paste(text3,namez[1],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),          text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)    
  
    i=1
  addText(rtf, paste(text3,namez[2],"=",round(n[[i]],r),  ",exposure in years=",round(x[[i]],0),          text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
 
  i=3
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),        text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  
  i=4
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ",exposure in years=",round(x[[i]],0),        text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
       addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years). ")
              addText(rtf, "The lower table depicts the number of patients at risk and the cumulative number of events in each stratum.\n\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                              ",status,   sep=" ") )
  done(rtf) 

}


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
id<-4
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, var="SEX",
              text1=namex2, 
              text2=title2 , 
              text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')
          
 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ GASTRDIS_STATUS, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.5, 1)))


```
  
## f_4.2-2.5 Kaplan Meier plot of patient infections and infestation candida infection free probability (All patients)


```{r  delete10, eval=FALSE,echo=TRUE}

title2 <- "Kaplan Meier plot of patient infections and infestation candida infection free probability (All patients)"
title1 <-'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x" ,var="x") {
  
  # Plot function  
  
  dt1 <- data
  
  myplot <- function() {
  
  d2 <- dt1
  
  #S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
   
  fit <- do.call(survfit,
                 list(formula = S ~ INFECT_STATUS, data = d2))  #eval(parse(text=var))
  
            print(addmargins(table(d2$INFECT_STATUS, d2$Treatment)))
 
    print(addmargins(table(d2$INFECT_STATUS, d2$e)))
  

  
    # n<-fit$numevents
    # x<-fit$exposure
    # h<-fit$numevents/fit$exposure
  
    gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better
    #readability 
    #legend.labs=c("NEVER","CURRENT","PRIOR","PRIOR/CURRENT"),
    legend.labs=c("Never","Current","Prior","Prior/current"),
    ylim = c(.75, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 7),
    palette = "jco",
        legend.title = "Infections and infestation"
    )
     

 
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
#  gg <- gg +  
#  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), #log(400),log(1000),log(2000) ),
 # labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
  cat("times that n risk and n of events are displayed\n")
 print(unique(gg$data.survtable$time)) 
 print(gg$data.survtable)
     cat("\n")     
 
}

 
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- Surv((dt1$TRTDURE), dt1$e)  
  
  ex<-survreg(S ~ INFECT_STATUS, dt1 , dist="exponential")
  exp(-coefficients(ex)[1])  # hazard treat
  exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz placebo
  exp(-coefficients(ex)[1] + -coefficients(ex)[3]) #haz placebo
  exp(-coefficients(ex)[1] + -coefficients(ex)[4]) #haz placebo
  require(rms)
  f2 <- npsurv(S ~ INFECT_STATUS, dt1)
      
  n <- f2$numevents
  x <- f2$exposure/365.25 # change to years rather than days
  tapply(dt1$TRTDURE, dt1$INFECT_STATUS,sum,na.rm=TRUE)
  h <- n/x #f2$numevents/f2$exposure
  namez <- gsub(".*=","",names(h))
  namez=c("Never","Current","Prior","Prior/current")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)  #8.5
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
  addText(rtf, " \n", bold=F, italic=F)
  
   r=7
  r2=3
  
   i=2
 
  addText(rtf, paste(text3,namez[1],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),    text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)    
  
    i=1
  addText(rtf, paste(text3,namez[2],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
 
  i=3
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  
  i=4
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
       addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years). ")
              addText(rtf, "The lower table depicts the number of patients at risk and the cumulative number of events in each stratum.\n\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                              ",status,   sep=" ") )
  done(rtf) 

}


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
id<-5
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, var="SEX",
              text1=namex2, 
              text2=title2 , 
              text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')
          

 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ INFECT_STATUS, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.5, 1)))




```
  
## f_4.2-2.6 Kaplan Meier plot of patient glucocorticosteroid use candida infection free probability (All patients)


```{r delete11, eval=FALSE,echo=TRUE}

title2 <- "Kaplan Meier plot of patient glucocorticosteroid use candida infection free probability (All patients)"
title1 <- 'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x" ,var="x") {
  
  # Plot function  
  
  dt1 <- data
  
  dt1$predictor <- dt1$H02A #!

  
  myplot <- function() {
  
      
  d2 <- dt1
  
  print(addmargins(table(d2$predictor, d2$Treatment)))
  
  
  #S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
   
  fit <- do.call(survfit,
                 list(formula = S ~ predictor, data = d2))  #eval(parse(text=var))
  
          print(addmargins(table(d2$H02A, d2$Treatment)))
 
    print(addmargins(table(d2$H02A, d2$e)))
  
  
    # n<-fit$numevents
    # x<-fit$exposure
    # h<-fit$numevents/fit$exposure
  
    gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better
    #readability 
    #legend.labs=c("NEVER","CURRENT","PRIOR","PRIOR/CURRENT"),
    legend.labs=c("Never","Current","Prior","Prior/current"),
    ylim = c(.75, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 7),
    palette = "jco",
        legend.title = "Glucocorticosteroid use"
    )
     

 
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
#  gg <- gg +  
#  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), #log(400),log(1000),log(2000) ),
 # labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
  cat("times that n risk and n of events are displayed\n")
 print(unique(gg$data.survtable$time)) 
 print(gg$data.survtable)
     cat("\n")     
}

 
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- Surv((dt1$TRTDURE), dt1$e)  
  
  ex<-survreg(S ~ predictor, dt1 , dist="exponential")
  exp(-coefficients(ex)[1])  # hazard treat
  exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz placebo
  exp(-coefficients(ex)[1] + -coefficients(ex)[3]) #haz placebo
  exp(-coefficients(ex)[1] + -coefficients(ex)[4]) #haz placebo
  require(rms)
  f2 <- npsurv(S ~ predictor, dt1)
      
  n <- f2$numevents
  x <- f2$exposure/365.25 # change to years rather than days
  tapply(dt1$TRTDURE, dt1$predictor,sum,na.rm=TRUE)
  h <- n/x #f2$numevents/f2$exposure
  namez <- gsub(".*=","",names(h))
  namez=c("Never","Current","Prior","Prior/current")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)  #8.5
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
  addText(rtf, " \n", bold=F, italic=F)
  
   r=7
  r2=3
  
   i=2
 
  addText(rtf, paste(text3,namez[1],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)    
  
    i=1
  addText(rtf, paste(text3,namez[2],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
 
  i=3
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  
  i=4
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
       addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years). ")
              addText(rtf, "The lower table depicts the number of patients at risk and the cumulative number of events in each stratum.\n\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                              ",status,   sep=" ") )
  done(rtf) 

}


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
id<-6
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, var="SEX",
              text1=namex2, 
              text2=title2 , 
              text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')
          


 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ H02A, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.5, 1)))

```

## f_4.2-2.7 Kaplan Meier plot of patient antibiotics use (ATC code: J01) candida infection free probability (All patients)


```{r delete12, eval=FALSE,echo=TRUE}

# new KMs request 17 Jan:

# 1.	Antibiotics (ATC code: J01)
# 2.	Proton pump inhibitors (ATC code: A02BC)
# 3.	Gastroesophageal reflux disease
# 4.	Diabetes mellitus (please combine here Diabetes mellitus AND Type 1 diabetes mellitus AND type 2 diabetes mellitus)

# Again for the categories current / prior / never. While it might make sense for diabetes to use only current / never or yes / no. I believe for those 4 we might also see tendencies.


title2 <- "Kaplan Meier plot of patient antibiotic use candida infection free probability (All patients)"
title1 <- 'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x" ,var="x") {
  
  # Plot function  
  
  dt1 <- data
  
  dt1$predictor <- dt1$J01
  
  myplot <- function() {
  
  d2 <- dt1
  
    print(addmargins(table(d2$predictor, d2$Treatment)))
 
    print(addmargins(table(d2$predictor, d2$e)))
  
  #S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
   
  fit <- do.call(survfit,
                 list(formula = S ~ predictor, data = d2))  #eval(parse(text=var))
  
    # n<-fit$numevents
    # x<-fit$exposure
    # h<-fit$numevents/fit$exposure
  
    gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better
    #readability 
    #legend.labs=c("NEVER","CURRENT","PRIOR","PRIOR/CURRENT"),
    legend.labs=c("Never","Prior","Prior/Current","Current"),
    ylim = c(.75, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 7),
    palette = "jco",
        legend.title = "antibiotic use"
    )
     

 
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
#  gg <- gg +  
#  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), #log(400),log(1000),log(2000) ),
 # labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
  cat("times that n risk and n of events are displayed\n")
 print(unique(gg$data.survtable$time)) 
 print(gg$data.survtable)
     cat("\n")     
}

 
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- Surv((dt1$TRTDURE), dt1$e)  
  
  ex<-survreg(S ~ predictor, dt1 , dist="exponential")
  exp(-coefficients(ex)[1])  # hazard prior
  exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz prior current
  exp(-coefficients(ex)[1] + -coefficients(ex)[3]) #haz current
  exp(-coefficients(ex)[1] + -coefficients(ex)[4]) #haz never
  require(rms)
  f2 <- npsurv(S ~ predictor, dt1)
      
  n <- f2$numevents
  x <- f2$exposure/365.25 # change to years rather than days
  tapply(dt1$TRTDURE, dt1$predictor,sum,na.rm=TRUE)
  h <- n/x #f2$numevents/f2$exposure
  namez <- gsub(".*=","",names(h))
  namez=c("Prior","Prior/current","Current","Never")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)  #8.5
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
  addText(rtf, " \n", bold=F, italic=F)
  
   r=7
  r2=3
  
   i=2
 
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)    
  
    i=1
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
 
  i=3
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  
  i=4
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
       addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years). ")
              addText(rtf, "The lower table depicts the number of patients at risk and the cumulative number of events in each stratum.\n\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                              ",status,   sep=" ") )
  done(rtf) 

}


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
id<-7
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, var="SEX",
              text1=namex2, 
              text2=title2 , 
              text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')

 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ J01, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.75, 1)))

```  

## f_4.2-2.8 Kaplan Meier plot of patient proton pump inhibitor use (ATC code: A02BC) candida infection free probability (All patients)


```{r delete13, eval=FALSE,echo=TRUE}

# new KMs request 17 Jan:

# 1.	Antibiotics (ATC code: J01)
# 2.	Proton pump inhibitors (ATC code: A02BC)
# 3.	Gastroesophageal reflux disease
# 4.	Diabetes mellitus (please combine here Diabetes mellitus AND Type 1 diabetes mellitus AND type 2 diabetes mellitus)

# Again for the categories current / prior / never. While it might make sense for diabetes to use only current / never or yes / no. I believe for those 4 we might also see tendencies.


title2 <- "Kaplan Meier plot of patient proton pump inhibitor use candida infection free probability (All patients)"
title1 <- 'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x" ,var="x") {
  
  # Plot function  
  
  dt1 <- data
  
  dt1$predictor <- dt1$A02BC
  
  myplot <- function() {
  
  d2 <- dt1
  
    print(addmargins(table(d2$predictor, d2$Treatment)))
 
    print(addmargins(table(d2$predictor, d2$e)))
  
  
  #S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
   
  fit <- do.call(survfit,
                 list(formula = S ~ predictor, data = d2))  #eval(parse(text=var))
  
    # n<-fit$numevents
    # x<-fit$exposure
    # h<-fit$numevents/fit$exposure
  
    gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better
    #readability 
    #legend.labs=c("NEVER","CURRENT","PRIOR","PRIOR/CURRENT"),
    legend.labs=c("Never","Prior","Prior/current","Current"),
    ylim = c(.75, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 7),
    palette = "jco",
        legend.title = "Proton pump inhibitor use"
    )
     

 
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
#  gg <- gg +  
#  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), #log(400),log(1000),log(2000) ),
 # labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
  cat("times that n risk and n of events are displayed\n")
 print(unique(gg$data.survtable$time)) 
 print(gg$data.survtable)
     cat("\n")     
}

 
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- Surv((dt1$TRTDURE), dt1$e)  
  
  ex<-survreg(S ~ predictor, dt1 , dist="exponential")
  exp(-coefficients(ex)[1])  # hazard treat
  exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz placebo
  exp(-coefficients(ex)[1] + -coefficients(ex)[3]) #haz placebo
  exp(-coefficients(ex)[1] + -coefficients(ex)[4]) #haz placebo
  require(rms)
  f2 <- npsurv(S ~ predictor, dt1)
      
  n <- f2$numevents
  x <- f2$exposure/365.25 # change to years rather than days
  tapply(dt1$TRTDURE, dt1$predictor,sum,na.rm=TRUE)
  h <- n/x #f2$numevents/f2$exposure
  namez <- gsub(".*=","",names(h))
  namez=c("Prior","Prior/current","Current","Never")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)  #8.5
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
  addText(rtf, " \n", bold=F, italic=F)
  
   r=7
  r2=3
  
   i=1
 
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)    
  
    i=2
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
 
  i=3
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  
  i=4
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
       addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years). ")
              addText(rtf, "The lower table depicts the number of patients at risk and the cumulative number of events in each stratum.\n\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                              ",status,   sep=" ") )
  done(rtf) 

}


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
id<-8
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, var="SEX",
              text1=namex2, 
              text2=title2 , 
              text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')


 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ A02BC, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.5, 1)))

```  


## f_4.2-2.9 Kaplan Meier plot of patient gastroesophageal reflux disease (ATC code: GRD_STATUS) candida infection free probability (All patients)


```{r delete14, eval=FALSE,echo=TRUE}

# new KMs request 17 Jan:

# 1.	Antibiotics (ATC code: J01)
# 2.	Proton pump inhibitors (ATC code: A02BC)
# 3.	Gastroesophageal reflux disease
# 4.	Diabetes mellitus (please combine here Diabetes mellitus AND Type 1 diabetes mellitus AND type 2 diabetes mellitus)

# Again for the categories current / prior / never. While it might make sense for diabetes to use only current / never or yes / no. I believe for those 4 we might also see tendencies.


title2 <- "Kaplan Meier plot of gastroesophageal reflux disease candida infection free probability (All patients)"
title1 <- 'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x" ,var="x") {
  
  # Plot function  
  
   dt1 <- data
  
  dt1$predictor <- dt1$GRD_STATUS
  
  myplot <- function() {
  
  d2 <- dt1
  
#thisnext tabulation matches the outputs provided(at leastfor treated)  
# medhis 4131
##Gastrooesophageal reflux disease 1 ( 0.7) 6 ( 0.2) 7 ( 0.3)
#current 4132
#Gastrooesophageal reflux disease 11 ( 8.1) 148 ( 5.6) 159 ( 5.7)
    
    print(addmargins(table(d2$predictor, d2$Treatment)))
 
    print(addmargins(table(d2$predictor, d2$e)))
  
  
  #S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
   # 3 observations deleted due to missingness
  fit <- do.call(survfit,
                 list(formula = S ~ GRD_STATUS, data = d2))  #eval(parse(text=var))
  
    # n<-fit$numevents
    # x<-fit$exposure
    # h<-fit$numevents/fit$exposure
  
    gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better
    #readability 
    #legend.labs=c("NEVER","CURRENT","PRIOR","PRIOR/CURRENT"),
    legend.labs=c("Current","Never","Prior"),
    #legend.labs=c("Never","Current","Prior","Prior/current"),
    ylim = c(.5, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 7),
    palette = "jco",
        legend.title = "Gastroesophageal reflux disease"  ###GRD_STATUS
    )
     

 
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
#  gg <- gg +  
#  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), #log(400),log(1000),log(2000) ),
 # labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
  cat("times that n risk and n of events are displayed\n")
 print(unique(gg$data.survtable$time)) 
 print(gg$data.survtable)
     cat("\n")     
}

 
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- Surv((dt1$TRTDURE), dt1$e)  
  
  ex<-survreg(S ~ predictor, dt1 , dist="exponential")
  exp(-coefficients(ex)[1])  # hazard treat
  exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz placebo
  exp(-coefficients(ex)[1] + -coefficients(ex)[3]) #haz placebo
 # exp(-coefficients(ex)[1] + -coefficients(ex)[4]) #haz placebo
  require(rms)
  f2 <- npsurv(S ~ predictor, dt1)
      
  n <- f2$numevents
  x <- f2$exposure/365.25 # change to years rather than days
  tapply(dt1$TRTDURE, dt1$predictor,sum,na.rm=TRUE)
  h <- n/x #f2$numevents/f2$exposure
  namez <- gsub(".*=","",names(h))
  #namez=c("Never","Current","Prior","Prior/current")
  namez=c("Current","Never","Prior")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)  #8.5
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
  addText(rtf, " \n", bold=F, italic=F)
  
   r=7
  r2=3
  
   i=1
 
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)    
  
    i=2
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
 
  i=3
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  
  i=4
 # addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
       addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years). ")
              addText(rtf, "The lower table depicts the number of patients at risk and the cumulative number of events in each stratum.\n\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                              ",status,   sep=" ") )
  done(rtf) 

}


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
id<-9
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, var="SEX",
              text1=namex2, 
              text2=title2 , 
              text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')


 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ GRD_STATUS, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.5, 1)))

  
``` 


## f_4.2-2.10 Kaplan Meier plot of patient diabetes status candida infection free probability (All patients)


```{r delete15, eval=FALSE,echo=TRUE}

# new KMs request 17 Jan:

# 4.	Diabetes mellitus (please combine here Diabetes mellitus AND Type 1 diabetes mellitus AND type 2 diabetes mellitus)

title2 <- "Kaplan Meier plot of diabetes status (both type I and II) candida infection free probability (All patients)"
title1 <- 'CAIN457 Post hoc analysis and regression modelling of risk factors for candida infection during secukinumab therapy'


surv.plot <- function(data, figname="x", projcode="x", text1="x", text2="x", text3="x", text4="x", text5="x", text6="x" ,var="x") {
  
  # Plot function  
  
   dt1 <- data
  
  dt1$predictor <- dt1$diab_status
  
  myplot <- function() {
  
  d2 <- dt1
  
#thisnext tabulation matches the outputs provided(at leastfor treated)  
# medhis 4131
##Gastrooesophageal reflux disease 1 ( 0.7) 6 ( 0.2) 7 ( 0.3)
#current 4132
#Gastrooesophageal reflux disease 11 ( 8.1) 148 ( 5.6) 159 ( 5.7)
    
    print(addmargins(table(d2$predictor, d2$Treatment)))
 
    print(addmargins(table(d2$predictor, d2$e)))
  
  
  #S <- Surv(log(d2$TRTDURE), d2$e)  #log base e time there are no zeros so I do not add 1

  S <- Surv(d2$TRTDURE, d2$e)  #natural log time there are no zeros so I do not add 1
  
   # 3 observations deleted due to missingness
  fit <- do.call(survfit,
                 list(formula = S ~ diab_status, data = d2))  #eval(parse(text=var))
  
    # n<-fit$numevents
    # x<-fit$exposure
    # h<-fit$numevents/fit$exposure
  
    gg<- ggsurvplot(fit,        
    pval = FALSE,             #displays p-value of log-rank test 
    conf.int = TRUE,          #plots a confidence interval for each curve
    ylab = "Candida infection \n free probability",
    xlab = "Time in days",
    #break.x.by = 10, 
    #break.time.by = 250,     # break X axis in time intervals by 100.
    ggtheme = theme_light(),  # customize theme with a grid for better
    #readability 
    #legend.labs=c("NEVER","CURRENT","PRIOR","PRIOR/CURRENT"),
    legend.labs=c("No diabetes","Diabetes"),
    #legend.labs=c("Never","Current","Prior","Prior/current"),
    ylim = c(.75, 1),
    risk.table = "nrisk_cumevents",#abs_pct",   # absolute number and percentage at risk
    #tables.theme = theme_survminer(font.main = 8),
    risk.table.y.text.col = FALSE,# colour risk table text annotations
    risk.table.y.text = TRUE,# show bars instead of names in legend of risk table.
    risk.table.height =.25,       
    fontsize = 3, #3
    ncensor.plot = FALSE#,      # plot the number of censored subjects at time t
    #surv.median.line = "hv"  
    ,tables.col = "strata", 
    tables.theme = theme_cleantable(font.main = 7),
    palette = "jco",
        legend.title = "Diabetes status"  ###GRD_STATUS
    )
     

 
       
  #gg <- gg +  scale_x_continuous(
  #breaks=c(log(5),log(10),log(25),log(50),log(100),log(200), log(400),log(1000),log(2000)), 
  #labels=c(5,     10,     25,     50,     100,     200,      400,     1000,     2000))       
  
#  gg <- gg +  
#  scale_x_continuous(breaks=c(0, log(3),log(10),log(25),log(50),log(100),log(200), #log(400),log(1000),log(2000) ),
 # labels=c(0, 3, 10, 25, 50, 100, 200, 400,1000,2000 ))
     
  
     gg$table <- gg$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=8)
    # axis.text.x = element_text(face="bold", color="#993333", 
    #                        size=14, angle=45)
    ) 

      

 print(gg)
     
  cat("times that n risk and n of events are displayed\n")
 print(unique(gg$data.survtable$time)) 
 print(gg$data.survtable)
     cat("\n")     
}

 
  #####################################

  ### RTF function calls plot function

  # run the simple analysis on untransforemed time
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- Surv((dt1$TRTDURE), dt1$e)  
  
  ex<-survreg(S ~ predictor, dt1 , dist="exponential")
  exp(-coefficients(ex)[1])  # hazard treat
  exp(-coefficients(ex)[1] + -coefficients(ex)[2]) #haz placebo
 # exp(-coefficients(ex)[1] + -coefficients(ex)[3]) #haz placebo
 # exp(-coefficients(ex)[1] + -coefficients(ex)[4]) #haz placebo
  require(rms)
  f2 <- npsurv(S ~ predictor, dt1)
      
  n <- f2$numevents
  x <- f2$exposure/365.25 # change to years rather than days
  tapply(dt1$TRTDURE, dt1$predictor,sum,na.rm=TRUE)
  h <- n/x #f2$numevents/f2$exposure
  namez <- gsub(".*=","",names(h))
  #namez=c("Never","Current","Prior","Prior/current")
  namez=c("No diabetes","Diabetes")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(rtf)
  rtf <- RTF(paste(path.export,figname, sep=''), font.size=11, omi=c(1,1,1,1), width=11, height=9)  #8.5
  addText(rtf, paste(projcode,"\n"), bold=F, italic=F) 
  addText(rtf, paste("Figure",text1, "(Page 1 of 1) \n"), bold=F, italic=F)
  addText(rtf, paste(text2,"\n"), bold=F, italic=F)  
  addText(rtf, " \n", bold=F, italic=F)
  
   r=7
  r2=3
  
   i=1
 
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)    
  
    i=2
  addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),      text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)
 
  # i=3
  # addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  
  # i=4
 # addText(rtf, paste(text3,namez[i],"=",round(n[[i]],r),  ", exposure in years=",round(x[[i]],0),     text4,round(h[[i]],r2), "\n",    sep=" "), bold=F, italic=F)  

  addPlot(rtf, plot.fun=myplot, width=8, height=4.5, res=rez) 
  #addText(rtf, "- The x-axis is presented using a log transformed scale. The y-axis is truncated to show the probability between 0.7 and 1.0\n")
 # addText(rtf, "- ")
    addText(rtf, "- The vertical lines on the curves denote censored patients. Coloured bands around each curve denote 95% confidence intervals.\n")
       addText(rtf, "- Hazard is the exponential distribution hazard rate estimate, the number of events divided by exposure (the sum of all failure and censoring times in years). ")
              addText(rtf, "The lower table depicts the number of patients at risk and the cumulative number of events in each stratum.\n\n")
      #addText(rtf, paste0("Source:",path.script,namex3,"           ",format(Sys.time(), "%b %d %Y, %H:%M:%S"), sep=" "))
  addText(rtf, paste0(path.short,namex3,"  ",format(Sys.time(), "%d%b%Y:%H:%M:%S"),
"                                                                              ",status,   sep=" ") )
  done(rtf) 

}


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  execute RTF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
id<-10
file.id <- paste0("4.2-2.",id)

namex1 <- paste0("f_",file.id,".rtf")
namex2 <- file.id
namex3 <- paste0("f_",file.id,".R ")  # space after R for better looking output 

#assign("mdata",get(paste0("d", id)))
assign("mdata",get(paste0("d")))

surv.plot(data=mdata, figname =namex1, projcode=title1, var="Diabetes",
              text1=namex2, 
              text2=title2 , 
              text3='Number of events',text4=', hazard=', 
              text5='Number of events', text6='')

 #quick check
 fit <- (do.call(survfit,
                 list(formula = S ~ diab_status, data = mdata)))
  print(ggsurvplot(fit,  ylim = c(.5, 1)))

``` 

\clearpage


## COMPUTING ENVIRONMENT

```{r computing envirnoment, echo=FALSE}

      options(width=70)
     # opts_knit$set(root.dir = wd.code)   ##THIS SETS YOUR WORKING DIRECTORY
      sessionInfo()
      print(getwd())
      stopTime<-proc.time()

```

This took `r (stopTime-startTime)[1][[1]]` seconds to execute.

```{r saving code and output, echo=TRUE, eval=TRUE, cache=FALSE}

```     

### Using knitr:purl to create R code
    
```{r knitr purl code}

    # move stangle R file to a folder in GPS
    # put this at bottom and give it the same name as the RMD file , replace any blanks with underscore
    # https://amywhiteheadresearch.wordpress.com/2014/11/12/copying-files-with-r/

    # https://stackoverflow.com/questions/21101573/need-the-filename-of-the-rmd-when-knitr-runs
    x <- knitr::current_input() 
    x
    rcode <-  gsub(' ','_', trimws(x))                       # replace blank with underscore, this is needed
    file.copy(rcode, path.script,  overwrite=TRUE)           # make a copy of the rcode in a folder of choice
    
    
    # saving the html to GPS!
    # https://stackoverflow.com/questions/28894515/rmarkdown-directing-output-file-into-a-directory
    # x1 <- gsub('\\.Rmd','\\.html', trimws(x)) 
    # rcode <-  gsub(' ','_', trimws(x1))   
    # knit: (function(inputFile, encoding) { 
    #   out_dir <- path.script;
    #   rmarkdown::render(inputFile,
    #                     encoding=encoding, 
    #                     output_file=file.path(dirname(x), out_dir, rcode)) })
    # 
    
    #http://felixfan.github.io/extract-r-code/
    #https://www.rdocumentation.org/packages/knitr/versions/1.19/topics/knit
    #y <- knitr::current_input()
    #knitr::purl(y, output = "test.R", documentation = 2)
 
    input  = knitr::current_input()  # filename of input document
    cat("checking purl..")
    input
    output = paste(tools::file_path_sans_ext(input), 'R', sep = '.')
    output
    #tools::file_path_sans_ext(input)
    knitr::purl(input,paste0(path.script,output),documentation=1,quiet=T)  ##where do we wish to save R code?
    getwd()
    
```     

### Using purl to create R code
    
```{r not used}

#    https://stackoverflow.com/questions/36868287/purl-within-knit-duplicate-label-error
#    This chunk automatically generates a text .R version of this script when     running within knitr.
#    input  = knitr::current_input()  # filename of input document
#    output = paste(tools::file_path_sans_ext(input), 'R', sep = '.')
#    knitr::purl(input,output,documentation=1,quiet=T)

```
    
### Using purl to create R code seperately for each chunk
    
```{r chunk output code,echo=TRUE, results='asis', eval=TRUE}

# https://stackoverflow.com/questions/25800604/how-to-purl-each-chunks-in-rmd-file-to-multiple-r-files-using-knitr

    library("knitr")
    # p <- purl("test.Rmd")
    p <-purl(knitr::current_input())
    read_chunk(p)
    chunks <- knitr:::knit_code$get()
    
    invisible(mapply(function(chunk, name) {
        writeLines(c(paste0("## ----",name,"----"), chunk), paste0(path.script,"",name,".R"))
    }, chunks, names(chunks)))
    unlink(p)                   # delete the original purl script
    knitr:::knit_code$restore() # remove chunks from current knitr session


 setwd(wd.code)





