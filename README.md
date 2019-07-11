In any war, it is always of value to one side to have good intelligence on the weapons resources of the other side. During the Second World War, for example, Allied military planners eagerly searched for ways to accurately estimate the Axis production of tanks, aircrafts and numerous other weapons platforms. In the specific case of German tanks, a very clever way to do that was based on using either the stamped serial numbers or gearbox markings on captured Mark I and Mark V tanks, respectively. This type of problem has far wider applications.

We can experimentally see how well the formula that was used works in practice with a Monte Carlo Simulation. That is, the program first randomly picks a value for N (integer between 100 and 1000) with a value for the the sample size as a percentage of N.

The program then generates, randomly, n different integers in the interval 1 to N, the maximum of those integers is then used in the estimation formula to estimate a value for N. This estimate can be compared to the actual value of N to determine how well the formula has performed. We investigate samples that are 2%, 5%, 10% and 20% of the population.



```

# Bayesian Proportional odds model , run once output saved

```{r , eval=FALSE}

# error when running stan
# had to backtrack as StanHeaders was ahead of rstan!
# from Rstudio help page
packageurl <- "http://cran.r-project.org/src/contrib/Archive/StanHeaders/StanHeaders_2.15.0-1.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

# information on k https://stat.ethz.ch/pipermail/r-help/2007-October/143569.html


  options(mc.cores=parallel::detectCores ())
  library(brms)

  # ~ takes about 10 minutes
  po1 <- brm(dlqi ~ s(avisit, k=3) + uas + vas + (avisit|usubjid), 
            data=both, family=cumulative("logit"), 
            threshold="flexible")
  
  # don't assume vas is linear
  po1a <- brm(dlqi ~ s(avisit, k=3) + uas + s(vas) + (avisit|usubjid), 
            data=both, family=cumulative("logit"), 
            threshold="flexible")
  


  # takes longer than 10 minutes
  # random intercept and slopes? - UAS and VAS are not random so not this?
  # https://stats.stackexchange.com/questions/62033/how-to-specify-uncorrelated-random-slopes-in-lmer-syntax
  po2 <- brm(dlqi ~ s(avisit, k=3) + uas + s(vas) + (uas+ vas+ avisit|usubjid), 
            data=both, family=cumulative("logit"), 
            threshold="flexible")

  # https://stackoverflow.com/questions/40056566/mgcv-how-to-set-number-and-or-locations-of-knots-for-splines
## fitting natural cubic splines with notes at quantiles
## adapt delta increased to 0.99 default 0.8 will be slower to sample
  po3 <- brm(dlqi ~ s(avisit, bs = 'cr', k = 3) + uas + s(vas, bs = 'cr', k = 5) + (avisit|usubjid), 
            data=both, family=cumulative("logit"), seed=123, control = list(adapt_delta = .99),
            threshold="flexible")
  


  


  o.all.loo <- LOO(po1,po1a, po3, compare=FALSE)
#        LOOIC     SE
#po1  16076.87 112.86
#po1a 15950.27 114.59
#po3  15964.34 114.51


######################################################################################
  ## priors investigation


        prior1 <- get_prior(dlqi ~ s(avisit, bs = 'cr', k = 3) + 
                              uas + s(vas, bs = 'cr', k = 5) + 
                              (avisit|usubjid),     
                            data=both, family=cumulative("logit"), 
            threshold="flexible" )

      #

        m <- tapply(both$dlqi, both$avisit, mean)[1][[1]]  # mean dlqi at baseline
        s <-  tapply(both$dlqi, both$avisit, sd)[1][[1]]  # sd dlqi at baseline

        test <- paste("normal(",as.character(m),
                ",",as.character(s), ")",
                sep="")
 

        # normally improper flat prior over the reals for fixed effects, see doc 
        prior1 <- c(
        prior_string(test, class = 'Intercept', coef = ''),  #'test' using data to inform prior
        prior_(~cauchy(0,1),  coef = ~savisit_1 ),  
        prior_(~cauchy(0,1),  coef = ~svas_1 )  ,
        prior_(~cauchy(0,1),  coef = ~uas1 )  ,
        prior_(~cauchy(0,1),  coef = ~uas2 )  ,
        prior_(~cauchy(0,1),  coef = ~uas3 )  ,
        prior_(~cauchy(0,1),  coef = ~uas4 )  ,
        prior_(~cauchy(0,1),  coef = ~uas5 )  ,
        prior_(~cauchy(0,1),  coef = ~uas6 )  ,
        prior_(~cauchy(0,1), class=~sd, group = ~usubjid, coef = ~avisit),
        prior_(~cauchy(0,1), class=~sd, group = ~usubjid, coef = ~Intercept)
        )
            
        prior1$prior[14:20] <-  "student_t(10, 0, 5)"    

          make_stancode(dlqi ~ s(avisit, bs = 'cr', k = 3) + uas + 
                          s(vas, bs = 'cr', k = 5) + (avisit|usubjid), 
            data=both, family=cumulative("logit"), prior=prior1,
            threshold="flexible")  
        
            prior.test.1 <- brm(dlqi ~ s(avisit, bs = 'cr', k = 3) + 
                                  uas + s(vas, bs = 'cr', k = 5) + (avisit|usubjid), 
            data=both, family=cumulative("logit"),  seed=123, prior=prior1,
            threshold="flexible")




#         
#         set_prior("normal(0,5)", class = "b", coef = "Intercept")
#         set_prior("normal(0,5)", class = "b", coef = "blast0.25")
#         set_prior("normal(0,5)", class = "b", coef = "enrich")
#         set_prior("normal(0,5)", class = "b", coef = "Intercept")
#         set_prior("normal(0,5)", class = "b", coef = "TrtU")
# 
   #https://rdrr.io/cran/brms/man/get_prior.html
       prior1$prior[14:20] <-  "student_t(10, 0, 5)"    
  
##################################################################################
# https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
        prior1 <- get_prior(dlqi ~ s(avisit, bs = 'cr', k = 3) + 
                              uas + s(vas, bs = 'cr', k = 5) + 
                              (avisit|usubjid),       data=both, family=cumulative("logit"), 
            threshold="flexible" )

#

        m <- tapply(both$dlqi, both$avisit, mean)[1][[1]]  # mean dlqi at baseline
        s <-  tapply(both$dlqi, both$avisit, sd)[1][[1]]  # sd dlqi at baseline

        test <- paste("normal(",as.character(m),
                ",",as.character(s), ")",
                sep="")
 

        # normally improper flat prior over the reals for fixed effects, see doc 

        prior1$prior[1:10] <-   "student_t(5, 0, 5)"     
        prior1$prior[14:20] <-  "student_t(5, 0, 5)" 
        prior1$prior[11] <-"lkj(4)"

         #check priors are indeed included in the model

          make_stancode(dlqi ~ s(avisit, bs = 'cr', k = 3) + uas + 
                          s(vas, bs = 'cr', k = 5) + (avisit|usubjid), 
            data=both, family=cumulative("logit"), prior=prior1,
            threshold="flexible")  
        
            prior.test.2 <- brm(dlqi ~ s(avisit, bs = 'cr', k = 3) + 
                                  uas + s(vas, bs = 'cr', k = 5) + (avisit|usubjid), 
            data=both, family=cumulative("logit"),  seed=123, prior=prior1,
            threshold="flexible")




#         
#         set_prior("normal(0,5)", class = "b", coef = "Intercept")
#         set_prior("normal(0,5)", class = "b", coef = "blast0.25")
#         set_prior("normal(0,5)", class = "b", coef = "enrich")
#         set_prior("normal(0,5)", class = "b", coef = "Intercept")
#         set_prior("normal(0,5)", class = "b", coef = "TrtU")
# 
   #https://rdrr.io/cran/brms/man/get_prior.html
  #     prior1$prior[14:20] <-  "student_t(10, 0, 5)"    
  
##################################################################################
    # frequentist run now andsave will look at the end of the program
    # subset as data to large
    rfoo <- subset(both, usubjid %in% sample(unique(usubjid ), 100))
    
    # brms will throw an error about not wnough levels if subjects sample are too small?
    newdata <- rfoo[order(rfoo$dlqi),] 
    newdata$dlqi2 <- cumsum(c(TRUE, newdata$dlqi[-1]!=newdata$dlqi[-length(newdata$dlqi)]))
    newdata$dlqi <- NULL
    names(newdata)[names(newdata)=="dlqi2"] <- "dlqi"
    newdata$dlqi <- newdata$dlqi-1 
    rfoo <- newdata
    newdata <- NULL


    # frequentist and bayesian comparison 

    require(ordinal)
    mm1 <- clmm(factor(dlqi) ~ avisit+ factor(uas) + vas + 
                  (avisit|usubjid), data = rfoo, link = "logit", 
                  threshold = "flexible")

    bb1 <- brm(dlqi ~ avisit + uas + vas + (avisit|usubjid), 
            data= rfoo, family=cumulative("logit"), 
            seed=123, control = list(adapt_delta = .99),
            threshold="flexible")
  
    # save the workspace
   save.image(file=paste(wd.data,"\\brms prop odds model.RData", sep=""))

  
 
```

### Load analysis

```{r}

  load(file=paste(wd.data,"\\brms prop odds model.RData", sep=""))
  library(brms)


```

### brms model output and plots

```{r}

  f <- po3 # selected model, default priors used,
  
  print(f)

  pp_check(f)  # shows dens_overlay plot by default
  pp_check(f, type = "error_hist", nsamples = 11)    # 11 draws
  pp_check(f, type = "scatter_avg", nsamples = 100)  # mean on x axis y observed
  #pp_check(f, type = "stat_2d")  ## took a long time so cancelled this
  #pp_check(f ,x="vas", type = "intervals") ## took a long time so cancelled this
  pp_check(f, type = "scatter", nsamples=2)
  
 # pp_check(f,x="avisit",type = "error_scatter_avg_vs_x", nsamples = 10)  #throwing error
 # pp_check(f,x="vas", type = "ribbon", nsamples = 20) #throwing error
  pp_check(f ,type = "error_scatter", nsamples = 6)




 # launch_shiny(f)


#  conditions <- data.frame(vas = c(50), avisit= c(9))
#  
#  plot(marginal_effects(f, effects = "uas", 
#                        conditions = conditions))
# 
# 
#  conditions <- data.frame(vas = c(50), avisit= c(0))
# 
#  plot(marginal_effects(f, effects = "uas", 
#                        conditions = conditions))

conditions <- data.frame(vas = c(50), avisit= 0, uas=3)

x1 <- (marginal_effects(f, method='predict', conditions=conditions))

plot(x1, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[1]] 
plot(x1, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[2]] 
plot(x1, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[3]] 

x2 <- (marginal_effects(f, method='fitted', conditions=conditions))

plot(x2, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[1]] 
plot(x2, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[2]] 
plot(x2, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[3]] 

cat("\nInclude all random effects\n")

x3 <- (marginal_effects(f, method='predict', conditions=conditions, re_formula=NULL ))

plot(x3, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[1]] 
plot(x3, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[2]] 
plot(x3, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[3]] 

cat("\nInclude no random effects\n")
x4 <- (marginal_effects(f, method='predict', conditions=conditions, re_formula=NA ))

plot(x4, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[1]] 
plot(x4, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[2]] 
plot(x4, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[3]] 

cat("\nFitted mean include all random effects\n")
x5 <- (marginal_effects(f, method='fitted', conditions=conditions, re_formula=NULL ))

plot(x5, plot = FALSE, point=TRUE, rug=T, jitter_width=.3, main="x")[[1]] 
plot(x5, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[2]] 
plot(x5, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[3]] 


par(mfrow=c(1,3))
cat("\nIndividual predictions include all random effects\n")
plot(x3, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[3]] ##indiv predictions all rand effects
cat("\nIndividual predictions no random effects\n")
plot(x4, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[3]] ##indiv predictions no rand effects
cat("\nMean predictions all random effects\n")
plot(x5, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[3]] ##mean all random effect
cat("\nMean predictions no random effects\n")
plot(x2, plot = FALSE, point=TRUE, rug=T, jitter_width=.3)[[3]] ##mean no random effects
par(mfrow=c(1,1))

```

#### Plots of predicted DLQI with UAS (Bayesian)
```{r}
  

 
      #both <- foo
      
      J <- sort(unique(both$avisit))
      k <- length(J) 
      
      
      for (i in 1:k)  {
        
        time=J[i] 
        var="uas"  # change to VAS
        
        # I look at effect of interest at time so other values than time dont matter
        conditions <- data.frame(vas = c(50), avisit= time, uas=3)
        
        x <- (marginal_effects(f, method='fitted', conditions=conditions))
#                              effects = "UAC"))
        
        #ribbon
        d <-(eval(parse(text=paste("x$",var, sep=""))))  # pull data FRAME
        hold <- unique(d$vas) 
        d <- d[,c("uas","estimate__","lower__","upper__")]
        
        d$uas <- as.numeric(levels(d$uas))[d$uas]

        #datapoints
        dx <- na.omit(foo[,c("avisit",var,"dlqi")])
        z <- dx[dx$avisit %in% time,]
        dx <- NULL
        n <- nrow(z)
        
        
        #so we have d and z 
        zz <- merge(d,z, all=T)
      
         px <- NULL
        p1a <- ggplot(z,aes(x=uas,y=dlqi) )+ 
          geom_point(colour="blue",  size=.3, #points from z data
                     position = position_jitter(w = .3, h = .2)) +
          geom_errorbar(data=d, 
                mapping=aes(y=estimate__ , ymin=lower__, ymax=upper__ ), 
                color = "red",
                width=0.2, lwd = .2) +
        
          EnvStats::stat_n_text(size = 3, y.pos =  32 , y.expand.factor=.0, 
                                angle = 0, hjust = 0.5, family = "mono", fontface = "bold")  +  
          scale_y_continuous(limits = c(-1,32)) +
          scale_x_continuous(breaks = c(unique(z$uas)),
                            labels = 
                             c(unique(z$uas)) )+
          
          theme(legend.position="none") +
          xlab("UAS total score") +
          ylab("DLQI total score")  
        
        px <-  p1a + 
        ggtitle(paste0("Plot of DLQI with UAS month "
                       ,time,", N=",n," subjects/data points\nPredicted 95% confidence  error bars for the predicted DLQI. Adj. to VAS=", hold))
        print(px)   
}
 
```

#### Now vas  (Bayesian)          

```{r}

      
      J <- sort(unique(both$avisit))
      k <- length(J) 
      
      
      for (i in 1:k)  {
        
        time=J[i] 
        var="vas"  # change to VAS
        
        # I look at effect of interest at time so other values than time dont matter
        conditions <- data.frame(vas = c(50), avisit= time, uas=3)
        
        x <- (marginal_effects(f, method='fitted', conditions=conditions))
#                              effects = "UAC"))
        
        #ribbon
        d <-(eval(parse(text=paste("x$",var, sep=""))))  # pull data FRAME
        hold <- unique(d$uas) 
        d <- d[,c("vas","estimate__","lower__","upper__")]
        
        #datapoints
        dx <- na.omit(foo[,c("avisit",var,"dlqi")])
        z <- dx[dx$avisit %in% time,]
        dx <- NULL
        n <- nrow(z)
        
        z$dlqi <- as.numeric(as.character(z$dlqi))
        #so we have d and z 
        zz <- merge(d,z, all=T)
      
         px <- NULL
        p1a <- ggplot(z,aes(x=vas,y=dlqi) )+ 
          geom_point(colour="blue",  size=.3, #points from z data
                     position = position_jitter(w = .3, h = .2)) +
  
     geom_line(data=d,     aes(x=vas,y=estimate__), colour="blue") + #pred from d data
            
    geom_ribbon(data=d, aes(x=vas, y=estimate__, ymin=lower__, ymax=upper__,
                                     alpha=0.2,fill="purple") )  +
    
       scale_y_continuous(limits = c(-1,32)) +
          #scale_x_continuous(breaks = c(unique(z$uas)),
                    #         labels = 
                      #         c(unique(z$uas)) )+
          
          theme(legend.position="none") +
          xlab("UAS total score") +
          ylab("DLQI total score")  
        
        px <-  p1a + 
        ggtitle(paste0("Plot of DLQI with UAS month "
                       ,time,", N=",n," subjects/data points\nPredicted 95% confidence  error bars for the predicted DLQI. Adj. to UAS=", hold))
        print(px)   
}

```

#### Print PO model regression tables

```{r}

cat("\nfinal model\n")
f
cat("\nsensitivity analysis to differnet priors\n")
prior.test.1
cat("\nsensitivity analysis to differnet priors\n")
prior.test.2


```

#### Print Bayesian PO model and frequentist PO model using 100 subjects selected randomly

```{r}

cat("\nBayesian\n")
bb1
cat("\nFrequentist\n")
mm1$coef


```
#### Frequentist, not working!

```{r, eval=FALSE}


    ## see earlier where ordinal package is used. 
 
    # uas 0 cumulative probs
    uas0 <- exp(mm1$alpha ) / (1 + exp(mm1$alpha ))   # intercepts only
   
    n=7; m = 28
    o <- array(NA, c(n,m))

    for( i in 0:6 ) {
    
       if (i %in% 0 ) {
            tmp <- (exp(mm1$alpha ) / (1 + exp(mm1$alpha )))
            o[i+1,] <- c(tmp[1],diff(tmp))
            
       } else { 
            a <- mm1$beta[paste0("factor(uas)",i)]
            tmp <- (exp(a + mm1$alpha) / (1 + exp(a + mm1$alpha)))
            o[i+1,] <- c(tmp[1],diff(tmp))
    }}
    
    t(o)
  
    1-t(o)

```

#### ordinal regression no random effects (rms)

```{r , eval=FALSE}

# ordinal regression no random effects rms
fx <- both
fx$uas <- factor(fx$uas)
fh<- lrm(factor(dlqi) ~   uas + vas ,data = fx, x=T)  # basic model
plot(foo$dlqi ~ foo$uas)  # better plots prior to this
L <- predict(fh, se.fit=TRUE)

plogis(with(L, linear.predictors + 1.96*cbind(-se.fit,se.fit)))
predict(fh, type="fitted.ind")[1:10,]   #gets Prob(better) and all others
d <- data.frame(uas=c(0,1,2,3,4,5,6), vas=c(0,0,0,0,0,0,0))
predict(fh, d, type="fitted")        # Prob(Y>=j) for new observation
f  # Prob(Y=j)
predict(fh, d, type='mean', codes=TRUE) # predicts mean(y) using codes 1,2,3
m <- Mean(fh, codes=TRUE)
lp <- predict(fh, d)
m(lp)
# Can use function m as an argument to plot.Design or nomogram to
# get predicted means instead of log odds or probabilities
# Don't use non.slopes argument to plot.Design for this
dd <- datadist(foo); options(datadist='dd')
m
plot(fh, vas=NA, fun=m, ylab='Predicted Mean')  ##

# look at fh coefficients,
plot(Predict(fh, uas, vas=0))
(Predict(fh, uas, vas=0))  # intercept, intercept+uas=1, intercept + uas=2....


summary(fh, uas=0)
plot(summary(fh, uas=0), log=TRUE)
plot(Predict(fh, vas)) 
ggplot(Predict(fh, vas, fun=m)) 
ggplot(Predict(fh, uas))
bg <- ggplot(Predict(fh, uas, fun=m) )
bg
 
# For any level of DLQI, the estimated odds that a 
# response is in the higher direction rather than the lower direction 
# is about x times the odds for lower UAS .
# eg x is 3.13 comparin UAS 1 to 0.
  
i <- predict(fh, d, type="fitted.ind")  
plot(i[,1]~c(1:7))  #prob of DLQI 1 when 



contrast(fh,list(uas=0),list(uas=5))

plot(Predict(fh, vas, fun=m, kint=10))  ##?
plot(Predict(fh, vas, fun=m, kint=2))  ##?

bar <- Mean(fh)
plot(Predict(fh, fun=bar), ylab='Predicted Mean', kint=10)

bar <- Mean(fh)
plot(Predict(fh, fun=bar), ylab='Predicted Mean', kint=1)

fit <- update(fh)   # make new reference value take effect
plot(Predict(fit, uas, ref.zero=FALSE, fun=m))
plot(Predict(fit, uas, ref.zero=TRUE , fun=m))


```

# manual probability calculations and cross checking results from rms functions

```{r}
     
    cat("\nShow understanding of the intercepts in PO model\n")
     # Model without random effects using all the data
    fx <- both
    fx$uas <- factor(fx$uas)
    fh<- lrm(factor(dlqi) ~   uas + vas ,data = fx, x=T)  # basic model
    fh
    
    cat("\nPredict membership of each DLQI level if all variables are 0, using functions\n")

    cat("\nprobability of membership when vas = 0 and uas =0\n")
    d <- data.frame(uas=c(0), vas=c(0))       
    cat("\n cumulative probabilities\n")
    predict(fh, d, type="fitted")           
    cat("\n individual probabilities\n")
    predict(fh, d, type="fitted.ind")         
    sum(predict(fh, d, type="fitted.ind") )   # checking probabilities sum to 1


    cat("\nPredict membership of each DLQI level if all variables are 0, manually\n")
    cat("\nCummulative probabilities\n")

    cat("\nJust select the intercepts from the all the model coefficients\n")
    intercepts <- fh$coef[grepl(">=", names(fh$coef))]

    cat("\nconvert log odds to cum probabilities for the intercepts of cumulative probabilities\n")
    exp(intercepts)/(1+exp(intercepts))       
  
    cat("\nIndividual probabilities\n")
    lo <- intercepts                               # grab intercepts
    p <- exp(lo)/(1+exp(lo))                              # create probabilities

    p0 <- 1 - p[1] # 1-(exp(3.3363064)/(1+exp(3.3363064)))# prob of membership to dlqi=0
    p30 <- p[30]                                          # prob of membership to dlqi=30
    p1.29 <- abs(diff(p[1:30]))                           # prob of membership to dlqi=1,2,3....29
    (prob <- c(p0, p1.29, p30))
    sum(prob)                                             # check my calculation sum to 1
    predict(fh, d, type="fitted.ind") 

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cat("\nNow do the same for the Bayesian model\n")
    cat("\nrms written in terms of p(y>=j), Bayesian (Py<=j)\n")

    cat("\nextract the intercept coefficients mcmcs\n")
   samples <- posterior_samples(f,  pars = "b")
   #names(samples)
    
   cat("\nlook at a couple of intercepts\n")
   logit <- samples[,"b_Intercept[30]"]
   pr <- 1-(exp(logit) / (1+exp(logit)))     # SUBTRACT FROM 1 HERE TO MATCH RMS PARAMETERIZATION
   summary(pr)
 
   logit <- samples[,"b_Intercept[7]"]
   pr <- 1-(exp(logit) / (1+exp(logit)))
   summary(pr)


   cat("\nLook at all intercepts, grep them all, so we have only intercept mcmc log odds estimates\n")
   intercepts <- samples[grepl("b_Intercept", names(samples))]
   logit <- intercepts
   cat("\ncheck this matches the model intercepts\n")
   apply(logit,2,mean)  


   cat("\nconvert to log odds to cumulative probability space, p(y>=j) \n")
   pr <- 1 - (exp(logit) / (1+exp(logit)))
   apply(pr,2,mean)                      # p(y>=j)
    
   cat("\nconvert cumulative probability to individual probability space\n") 
   cat("\nindividual probability of membership of dlqi=j \n") 

   p0 <- 1 - pr[1]    # prob of membership to dlqi=0
   p30 <- pr[30]      # prob of membership to dlqi=30
   p1.29 <- t(abs(apply(pr[1:30], 1,function(x) diff(x) )))  # individual probability 1:29

   cat("\nbind individual probabilities together\n") 
   prob <- cbind(p0, p1.29, p30)
   cat("\nCheck probabilties sum to 1\n") 
   head(apply(prob,1,sum)) ## it worked all sum to 1

   cat("\nnow we have the individual probability of membership of dlqi=j \n") 
   apply(prob,2,mean)        # individual probabilities when all covariate are 0
   sum(apply(prob,2,mean))   # sums to 1 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 cat("\nBayesian predict another way to look at the predictions?\n")
 cat("\nhttps://github.com/paul-buerkner/brms/issues/68\n")
 cat("\nhttps://github.com/paul-buerkner/brms/issues/82\n")

#                      uas = 3, vas=50, usubjid=NA)
#predict(f, newdata = newdata)
  cat("\nHow does this compare to previous section individual probabilities?\n")
  conditions <- data.frame(vas = c(0), avisit= 0, uas=0)
  cat("\nno random effects\n")
  list_of_data <- marginal_effects(f, method='fitted', conditions=conditions, re_formula=NA)
  head(list_of_data$uas)
  x <- predict(f, newdata = list_of_data$uas, summary = FALSE)
  head(x)
  cat("\nLook at the the first column that is UAS=0\n")
  table(x[,1])/length(x[,1])
 





  list_of_data <- marginal_effects(f, conditions= conditions)
  x <- predict(f, newdata = list_of_data$uas, summary = FALSE)
  table(x[,1])/length(x[,1])

```

#### dumping code that needs explanation

```{r, eval=FALSE}

# prob of identifying as DLQI 0 given VAS 0,1,2,3,4,5,6
#https://mvuorre.github.io/post/2017/bayes-factors-with-brms/

  x<-seq(1:6)
  est <- exp(0.57*x +-5.05)
  p <- est/(1+est);p
  
  
  
  # g <- unlist(mcmc[1:4,1,1])
  # 
  # g <- as.mcmc("b_Intercept[1]", pars = NA, exact_match = FALSE, combine_chains = TRUE, 
  #              inc_warmup = FALSE)
  # 
  # mcmc <- as.mcmc(ot1)
  # 
  # g <- unlist(mcmc[1,1:1782,1])
  # int <- (as.mcmc(ot1, pars = "b_Intercept", combine_chains = TRUE))
  # uas <-
  #   
  #   
    
    
   samples <- posterior_samples(f,  pars = "b")
   names(samples)
    
   logit0v1 <- samples[,"b_Intercept[1]"]
   uas <- samples[,"b_uas1"]
    
   mean(uas)
   mean(logit0v1)
   
   #################################################################################
   
   #what's the prob vas 1,2,3,4,5,6 is DLQI 0 V DLQI 1-30?
  
   cprob <- rep(NA, 6)
   
   for(j in 1:6) {
     
       br <-  uas*j+logit0v1      
       est <- exp(br)
       q1 <-  est/(1+est)
       cprob[j] <- mean(q1)
       
   }
   
   cprob
   diff(cprob)
   
#################################################################################
   
   #what's the prob vas0  is DLQI 0 V DLQI 1-30 etc ?
   
   intercepts <- posterior_samples(ot1,  pars = "b_Intercep")
   
   uas <- posterior_samples(ot1,  pars = "b_UAS")
   
   cprob <- rep(NA, 30)
   
   i=6
   for(j in 1:30) {
     
     br <-  uas*i+intercepts[,j]      
     est <- exp(br)
     q1 <-  est/(1+est)
     q1 <- as.vector(unlist(q1))
     cprob[j] <- mean(q1)
     
   }
   
   cprob
   p <- diff(cprob)
   p
   
   #################################################################################
   
   q3 <-  as.mcmc(f, pars = NA, exact_match = FALSE,
  combine_chains = TRUE, inc_warmup = FALSE)

 
```      

\clearpage

# COMPUTING ENVIRONMENT

```{r, echo=FALSE}

options(width=70)
opts_knit$set(root.dir = wd.code)   ##THIS SETS YOUR WORKING DIRECTORY
sessionInfo()
print(getwd())
stopTime<-proc.time()

```

This took `r (stopTime-startTime)[1][[1]]` seconds to execute.

```{r echo=FALSE}

# move stangle R file to a folder in GPS
# put this at bottom and give it the same name as the RMD file , replace any blanks with underscore
# https://amywhiteheadresearch.wordpress.com/2014/11/12/copying-files-with-r/

rcode <-  gsub(' ','_', trimws(namex3))                  # replace blank with underscore, this is needed
file.copy(rcode, path.script,  overwrite=TRUE)           # make a copy of the rcode in a folder of choice

```
