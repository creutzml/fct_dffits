#This is the Code for |DFBETAS0| Simulation and does 1 iteration


##Can manually change things like BasisNum, the B0fun, B1fun, 
##B1funout, functions, and length
library(fda)
library(MASS)
library(corpcor) 

sin_a=function()sample(c(-3,-2,-1,0,1,2,3),1)
sin_k=function()sample(c(-300,-200,-100,100,200,300),1)
sin_d=function()sample(c(-100,-50,0,50,100),1)
sin_c=function()sample(c(-3,-2,-1,0,1,2,3),1)

cos_a=function()sample(c(-3,-2,-1,0,1,2,3),1)
cos_k=function()sample(c(-300,-200,-100,100,200,300),1)
cos_d=function()sample(c(-100,-50,0,50,100),1)
cos_c=function()sample(c(-3,-2,-1,0,1,2,3),1)

#used to generate X data
mu2t<-function(t){ 
  1/12*t + (sin_a()*sin(1/sin_k()*(t-sin_d())) + 
              sin_c())*(cos_a()*cos(1/cos_k()*(t-cos_d()))+cos_c())
}


OUonFloodEM<-function(sigma=1,drift=1,times,endtime,Data){
  # want this to be around .5 for most situations to keep it stable
  dt<-endtime/times; 
  #number of flood events I have
  objects<-ncol(Data) 
  #this is randn from matlab code
  normalvec<-matrix(rnorm(objects*times,sd=1),
                    nrow=objects,ncol=times); 
  #null matrix of 0s
  ouvalue<-matrix(0:0,nrow=objects,ncol=times) 
  
  # Defining the cluster sizes and the true clustering structure: 
  # No clusters here so just 1
  timevec<-seq(0,endtime,length=times)
  
  # loop to create approximation of Ornstein-Uhlenbeck process	
  # Using Euler-Maruyama approximation
  ouvalue[1:objects,1]=sigma*normalvec[,1] #So we don't start at 0.
  for (timeint in 2:times) {
    ouvalue[1:objects,timeint]<- ouvalue[1:objects,timeint-1] - 
      drift*ouvalue[1:objects,timeint-1]*dt + 
      sigma*normalvec[1:objects,timeint]*sqrt(dt)
    # normalvec[1:clussizes[1],timeint] n random values from 
    # normal(0, sigma)
  }
  # makes sure it is the same size as the vector it is being applied
  samesizeOU<-matrix(0:0,nrow=nrow(Data),ncol=objects) 
  for(i in 1:objects){
    samesizeOU[,i]<-approx(x = 1:ncol(ouvalue),
                           y = ouvalue[i,], 
                           n = nrow(Data))$y
  }
  # Defining and naming the observed curves
  obscurves<-Data+samesizeOU;
  out<-list(1:nrow(Data),obscurves,samesizeOU)
  names(out)<-c("timevec","obscurves","samesizeOU")
  return(out)
}

#Need this to be established for the above function to run
B0fun<-function(t) cos(t/200)+2
B0fun1<-function(t) t/5000
B1fun<-function(t) sin(t/200)+2
B1fun1<-function(t) t/500
B1funout<-function(t,lambda) lambda*sin(t/200)+2
B1funout1<-function(t,lambda) (lambda-1) + t/500
B1funout2<-function(t,lambda) (lambda-1)*(t > 500) + t/500
length = 1000
t=1:length #This is the length

#Function to return some measures from the concurrent model
BetaOnlyBSpline<-function(xData,yData,BasisNum=15, 
                          BetaOnly=FALSE,lambda=10^(-1)) {
  n<-nrow(xData)
  nPred<-ncol(xData) #Number of predictors (ie flood events)
  
  gaittime <- seq(1:n)
  gaitrange <- c(0,n)
  gaitfine = seq(0,n,.2)
  
  mygaitExp <- array(NA, dim = c(n,nPred,2))
  mygaitExp[1:n, ,] <- seq(1:n)
  
  for(i in 1:nPred){
    mygaitExp[,i, 1] <- xData[,i]
    mygaitExp[,i, 2] <- yData[,i]
  }
  
  #original 20 leaving norder at 6.
  gaitbasis <- create.bspline.basis(gaitrange, 
                                    nbasis=BasisNum,
                                    norder=6) 
  D2fdPar = fdPar(gaitbasis, lambda=lambda)
  gaitSmooth = smooth.basis(gaittime, mygaitExp, D2fdPar)
  betalist  = list(const=D2fdPar, cong=D2fdPar)
  gaitfd = gaitSmooth$fd
  
  names(gaitfd$fdnames) = c("Normalized time", "Event", "Height")
  gaitfd$fdnames[[3]] = c("Cong", "Cedar")
  
  congfd  = gaitfd[,1] #Predictor Stuff
  cedfd = gaitfd[,2] # Response Stuff
  
  xfdlist   = list(const=rep(1,nPred), cong=congfd)
  
  gaitRegress= fRegress(cedfd, xfdlist, betalist)
  betaestlist = gaitRegress$betaestlist
  cedIntercept = predict(betaestlist$const$fd, gaitfine) #B0
  congCoef = predict(betaestlist$cong$fd, gaitfine) #Slope Term B1
  
  Intercept<-approx(x=1:length(cedIntercept), y=cedIntercept, n = n)$y
  Slope<-approx(x=1:length(congCoef), y=congCoef , n = n)$y
  
  if(BetaOnly==FALSE){
    #$fd #THIS IS WHERE THE ERROR OCCURS unless using fda version 5.1.5.1
    cedhatfd = gaitRegress$yhatfd 
    #1 to 2000 compare this to yHat.t.mat
    cedhatmat = eval.fd(gaittime, cedhatfd)  
    #The 2 represents cedar (response)
    resmat. = mygaitExp[,,2] - cedhatmat 
    SigmaE = cov(t(resmat.))
    
    cedfinemat   = eval.fd(gaitfine, cedfd)
    cedmeanvec   = eval.fd(gaitfine, mean.fd(cedfd))
    cedhatfinemat= eval.fd(gaitfine, cedhatfd)
    #finer grid than the . needed for computation below
    resmat        = cedfinemat - cedhatfinemat 
    #Says the number of flood events
    ncurve        = dim(mygaitExp)[2] 
    resmat0 = cedfinemat - cedmeanvec %*% matrix(1,1,ncurve)
    SSE0 = apply((resmat0)^2, 1, sum)
    SSE1 = apply(resmat^2, 1, sum)
    ced.R2 = (SSE0-SSE1)/SSE0
    #And this is what I had plotted before exactly.
    ylim2=c(0, max(congCoef, ced.R2)) 
    gaitbasismat = eval.basis(gaitfine, gaitbasis)
    #Needed to get error
    y2cMap = gaitSmooth$y2cMap 
    #$fd #
    cedhatfd = gaitRegress$yhatfd 
    #1 to 2000 compare this to yHat.t.mat
    cedhatmat = eval.fd(gaittime, cedhatfd)  
    #The 2 represents cedar (response)
    resmat. = mygaitExp[,,2] - cedhatmat 
    SigmaE = cov(t(resmat.))
    
    fRegressList1 = fRegress(cedfd, xfdlist, betalist,
                             y2cMap=y2cMap, SigmaE=SigmaE)
    
    fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
    betastderrlist = fRegressList2$betastderrlist
    
    error.B0<-predict(betastderrlist[[1]], gaitfine)
    error.B1<-predict(betastderrlist[[2]], gaitfine)
    #This gives the standard errors for B0t and B1t to use for DFBETAS
    stderr.B0<-approx(x=1:length(error.B0), y=error.B0, n = n)$y
    stderr.B1<-approx(x=1:length(error.B1), y=error.B1, n = n)$y
    
    
    #Adjust to be nRow Ncol
    meanHat_i<-c()
    #adds the column of 1's
    data1<-cbind(rep(1,nPred),t(xData)) 
    Hatiit<-matrix(NA,nrow=nPred,ncol=n)
    for(t in 2:(n+1)){
      X<-data1[,c(1,t)]
      Hat<-X%*%solve(t(X)%*%X)%*%t(X)
      #Each row is an event. The corresponding column is the Hii of 
      # that event at each timepoint
      Hatiit[,t-1]<-diag(Hat) 
    }
    #will return an average Hii for each of the 10 sampled curve
    meanHat_i<-c(meanHat_i,apply(Hatiit,1,mean))
    
    out<-list(Intercept,Slope,stderr.B0,stderr.B1,resmat.,
              Hatiit,meanHat_i,cedhatmat)
    names(out)<-c("B0t","B1t","stderr.B0t","stderr.B1t","Residuals",
                  "AdjHatMatrix","MeanHat","FullyHat")
  }
  else {out<-list(Intercept,Slope)
  names(out)=c("B0t","B1t")
  }
  return(out)
}



##Calculates dffits for each t and each observation
BestDFFITS_fun<-function(Xtmat,Ytmat,BasisNum = 15){
  FullResults<-BetaOnlyBSpline(Xtmat,Ytmat,BasisNum=BasisNum)
  #Then use the above function to get Betas with things left out.
  B0_i_<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  B1_i_<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  stderrB1_i_2<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  stderrB0_i_2<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  BasisToUse<-c()
  for(i in 1:ncol(Xtmat)){
    #calculate B0t and B1t without each i.
    tempBetas<-BetaOnlyBSpline(Xtmat[,-i],Ytmat[,-i],BasisNum=BasisNum) 
    B0_i_[,i]<-tempBetas$B0t
    B1_i_[,i]<-tempBetas$B1t
    stderrB0_i_2[,i]<-tempBetas$stderr.B0t
    stderrB1_i_2[,i]<-tempBetas$stderr.B1t
  }
  yHat_mini<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  for(i in 1:ncol(Xtmat)){
    yHat_mini[,i]<-B0_i_[,i]+B1_i_[,i]*Xtmat[,i]
  }
  # DFFIT<-FullResults$FullyHat-yHat_mini
  MSE_mini_t<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  yHatsNoi<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat)-1)
  n=1:ncol(Xtmat)
  for(i in n){
    k=1
    #go through all observations except the ith
    for(j in n[-i]){ 
      #Gets yHat for all N-1 observations with i removed
      yHatsNoi[,k]<-B0_i_[,i]+B1_i_[,i]*Xtmat[,j] 
      k=k+1
    }
    #Get an MSE across the other events when event i is left out
    MSE_mini_t[,i]<-apply((Ytmat[,-i]-yHatsNoi)^2,1,mean) 
  }
  altDFFITS<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  for(t in 1:nrow(Xtmat)){ 
    #at each timepoint calculate equation from pg 396/401 of kutner
    altDFFITS[t,]<- (
      FullResults$Residuals[t,]/
        sqrt(MSE_mini_t[t,]*(1-FullResults$AdjHatMatrix[,t]))
      )*sqrt(FullResults$AdjHatMatrix[,t]/
               (1-FullResults$AdjHatMatrix[,t]))
  }
  out<-list(altDFFITS)
  names(out)<-c("altdffits")
  return(out)
}

#Generates Functional Data with an outlier
GenerateFunctionalDataOut1<-function(N = 100,length=1000,
                                     lambda=1, n_inf = 1) {
  
  NewXout<-matrix(NA, nrow = length,ncol = N)
  t<-1:length
  for(i in 1:N){
    NewXout[,i]<-mu2t(t)
  }
  #Runs the functions created at the top of the page 
  B0t<-B0fun(t)
  B1t<-B1fun(t)
  B1funout<-B1funout(t,lambda)
  #add an outlier
  outlier_obs<-sample(1:N, n_inf)
  holdY<-NewXout ##just sets the dimensions
  for(i in 1:N){
    if(!(i %in% outlier_obs)){
      holdY[,i]<-NewXout[,i]*B1t+B0t
    }
    else holdY[,i]<-NewXout[,i]*B1funout+B0t
  }
  
  NewYout<-holdY #just gives it same dimensions before editing
  run2<-OUonFloodEM(sigma=15,drift=.5,times= 2*length, 
                    endtime = length,Data=as.matrix(holdY))
  for(i in 1:N){
    NewYout[,i]<-supsmu(x=run2$timevec,y=run2$obscurves[,i])$y
  }
  out<-list(NewXout, NewYout,outlier_obs)
  names(out)<-c("Xdata","Ydata","OutlierNumber")
  return(out)
}

#Generates Functional Data with an outlier
GenerateFunctionalDataOut2<-function(N = 100,length=1000,
                                     lambda=1, n_inf = 1) {
  NewXout<-matrix(NA, nrow = length,ncol = N)
  t<-1:length
  for(i in 1:N){
    NewXout[,i]<-mu2t(t)
  }
  
  #Runs the functions created at the top of the page 
  B0t<-B0fun1(t)
  B1t<-B1fun1(t)
  B1funout<-B1funout1(t,lambda)
  #add an outlier
  outlier_obs<-sample(1:N, n_inf)
  ##just sets the dimensions
  holdY<-NewXout 
  for(i in 1:N){
    if(!(i %in% outlier_obs)){
      holdY[,i]<-NewXout[,i]*B1t+B0t
    }
    else holdY[,i]<-NewXout[,i]*B1funout+B0t
  }
  
  #just gives it same dimensions before editing
  NewYout<-holdY 
  run2<-OUonFloodEM(sigma=15,drift=.5,times= 2*length, 
                    endtime = length,Data=as.matrix(holdY))
  for(i in 1:N){
    NewYout[,i]<-supsmu(x=run2$timevec,y=run2$obscurves[,i])$y
  }
  out<-list(NewXout, NewYout,outlier_obs)
  names(out)<-c("Xdata","Ydata","OutlierNumber")
  return(out)
}


#Generates Functional Data with an outlier
GenerateFunctionalDataOut3<-function(N = 100,length=1000,
                                     lambda=1, n_inf = 1) {
  NewXout<-matrix(NA, nrow = length,ncol = N)
  t<-1:length
  for(i in 1:N){
    NewXout[,i]<-mu2t(t)
  }
  #Runs the functions created at the top of the page 
  B0t<-B0fun1(t)
  B1t<-B1fun1(t)
  B1funout<-B1funout2(t,lambda)
  #add an outlier
  outlier_obs<-sample(1:N, n_inf)
  holdY<-NewXout ##just sets the dimensions
  for(i in 1:N){
    if(!(i %in% outlier_obs)){
      holdY[,i]<-NewXout[,i]*B1t+B0t
    }
    else holdY[,i]<-NewXout[,i]*B1funout+B0t
  }
  
  #just gives it same dimensions before editing
  NewYout<-holdY 
  run2<-OUonFloodEM(sigma=15,drift=.5,times= 2*length, 
                    endtime = length,Data=as.matrix(holdY))
  for(i in 1:N){
    NewYout[,i]<-supsmu(x=run2$timevec,y=run2$obscurves[,i])$y
  }
  out<-list(NewXout, NewYout,outlier_obs)
  names(out)<-c("Xdata","Ydata","OutlierNumber")
  return(out)
}


#generate data before using this function. 
# This function bootstraps given data
#Give it the outliar number from the data
# (MLC 1/2/24): changed function so that alpha in Pittman is 0.5, one
# that performs well for small sample sizes. Then, alpha will 
# correspond to the percentile, rather than their parameter
absdffit_bootsamedata<-function(alpha, 
                                xdata = sim_data$Xdata, 
                                outlier = sim_data$OutlierNumber,
                                ydata = sim_data$Ydata, 
                                nboot = 100) {
  
  #Say the number of observations we have
  N = ncol(xdata) 
  
  #Now run the concurrent model and calculate mean |dffits| for each
  # of the observation
  #pointwise DFFITS for each observation
  raw_dffits<-BestDFFITS_fun(Xtmat = xdata,
                             Ytmat = ydata,
                             BasisNum = 15) 
  mean_abs_dffits<-apply(abs(raw_dffits$altdffits),2,mean)
  adjusted_obs_absdffits<-mean_abs_dffits[outlier]
  theta<-as.vector((1/mean_abs_dffits)^0/ # this is alpha_b
                     sum((1/mean_abs_dffits)^0))
  #Need to set up the bootstrap next SEE OrdinaryRegressionDistribution 
  # for Details
  #Set blank objects here
  boot_meanabsdffits<-c()
  for(i in 1:nboot){
    #resamples the needed data with the given probability 
    index_samp<-sample(1:N,N,replace=TRUE,prob = theta) 
    
    # Error with negative eigenvalues on fRegress routine
    temp_dffits <- tryCatch(BestDFFITS_fun(Xtmat = xdata[,index_samp],
                                           Ytmat = ydata[,index_samp],
                                           BasisNum = 15)$altdffits, 
                            error = function(e) e
    )
    
    # Skip if an error was caught
    if(inherits(temp_dffits, "error")) next
    
    # Otherwise save result
    boot_meanabsdffits<-c(boot_meanabsdffits,
                          apply(abs(temp_dffits),2,mean))
  }
  #returns percentiles from the N*nboot |dffits| calculated
  boot_percentiles<-quantile(boot_meanabsdffits,1-alpha)
  percentile <- ecdf(boot_meanabsdffits)
  adjusted_point_percentile<- percentile(mean_abs_dffits[outlier])
  out<-list(boot_percentiles,
            adjusted_point_percentile,
            mean_abs_dffits, 
            raw_dffits)
  #And returns the percentile of this list that the observed |dffits| 
  #falls into
  names(out)<-c("BootPercentiles",
                "AdjPercentile",
                "AllObsMeanAbsDFFITS",
                "RawDFFITS") 
  return(out)
}
# 
# ##Can adjust this so I input more into the function like alpha, lambda 
# ## etc
# one_fulltable_bootstrap<-function(N = 100, nboot = 100, alpha, 
#                                   lambda, length = 1000, 
#                                   data_n_inf = 1) {
#   
#   #Will store the 95%ile from the bootstrap results
#   results_output95<-matrix(NA,nrow = length(lambda), 
#                            ncol = length(alpha)) 
#   #Will store the percentile that the adjusted observation falls into 
#   results_outputpercentile<-matrix(NA,nrow = length(lambda), 
#                                    ncol = length(alpha)) 
#   #Will store the percentile that the adjusted observation falls into 
#   results_numberobs_over95<-matrix(NA,nrow = length(lambda), 
#                                    ncol = length(alpha)) 
# 
# 
#   colnames(results_output95)<-as.character(alpha)
#   rownames(results_output95)<-as.character(lambda)
#   colnames(results_outputpercentile)<-as.character(alpha)
#   rownames(results_outputpercentile)<-as.character(lambda)
#   colnames(results_numberobs_over95)<-as.character(alpha)
#   rownames(results_numberobs_over95)<-as.character(lambda)
#   
#   ##Below will do 1 iteration of each cell
#   i = 1 #starts filling in column 1
#   for(l in lambda){
#     #Generate new data for each lambda. then keep same for different
#     #alpha functions are defined at the top
#     B0fun<-function(t) cos(t/200)+2
#     B1fun<-function(t) sin(t/200)+2
#     B1funout<-function(t,lambda) l*sin(t/200)+2
#     t=1:length
#     sim_data<-GenerateFunctionalDataOut(
#       N=N,length=length,lambda=l, n_inf = data_n_inf
#     )
#     j = 1 #start on row 1 each time
#     for(a in alpha){
#       tempout<-absdffit_bootsamedata(
#         xdata = sim_data$Xdata, 
#         ydata = sim_data$Ydata, 
#         outlier = sim_data$OutlierNumber, 
#         alpha = a,  
#         nboot = nboot
#       )
#       results_output95[i,j]<- tempout$Boot95
#       results_outputpercentile[i,j]<- tempout$AdjPercentile
#       results_numberobs_over95[i,j] <-sum(
#         ifelse(tempout$AllObsMeanAbsDFFITS>=tempout$Boot95,1,0)
#       )
#       #results_numberobs_over95[i,j]   count number from tempout$raw
#       # that are over tempout$boot95
#       j=j+1 #move to the next column which is a new Alpha
#     }
#     i=i+1 #moves to the next column 
#   }
#   
#   out<-list(results_output95,
#             results_outputpercentile,
#             results_numberobs_over95)
#   names(out)<-c("the95percentile", 
#                 "outlierpercentile",
#                 "numberobs_over95")
#   return(out)
# }

##Run that function and repeat "iterations" time and get the mean of each cell
# alpha<-c(0,.1,.3,.5,.75)
# lambda<-c(.5,.75,.9,1,1.1,1.25,1.5,1.75,2)
# s_time <- Sys.time()
 # hold<-one_fulltable_bootstrap(N=10, nboot = 10, alpha = .1, lambda = 1.5, length = 1000, data_n_inf = 2)
# e_time <- Sys.time()
# paste0(round(as.numeric(difftime(time1 = e_time, time2 = s_time, units = "secs")), 3), " Seconds")
   #hold$the95percentile
#hold$outlierpercentile
#  Adjobs_above_95<-ifelse(hold$outlierpercentile>=0.95,1,0)
# 
# #Want to save avg_95 and avg_percentile
#   filename95=paste("~/InfluenceMeasures/data/New_avg_95_",iter,".RData",sep="")
#   x<-hold$the95percentile
#   save(x, file=filename95)
# 
#   filenamepercentile=paste("~/InfluenceMeasures/data/New_avg_adjpercentile_",iter,".RData",sep="")
#   y<-hold$outlierpercentile
#   save(y, file=filenamepercentile)
#   
#   filenameOver95=paste("~/InfluenceMeasures/data/New_avg_numberobs_over95_",iter,".RData",sep="")
#   z<-hold$numberobs_over95
#   save(z, file=filenameOver95)
#   
#   filenameAdjpercentileover95=paste("~/InfluenceMeasures/data/New_Above_95_",iter,".RData",sep="")
#   save(Adjobs_above_95, file=filenameAdjpercentileover95)




### Figure 4.2 in Chapter 4, Creutzinger (2024)
### Create an example plot of the three models
#####################################################################
newData1 <- GenerateFunctionalDataOut1(lambda = 2, N = 30)
newData2 <- GenerateFunctionalDataOut2(lambda = 2, N = 30)
newData3 <- GenerateFunctionalDataOut3(lambda = 2, N = 30)
# par(mfrow = c(1,3))
# with(test_data1, {
#   matplot(Ydata[,-OutlierNumber], type = "l", col = "gray",
#           main = "Model 1", xlab = "Sampling Point (t)",
#           ylab = "Y(t)")
#   matlines(Ydata[,OutlierNumber], col = "red")
# })
# with(test_data2, {
#   matplot(Ydata[,-OutlierNumber], type = "l", col = "gray",
#           main = "Model 2", xlab = "Sampling Point (t)",
#           ylab = "Y(t)")
#   matlines(Ydata[,OutlierNumber], col = "red")
# })
# with(test_data3, {
#   matplot(Ydata[,-OutlierNumber], type = "l", col = "gray",
#           main = "Model 3", xlab = "Sampling Point (t)",
#           ylab = "Y(t)")
#   matlines(Ydata[,OutlierNumber], col = "red")
# })
# 
# matplot(test_data1$Xdata, type = "l")

# Plot of the three models with influential marked as well
y_gg_plot_data1 <- newData1$Ydata %>%
  as.data.frame() %>%
  dplyr::mutate(sp = 1:1000) %>%
  tidyr::pivot_longer(-sp,
                      names_to = "Observation",
                      values_to = "Y(t)") %>%
  dplyr::mutate(
    inf_ob = ifelse(
      as.numeric(substring(Observation, 2))==newData1$OutlierNumber,
      "Influential",
      "Not"),
    Model = 1)

y_gg_plot_data2 <- newData2$Ydata %>%
  as.data.frame() %>%
  dplyr::mutate(sp = 1:1000) %>%
  tidyr::pivot_longer(-sp,
                      names_to = "Observation",
                      values_to = "Y(t)") %>%
  dplyr::mutate(
    inf_ob = ifelse(
      as.numeric(substring(Observation, 2))==newData2$OutlierNumber,
      "Influential",
      "Not"),
    Model = 2)

y_gg_plot_data3 <- newData3$Ydata %>%
  as.data.frame() %>%
  dplyr::mutate(sp = 1:1000) %>%
  tidyr::pivot_longer(-sp,
                      names_to = "Observation",
                      values_to = "Y(t)") %>%
  dplyr::mutate(
    inf_ob = ifelse(
      as.numeric(substring(Observation, 2))==newData3$OutlierNumber,
      "Influential",
      "Not"),
    Model = 3)

y_gg_plot_data <- y_gg_plot_data1 %>%
  dplyr::bind_rows(y_gg_plot_data2) %>%
  dplyr::bind_rows(y_gg_plot_data3) %>%
  dplyr::mutate(inf_ob = factor(inf_ob,
                                levels = c("Not", "Influential")))

# Plot
y_data_plot <- ggplot() +
  geom_line(aes(x = sp, y = `Y(t)`, group = Observation),
            color = "lightgray", alpha = 0.5,
            data = y_gg_plot_data) +
  geom_line(aes(x = sp, y = `Y(t)`, group = Observation),
            color = "#E69F00", size = 1.1,
            linetype = "longdash",
            data = y_gg_plot_data %>%
              filter(inf_ob == "Influential")) +
  facet_grid(cols = vars(Model),
             labeller = label_bquote(cols = "Model "*.(Model))) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = unit(c(.5, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(-25, 200)) +
  labs(x = "Sampling Point (t)")


## Next, plot the dffits from those same models
fReg_dffits1 <- fRegress_concurrent(y_mat = newData1$Ydata, 
                                    x_array = newData1$Xdata)$dffits_mat

dffits_gg_plot_data1 <- fReg_dffits1 %>%
  as.data.frame() %>%
  dplyr::mutate(sp = 1:1000) %>%
  tidyr::pivot_longer(-sp,
                      names_to = "Observation",
                      values_to = "Y(t)") %>%
  dplyr::mutate(
    inf_ob = ifelse(
      as.numeric(substring(Observation, 2))==newData1$OutlierNumber,
      "Influential",
      "Not"),
    Model = 1)

fReg_dffits2 <- fRegress_concurrent(y_mat = newData2$Ydata, 
                                    x_array = newData2$Xdata)$dffits_mat

dffits_gg_plot_data2 <- fReg_dffits2 %>%
  as.data.frame() %>%
  dplyr::mutate(sp = 1:1000) %>%
  tidyr::pivot_longer(-sp,
                      names_to = "Observation",
                      values_to = "Y(t)") %>%
  dplyr::mutate(
    inf_ob = ifelse(
      as.numeric(substring(Observation, 2))==newData2$OutlierNumber,
      "Influential",
      "Not"),
    Model = 2)

fReg_dffits3 <- fRegress_concurrent(y_mat = newData3$Ydata, 
                                    x_array = newData3$Xdata)$dffits_mat

dffits_gg_plot_data3 <- fReg_dffits3 %>%
  as.data.frame() %>%
  dplyr::mutate(sp = 1:1000) %>%
  tidyr::pivot_longer(-sp,
                      names_to = "Observation",
                      values_to = "Y(t)") %>%
  dplyr::mutate(
    inf_ob = ifelse(
      as.numeric(substring(Observation, 2))==newData3$OutlierNumber,
      "Influential",
      "Not"),
    Model = 3)

dffits_gg_plot_data <- dffits_gg_plot_data1 %>%
  dplyr::bind_rows(dffits_gg_plot_data2) %>%
  dplyr::bind_rows(dffits_gg_plot_data3) %>%
  dplyr::mutate(inf_ob = factor(inf_ob,
                                levels = c("Not", "Influential")))

# Plot
dffits_data_plot <- ggplot() +
  geom_line(aes(x = sp, y = `Y(t)`, group = Observation),
            color = "lightgray", alpha = 0.5,
            data = dffits_gg_plot_data) +
  geom_line(aes(x = sp, y = `Y(t)`, group = Observation),
            color = "#E69F00", size = 1.1,
            linetype = "longdash",
            data = dffits_gg_plot_data %>%
              filter(inf_ob == "Influential")) +
  facet_grid(cols = vars(Model),
             labeller = label_bquote(cols = "Model "*.(Model))) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = unit(c(.5, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Sampling Point (t)", 
       y = "DFFITS(t)")

library(cowplot)
plot_grid(y_data_plot, dffits_data_plot,
          align = "v", nrow = 2, labels = c("(a)", "(b)"))

# # Plot of X(t) sample
# x_gg_plot_data <- test_data1$Xdata %>%
#   as.data.frame() %>%
#   dplyr::mutate(sp = 1:1000) %>%
#   tidyr::pivot_longer(-sp, 
#                       names_to = "Observation", 
#                       values_to = "X(t)")
# 
# ggplot(data = x_gg_plot_data) +
#   geom_line(aes(x = sp, y = `X(t)`, group = Observation)) +
#   scale_x_continuous(expand = c(0,0)) +
#   theme_bw() + 
#   theme(text = element_text(size = 16), 
#         plot.margin = unit(c(.5, 1, 1, 1), "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   labs(x = "Sampling Point (t)")
