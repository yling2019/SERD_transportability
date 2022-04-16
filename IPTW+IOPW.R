### load libraries
library("survival")
library("survminer")
library("mice")

### logit function
logit <- function(x){
  return(log(x/(1-x)))
}

### function for Tipton index
### reference: Tipton, Elizabeth. "How generalizable is your experiment? An index for comparing 
### experimental samples and populations." Journal of Educational and Behavioral Statistics 39.6 (2014): 478-501.
Bindex <- function(dat1B,dat2B) {
  # bandwidth
  h = function(x){
    n = length(x) 
    return((4*sqrt(var(x))^5/(3*n))^(1/5))
  }#kernel density
  
  kg = function(x,data){
    hb = h(data) #bin width
    k = r = length(x)
    for(i in 1:k) r[i] = mean(dnorm((x[i]-data)/hb))/hb 
    return(r)
  }
  ##B index calculation
  return( as.numeric(integrate(function(x) sqrt(kg(x,dat1B)*kg(x,dat2B)),-Inf,Inf)$value))
}


# dataframe d contains the following columns:
# S: 1 for second-line patients and 0 first-line
# covariates: as outlined in the main text (some can be missing values)
# treatment: 1 for fulvestrant plus palbociclib and 0 letrozole plus palbociclib
# event: as defined in the main text
# timetoEvent: as defined in the main text 

### Step 1: IPTW ###
# determine m, the number of datasets to impute
d_iptw = d[which(d$S==1),]
m=floor((1-as.numeric(rownames(missing_pattern(d_iptw))[1])/nrow(d_iptw))*100)

# implement multiple imputation (MI) (MI-passive INT-within)
mi.output<-as.data.frame(matrix(data=NA, nrow=m, ncol=5))
colnames(mi.output)<-c("m", "coef", "coef.exp", "robust.se", "Schoenfeld")
mi.output$m<-1:m

imputed_Data <- mice(d_iptw[,which(colnames(d_iptw) %in% c(PSvariables, "treatment","event", "timeToEvent"))], m=m, maxit = 5, seed = 500, print=FALSE, remove_collinear = FALSE)
imp.data <- complete(imputed_Data, "long")

# MI-passive: conduct IPTW within each imputed dataset
for (i in 1:m){
  data<-imp.data[which(imp.data$.imp==i),]

  # fit PS model
  m_ps <- glm(as.formula(paste("treatment ~", paste(PSvariables, collapse="+"))), data=data, family = binomial)

  # estimate stablized weights
  data$PS<-m_ps$fitted
  prev<-length(which(data$treatment==1))/nrow(data)
  data$weight_iptw <- ifelse(data$treatment==1, prev/data$PS, (1-prev)/(1-data$PS))

  # fit Cox model with weights
  model<-coxph(formula = Surv(timeToEvent, event) ~ treatment, data = data, weights=data$weight_iptw, ties = "breslow", robust = TRUE)
  mi.output[i,2:4]<-summary(model)$coefficients[c(1,2,4)]

  # check Cox model assumptions
  test.ph <- cox.zph(model) 
  mi.output$Schoenfeld[i]<-test.ph$table[2,"p"]
}

# INT-within: pool treatment estimates and SE according to Rubin's Rule
mean.list<-mi.output[,2]
overall.mean<-mean(mean.list)
robust.overall.se<-sqrt((mean(mi.output$robust.se^2) + var(mean.list)*(1+1/m)))

estimate = overall.mean
SE = robust.overall.se
HR = exp(estimate)
LCL = exp(estimate - qnorm(0.975) * SE)
UCL = exp(estimate + qnorm(0.975) * SE)

### Step 2: IOPW ###
# determine m, the number of datasets to impute
m=floor((1-as.numeric(rownames(missing_pattern(d))[1])/nrow(d))*100)

# implement multiple imputation (MI) (MI-passive INT-within)
mi.output<-as.data.frame(matrix(data=NA, nrow=m, ncol=6))
colnames(mi.output)<-c("m", "coef", "coef.exp", "robust.se", "Tipton", "Schoenfeld")
mi.output$m<-1:m

imputed_Data <- mice(d[,which(colnames(d) %in% c(PSvariables, "treatment","event", "timeToEvent"))], m=m, maxit = 5, seed = 500, print=FALSE, remove_collinear = FALSE)
imp.data <- complete(imputed_Data, "long")

# MI-passive: conduct IOPW within each imputed dataset
for (i in 1:m){
    data=imp.data[imp.data$.imp == i,]
    data$index<-data$.id
    
    # get weights from IPTW
    d_iptw<-data[which(data$S==1),]
    m_ps1 <- glm(as.formula(paste("treatment ~", paste(PSvariables, collapse="+"))), data=d_iptw, family = binomial)
    d_iptw$PS<-m_ps1$fitted
    prev1<-length(which(d_iptw$treatment==1))/nrow(d_iptw)
    weight.iptw <- ifelse(d_iptw$treatment==1, prev1/d_iptw$PS, (1-prev1)/(1-d_iptw$PS))

    # get additional weights for IOPW
    m_ps2<-glm(as.formula(paste("S ~", paste(PSvariables, collapse="+"), sep="")), data=data, family=binomial)
    ps2 <- m_ps2$fitted
    weight.iopw = (1-ps2)/ps2 * length(which(data$S==1))/length(which(data$S==0))
    final_weight <-weight.iptw*data$weight.iopw[which(data$S==1)]

    # calculate Tipton index
    ps_logit<- logit(ps2)
    Bindex_v = Bindex(ps_logit[which(data$S==1)], ps_logit[which(data$S==0)])
    mi.output$Tipton[i]<-Bindex_v

    # fit Cox model with weights
    model<-coxph(formula = Surv(timeToEvent, event) ~ treatment, data = data[which(data$S==1),], weights=final_weight, ties = "breslow", robust = TRUE)
    mi.output[i,2:4]<-summary(model)$coefficients[c(1,2,4)]

    # check Cox model assumptions
    test.ph <- cox.zph(model) 
    mi.output$Schoenfeld[i]<-test.ph$table[2,"p"]
}

# INT-within: pool treatment estimates and SE according to Rubin's Rule
mean.list<-mi.output[,2]
overall.mean<-mean(mean.list)
robust.overall.se<-sqrt((mean(mi.output$robust.se^2) + var(mean.list)*(1+1/m)))

estimate = overall.mean
SE = robust.overall.se
HR = exp(estimate)
LCL = exp(estimate - qnorm(0.975) * SE)
UCL = exp(estimate + qnorm(0.975) * SE)