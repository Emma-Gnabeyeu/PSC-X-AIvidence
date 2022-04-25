library("tableone")
library("MatchIt")
library("Matching")

#1-Reading the data
#First download the lalonde.csv file available on Internet, link in the Report page 29
lalonde_data<-read.table("lalonde.csv",sep=",",dec=".",header=TRUE)
View(lalonde_data)

#collecting the variables
age<-lalonde_data$age
educ<-lalonde_data$educ
black<-lalonde_data$black
hisp<-lalonde_data$hisp
married<-lalonde_data$married
nodegr<-lalonde_data$nodegr
re74<-lalonde_data$re74
re75<-lalonde_data$re75
re78<-lalonde_data$re78
treat<-lalonde_data$treat

mydata<-cbind(age,educ,black,hisp,married,nodegr,re74,re75,re78,treat)
mydata<-data.frame(mydata)

xvars<-c("age","educ","black","hisp","married","nodegr","re74","re75")

table1<-CreateTableOne(vars=xvars, strata="treat",data=mydata, test=FALSE, smd=TRUE)
print(table1, smd=TRUE)
y_trt<-mydata$re78[mydata$treat==1]
y_con<-mydata$re78[mydata$treat==0]

diff_earn<-mean(y_trt)-mean(y_con)
diff_earn

#2-Matching on Covariates
#Perform the greedy matching and assess balance
greedymatch<-Match(Tr=treat,M=1,X=mydata[xvars])
matched0<-mydata[unlist(greedymatch[c("index.treated","index.control")]),]
matchedtab0<-CreateTableOne(vars=xvars, strata="treat", data=matched0, test=FALSE)
print(matchedtab0, smd=TRUE)

y_trt<-matched0$re78[matched0$treat==1]
y_con<-matched0$re78[matched0$treat==0]
diff_earn<-mean(y_trt)-mean(y_con)
diff_earn


psmodel<-glm(treat~age+educ+black+hisp+married+nodegr+re74+re75, family=binomial(),data=mydata)
pscore<-psmodel$fitted.values

summary(pscore)

#3- Matching on propensity score
set.seed(931139)

psmatch<-Match(Tr=treat,M=1,X=pscore,replace=FALSE)
matched=mydata[unlist(psmatch[c("index.treated","index.control")]),]

matchedtab1<-CreateTableOne(vars=xvars, strata="treat",data=matched, test=FALSE, smd=TRUE)
print(matchedtab1, smd=TRUE)

y_trt<-matched$re78[matched$treat==1]
y_con<-matched$re78[matched$treat==0]
diff_earn<-mean(y_trt)-mean(y_con)
diff_earn


set.seed(931139)
psmatch2<-Match(Tr=treat,M=1,X=pscore,replace=FALSE,caliper=0.1)
matched2=mydata[unlist(psmatch2[c("index.treated","index.control")]),]

matchedtab2<-CreateTableOne(vars=xvars, strata="treat",data=matched2, test=FALSE, smd=TRUE)
print(matchedtab2, smd=TRUE)

y_treat<-matched2$re78[matched2$treat==1]
y_con<-matched2$re78[matched2$treat==0]

diff_y<-y_treat-y_con
mean(diff_y)

t.test(diff_y)


#4- Inverted Probability of Treatment Weighting approach

weight=ifelse(treat==1,1/pscore,1/(1-pscore))
summary(weight)

weighteddata<-svydesign(ids=~1, data=mydata, weights=~weight)
weightedtable<-svyCreateTableOne(vars = xvars, strata = "treat", data =weighteddata, test =FALSE)

print(weightedtable, smd=TRUE)

#Truncate weights
#To remove excessive weights that might skew the results
truncweight<-replace(weight,weight>10,10)

#get causal risk difference
glm.obj<-glm(re78~treat, weights=truncweight, family=gaussian(link=identity))
coef(glm.obj)

 
weightmodel<-ipwpoint(exposure=treat, family="binomial", link="logit",
denominator=~age+educ+black+hisp+married+nodegr+re74+re75, data=mydata, trunc=.01)

#numeric summary of weights
summary(weightmodel$weights.trun)

#plot of weights
#a density kind of plot
ipwplot(weights=weightmodel$weights.trun, logscale=FALSE,
main = "weights", xlim=c(0,13))

#Fit the MSM (for estimating the risk difference)
#recall that svydesign() affects weights to data
mydata$wt<-weightmodel$weights.trun
msm<-(svyglm(re78~treat, design=svydesign(~1, weights=~wt, data=mydata)))
coef(msm)
confint(msm)