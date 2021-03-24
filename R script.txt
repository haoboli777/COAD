##############################################################################################################################
setwd("C:\\Users\\dell\\Desktop\\COAD")
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
library(glmnet)
library(timeROC)
library(survivalROC)
library(pheatmap)
library(rms)

#Sankey diagram
rt=read.table("network.txt",sep = "\t",header = T)                   
cox=read.table("multi.Cox.txt",sep="\t",header=T,row.names=1)
protectGene=row.names(cox[cox$HR<1,])
riskType=ifelse(rt$lncRNA%in%protectGene,"Protect","Risk")
newData=cbind(rt[,c(1,2)],riskType)
corLodes=to_lodes_form(newData, axes = 1:3, id = "Cohort")
pdf(file="ggalluvial.pdf",width=7,height=6)
mycol <- rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),5)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 geom_flow(width = 1/8,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 1/10) +
	 scale_fill_manual(values = mycol) +
	 geom_text(stat = "stratum", size = 2.4,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + 
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()

# Univariate cox regression
coxPfilter=0.05 
rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)     
rt$futime=rt$futime/365
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
	if(sd(rt[,i])<0.001){next}
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<coxPfilter){
	        sigGenes=c(sigGenes,i)
			outTab=rbind(outTab,
			             cbind(id=i,
			             HR=coxSummary$conf.int[,"exp(coef)"],
			             HR.95L=coxSummary$conf.int[,"lower .95"],
			             HR.95H=coxSummary$conf.int[,"upper .95"],
			             pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			             )
	}
}
write.table(outTab,file="uni.Cox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uni.SigExp.txt",sep="\t",row.names=F,quote=F)

# LASSO regression
rt=read.table("uni.SigExp.txt",header=T,sep="\t",row.names=1,check.names=F)       
rt$futime[rt$futime<=0]=0.003
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)

# Multivariate cox regression
rt=read.table("lasso.SigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")   
multiCoxSum=summary(multiCox)
outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multi.Cox.txt",sep="\t",row.names=F,quote=F)
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
riskOut=cbind(rt[,outCol],riskScore,risk)
riskOut=cbind(id=rownames(riskOut),riskOut)
write.table(riskOut,file="geneRisk.txt",sep="\t",quote=F,row.names=F)
pdf(file="multi.forest.pdf",width = 10,height = 6,onefile = FALSE)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

# K-M survival analysis
inputFile="geneRisk.txt"
survFile="survival.pdf" 
rt=read.table(inputFile,header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=5,
		           risk.table=TRUE,
		           legend.labs=c("High risk", "Low risk"),
		           legend.title="Risk",
		           xlab="Time(years)",
		           break.time.by = 1,
		           risk.table.title="",
		           palette=c("red", "blue"),
		           risk.table.height=.25)
pdf(file=survFile,onefile = FALSE,width = 12.5,height =9)
print(surPlot)
dev.off()

# ROC curve
rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =1, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

rocPlot(inputFile="geneRisk.txt",outPdf="1rocTrain.pdf")

#riskScore
rt=read.table("generisk.txt",sep="\t",header=T,row.names=1,check.names=F)       
rt=rt[order(rt$riskScore),]                                     
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file="riskScore.pdf",width = 10,height = 3.5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
     rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
dev.off()

#survStat
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStat.pdf",width = 10,height = 3.5)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

#heatmap
rt1=log2(rt[c(3:(ncol(rt)-2))]+0.01)
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmap.pdf",width = 10,height = 3)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         color = colorRampPalette(c("green", "black", "red"))(50) )
dev.off()

#Univariate independent prognostic analysis
rt=read.table("indepInput.txt",header=T,sep="\t",check.names=F,row.names=1)
uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	 coxSummary = summary(cox)
	 uniTab=rbind(uniTab,
	              cbind(id=i,
	              	    B=coxSummary$coefficients[,"coef"],
	                    SE=coxSummary$coefficients[,"se(coef)"],
	                    z=coxSummary$coefficients[,"z"],
	                    HR=coxSummary$conf.int[,"exp(coef)"],
	                    HR.95L=coxSummary$conf.int[,"lower .95"],
	                    HR.95H=coxSummary$conf.int[,"upper .95"],
	                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
	              )
}
write.table(uniTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

#multivariate independent prognostic analysis
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
	         B=multiCoxSum$coefficients[,"coef"],
	         SE=multiCoxSum$coefficients[,"se(coef)"],
	         z=multiCoxSum$coefficients[,"z"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
bioForest=function(coxFile=null,forestCol=null,forestFile=null){
		rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
		gene <- rownames(rt)
		hr <- sprintf("%.3f",rt$"HR")
		hrLow  <- sprintf("%.3f",rt$"HR.95L")
		hrHigh <- sprintf("%.3f",rt$"HR.95H")
		Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
		pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		pdf(file=forestFile, width = 6,height = 4.2)
		n <- nrow(rt)
		nRow <- n+1
		ylim <- c(1,nRow)
		layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		xlim = c(0,3)
		par(mar=c(4,2.5,2,1))
		plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
		text.cex=0.8
		text(0,n:1,gene,adj=0,cex=text.cex)
		text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
		text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
		par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
		xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
		plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
		arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
		abline(v=1,col="black",lty=2,lwd=2)
		boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
		points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
		axis(1)
		dev.off()
}
bioForest(coxFile="uniCox.txt", forestCol="green", forestFile="uniForest.pdf")
bioForest(coxFile="multiCox.txt", forestCol="red", forestFile="multiForest.pdf")

#normogram
riskFile="indepInput.txt"     
outFile="Nomogram.pdf"      
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)       
rt=risk[,1:(ncol(risk)-0)]
rt[,"futime"]=rt[,"futime"]/365  
dd <- datadist(rt)
options(datadist="dd")
f <- cph(Surv(futime, fustat) ~ age		+stage 		+riskScore
, x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), 
    lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.5, 0.3,0.1,0.01))  
pdf(file=outFile,height=6,width=9)
plot(nom)
dev.off()

#calibration curve
time=5  
f <- cph(Surv(futime, fustat) ~ age		+stage		+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=110, B=1000)
pdf(file="5calibration.pdf",height=6,width=8)
plot(cal,xlab="Nomogram-Predicted Probability of 5-Year OS",ylab="Actual 5-Year OS(proportion)",col="red",sub=F)
dev.off()

#C-INDEX
fmla1 <- cph(Surv(futime, fustat) ~ age		+stage		+riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
cox2 <- coxph(Surv(futime, fustat) ~age		+stage		+riskScore ,data = rt)
summary(cox2)

#~!!!