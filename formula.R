
downloadNHANES <- function(fileprefix){
  print (fileprefix)
  outdf <- data.frame(NULL)
  for (j in 1:length(letters)){
    urlstring <- paste('https://wwwn.cdc.gov/nchs/nhanes/',yrs[j],'/',fileprefix,letters[j],'.XPT', sep='')
    download.file(urlstring, tf <- tempfile(), mode="wb")
    tmpframe <- foreign::read.xport(tf)
    outdf <- bind_rows(outdf, tmpframe)
  }
  return(outdf)
}
######################################################???????غ???
find_na <- function(all){temp <- all
for(i in colnames(temp)){temp[,i] <- ifelse(is.na(temp[,i]),1,0)

}
temp <- apply(temp,2,sum)
return(temp)
}
plotrcs <- function(module,newdata,title,family='gaoss',x.axis.min=1.5,x.axis.max=9,y.axis.min=0,y.axis.max=5,main_var,lab.x='%UPF',lab.y="OR (95%CI)",pwald,pwald.y=1.003,pwald.x=1,margin=c(5,5,5,5))
{
  preds <- as.data.frame(predict(module,newdata = newdata,type = 'response'))
  if (family=='bio'){
    preds <- mutate(preds,upper=response+1.96*SE,lower=response-1.96*SE,)%>%select(.,response,upper,lower)%>%sapply(., function(x) exp(x))%>%as.data.frame(.)}
  else {
    preds <- mutate(preds,upper=response+1.96*SE,lower=response-1.96*SE,)%>%select(.,response,upper,lower)}
  preds <- cbind(newdata,preds)
  ggplot()+geom_line(data=preds, aes_string(main_var,'response'),linetype=1,size=1,alpha = 0.9,colour="black")+
    geom_ribbon(data=preds, aes_string(main_var,ymin = 'lower', ymax = 'upper'),alpha = 0.3,fill="grey40")+
    theme_classic()+  
    labs(title = title, x=lab.x, y=lab.y)+xlim(x.axis.min,x.axis.max)+ylim(y.axis.min,y.axis.max)+
    theme(plot.title=element_text(size=20, 
                                  face="bold", 
                                  hjust=0.5, #????λ?ã????м?
                                  lineheight=1.2)
          ,axis.title.x = element_text(size = 15,
          ),
          axis.title.y=element_text(size=15,
                                    color = "black"),
          axis.text.x =element_text(size = 10,
                                    colour = 'black',
                                    face="bold"
          ), 
          axis.text.y =element_text(size = 10,
                                    colour = 'black',
                                    face="bold"),
          plot.margin=unit(margin,'cm'))+
    annotate("text",x=pwald.x,y=pwald.y
             ,label=paste0('P-nonlinear=',pwald)
             ,size=5.5,alpha=0.9,colour="black")}
creat_ka <- function(varies,design){formul <- sapply(varies,function(x) as.formula(paste0('~Diigr+',x)))
ka <- lapply(formul,function(x) count_Diigr(x,design))
ka <- lapply(varies, function(x) dcast(ka[[x]],as.formula(paste0(x,'~Diigr')),value.var='rio'))
names(ka) <- varies
orig <- data.frame(matrix(ncol = 4))
colnames(orig) <- 1:4
for (i in varies) {rownames(ka[[i]]) <- paste0(i,'  ',ka[[i]][,1])
orig <- rbind(orig,ka[[i]][,-1])
}
return(na.omit(orig))}
glp <- function(a,coa,ma,design,family='gaussian',level=0.95,test=1,interact=''){
  if (interact==''){
    formu <-sapply(paste0(ma,'~',a,'+',coa),function(x) as.formula(x))
  }
  else {intra <- paste0("+",a,':',interact)
  form <- paste0(ma,'~',a,'+',coa)
  formu <- sapply(paste0(form,intra),function(x) as.formula(x))}
  names(formu) <- a  
  fit <- lapply(formu, function(x) svyglm(x,family=family,design=design))  
  tes=switch (test,a,paste0(a,':',interact))
  tes=sapply(paste0('~',tes),function(x) as.formula(x))
  names(tes) <- a
  Pwald <- lapply(a, function(x){regTermTest(fit[[x]], tes[[x]])})
  names(Pwald) <- a
  Pwald <- lapply(a, function(x) Pwald[[x]][['p']])
  names(Pwald) <- a
  coe <- lapply(fit, function(x) coef(summary(x)))
  coef <- lapply(a, function(x) coe[[x]][grepl(paste0(x,'*'),rownames(coe[[x]]))][1:2])
  names(coef) <- a
  mult <- qnorm((1+level)/2)
  if(family=='gaussian'){upper <- lapply(a, function(x) coef[[x]][1]+coef[[x]][2]*mult)
  names(upper) <- a
  lower <- lapply(a, function(x) coef[[x]][1]-coef[[x]][2]*mult)
  names(lower) <- a
  CI <- lapply(a, function(x) paste0(signif(coef[[x]][1],2),' (',signif(lower[[x]],2),', ',signif(upper[[x]],2),')'))}
  else{upper <- lapply(a, function(x) exp(coef[[x]][1]+coef[[x]][2]*mult))
  names(upper) <- a
  lower <- lapply(a, function(x) exp(coef[[x]][1]-coef[[x]][2]*mult))
  names(lower) <- a
  CI <- lapply(a, function(x) paste0(signif(exp(coef[[x]][1]),2),' (',signif(lower[[x]],2),', ',signif(upper[[x]],2),')'))}
  names(CI) <- a
  allvalues <- list(unlist(Pwald),CI)
  names(allvalues) <- c('p','CI')
  return(allvalues)
}
glpr <- function(a,coa,ma,design,family='gaussian',level=0.95,interact=''){if (interact==''){
  formu <-sapply(paste0(ma,'~',a,'+',coa),function(x) as.formula(x))
}
  else {intra <- paste0("+",a,':',interact)
  form <- paste0(ma,'~',a,'+',coa)
  formu <- sapply(paste0(form,intra),function(x) as.formula(x))}
  names(formu) <- a  
  fit <- lapply(formu, function(x) svyglm(x,family=family,design=design))  
  summ <- lapply(fit, function(x) summary(x))
  Pv <- lapply(a, function(x) summ[[x]][["coefficients"]][2:4,'Pr(>|t|)'])
  coe <- lapply(fit, function(x) coef(summary(x)))
  coef <- lapply(a, function(x) coe[[x]][grepl(paste0(x,'*'),rownames(coe[[x]])),1:2])
  names(coef) <- a
  mult <- qnorm((1+level)/2)
  if(family=='gaussian'){coef <- lapply(a, function(x) mutate(as.data.frame(coef[[x]]),upper=round(Estimate+coef[[x]][,2]*mult,2),lower=round(Estimate-coef[[x]][,2]*mult,2)))
  names(coef) <- a
  coef <- lapply(a, function(x) mutate(as.data.frame(coef[[x]]),Estimate=round(Estimate,2)))}
  else{coef <- lapply(a, function(x) mutate(as.data.frame(coef[[x]]),upper=round(exp(Estimate+coef[[x]][,2]*mult),2),lower=round(exp(Estimate-coef[[x]][,2]*mult),2)))
  names(coef) <- a
  coef <- lapply(a, function(x) mutate(as.data.frame(coef[[x]]),Estimate=round(exp(Estimate),2)))}
  names(coef) <- a
  coef <- lapply(a, function(x) paste0(coef[[x]][,1],' (',coef[[x]][,4],' ,',coef[[x]][,3],')'))
  coef <- sapply(coef, '[',1:3)
  coef <- as.matrix(coef)
  dim(coef) <- c(length(a)*3,1)
  allvalues <- data.frame(unlist(Pv),coef)
  names(allvalues) <- paste0(c('p','CI'),'-',ma)
  return(allvalues)}
creat_newtab <- function(tab1,tab_isn,ins_col=2:5){
  if(!(inherits(tab1,'data.frame')&inherits(tab_isn,'data.frame')))
  {stop('error: tab and tab_isn must be dataframe')}
  if(!(nrow(tab1)==nrow(tab_isn))&ncol(tab1)==ncol(tab_isn))
  {stop('error :tab and tab_isn must have same rows and cols')}
  tab=tab1
  index=grep(') $',tab[,ins_col[1]])
  for (i in ins_col) {insd1=sapply(strsplit(grep(') $',tab[,i],value = T),'(',fixed = T), '[',1)
  insd2=sapply(strsplit(grep(') $',tab_isn[,i],value = T),'(',fixed = T), '[',2)
  tab[,i][index]=paste0(insd1,'(',insd2)
  }
  return(tab)}

meff <- function(x,method='spearman'){cor_data <- cor(x,method =method)
eign <- eigen(cor_data)[["values"]]-sapply(eigen(cor_data)[["values"]],floor)+sapply(eigen(cor_data)[["values"]], function(x) ifelse(x>=1,1,0))
return(sum(eign))}
Pmeff <- function(P,meff){p.ad <- sort(P)
M <- length(P)
p.fdr <- vector(length = length(P))
names(p.fdr) <- names(p.ad)
for(i in p.ad) {
  
  
  p.fdr[which(p.ad==i)]=i*meff*(M-1)/((M-1)+(which(p.ad==i)-1)*(meff-1))
  
  
  
  
}
for (i in 2:length(P)) {
  p.fdr[i-1]=ifelse(p.fdr[i]>p.fdr[i-1],p.fdr[i-1],p.fdr[i])  
}
p.fdr <- p.fdr[names(P)]
return(p.fdr)}

modulefc <- function(module){va=cbind(as.data.frame(exp(coef(module))),as.data.frame(exp(confint(module))))
modulep <- function(module){test=summary(module)
return(test$coefficients[,'Pr(>|t|)'])}
va <- round(va,3)
va$re <- paste0(va[,1],' (',va[,2],', ',va[,3],')')
return(va)}
modulep <- function(module){test <- summary(module)
return(test$coefficients[,'Pr(>|t|)'])}
modulefc_gasso <- function(module){va=cbind(as.data.frame(coef(module)),as.data.frame(confint(module)))
va <- round(va,3)
va$re <- paste0(va[,1],' (',va[,2],', ',va[,3],')')
return(va)}
creatforest <- function(module,label,ref.set=1){OR <- modulefc(module)[2:4,]
names(OR) <- c('OR','lower','upper','range')
OR <- add_row(OR=1.00,lower=1.00,upper=1.00,range='[Reference]',.before = ref.set,.data = OR)
OR$group <- paste0('Quartile',1:4)
OR <- add_row(group=label,.data = OR,.before = 1)}
creatforest_gasso <- function(module,label,ref.set=1){OR <- modulefc_gasso(module)[2:4,]
names(OR) <- c('OR','lower','upper','range')
OR <- add_row(OR=0.00,lower=0.00,upper=0.00,range='[Reference]',.before = ref.set,.data = OR)
OR$group <- paste0('Quartile',1:4)
OR <- add_row(group=label,.data = OR,.before = 1)}