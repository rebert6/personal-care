if(!require("survey")) install.packages("survey")
if(!require("dplyr")) install.packages("dplyr")
if(!require("rms")) install.packages("rms")
if(!require("tableone")) install.packages("tableone")
if(!require("lsmeans")) install.packages("lsmeans")
if(!require('ggcor')) install.packages("ggcor")####若提示版本问题，请GitHub下载安装https://github.com/mj163163/ggcor-1
if(!require('cowplot')) install.packages("cowplot")
if(!require('corrplot')) install.packages("corrplot")
if(!require('Hmisc')) install.packages("Hmisc")
if(!require('cluster')) install.packages("cluster")
if(!require('reshape2')) install.packages("reshape2")
if(!require('factoextra')) install.packages("factoextra")
if(!require('gridExtra')) install.packages("gridExtra")
if(!require('Rtsne')) install.packages("Rtsne")
if(!require('scatterplot3d')) install.packages("scatterplot3d")
if(!require('fmsb')) install.packages("fmsb")
if(!require('bkmr')) install.packages("bkmr")
options(survey.lonely.psu='adjust')
source('formula.R')
######################耳毒性药物
du <- c("h00026","h00031","d03456","d03455","d03425" ,"d03346", "d04269" ,"d05644" ,"h00003", "d06872" ,"d03319",
        "a11122" ,"d03423" ,"d03342" ,"d03459" ,"d03428" ,"h00027" ,"d03431", "d03682" ,"d03297" ,"d03434" ,"d03289",
        "d04333", "d04766", "d00170" ,"d03458" ,"d03457" ,"d03426", "d03430" ,"d03472", "d03435" ,"d03469", "d03470",
        "d03424" ,"d04497" ,"d03429" ,"d03448" ,"d03468", "d03432", "d03291" ,"d00035" ,"d00070", "d00014", "c00061",
        "d00195", "d00185", "d07868" ,"d00179")
######################数据下载
yrs <- c('2005-2006','2007-2008','2009-2010','2011-2012','2013-2014','2015-2016')
letters <- c('_D','_E','_F','_G','_H','_I')
DEMO <- downloadNHANES('DEMO')
DEMO <- DEMO[,c("SEQN",'SDDSRVYR',"RIAGENDR","RIDAGEYR","SDMVSTRA","DMDEDUC2","RIDRETH1","SDMVPSU",'INDFMPIR')]
######
BPQ <- downloadNHANES('BPQ')
BPQ <- BPQ[,c("SEQN","BPQ020")]
#####
DIQ <- downloadNHANES('DIQ')
DIQ <- DIQ[,c("SEQN","DIQ010")]
######
SMQ <- downloadNHANES('SMQ')
SMQ <- SMQ[,c("SEQN","SMQ020","SMQ040")]
SMQ <- mutate(SMQ,smoke=case_when(SMQ020==2~1,SMQ020==1&(SMQ040==1|SMQ040==2)~3,SMQ020==1&SMQ040==3~2))%>%select(.,SEQN,smoke)
#####
BMX <- downloadNHANES('BMX')
BMX <- BMX[,c("SEQN","BMXBMI")]
#####
RXQ <- downloadNHANES('RXQ_RX')
RXQ <- RXQ[,c("SEQN","RXDUSE","RXDDRGID")]
RXQ <- mutate(RXQ,RXDUSE=ifelse(RXDUSE==7|RXDUSE==9,NA,RXDUSE))
RXQ <- mutate(RXQ,medicine=case_when(RXDDRGID%in%du~1,!RXDDRGID%in%du~0))%>%
  mutate(.,medicine=ifelse(is.na(RXDUSE),NA,medicine))%>%select(.,SEQN,medicine)
undup <- as.data.frame(tapply(RXQ$medicine,RXQ$SEQN, function(x) mean(x,na.rm = T)))
undup$SEQN <- rownames(undup)
rownames(undup) <- 1:nrow(undup)
colnames(undup)[1] <-"medicine"
undup <- mutate(undup,medicine=ifelse(medicine==0,0,1))
RXQ <- undup
RXQ$SEQN <- as.numeric(RXQ$SEQN)
rm(undup)
####
MCQ <- downloadNHANES('MCQ')
MCQ <- MCQ[,c("SEQN","MCQ160B","MCQ160F")]
colnames(MCQ) <- c("SEQN","HartF","stroke")
####
yrs <- c('2005-2006','2007-2008','2009-2010','2011-2012','2013-2014','2015-2016')
letters <- c('_D','_E','_F','_G','_H','_I')
UCR <- downloadNHANES('ALB_CR')
UCR <- UCR[,c('SEQN','URXUCR')]
yrs <- c('2005-2006','2007-2008','2009-2010')
letters <- c('_D','_E','_F')
PP0510 <- downloadNHANES('PP')
PP0510 <- PP0510[,c("SEQN",'WTSB2YR','URX14D','URD14DLC','URXDCB','URDDCBLC')]
EP0510 <- downloadNHANES('EPH')
EP0510 <- EP0510[,c("SEQN",'URXBPH','URDBPHLC','URXBP3','URDBP3LC','URXTRS','URDTRSLC','URXMPB','URDMPBLC','URXPPB','URDPPBLC')]
EPH0510 <- left_join(PP0510,EP0510,by='SEQN')
yrs <- c('2011-2012')
letters <- c('_G')
PP1112 <- downloadNHANES('PP')
PP1112 <- PP1112[,c("SEQN",'WTSA2YR','URX14D','URD14DLC','URXDCB','URDDCBLC')]
names(PP1112)[2] <- 'WTSB2YR'
EP1112 <- downloadNHANES('EPH')
EP1112 <- EP1112[,c("SEQN",'URXBPH','URDBPHLC','URXBP3','URDBP3LC','URXTRS','URDTRSLC','URXMPB','URDMPBLC','URXPPB','URDPPBLC')]
EPH1112 <- left_join(PP1112,EP1112,by='SEQN')
yrs <- c('2013-2014','2015-2016')
letters <- c('_H','_I')
EPH1316 <- downloadNHANES('EPHPP')
EPH1316 <- EPH1316[,c("SEQN",'WTSB2YR','URX14D','URD14DLC','URXDCB','URDDCBLC','URXBPH','URDBPHLC','URXBP3','URDBP3LC','URXTRS','URDTRSLC','URXMPB','URDMPBLC','URXPPB','URDPPBLC')]
EPH <- rbind(EPH0510,EPH1112,EPH1316)
EPH_ID <-select(EPH,ends_with('LC'))
EPH_ID <- na.omit(EPH_ID)
rio <- sapply(apply(EPH_ID,2,sum), function(x) 1-(x/nrow(EPH_ID)))
EPHun <- select(EPH,-ends_with('LC'))
EPH <- left_join(EPHun,UCR,by='SEQN')
EPH[,3:9] <- sapply(EPH[,3:9],function(x)  round(x/EPH[,10]*100,2))
EPH <- EPH[,-10]
write.csv(EPH,'EPH.csv')
########听力
yrs <- c('2005-2006','2007-2008','2009-2010','2011-2012','2015-2016')
letters <- c('_D','_E','_F','_G','_I')
AUX <- downloadNHANES('AUX')
AUX0 <- select(AUX,SEQN,AUXOTSPL,AUXROTSP,AUXTCOMR,AUXTCOML,AUATYMTL,AUATYMTR,AUXTMEPR,AUXTMEPL,starts_with('AUXU')) %>%select(.,-AUXU1K2R,-AUXU1K2L)%>%
  filter(.,AUXROTSP==1&AUXTCOMR>0.3&AUXOTSPL==1&AUXTCOML>0.3) %>% 
  mutate_at(vars(starts_with('AUXU')),~ifelse(.==666,120,.)) %>% 
  mutate_at(vars(starts_with('AUXU')),~ifelse(.==888|.==777,NA,.))%>% 
  mutate_at(vars(starts_with('AUXTC')),~ifelse(.==888|.==777,NA,.))%>% 
  mutate(.,LFR=apply(select(.,AUXU1K1R,AUXU500R,AUXU2KR),1,mean),
         LFL=apply(select(.,AUXU1K1L,AUXU500L,AUXU2KL),1,mean),
         HFR=apply(select(.,AUXU3KR,AUXU4KR,AUXU6KR,AUXU8KR),1,mean),
         HFL=apply(select(.,AUXU3KL,AUXU4KL,AUXU6KL,AUXU8KL),1,mean)) %>%
  mutate(noise_R=ifelse(AUXU3KR<=AUXU4KR&AUXU4KR>=AUXU6KR&AUXU4KR>25,1,0),
         noise_L=ifelse(AUXU3KL<=AUXU4KL&AUXU4KL>=AUXU6KL&AUXU4KL>25,1,0)) %>% 
  mutate(.,LF=ifelse(is.na(LFR),ifelse(is.na(LFL),NA,LFL),ifelse(is.na(LFL),LFR,ifelse(LFR>LFL,LFR,LFL))),
         HF=ifelse(is.na(HFR),ifelse(is.na(HFL),NA,HFL),ifelse(is.na(HFL),HFR,ifelse(HFR>HFL,HFR,HFL))),
         noise_4k=ifelse(is.na(noise_R),ifelse(is.na(noise_L),NA,noise_L),ifelse(is.na(noise_L),noise_R,ifelse(noise_L==1|noise_R==1,1,0)))) %>% 
  mutate(.,LFHL=ifelse(LF>25,1,0),
         HFHL=ifelse(HF>25,1,0))
AUX1 <- select(AUX0,SEQN,LF,HF,LFR,HFR,LFL,HFL,LFHL,HFHL,noise_4k,noise_R,noise_L)
####################数据合并
all <- left_join(EPH,DEMO,by='SEQN')%>% left_join(.,BPQ,by='SEQN') %>% left_join(.,DIQ,by='SEQN') %>% left_join(.,MCQ,by='SEQN') %>% 
  left_join(.,RXQ,by='SEQN') %>% left_join(.,SMQ,by='SEQN') %>% left_join(.,BMX,by='SEQN') %>% 
  left_join(.,dplyr::select(AUX1,SEQN,LF,HF),by='SEQN')%>% left_join(.,UCR,by='SEQN')
all <- filter(all,!is.na(all$WTSB2YR))
all <- mutate(all,WTMEC12YR=WTSB2YR/6)%>%
  mutate(.,DIQ010=ifelse(DIQ010==3,1,ifelse(DIQ010==9|DIQ010==7,2,DIQ010)))%>%
  mutate(.,BPQ020=ifelse(BPQ020==9,2,BPQ020))%>%
  mutate(.,HartF=ifelse(HartF==9|HartF==7,2,HartF))%>%
  mutate(.,stroke=ifelse(stroke==9|stroke==7,2,stroke))%>%
  mutate(.,DMDEDUC2=ifelse(DMDEDUC2==9|DMDEDUC2==7,NA,DMDEDUC2))%>%
  mutate(.,Gender = factor(RIAGENDR, levels = c(1,2),labels=c("Men", "Women")),
         Age.Group = cut(RIDAGEYR, breaks=c(-Inf,19,44,59,Inf),labels=c(NA, NA,"45-59","60 and over")),
         BMI=cut(BMXBMI,breaks = c(-Inf,24.99,29.99,Inf),labels = c(1,2,3)),
         race=case_when(RIDRETH1==3~0,RIDRETH1==4~-1,RIDRETH1==1|RIDRETH1==2|RIDRETH1==5~1),
         PIR=case_when(INDFMPIR<=1~1,INDFMPIR>1&INDFMPIR<=3~2,INDFMPIR>3~3)
  )%>%
  mutate(.,education=case_when(DMDEDUC2==1|DMDEDUC2==2~1,DMDEDUC2==3~2,DMDEDUC2==4~3,DMDEDUC2==5~4))
fac <- c('race','Age.Group','education','RIDRETH1','DMDEDUC2','Gender','medicine','smoke','BPQ020','DIQ010','stroke','HartF','BMI','PIR')
all[,fac] <- lapply(all[,fac], factor)
##################数据筛选建模
hea<- filter(all,RIDAGEYR>=45) %>% na.omit(.)
heamodel <- mutate(all,inAnalysis= SEQN%in%hea$SEQN)
heaM <- svydesign(data=heamodel, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC12YR, nest=TRUE)
heamod <- subset(heaM, inAnalysis)
##################流程图中的数字
DEMO_be <- filter(DEMO,RIDAGEYR>=45)
EPH_be <- na.omit(EPH)
ce_be <- left_join(EPH_be,DEMO_be,by='SEQN')
ce1 <- filter(ce_be,!is.na(RIDAGEYR))
print(length(ce1))
ting_be <- select(AUX,SEQN,AUXOTSPL,AUXROTSP,AUXTCOMR,AUXTCOML,starts_with('AUXU')) %>%select(.,-AUXU1K2R,-AUXU1K2L)
ting <- left_join(ce1,ting_be,by='SEQN')
ce2 <- filter(ting,AUXROTSP==1&AUXTCOMR>0.3&AUXOTSPL==1&AUXTCOML>0.3)
print(length(ce2))
ce2_be <- mutate_at(ce2,vars(starts_with('AUXU')),~ifelse(.==666,120,.)) %>% 
  mutate_at(vars(starts_with('AUXU')),~ifelse(.==888|.==777,NA,.))%>% 
  mutate_at(vars(starts_with('AUXTC')),~ifelse(.==888|.==777,NA,.))%>% 
  mutate(.,LFR=apply(select(.,AUXU1K1R,AUXU500R,AUXU2KR),1,mean),
         LFL=apply(select(.,AUXU1K1L,AUXU500L,AUXU2KL),1,mean),
         HFR=apply(select(.,AUXU3KR,AUXU4KR,AUXU6KR,AUXU8KR),1,mean),
         HFL=apply(select(.,AUXU3KL,AUXU4KL,AUXU6KL,AUXU8KL),1,mean)) %>% 
  mutate(.,LF=ifelse(is.na(LFR),ifelse(is.na(LFL),NA,LFL),ifelse(is.na(LFL),LFR,ifelse(LFR>LFL,LFR,LFL))),
         HF=ifelse(is.na(HFR),ifelse(is.na(HFL),NA,HFL),ifelse(is.na(HFL),HFR,ifelse(HFR>HFL,HFR,HFL)))) %>% 
  mutate(.,LFHL=ifelse(LF>25,1,0),
         HFHL=ifelse(HF>25,1,0)) 
ce3 <- filter(ce2_be,(!is.na(LF))&(!is.na(HF)))
print(length(ce3))
ce3_be <- left_join(ce3,BPQ,by='SEQN') %>% left_join(.,DIQ,by='SEQN') %>% left_join(.,MCQ,by='SEQN') %>% 
  left_join(.,RXQ,by='SEQN') %>% left_join(.,SMQ,by='SEQN') %>% left_join(.,BMX,by='SEQN') %>%left_join(.,UCR,by='SEQN')
find_na(ce3_be)
#####################相关性分析和有效检验次数
Me <- meff(log(hea[,3:9]))
expo <-as.matrix(scale(log(hea[,3:9])))
colnames(expo) <- c('2,4-二氯苯酚','2,5-二氯苯酚','双酚A','二苯甲酮-3','三氯生','对羟基苯甲酸甲脂','对羟基苯甲酸丙脂')
jpeg('图2.jpg',width = 4800, height = 4800)
quickcor(expo, cor.test = T,method='spearman')+
  geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
  geom_mark(data = get_data(type = "upper", show.diag = FALSE), size = 5)
dev.off()
###################人口统计信息
table <- svyCreateTableOne(vars = c('HF','LF','race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF','Gender','BMXBMI','Age.Group','BMI','PIR','URXUCR','RIDAGEYR'),
                           data = heamod,factorVars = c('race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF','Age.Group','BMI','PIR'))
table_1 <- as.data.frame(print(table, showAllLevels = TRUE))
tablen <- CreateTableOne(vars = c('HF','LF','race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF','Gender','BMXBMI','Age.Group','BMI','PIR','URXUCR','RIDAGEYR'),
                         data = hea)
tablen_1 <- as.data.frame(print(tablen, showAllLevels = TRUE))
table_tot <- creat_newtab(tab1 = tablen_1,tab_isn = table_1,ins_col = 2)
table_tot[c(2,3,35,36),] <- table_1[c(2,3,35,36),]
write.csv(table_tot,'表1.csv')
#################单个化学物质
aa <- c('URX14D','URXDCB','URXBPH','URXBP3','URXTRS','URXMPB','URXPPB')
logaa <- paste0('log(',aa,')')
logrcs <- paste0('rcs(log(',aa,'), 3)')
mian_test <- c('log(HF+12)','log(LF+11)')
coa <- 'Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010'
####加权线性回归
fitHF_vif <- svyglm(log(HF+12)~log(URXTRS)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)

vif(fitHF_vif)######方差膨胀因子
all_line_results <- lapply(mian_test,function(x) glp(a=logaa,ma=x,coa = coa,
                                                     design =heamod))
names(all_line_results) <- mian_test
P_line_HF <- Pmeff(all_line_results[["log(HF+12)"]][["p"]],Me)
CI_HF <- unlist(all_line_results[["log(HF+12)"]][["CI"]])
P_line_LF<- Pmeff(all_line_results[["log(LF+11)"]][["p"]],Me)
CI_LF <- unlist(all_line_results[["log(LF+11)"]][["CI"]])
all_line_table <- data.frame(CI_HF,round(P_line_HF,3),CI_LF,round(P_line_LF,3))
write.csv(all_line_table,'all_line_table.csv')
########性别交互作用
all_line_results_in_gender <- lapply(mian_test,function(x) glp(a=logaa,ma=x,coa = coa,interact = 'Gender',test = 2,
                                                               #family = "quasibinomial",
                                                               design =heamod))
names(all_line_results_in_gender) <- mian_test
all_rcs_results_in_gender <- lapply(mian_test,function(x) glp(a=logrcs,ma=x,coa = coa,interact = 'Gender',test = 2,
                                                              #family = "quasibinomial",
                                                              design =heamod))
names(all_rcs_results_in_gender) <- mian_test

P_line_HF_gender <- signif(all_line_results_in_gender[["log(HF+12)"]][["p"]],2)
P_line_LF_gender <- signif(all_line_results_in_gender[["log(LF+11)"]][["p"]],2)
P_rcs_HF_gender <- signif(all_rcs_results_in_gender[["log(HF+12)"]][["p"]],2)
P_rcs_LF_gender <- signif(all_rcs_results_in_gender[["log(LF+11)"]][["p"]],2)

inter_gender <- data.frame(P_line_HF_gender,P_rcs_HF_gender,P_line_LF_gender,P_rcs_LF_gender)
colnames(inter_gender) <- c('高频线性','高频非线性','低频线性','低频非线性')
rownames(inter_gender) <- c('2,4-二氯苯酚','2,5-二氯苯酚','双酚A','二苯甲酮-3','三氯生','对羟基苯甲酸甲脂','对羟基苯甲酸丙脂')
write.csv(inter_gender,'性别交互.csv',fileEncoding = "GB18030")
####限制性立方样条
all_rcs_results <- lapply(mian_test,function(x) glp(a=logrcs,ma=x,coa = coa,
                                                    #family = "quasibinomial",
                                                    design =heamod))
names(all_rcs_results) <- mian_test

P_HF_rcs_adj <-Pmeff(all_rcs_results[["log(HF+12)"]][["p"]],Me)
#高频听力阈值限制性立方样条调整后P值
P_HF_rcs_adj
P_LF_rcs_adj <-Pmeff(all_rcs_results[["log(LF+11)"]][["p"]],Me)
#低频听力阈值限制性立方样条调整后P值
P_LF_rcs_adj
##模型对比--赤池信息准则（AIC）
fitHF_TCS <- svyglm(log(HF+12)~rcs(log(URXTRS),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitHF_TCS_line <- svyglm(log(HF+12)~log(URXTRS)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                         ,design=heamod)
AIC(fitHF_TCS_line,fitHF_TCS)#####AIC值
###画图
fitHF_14D <- svyglm(log(HF+12)~rcs(log(URX14D),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitHF_BPH <- svyglm(log(HF+12)~rcs(log(URXBPH),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitHF_BP3 <- svyglm(log(HF+12)~rcs(log(URXBP3),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitHF_MPB <- svyglm(log(HF+12)~rcs(log(URXMPB),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitHF_PPB <- svyglm(log(HF+12)~rcs(log(URXPPB),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitHF_DCB <- svyglm(log(HF+12)~rcs(log(URXDCB),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
#####
newfc <- c('race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF')
newHL <- expand.grid(URXTRS=unique(hea$URXTRS),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P1 <- plotrcs(fitHF_TCS,newHL,main_var = 'log(URXTRS)',x.axis.min  =-1.8 ,x.axis.max = 9,y.axis.min = 3.3,y.axis.max = 4.1,title = '三氯生与高频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='浓度（g/g肌酐）',pwald =round(P_HF_rcs_adj[5],3), pwald.x = 5,pwald.y =4.1)

newHL <- expand.grid(URX14D=unique(hea$URX14D),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P2 <- plotrcs(fitHF_14D,newHL,main_var = 'log(URX14D)',x.axis.min  =-3.5 ,x.axis.max = 11,y.axis.min = 3,y.axis.max = 4.1,title = '2,4二氯苯酚与高频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='浓度（g/g肌酐）ֵ',pwald =round(P_HF_rcs_adj[1],3), pwald.x = 5,pwald.y =4.1)

newHL <- expand.grid(URXDCB=unique(hea$URXDCB),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P3 <- plotrcs(fitHF_DCB,newHL,main_var = 'log(URXDCB)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 3,y.axis.max = 4.1,title = '2,5二氯苯酚与高频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_HF_rcs_adj[1],3), pwald.x = 4.5,pwald.y =4.1)

newHL <- expand.grid(URXBPH=unique(hea$URXBPH),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P4 <- plotrcs(fitHF_BPH,newHL,main_var = 'log(URXBPH)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 3,y.axis.max = 4.1,title = '双份A与高频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_HF_rcs_adj[3],3), pwald.x =3,pwald.y =4.1)

newHL <- expand.grid(URXBP3=unique(hea$URXBP3),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P5 <- plotrcs(fitHF_BP3,newHL,main_var = 'log(URXBP3)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 3,y.axis.max = 4.1,title = '二苯甲酮-3与高频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_HF_rcs_adj[4],3), pwald.x = 3,pwald.y =4.1)

newHL <- expand.grid(URXMPB=unique(hea$URXMPB),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P6 <- plotrcs(fitHF_MPB,newHL,main_var = 'log(URXMPB)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 3,y.axis.max = 4.1,title = '对羟基苯甲酸甲酯与高频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_HF_rcs_adj[6],3), pwald.x = 3,pwald.y =4.1)

newHL <- expand.grid(URXPPB=unique(hea$URXPPB),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P7 <- plotrcs(fitHF_PPB,newHL,main_var = 'log(URXPPB)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 3,y.axis.max = 4.1,title = '对羟基苯甲酸丙酯与高频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_HF_rcs_adj[7],3), pwald.x = 3,pwald.y =4.1)

#######
jpeg(file='高频非线性.jpg',width=15000,height=8500,res=600)
ggdraw() +
  draw_plot(P1, x = 0, y = .5, width = .3, height = .5) +
  draw_plot(P2, x = .25, y = .5, width = .3, height = .5) +
  draw_plot(P3, x = .5, y = .5, width = .3, height = .5) +
  draw_plot(P4, x = .75, y = .5, width = .3, height = .5) +
  draw_plot(P5, x = 0, y = .0, width = .3, height = .5) +
  draw_plot(P6, x = .25, y = .0, width = .3, height = .5) +
  draw_plot(P7, x = .5, y = .0, width = .3, height = .5) +
  draw_plot_label(label = c("A","B","C",'D','E','F','G'), size = 22,
                  x = c(0.05,0.3,0.55,.8,0.05,.3,.55), y = c(0.9,0.9,0.9,0.9,0.4,0.4,0.4))
dev.off()



######低频
fitLF_TCS <- svyglm(log(LF+11)~rcs(log(URXTRS),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitLF_14D <- svyglm(log(LF+11)~rcs(log(URX14D),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitLF_BPH <- svyglm(log(LF+11)~rcs(log(URXBPH),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitLF_BP3 <- svyglm(log(LF+11)~rcs(log(URXBP3),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitLF_MPB <- svyglm(log(LF+11)~rcs(log(URXMPB),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitLF_PPB <- svyglm(log(LF+11)~rcs(log(URXPPB),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)
fitLF_DCB <- svyglm(log(LF+11)~rcs(log(URXDCB),3)+Gender+RIDAGEYR+education+URXUCR+INDFMPIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    ,design=heamod)



#####
newfc <- c('race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF')
newHL <- expand.grid(URXTRS=unique(hea$URXTRS),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P1 <- plotrcs(fitLF_TCS,newHL,main_var = 'log(URXTRS)',x.axis.min  =-1.8 ,x.axis.max = 9,y.axis.min = 2.5,y.axis.max = 3.8,title = '三氯生与低频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='浓度（g/g肌酐）',pwald =round(P_LF_rcs_adj[5],3), pwald.x = 5,pwald.y =3.8)

newHL <- expand.grid(URX14D=unique(hea$URX14D),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P2 <- plotrcs(fitLF_14D,newHL,main_var = 'log(URX14D)',x.axis.min  =-3.5 ,x.axis.max = 11,y.axis.min = 2.5,y.axis.max = 3.8,title = '2,4二氯苯酚与低频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='浓度（g/g肌酐）ֵ',pwald =round(P_LF_rcs_adj[1],3), pwald.x = 5,pwald.y =3.8)

newHL <- expand.grid(URXDCB=unique(hea$URXDCB),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P3 <- plotrcs(fitLF_DCB,newHL,main_var = 'log(URXDCB)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 2.5,y.axis.max = 3.8,title = '2,5二氯苯酚与低频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_LF_rcs_adj[1],3), pwald.x = 4.5,pwald.y =3.8)

newHL <- expand.grid(URXBPH=unique(hea$URXBPH),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P4 <- plotrcs(fitLF_BPH,newHL,main_var = 'log(URXBPH)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 2.5,y.axis.max = 3.8,title = '双份A与低频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_LF_rcs_adj[3],3), pwald.x =3,pwald.y =3.8)

newHL <- expand.grid(URXBP3=unique(hea$URXBP3),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P5 <- plotrcs(fitLF_BP3,newHL,main_var = 'log(URXBP3)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 2.5,y.axis.max = 3.8,title = '二苯甲酮-3与低频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_LF_rcs_adj[4],3), pwald.x = 3,pwald.y =3.8)

newHL <- expand.grid(URXMPB=unique(hea$URXMPB),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P6 <- plotrcs(fitLF_MPB,newHL,main_var = 'log(URXMPB)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 2.5,y.axis.max = 3.8,title = '对羟基苯甲酸甲酯与低频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_LF_rcs_adj[6],3), pwald.x = 3,pwald.y =3.8)

newHL <- expand.grid(URXPPB=unique(hea$URXPPB),race=0,URXUCR=113,Gender='Women',INDFMPIR=2.7,education=4,medicine=0,BMXBMI=10,smoke=1,stroke=2,HartF=2,RIDAGEYR=60,BPQ020=2,DIQ010=2)
newHL[,newfc] <- lapply(newHL[,newfc],factor)
P7 <- plotrcs(fitLF_PPB,newHL,main_var = 'log(URXPPB)',x.axis.min  =-4.5 ,x.axis.max = 8,y.axis.min = 2.5,y.axis.max = 3.8,title = '对羟基苯甲酸丙酯与低频听力阈值',lab.y = '听力阈值（ln（dB））',lab.x ='ֵ浓度（g/g肌酐）',pwald =round(P_LF_rcs_adj[7],3), pwald.x = 3,pwald.y =3.8)



jpeg(file='低频非线性.jpg',width=15000,height=8500,res=600)
ggdraw() +
  draw_plot(P1, x = 0, y = .5, width = .3, height = .5) +
  draw_plot(P2, x = .25, y = .5, width = .3, height = .5) +
  draw_plot(P3, x = .5, y = .5, width = .3, height = .5) +
  draw_plot(P4, x = .75, y = .5, width = .3, height = .5) +
  draw_plot(P5, x = 0, y = .0, width = .3, height = .5) +
  draw_plot(P6, x = .25, y = .0, width = .3, height = .5) +
  draw_plot(P7, x = .5, y = .0, width = .3, height = .5) +
  draw_plot_label(label = c("A","B","C",'D','E','F','G'), size = 22,
                  x = c(0.05,0.3,0.55,.8,0.05,.3,.55), y = c(0.9,0.9,0.9,0.9,0.4,0.4,0.4))
dev.off()
########bkmr
YH <- log(hea$HF+12)
YL <- log(hea$LF+11)
expo <-as.matrix(scale(log(hea[,3:9])))
set.seed(123)
cora <- sapply(hea[,c('RIDAGEYR','PIR','Gender','race','URXUCR','medicine','BMXBMI','smoke','BPQ020','DIQ010','stroke','HartF')],as.numeric)
fitkmH <- kmbayes(YH,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5))
PIPHF <- ExtractPIPs(fitkmH)
PIPHF[,c(3,4)] <- round(PIPHF[,c(3,4)],2)
PIPHF$variable <- c('2,5-二氯苯酚','2,4-二氯苯酚','双酚A','二苯甲酮-3','三氯生','对羟基苯甲酸甲脂','对羟基苯甲酸丙脂')
#######
risks.overall <- OverallRiskSummaries(fit = fitkmH,q.fixed = 0.25)
Ptot1 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+ylim(-0.15,0.03)+xlim(0.23,0.77)+
  geom_pointrange()+annotation_custom(grob = tableGrob(PIPHF, rows = NULL),
                                      xmin = 0.368,xmax = 0.369,
                                      ymin = -0.105,ymax = -0.1055)+ylab('Estimation')+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"),axis.text.y = element_text(size=15,
                                                                               color = "black"),
        axis.text.x = element_text(size=15,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmH,q.fixed = 0.25)
pred.resp.univar <- mutate_at(pred.resp.univar,vars('variable'),
                              ~case_when(.=='URX14D'~'2,5-二氯苯酚',
                                         .=='URXDCB'~'2,4-二氯苯酚',
                                         .=='URXBPH'~'双酚A',
                                         .=='URXBP3'~'二苯甲酮-3',
                                         .=='URXTRS'~'三氯生',
                                         .=='URXMPB'~'对羟基苯甲酸甲脂',
                                         .=='URXPPB'~'对羟基苯甲酸丙脂'))
Ptot2 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"),axis.text.y = element_text(size=15,
                                                                                                    color = "black"),
                              axis.text.x = element_text(size=15,
                                                         color = "black"))
#####
fitkmL <- kmbayes(YL,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5),control.params=list(a.p0=10,b.p0=10))
PIPLF <- ExtractPIPs(fitkmL)
PIPLF[,c(3,4)] <- round(PIPLF[,c(3,4)],2)
PIPLF$variable <- c('2,5-二氯苯酚','2,4-二氯苯酚','双酚A','二苯甲酮-3','三氯生','对羟基苯甲酸甲脂','对羟基苯甲酸丙脂')
###????
risks.overall <- OverallRiskSummaries(fit = fitkmL,q.fixed = 0.25)
Ptot3 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+ylim(-0.15,0.03)+xlim(0.23,0.77)+
  geom_pointrange()+annotation_custom(grob = tableGrob(PIPLF, rows = NULL),
                                      xmin = 0.368,xmax = 0.369,
                                      ymin = -0.105,ymax = -0.1055)+ylab('Estimation')+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"),axis.text.y = element_text(size=15,
                                                                               color = "black"),
        axis.text.x = element_text(size=15,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmL)
pred.resp.univar <- mutate_at(pred.resp.univar,vars('variable'),
                              ~case_when(.=='URX14D'~'2,5-二氯苯酚',
                                         .=='URXDCB'~'2,4-二氯苯酚',
                                         .=='URXBPH'~'双酚A',
                                         .=='URXBP3'~'二苯甲酮-3',
                                         .=='URXTRS'~'三氯生',
                                         .=='URXMPB'~'对羟基苯甲酸甲脂',
                                         .=='URXPPB'~'对羟基苯甲酸丙脂'))
Ptot4 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"),axis.text.y = element_text(size=15,
                                                                                                    color = "black"),
                              axis.text.x = element_text(size=15,
                                                         color = "black"))
jpeg(file='图5.jpg',width=15000,height=8500,res=600)
ggdraw() +
  draw_plot(Ptot1, x = 0, y = .5, width = .35, height = .48) +
  draw_plot(Ptot2, x = 0, y = .0, width = .35, height = .48) +
  draw_plot(Ptot3, x = .35, y = .5, width = 0.35, height = 0.48) +
  draw_plot(Ptot4, x = .35, y = .0, width = 0.35, height = 0.48) + 
  draw_plot_label(label = c("A","B","D",'D'), size = 22,
                  x = c(0,0.35,0,0.35), y = c(1,1,0.5,0.5))
dev.off()
#######################################################
data <- hea[,c(1,3:9)]
rownames(data) <- data[,1]
data <- data[,-1]
data <- scale(log(data))
#################
par(mfrow=c(1,2))
fviz_nbclust(data, cluster::pam, method = "gap_stat",k.max = 10) +
  geom_vline(xintercept = 4, linetype = 2)######图6左
fviz_nbclust(data, cluster::pam, method = "wss",k.max = 10) +
  geom_vline(xintercept = 4, linetype = 2)#####图6右
########
cl <- pam(data,k=4)
cluster <- cl[["clustering"]]
data <- as.data.frame(data)
data$cluster <- cluster
data$cluster <- factor(data$cluster,order=F)
####
set.seed(123)
tsne <- Rtsne(data,dims=3,check_duplicates = FALSE)
tsnedata <- tsne$Y
colnames(tsnedata) <- c('tsne-x','tsne-y','tsne-z')
colors_border=c( rgb(1,0,0,0.9), rgb(0,0,0.55,0.9) ,rgb(0,1,0,0.9), rgb(1 ,1, 0,0.9) )
colors_in=c( rgb(1,0,0,0.5), rgb(0,0,0.55,0.3) , rgb(0,1,0,0.6) ,rgb(1,1, 0,0.5))
color.lib <- colors_border
color <- color.lib[as.numeric(data$cluster)]
radar <- data[,-9]
radar_ccol <- names(radar)[-8]
mean_ra <- as.data.frame(sapply(radar_ccol, function(x) tapply(radar[,x],radar[,'cluster'],mean)))
rownames(mean_ra) <- paste0('集群',1:4)
mean_ra <- rbind(rep(-1.05,7),mean_ra)
mean_ra <- rbind(rep(1.8,7),mean_ra)
colnames(mean_ra) <-  c('2,5-二氯苯酚','2,4-二氯苯酚','双酚A','二苯甲酮-3','三氯生','对羟基苯甲酸甲脂','对羟基苯甲酸丙脂')
jpeg('图7.jpg',width=15000,height=8500,res=600)
par(mfrow=c(1,2),cex.main=2.5)
scatterplot3d(tsnedata,color = color,angle = 45,pch=16,
              cex.symbols  = 2,main = 'A                                                                      ')
radarchart(mean_ra, axistype=1 , 
           pcol= colors_border, pfcol=colors_in , plwd=4 , plty=1,
           
           cglcol="grey60", cglty=1, axislabcol="grey", cglwd=0.8,
           
           vlcex=2 ,caxislabels = c(-1, -0.5,0, 0.5, 1.8)
           ,title ='B                                                             ')
legend("bottomleft",legend = levels(data$cluster),col = color.lib,pch = 16,
       inset = -0.04,xpd = T,horiz = T,cex = 2.5,title = '集群',text.col='black')
dev.off()
############
data$SEQN <- as.numeric(rownames(data))
all <- left_join(all,data[,c(8,9)],by='SEQN')
all$cluster <- relevel(all$cluster,ref = 3)
###????????չʾ
hea<- filter(all,RIDAGEYR>=45) %>% na.omit(.)
heamodel <- mutate(all,inAnalysis= SEQN%in%hea$SEQN)
heaM <- svydesign(data=heamodel, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC12YR, nest=TRUE)
heamod <- subset(heaM, inAnalysis)
fitHF <- svyglm(log(HF+12)~cluster+Gender+RIDAGEYR+education+PIR+race+BMI+medicine+URXUCR+smoke+stroke+HartF+BPQ020+DIQ010
                ,design=heamod)
CIHF <- signif(confint(fitHF)[2:4,],2)
betaHF <- signif(coef(fitHF)[2:4],2)
betaCIHF <- paste0(betaHF,' (',CIHF[,1],', ',CIHF[,2],')')
fitLF <- svyglm(log(LF+11)~cluster+Gender+RIDAGEYR+education+URXUCR+PIR+race+medicine+BMI+smoke+stroke+HartF+BPQ020+DIQ010
                ,design=heamod)
betaLF <- summary(fitLF)
CILF <- signif(confint(fitLF)[2:4,],2)
betaLF <- signif(coef(fitLF)[2:4],2)
betaCILF <- paste0(betaLF,' (',CILF[,1],', ',CILF[,2],')')
CIdata <- data.frame(betaCIHF,betaCILF)
CIdata <- add_row(.data = CIdata,.before = 3,betaCIHF='Reference',betaCILF='Reference')
rownames(CIdata) <- paste0('Cluster',1:4)
colnames(CIdata) <- c('HF ?? (95% CI)','LF ?? (95% CI)')
PHF <- signif(summary(fitHF)$coefficients[,4][2:4],3)
PHF <- c(PHF[1:2],'Reference',PHF[3])
PHF[c(1,2,4)] <- paste0('P=',PHF[c(1,2,4)])
PLF <-  signif(summary(fitLF)$coefficients[,4][2:4],3)
PLF <- c(PLF[1:2],'Reference',PLF[3])
PLF[c(1,2,4)] <- paste0('P=',PLF[c(1,2,4)])
HF_GLM <- as.data.frame(lsmeans(fitHF,"cluster",rg.limit=90000))[,c(1,2,5,6)]
LF_GLm <- as.data.frame(lsmeans(fitLF,'cluster',rg.limit=90000))[,c(1,2,5,6)]
clu_names <- c("cluster","lsmean","LCL","UCL")
HF_GLM[,2:4] <- exp(HF_GLM[,2:4])-12
LF_GLm[,2:4] <- exp(LF_GLm[,2:4])-11
names(HF_GLM) <- clu_names
names(LF_GLm) <- clu_names
HF_GLM$cluster <- paste0('cluster',HF_GLM$cluster)
LF_GLm$cluster <- paste0('cluster',LF_GLm$cluster)
p_glmhf <- ggplot(HF_GLM, aes(lsmean,cluster))+ geom_point(size=3.6)+
  geom_errorbarh(aes(xmax =LCL, xmin = UCL), height = 0.4)+labs(x='高频听力阈值最小二乘几何均值 (dB)')+xlim(10,40)+
  theme(
    axis.title.x = element_text(size = 20,color = "black",face="bold"),
    axis.title.y=element_text(size=20,color = "black",face="bold"),
    axis.text.x =element_text(size = 15,colour = 'black',face="bold"), 
    axis.text.y =element_text(size = 15,colour = 'black',face="bold"))+ annotate("text",x=40,y=1:4,label=PHF,size=5.5,alpha=0.9,color="black")+coord_flip()
p_glmlf <- ggplot(LF_GLm, aes(lsmean,cluster))+ geom_point(size=3.6)+
  geom_errorbarh(aes(xmax =LCL, xmin = UCL), height = 0.4)+labs(x='低频听力阈值最小二乘几何均值 (dB)')+xlim(10,40)+
  theme(
    axis.title.x = element_text(size = 20,color = "black",face="bold"),
    axis.title.y=element_text(size=20,color = "black",face="bold"),
    axis.text.x =element_text(size = 15,colour = 'black',face="bold"), 
    axis.text.y =element_text(size = 15,colour = 'black',face="bold"))+ annotate("text",x=25,y=1:4,label=PLF,size=5.5,alpha=0.9,color="black")+coord_flip()
jpeg(file='图8ֵ.jpg',width=15000,height=8500,res=600)
ggdraw( xlim = c(0, 1.2)) +
  draw_plot(p_glmhf, x = 0, y = .2, width = .5, height = .6) +
  draw_plot(p_glmlf, x = 0.6, y = .2, width = .5, height = .6) +
  draw_plot_label(label = c("A","B"), size = 28,
                  x = c(0,0.6), y = c(0.9,0.9))

dev.off()
######
table_cluster <- svyCreateTableOne(vars = c('race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF','Gender','Age.Group','BMI','URXUCR'),
                                   data = heamod,factorVars = c('race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF','Gender','Age.Group','BMI'),strata = 'cluster')
tab1 <- as.data.frame(print(table_cluster, showAllLevels = TRUE))
table_clustern <- CreateTableOne(vars = c('race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF','Gender','Age.Group','BMI','URXUCR'),
                                 data = hea,strata = 'cluster')
tab2 <- as.data.frame(print(table_clustern, showAllLevels = TRUE))

clustertab <- creat_newtab(tab1 = tab2,tab_isn = tab1)
clustertab[29,] <- tab1[29,]
write.csv(clustertab,'表4.csv')
########################敏感性分析
all <- left_join(EPH,DEMO,by='SEQN')%>% left_join(.,BPQ,by='SEQN') %>% left_join(.,DIQ,by='SEQN') %>% left_join(.,MCQ,by='SEQN') %>% 
  left_join(.,RXQ,by='SEQN') %>% left_join(.,SMQ,by='SEQN') %>% left_join(.,BMX,by='SEQN') %>% 
  left_join(.,select(AUX1,SEQN,LFR,HFR),by='SEQN')%>% left_join(.,UCR,by='SEQN')
all <- filter(all,!is.na(all$WTSB2YR))
aa <- c('URX14D','URXDCB','URXBPH','URXBP3','URXTRS','URXMPB','URXPPB')
logaa <- paste0('log(',aa,')')
ar <- paste0(aa,'gr')
all <- mutate(all,WTMEC12YR=WTSB2YR/6)%>%
  mutate(.,DIQ010=ifelse(DIQ010==3,1,ifelse(DIQ010==9|DIQ010==7,2,DIQ010)))%>%
  mutate(.,BPQ020=ifelse(BPQ020==9,2,BPQ020))%>%
  mutate(.,HartF=ifelse(HartF==9|HartF==7,2,HartF))%>%
  mutate(.,stroke=ifelse(stroke==9|stroke==7,2,stroke))%>%
  mutate(.,DMDEDUC2=ifelse(DMDEDUC2==9|DMDEDUC2==7,NA,DMDEDUC2))%>%
  mutate(.,Gender = factor(RIAGENDR, levels = c(1,2),labels=c("Men", "Women")),
         Age.Group = cut(RIDAGEYR, breaks=c(-Inf,19,44,59,Inf),labels=c(NA, NA,"45-59","60 and over")),
         BMI=cut(BMXBMI,breaks = c(-Inf,24.99,29.99,Inf),labels = c(1,2,3)),
         race=case_when(RIDRETH1==3~0,RIDRETH1==4~-1,RIDRETH1==1|RIDRETH1==2|RIDRETH1==5~1),
         PIR=case_when(INDFMPIR<=1~1,INDFMPIR>1&INDFMPIR<=3~2,INDFMPIR>3~3)
  )%>%
  mutate(.,education=case_when(DMDEDUC2==1|DMDEDUC2==2~1,DMDEDUC2==3~2,DMDEDUC2==4~3,DMDEDUC2==5~4))
fac <- c('race','Age.Group','education','RIDRETH1','DMDEDUC2','Gender','medicine','smoke','BPQ020','DIQ010','stroke','HartF','BMI','PIR')
all[,fac] <- lapply(all[,fac], factor)
##############
hea<- filter(all,RIDAGEYR>=45) %>% na.omit(.)
all <- mutate(all,URX14Dgr=ntile(URX14D,4),
              URXDCBgr=ntile(URXDCB,4),
              URXBPHgr=ntile(URXBPH,4),
              URXBP3gr=ntile(URXBP3,4),
              URXTRSgr=ntile(URXTRS,4),
              URXMPBgr=ntile(URXMPB,4),
              URXPPBgr=ntile(URXPPB,4))
all[,ar] <- lapply(all[,ar],factor)
heamodel <- mutate(all,inAnalysis= SEQN%in%hea$SEQN)
heaM <- svydesign(data=heamodel, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC12YR, nest=TRUE)
heamod <- subset(heaM, inAnalysis)
##???ݴ???
YHR <- log(hea$HFR+12)
YLR <- log(hea$LFR+11)
expo <-as.matrix(scale(log(hea[,3:9])))
set.seed(123)
cora <- sapply(hea[,c('RIDAGEYR','PIR','Gender','race','URXUCR','medicine','BMXBMI','smoke','BPQ020','DIQ010','stroke','HartF')],as.numeric)
#####
fitkmH2 <- kmbayes(YHR,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5))
PIPHF2 <- ExtractPIPs(fitkmH2)
PIPHF2[,c(3,4)] <- round(PIPHF2[,c(3,4)],2)
PIPHF2$variable <- gsub('URX','',PIPHF2$variable)
#####
risks.overall <- OverallRiskSummaries(fit = fitkmH2,q.fixed = 0.25)
Ptot1 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+ylim(-0.01,0.007)+xlim(0.23,0.77)+
  geom_pointrange()+annotation_custom(grob = tableGrob(PIPHF2, rows = NULL),
                                      xmin = 0.368,xmax = 0.369,
                                      ymin = -0.105,ymax = -0.1055)+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmH2,q.fixed = 0.25)
Ptot2 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"))

#####
pred.resp.bivar <- PredictorResponseBivar(fit = fitkmH2, min.plot.dist = 1)
pred.resp.bivar.levels <- PredictorResponseBivarLevels(
  pred.resp.df = pred.resp.bivar,Z=expo,qs = c(0.1, 0.5, 0.9))
jpeg(file='?̶?.jpg',width=15000,height=8500,res=600)
ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")
dev.off()
#####
tab_model
fitkmL2 <- kmbayes(YLR,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5))
PIPLF2 <- ExtractPIPs(fitkmL2)
PIPLF2[,c(3,4)] <- round(PIPLF2[,c(3,4)],2)
PIPLF2$variable <- gsub('URX','',PIPLF2$variable)
###????
risks.overall <- OverallRiskSummaries(fit = fitkmL2,q.fixed = 0.25)
Ptot3 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+ylim(-0.005,0.003)+xlim(0.23,0.77)+
  geom_pointrange()+annotation_custom(grob = tableGrob(PIPLF2, rows = NULL),
                                      xmin = 0.368,xmax = 0.369,
                                      ymin = -0.105,ymax = -0.1055)+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmL2)
Ptot4 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"))
####
pred.resp.bivar <- PredictorResponseBivar(fit = fitkmL2, min.plot.dist = 1)
pred.resp.bivar.levels <- PredictorResponseBivarLevels(
  pred.resp.df = pred.resp.bivar,Z=expo,qs = c(0.1, 0.5, 0.9))
jpeg(file='?̶?2.jpg',width=15000,height=8500,res=600)
ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")
dev.off()

jpeg(file='图.jpg',width=15000,height=8500,res=600)
ggdraw() +
  draw_plot(Ptot1, x = 0, y = .5, width = .35, height = .48) +
  draw_plot(Ptot2, x = 0, y = .0, width = .35, height = .48) +
  draw_plot(Ptot3, x = .5, y = .5, width = 0.35, height = 0.48) +
  draw_plot(Ptot4, x = .5, y = .0, width = 0.35, height = 0.48) + 
  draw_plot_label(label = c("A","B","C",'D'), size = 22,
                  x = c(0,0.5,0,0.5), y = c(1,1,0.5,0.5))
dev.off()
all <- left_join(EPHun,DEMO,by='SEQN')%>% left_join(.,BPQ,by='SEQN') %>% left_join(.,DIQ,by='SEQN') %>% left_join(.,MCQ,by='SEQN') %>% 
  left_join(.,RXQ,by='SEQN') %>% left_join(.,SMQ,by='SEQN') %>% left_join(.,BMX,by='SEQN') %>% 
  left_join(.,select(AUX1,SEQN,LF,HF),by='SEQN')%>% left_join(.,UCR,by='SEQN')
all <- filter(all,!is.na(all$WTSB2YR))
all <- mutate(all,WTMEC12YR=WTSB2YR/6)%>%
  mutate(.,DIQ010=ifelse(DIQ010==3,1,ifelse(DIQ010==9|DIQ010==7,2,DIQ010)))%>%
  mutate(.,BPQ020=ifelse(BPQ020==9,2,BPQ020))%>%
  mutate(.,HartF=ifelse(HartF==9|HartF==7,2,HartF))%>%
  mutate(.,stroke=ifelse(stroke==9|stroke==7,2,stroke))%>%
  mutate(.,DMDEDUC2=ifelse(DMDEDUC2==9|DMDEDUC2==7,NA,DMDEDUC2))%>%
  mutate(.,Gender = factor(RIAGENDR, levels = c(1,2),labels=c("Men", "Women")),
         Age.Group = cut(RIDAGEYR, breaks=c(-Inf,19,44,59,Inf),labels=c(NA, NA,"45-59","60 and over")),
         BMI=cut(BMXBMI,breaks = c(-Inf,24.99,29.99,Inf),labels = c(1,2,3)),
         race=case_when(RIDRETH1==3~0,RIDRETH1==4~-1,RIDRETH1==1|RIDRETH1==2|RIDRETH1==5~1),
         PIR=case_when(INDFMPIR<=1~1,INDFMPIR>1&INDFMPIR<=3~2,INDFMPIR>3~3)
  )%>%
  mutate(.,education=case_when(DMDEDUC2==1|DMDEDUC2==2~1,DMDEDUC2==3~2,DMDEDUC2==4~3,DMDEDUC2==5~4))
fac <- c('race','Age.Group','education','RIDRETH1','DMDEDUC2','Gender','medicine','smoke','BPQ020','DIQ010','stroke','HartF','BMI','PIR')
all[,fac] <- lapply(all[,fac], factor)
##################????
hea<- filter(all,RIDAGEYR>=45) %>% na.omit(.)
all <- mutate(all,URX14Dgr=ntile(URX14D,4),
              URXDCBgr=ntile(URXDCB,4),
              URXBPHgr=ntile(URXBPH,4),
              URXBP3gr=ntile(URXBP3,4),
              URXTRSgr=ntile(URXTRS,4),
              URXMPBgr=ntile(URXMPB,4),
              URXPPBgr=ntile(URXPPB,4))
all[,ar] <- lapply(all[,ar],factor)
heamodel <- mutate(all,inAnalysis= SEQN%in%hea$SEQN)
heaM <- svydesign(data=heamodel, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC12YR, nest=TRUE)
heamod <- subset(heaM, inAnalysis)
Me <- meff(log(hea[,3:9]))

all_line_results_un <- lapply(mian_test,function(x) glp(a=logaa,ma=x,coa = coa,
                                                        #family = "quasibinomial",
                                                        design =heamod))
names(all_line_results) <- mian_test
all_rcs_results_un <- lapply(mian_test,function(x) glp(a=logaa,ma=x,coa = coa,
                                                       #family = "quasibinomial",
                                                       design =heamod))
names(all_rcs_results_un) <- mian_test
########
##???ݴ???
YH <- log(hea$HF+12)
YL <- log(hea$LF+11)
expo <-as.matrix(scale(log(hea[,3:9])))
set.seed(123)
cora <- sapply(hea[,c('RIDAGEYR','PIR','Gender','race','URXUCR','medicine','BMXBMI','smoke','BPQ020','DIQ010','stroke','HartF')],as.numeric)
#####??ģ
###
fitkmH3 <- kmbayes(YH,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5))
PIPHF3 <- ExtractPIPs(fitkmH3)
PIPHF3[,c(3,4)] <- round(PIPHF3[,c(3,4)],2)
PIPHF3$variable <- c('2,5-dichlorophenol','2,4-dichlorophenol','BPA','benzophenone-3','TCS','Mep','PrP')
###
risks.overall <- OverallRiskSummaries(fit = fitkmH3,q.fixed = 0.25)
Ptot1 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+
  geom_pointrange()+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmH3,q.fixed = 0.25)
Ptot2 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"))
######
fitkmL3 <- kmbayes(YL,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5))
PIPLF3 <- ExtractPIPs(fitkmL3)
PIPLF3[,c(3,4)] <- round(PIPLF3[,c(3,4)],2)
PIPLF3$variable <- c('2,5-dichlorophenol','2,4-dichlorophenol','BPA','benzophenone-3','TCS','Mep','PrP')
######
risks.overall <- OverallRiskSummaries(fit = fitkmL3,q.fixed = 0.25)
Ptot3 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+
  geom_pointrange()+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmL3)
Ptot4 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"))
####
jpeg(file='图.jpg',width=15000,height=8500,res=600)
ggdraw() +
  draw_plot(Ptot1, x = 0, y = .5, width = .35, height = .48) +
  draw_plot(Ptot2, x = 0, y = .0, width = .35, height = .48) +
  draw_plot(Ptot3, x = .5, y = .5, width = 0.35, height = 0.48) +
  draw_plot(Ptot4, x = .5, y = .0, width = 0.35, height = 0.48) + 
  draw_plot_label(label = c("A","B","C",'D'), size = 22,
                  x = c(0,0.5,0,0.5), y = c(1,1,0.5,0.5))
dev.off()
####
hea_4k<- filter(all,RIDAGEYR>=45&noise_4k==0) %>% na.omit(.)
YH <- log(hea_4k$HF+12)
YL <- log(hea_4k$LF+11)
expo <-as.matrix(scale(log(hea_4k[,3:9])))
set.seed(123)
cora <- sapply(hea_4k[,c('RIDAGEYR','PIR','Gender','race','URXUCR','medicine','BMXBMI','smoke','BPQ020','DIQ010','stroke','HartF')],as.numeric)
fitkmL_4k <- kmbayes(YL,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5))
PIPLF_4k <- ExtractPIPs(fitkmL_4k)
PIPLF_4k[,c(3,4)] <- round(PIPLF3[,c(3,4)],2)
PIPLF_4k$variable <- c('2,5-dichlorophenol','2,4-dichlorophenol','BPA','benzophenone-3','TCS','Mep','PrP')
######
risks.overall <- OverallRiskSummaries(fit = fitkmL_4k,q.fixed = 0.25)
Ptot5 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+
  geom_pointrange()+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmL_4k)
Ptot6 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"))
fitkmH_4k <- kmbayes(YH,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5))
PIPHF_4k <- ExtractPIPs(fitkmH_4k)
PIPHF_4k[,c(3,4)] <- round(PIPHF_4k[,c(3,4)],2)
PIPHF_4k$variable <- c('2,5-dichlorophenol','2,4-dichlorophenol','BPA','benzophenone-3','TCS','Mep','PrP')
risks.overall <- OverallRiskSummaries(fit = fitkmH_4k,q.fixed = 0.25)
Ptot7 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+
  geom_pointrange()+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmH_4k)
Ptot8 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"))
######
hea_15<- filter(all,RIDAGEYR>=45&SDDSRVYR!=9) %>% na.omit(.)
YH <- log(hea_15$HF+12)
YL <- log(hea_15$LF+11)
expo <-as.matrix(scale(log(hea_15[,3:9])))
set.seed(123)
cora <- sapply(hea_15[,c('RIDAGEYR','PIR','Gender','race','URXUCR','medicine','BMXBMI','smoke','BPQ020','DIQ010','stroke','HartF')],as.numeric)
fitkmL_15 <- kmbayes(YL,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5))
PIPLF_15 <- ExtractPIPs(fitkmL_15)
PIPLF_15[,c(3,4)] <- round(PIPLF_15[,c(3,4)],2)
PIPLF_15$variable <- c('2,5-dichlorophenol','2,4-dichlorophenol','BPA','benzophenone-3','TCS','Mep','PrP')
risks.overall <- OverallRiskSummaries(fit = fitkmL_15,q.fixed = 0.25)
Ptot9 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+
  geom_pointrange()+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmL_15)
Ptot10 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"))
fitkmH_15 <- kmbayes(YH,Z=expo,X=cora,iter = 20000,varsel = T,groups = c(1,1,2,3,4,5,5))
PIPHF_15 <- ExtractPIPs(fitkmH_15)
PIPHF_15[,c(3,4)] <- round(PIPHF_15[,c(3,4)],2)
PIPHF_15$variable <- c('2,5-dichlorophenol','2,4-dichlorophenol','BPA','benzophenone-3','TCS','Mep','PrP')
risks.overall <- OverallRiskSummaries(fit = fitkmH_15,q.fixed = 0.25)
Ptot11 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + geom_hline(yintercept=0,color='red',linetype="dashed")+
  geom_pointrange()+
  theme(axis.title.y=element_text(size=20,
                                  color = "black"),
        axis.title.x =element_text(size=20,
                                   color = "black"))
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmH_15)
Ptot12 <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")+xlab('')+theme(axis.title.y=element_text(size=20,
                                                        color = "black"))


#########第三部分
yrs <- c('2011-2012','2013-2014','2015-2016')
letters <- c('_G','_H','_I')
DEMO <- downloadNHANES('DEMO')
DEMO <- DEMO[,c("SEQN",'SDDSRVYR',"RIAGENDR","RIDAGEYR","SDMVSTRA","DMDEDUC2","RIDRETH1","SDMVPSU",'INDFMPIR')]
######
BPQ <- downloadNHANES('BPQ')
BPQ <- BPQ[,c("SEQN","BPQ020")]
#####
DIQ <- downloadNHANES('DIQ')
DIQ <- DIQ[,c("SEQN","DIQ010")]
######
SMQ <- downloadNHANES('SMQ')
SMQ <- SMQ[,c("SEQN","SMQ020","SMQ040")]
SMQ <- mutate(SMQ,smoke=case_when(SMQ020==2~1,SMQ020==1&(SMQ040==1|SMQ040==2)~3,SMQ020==1&SMQ040==3~2))%>%select(.,SEQN,smoke)
#####
BMX <- downloadNHANES('BMX')
BMX <- BMX[,c("SEQN","BMXBMI")]
#####
RXQ <- downloadNHANES('RXQ_RX')
RXQ <- RXQ[,c("SEQN","RXDUSE","RXDDRGID")]
RXQ <- mutate(RXQ,RXDUSE=ifelse(RXDUSE==7|RXDUSE==9,NA,RXDUSE))
RXQ <- mutate(RXQ,medicine=case_when(RXDDRGID%in%du~1,!RXDDRGID%in%du~0))%>%
  mutate(.,medicine=ifelse(is.na(RXDUSE),NA,medicine))%>%select(.,SEQN,medicine)
undup <- as.data.frame(tapply(RXQ$medicine,RXQ$SEQN, function(x) mean(x,na.rm = T)))
undup$SEQN <- rownames(undup)
rownames(undup) <- 1:nrow(undup)
colnames(undup)[1] <-"medicine"
undup <- mutate(undup,medicine=ifelse(medicine==0,0,1))
RXQ <- undup
RXQ$SEQN <- as.numeric(RXQ$SEQN)
rm(undup)
####
MCQ <- downloadNHANES('MCQ')
MCQ <- MCQ[,c("SEQN","MCQ160B","MCQ160F")]
colnames(MCQ) <- c("SEQN","HartF","stroke")
#####
######
ORA <- downloadNHANES('OHQ')
ORA <- ORA[,c("SEQN",'OHQ845')]
####
UCR <- downloadNHANES('ALB_CR')
UCR <- UCR[,c('SEQN','URXUCR')]
####
yrs <- c('2011-2012')
letters <- c('_G')
PP1112 <- downloadNHANES('PP')
PP1112 <- PP1112[,c("SEQN",'WTSA2YR','URX14D','URD14DLC','URXDCB','URDDCBLC')]
names(PP1112)[2] <- 'WTSB2YR'
EP1112 <- downloadNHANES('EPH')
EP1112 <- EP1112[,c("SEQN",'URXBPH','URDBPHLC','URXBP3','URDBP3LC','URXTRS','URDTRSLC','URXMPB','URDMPBLC','URXPPB','URDPPBLC')]
EPH1112 <- left_join(PP1112,EP1112,by='SEQN')
yrs <- c('2013-2014','2015-2016')
letters <- c('_H','_I')
EPH1316 <- downloadNHANES('EPHPP')
EPH1316 <- EPH1316[,c("SEQN",'WTSB2YR','URX14D','URD14DLC','URXDCB','URDDCBLC','URXBPH','URDBPHLC','URXBP3','URDBP3LC','URXTRS','URDTRSLC','URXMPB','URDMPBLC','URXPPB','URDPPBLC')]
EPH <- rbind(EPH1112,EPH1316)
EPH_ID <-select(EPH,ends_with('LC'))
EPH_ID <- na.omit(EPH_ID)
rio <- sapply(apply(EPH_ID,2,sum), function(x) 1-(x/nrow(EPH_ID)))
EPHun <- select(EPH,-ends_with('LC'))
EPH <- left_join(EPHun,UCR,by='SEQN')
EPH[,3:9] <- sapply(EPH[,3:9],function(x)  round(x/EPH[,10]*100,2))
EPH <- EPH[,-10]
write.csv(EPH,'EPH.csv')
########听力
yrs <- c('2011-2012','2015-2016')
letters <- c('_G','_I')
AUX <- downloadNHANES('AUX')
AUX0 <- select(AUX,SEQN,AUXOTSPL,AUXROTSP,AUXTCOMR,AUXTCOML,AUATYMTL,AUATYMTR,AUXTMEPR,AUXTMEPL,starts_with('AUXU')) %>%select(.,-AUXU1K2R,-AUXU1K2L)%>%
  filter(.,AUXROTSP==1&AUXTCOMR>0.3&AUXOTSPL==1&AUXTCOML>0.3) %>% 
  mutate_at(vars(starts_with('AUXU')),~ifelse(.==666,120,.)) %>% 
  mutate_at(vars(starts_with('AUXU')),~ifelse(.==888|.==777,NA,.))%>% 
  mutate_at(vars(starts_with('AUXTC')),~ifelse(.==888|.==777,NA,.))%>% 
  mutate(.,LFR=apply(select(.,AUXU1K1R,AUXU500R,AUXU2KR),1,mean),
         LFL=apply(select(.,AUXU1K1L,AUXU500L,AUXU2KL),1,mean),
         HFR=apply(select(.,AUXU3KR,AUXU4KR,AUXU6KR,AUXU8KR),1,mean),
         HFL=apply(select(.,AUXU3KL,AUXU4KL,AUXU6KL,AUXU8KL),1,mean)) %>%
  mutate(noise_R=ifelse(AUXU3KR<=AUXU4KR&AUXU4KR>=AUXU6KR&AUXU4KR>25,1,0),
         noise_L=ifelse(AUXU3KL<=AUXU4KL&AUXU4KL>=AUXU6KL&AUXU4KL>25,1,0)) %>% 
  mutate(.,LF=ifelse(is.na(LFR),ifelse(is.na(LFL),NA,LFL),ifelse(is.na(LFL),LFR,ifelse(LFR>LFL,LFR,LFL))),
         HF=ifelse(is.na(HFR),ifelse(is.na(HFL),NA,HFL),ifelse(is.na(HFL),HFR,ifelse(HFR>HFL,HFR,HFL))),
         noise_4k=ifelse(is.na(noise_R),ifelse(is.na(noise_L),NA,noise_L),ifelse(is.na(noise_L),noise_R,ifelse(noise_L==1|noise_R==1,1,0)))) %>% 
  mutate(.,LFHL=ifelse(LF>25,1,0),
         HFHL=ifelse(HF>25,1,0))
AUX1 <- select(AUX0,SEQN,LF,HF,LFR,HFR,LFL,HFL,LFHL,HFHL,noise_4k,noise_R,noise_L)
#####
all <- left_join(DEMO,EPH,by='SEQN')%>% left_join(.,BPQ,by='SEQN') %>% left_join(.,DIQ,by='SEQN') %>% left_join(.,MCQ,by='SEQN') %>% 
  left_join(.,RXQ,by='SEQN') %>% left_join(.,SMQ,by='SEQN') %>% left_join(.,BMX,by='SEQN') %>% 
  left_join(.,dplyr::select(AUX1,SEQN,LF,HF,HFHL,LFHL),by='SEQN')%>% left_join(.,UCR,by='SEQN') %>% left_join(.,ORA,by='SEQN')

all <- filter(all,!is.na(all$WTSB2YR))
all <- mutate(all,URXTRSgr=ntile(URXTRS,4))
all <- mutate(all,WTMEC6YR=WTSB2YR/3)%>%
  mutate(.,DIQ010=ifelse(DIQ010==3,1,ifelse(DIQ010==9|DIQ010==7,2,DIQ010)))%>%
  mutate(.,BPQ020=ifelse(BPQ020==9,2,BPQ020))%>%
  mutate(.,HartF=ifelse(HartF==9|HartF==7,2,HartF))%>%
  mutate(.,stroke=ifelse(stroke==9|stroke==7,2,stroke))%>%
  mutate(.,DMDEDUC2=ifelse(DMDEDUC2==9|DMDEDUC2==7,NA,DMDEDUC2))%>%
  mutate(.,Gender = factor(RIAGENDR, levels = c(1,2),labels=c("Men", "Women")),
         Age.Group = cut(RIDAGEYR, breaks=c(-Inf,19,44,59,Inf),labels=c(NA, NA,"45-59","60 and over")),
         BMI=cut(BMXBMI,breaks = c(-Inf,24.99,29.99,Inf),labels = c(1,2,3)),
         race=case_when(RIDRETH1==3~0,RIDRETH1==4~-1,RIDRETH1==1|RIDRETH1==2|RIDRETH1==5~1),
         PIR=case_when(INDFMPIR<=1~1,INDFMPIR>1&INDFMPIR<=3~2,INDFMPIR>3~3),
         OHQ845=ifelse(OHQ845>=4,1,0),
         URXTRSgr2=ifelse(URXTRSgr==1,0,1))%>%
  mutate(.,education=case_when(DMDEDUC2==1|DMDEDUC2==2~1,DMDEDUC2==3~2,DMDEDUC2==4~3,DMDEDUC2==5~4))
  
fac <- c('URXTRSgr','race','Age.Group','education','RIDRETH1','DMDEDUC2','Gender','medicine','smoke','BPQ020','DIQ010','stroke','HartF','BMI','PIR','OHQ845')
all[,fac] <- lapply(all[,fac], factor)
#####
hea<- filter(all,RIDAGEYR>=45)%>% na.omit(.)
heamodel <- mutate(all,inAnalysis= SEQN%in%hea$SEQN)
heaM <- svydesign(data=heamodel, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC6YR, nest=TRUE)
heamod <- subset(heaM, inAnalysis)
#######
table <- svyCreateTableOne(vars = c('HF','LF','race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF','Gender','OHQ845','BMXBMI','Age.Group','BMI','PIR','URXUCR','RIDAGEYR'),
                           data = heamod,factorVars = c('race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF','Age.Group','BMI','PIR','OHQ845'))
table_1 <- as.data.frame(print(table, showAllLevels = TRUE))
tablen <- CreateTableOne(vars = c('HF','LF','race','education','medicine','smoke','BPQ020','DIQ010','stroke','HartF','Gender','OHQ845','BMXBMI','Age.Group','BMI','PIR','URXUCR','RIDAGEYR'),
                         data = hea)
tablen_1 <- as.data.frame(print(tablen, showAllLevels = TRUE))
table_tot <- creat_newtab(tab1 = tablen_1,tab_isn = table_1,ins_col = 2)
table_tot[c(2,3,35,36),] <- table_1[c(2,3,35,36),]
write.csv(table_tot,'第三部分人口特征.csv')
#####
DEMO_be <- filter(DEMO,RIDAGEYR>=45)
EPH_be <- na.omit(EPH)
ce_be <- left_join(EPH_be,DEMO_be,by='SEQN')
ce1 <- filter(ce_be,!is.na(RIDAGEYR))
print(length(ce1))
ting_be <- select(AUX,SEQN,AUXOTSPL,AUXROTSP,AUXTCOMR,AUXTCOML,starts_with('AUXU')) %>%select(.,-AUXU1K2R,-AUXU1K2L)
ting <- left_join(ce1,ting_be,by='SEQN')
ce2 <- filter(ting,AUXROTSP==1&AUXTCOMR>0.3&AUXOTSPL==1&AUXTCOML>0.3)
print(length(ce2))
ce2_be <- mutate_at(ce2,vars(starts_with('AUXU')),~ifelse(.==666,120,.)) %>% 
  mutate_at(vars(starts_with('AUXU')),~ifelse(.==888|.==777,NA,.))%>% 
  mutate_at(vars(starts_with('AUXTC')),~ifelse(.==888|.==777,NA,.))%>% 
  mutate(.,LFR=apply(select(.,AUXU1K1R,AUXU500R,AUXU2KR),1,mean),
         LFL=apply(select(.,AUXU1K1L,AUXU500L,AUXU2KL),1,mean),
         HFR=apply(select(.,AUXU3KR,AUXU4KR,AUXU6KR,AUXU8KR),1,mean),
         HFL=apply(select(.,AUXU3KL,AUXU4KL,AUXU6KL,AUXU8KL),1,mean)) %>% 
  mutate(.,LF=ifelse(is.na(LFR),ifelse(is.na(LFL),NA,LFL),ifelse(is.na(LFL),LFR,ifelse(LFR>LFL,LFR,LFL))),
         HF=ifelse(is.na(HFR),ifelse(is.na(HFL),NA,HFL),ifelse(is.na(HFL),HFR,ifelse(HFR>HFL,HFR,HFL)))) %>% 
  mutate(.,LFHL=ifelse(LF>25,1,0),
         HFHL=ifelse(HF>25,1,0)) 
ce3 <- filter(ce2_be,(!is.na(LF))&(!is.na(HF)))
ce4 <- left_join(ce3,ORA,by='SEQN')
print(length(ce3))
ce3_be <- left_join(ce3,BPQ,by='SEQN') %>% left_join(.,DIQ,by='SEQN') %>% left_join(.,MCQ,by='SEQN') %>% 
  left_join(.,RXQ,by='SEQN') %>% left_join(.,SMQ,by='SEQN') %>% left_join(.,BMX,by='SEQN') %>%left_join(.,UCR,by='SEQN')
find_na(ce3_be)
####
fitLF <- svyglm(log(LF+11)~OHQ845+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                    #,family = "quasibinomial"
                ,design=heamod)
summary(fitLF)
PLF <- c('Reference',signif(summary(fitLF)[["coefficients"]][,'Pr(>|t|)'][2],2))
fitHF <- svyglm(log(HF+12)~OHQ845+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                #,family = "quasibinomial"
                ,design=heamod)
summary(fitHF)
PHF <- c('Reference',signif(summary(fitHF)[["coefficients"]][,'Pr(>|t|)'][2],2))


HF_GLM <- as.data.frame(lsmeans(fitHF,"OHQ845",rg.limit=90000))[,c(1,2,5,6)]
LF_GLm <- as.data.frame(lsmeans(fitLF,'OHQ845',rg.limit=90000))[,c(1,2,5,6)]
clu_names <- c("ora_health","lsmean","LCL","UCL")
HF_GLM[,2:4] <- exp(HF_GLM[,2:4])-12
LF_GLm[,2:4] <- exp(LF_GLm[,2:4])-11
names(HF_GLM) <- clu_names
names(LF_GLm) <- clu_names
HF_GLM$ora_health <- c('好','一般/差')
LF_GLm$ora_health <- c('好','一般/差')
p_glmhf <- ggplot(HF_GLM, aes(lsmean,ora_health))+ geom_point(size=3.6)+
  geom_errorbarh(aes(xmax =LCL, xmin = UCL), height = 0.4)+labs(x='高频听力阈值最小二乘几何均值 (dB)',y='口腔健康')+xlim(10,40)+
  theme(
    axis.title.x = element_text(size = 20,color = "black",face="bold"),
    axis.title.y=element_text(size=20,color = "black",face="bold"),
    axis.text.x =element_text(size = 15,colour = 'black',face="bold"), 
    axis.text.y =element_text(size = 15,colour = 'black',face="bold"))+ annotate("text",x=40,y=1:2,label=PHF,size=5.5,alpha=0.9,color="black")+coord_flip()
p_glmlf <- ggplot(LF_GLm, aes(lsmean,ora_health))+ geom_point(size=3.6)+
  geom_errorbarh(aes(xmax =LCL, xmin = UCL), height = 0.4)+labs(x='低频听力阈值最小二乘几何均值 (dB)',y='口腔健康')+xlim(10,40)+
  theme(
    axis.title.x = element_text(size = 20,color = "black",face="bold"),
    axis.title.y=element_text(size=20,color = "black",face="bold"),
    axis.text.x =element_text(size = 15,colour = 'black',face="bold"), 
    axis.text.y =element_text(size = 15,colour = 'black',face="bold"))+ annotate("text",x=25,y=1:2,label=PLF,size=5.5,alpha=0.9,color="black")+coord_flip()

hea_tcs1<- filter(all,RIDAGEYR>=45&URXTRSgr==1)%>% na.omit(.)
hea_tcs1model <- mutate(all,inAnalysis= SEQN%in%hea_tcs1$SEQN)
hea_tcs1M <- svydesign(data=hea_tcs1model, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC6YR, nest=TRUE)
hea_tcs1mod <- subset(hea_tcs1M, inAnalysis)
####
fitLF_tcs1 <- svyglm(log(LF+11)~OHQ845+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                     #,family = "quasibinomial"
                     ,design=hea_tcs1mod)
summary(fitLF_tcs1)
PLF_tcs1 <- c('Reference',signif(summary(fitLF_tcs1)[["coefficients"]][,'Pr(>|t|)'][2],2))
fitHF_tcs1 <- svyglm(log(HF+12)~OHQ845+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                     #,family = "quasibinomial"
                     ,design=hea_tcs1mod)
summary(fitHF_tcs1)

PHF_tcs1 <- c('Reference',signif(summary(fitHF_tcs1)[["coefficients"]][,'Pr(>|t|)'][2],2))
######
HF_GLM <- as.data.frame(lsmeans(fitHF_tcs1,"OHQ845",rg.limit=90000))[,c(1,2,5,6)]
LF_GLm <- as.data.frame(lsmeans(fitLF_tcs1,'OHQ845',rg.limit=90000))[,c(1,2,5,6)]
clu_names <- c("ora_health","lsmean","LCL","UCL")
HF_GLM[,2:4] <- exp(HF_GLM[,2:4])-12
LF_GLm[,2:4] <- exp(LF_GLm[,2:4])-11
names(HF_GLM) <- clu_names
names(LF_GLm) <- clu_names
HF_GLM$ora_health <- c('好','一般/差')
LF_GLm$ora_health <- c('好','一般/差')
p_glmhf_tcs1 <- ggplot(HF_GLM, aes(lsmean,ora_health))+ geom_point(size=3.6)+
  geom_errorbarh(aes(xmax =LCL, xmin = UCL), height = 0.4)+labs(x='高频听力阈值最小二乘几何均值 (dB)',y='口腔健康')+xlim(5,50)+
  theme(
    axis.title.x = element_text(size = 20,color = "black",face="bold"),
    axis.title.y=element_text(size=20,color = "black",face="bold"),
    axis.text.x =element_text(size = 15,colour = 'black',face="bold"), 
    axis.text.y =element_text(size = 15,colour = 'black',face="bold"))+ annotate("text",x=50,y=1:2,label=PHF_tcs1,size=5.5,alpha=0.9,color="black")+coord_flip()
p_glmlf_tcs1 <- ggplot(LF_GLm, aes(lsmean,ora_health))+ geom_point(size=3.6)+
  geom_errorbarh(aes(xmax =LCL, xmin = UCL), height = 0.4)+labs(x='低频听力阈值最小二乘几何均值 (dB)',y='口腔健康')+xlim(5,40)+
  theme(
    axis.title.x = element_text(size = 20,color = "black",face="bold"),
    axis.title.y=element_text(size=20,color = "black",face="bold"),
    axis.text.x =element_text(size = 15,colour = 'black',face="bold"), 
    axis.text.y =element_text(size = 15,colour = 'black',face="bold"))+ annotate("text",x=25,y=1:2,label=PLF_tcs1,size=5.5,alpha=0.9,color="black")+coord_flip()
jpeg(file='口腔与听力及检测限下.jpg',width=15000,height=8500,res=600)
ggdraw() +
  draw_plot(p_glmhf, x = 0, y = .5, width = .35, height = .48) +
  draw_plot(p_glmlf, x = 0, y = .0, width = .35, height = .48) +
  draw_plot(p_glmhf_tcs1, x = .5, y = .5, width = 0.35, height = 0.48) +
  draw_plot(p_glmlf_tcs1, x = .5, y = .0, width = 0.35, height = 0.48) + 
  draw_plot_label(label = c("A","B","C",'D'), size = 22,
                  x = c(0,0.5,0,0.5), y = c(1,1,0.5,0.5))
dev.off()



jpeg(file='第三部分图ֵ1.jpg',width=15000,height=8500,res=600)
ggdraw( xlim = c(0, 1.2)) +
  draw_plot(p_glmhf, x = 0, y = .2, width = .5, height = .6) +
  draw_plot(p_glmlf, x = 0.6, y = .2, width = .5, height = .6) +
  draw_plot_label(label = c("A","B"), size = 28,
                  x = c(0,0.6), y = c(0.9,0.9))

dev.off()


#######
fitHF_ora_tcs <- svyglm(OHQ845~URXTRSgr+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+URXUCR+stroke+HartF+BPQ020+DIQ010
                ,family = "quasibinomial"
                ,design=heamod)
summary(fitHF_ora_tcs)

ora_tcs1 <- creatforest(fitHF_ora_tcs,'三氯生')
ora_tcs1 <- add_row(range='Odds Ratio (95% CI)',.data = ora_tcs1,.before = 1)
is.sum <- is.na(ora_tcs1$upper)
library(forestplot)
jpeg(file='三氯生-口腔健康.jpg',width=15000,height=8500,res=600)
forestplot(labeltext = as.matrix(ora_tcs1[,c(5,4)]),is.summary = is.sum,ora_tcs1$OR,ora_tcs1$lower,ora_tcs1$upper,zero = 1
           ,xticks =c(0,1,2),graph.pos = 3,boxsize = 0.2
           ,lwd.xaxis = 2,lwd.zero = 2,lwd.ci = 2,title="",graphwidth = unit(.4,"npc"),
           col = fpColors(all.elements = 'black'),txt_gp=fpTxtGp(label=gpar(cex=1.5),ticks=gpar(cex=1.5),xlab=gpar(cex = 1.5), title=gpar(cex = 2.0)))
dev.off()


#######
hea_ora1<- filter(all,RIDAGEYR>=45&OHQ845==0)%>% na.omit(.)
hea_ora1model <- mutate(all,inAnalysis= SEQN%in%hea_ora1$SEQN)
hea_ora1M <- svydesign(data=hea_ora1model, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC6YR, nest=TRUE)
hea_ora1mod <- subset(hea_ora1M, inAnalysis)
####
fitLF_ora1 <- svyglm(log(LF+11)~URXTRSgr+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+URXUCR+stroke+HartF+BPQ020+DIQ010
                #,family = "quasibinomial"
                ,design=hea_ora1mod)
summary(fitLF_ora1)
PLF_ora1 <- c('Reference',signif(summary(fitLF_ora1)[["coefficients"]][,'Pr(>|t|)'][2:4],2))
fitHF_ora1 <- svyglm(log(HF+12)~URXTRSgr+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+URXUCR+stroke+HartF+BPQ020+DIQ010
                #,family = "quasibinomial"
                ,design=hea_ora1mod)
summary(fitHF_ora1)
PHF_ora1 <- c('Reference',signif(summary(fitHF_ora1)[["coefficients"]][,'Pr(>|t|)'][2:4],2))
#####
#######
#######
hea_ora2<- filter(all,RIDAGEYR>=45&OHQ845==1)%>% na.omit(.)
hea_ora2model <- mutate(all,inAnalysis= SEQN%in%hea_ora2$SEQN)
hea_ora2M <- svydesign(data=hea_ora2model, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC6YR, nest=TRUE)
hea_ora2mod <- subset(hea_ora2M, inAnalysis)
####
fitLF_ora2 <- svyglm(log(LF+11)~URXTRSgr+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+URXUCR+stroke+HartF+BPQ020+DIQ010
                #,family = "quasibinomial"
                ,design=hea_ora2mod)
summary(fitLF_ora2)
PLF_ora2 <- c('Reference',signif(summary(fitLF_ora2)[["coefficients"]][,'Pr(>|t|)'][2:4],2))
fitHF_ora2 <- svyglm(log(HF+12)~URXTRSgr+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+URXUCR+stroke+HartF+BPQ020+DIQ010
                #,family = "quasibinomial"
                ,design=hea_ora2mod)
summary(fitHF_ora2)
PHF_ora2 <- c('Reference',signif(summary(fitHF_ora2)[["coefficients"]][,'Pr(>|t|)'][2:4],2))
######
HF_log_gr_model1 <- creatforest_gasso(fitHF_ora1,'口腔健康良好')
LF_log_gr_model1 <- creatforest_gasso(fitHF_ora2,'口腔健康一般/良好')
log_gr_model1 <- rbind(HF_log_gr_model1,LF_log_gr_model1)
is.sum <- is.na(log_gr_model1$upper)
jpeg(file='口腔分层.jpg',width=15000,height=8500,res=600)
forestplot(labeltext = as.matrix(log_gr_model1[,c(5,4)]),is.summary = is.sum,log_gr_model1$OR,log_gr_model1$lower,log_gr_model1$upper,zero = 0
           ,xticks =c(-0.25,0,0.25),graph.pos = 3,boxsize = 0.2
           ,lwd.xaxis = 2,lwd.zero = 2,lwd.ci = 2,title="",graphwidth = unit(.4,"npc"),
           col = fpColors(all.elements = 'black'),txt_gp=fpTxtGp(label=gpar(cex=1.5),ticks=gpar(cex=1.5),xlab=gpar(cex = 1.5), title=gpar(cex = 2.0)))
dev.off()
######
hea_tcs1<- filter(all,RIDAGEYR>=45&URXTRSgr==1)%>% na.omit(.)
hea_tcs1model <- mutate(all,inAnalysis= SEQN%in%hea_tcs1$SEQN)
hea_tcs1M <- svydesign(data=hea_tcs1model, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC6YR, nest=TRUE)
hea_tcs1mod <- subset(hea_tcs1M, inAnalysis)
####
fitLF_tcs1 <- svyglm(log(LF+11)~OHQ845+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                     #,family = "quasibinomial"
                     ,design=hea_tcs1mod)
summary(fitLF_tcs1)
PLF_tcs1 <- c('Reference',signif(summary(fitLF_tcs1)[["coefficients"]][,'Pr(>|t|)'][2],2))
fitHF_tcs1 <- svyglm(log(HF+12)~OHQ845+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+stroke+HartF+BPQ020+DIQ010
                     #,family = "quasibinomial"
                     ,design=hea_tcs1mod)
summary(fitHF_tcs1)

PHF_tcs1 <- c('Reference',signif(summary(fitHF_tcs1)[["coefficients"]][,'Pr(>|t|)'][2],2))
######
HF_GLM <- as.data.frame(lsmeans(fitHF_tcs1,"OHQ845",rg.limit=90000))[,c(1,2,5,6)]
LF_GLm <- as.data.frame(lsmeans(fitLF_tcs1,'OHQ845',rg.limit=90000))[,c(1,2,5,6)]
clu_names <- c("ora_health","lsmean","LCL","UCL")
HF_GLM[,2:4] <- exp(HF_GLM[,2:4])-12
LF_GLm[,2:4] <- exp(LF_GLm[,2:4])-11
names(HF_GLM) <- clu_names
names(LF_GLm) <- clu_names
HF_GLM$ora_health <- c('好','一般/差')
LF_GLm$ora_health <- c('好','一般/差')
p_glmhf_tcs1 <- ggplot(HF_GLM, aes(lsmean,ora_health))+ geom_point(size=3.6)+
  geom_errorbarh(aes(xmax =LCL, xmin = UCL), height = 0.4)+labs(x='高频听力阈值最小二乘几何均值 (dB)',y='口腔健康')+xlim(5,50)+
theme(
  axis.title.x = element_text(size = 20,color = "black",face="bold"),
  axis.title.y=element_text(size=20,color = "black",face="bold"),
  axis.text.x =element_text(size = 15,colour = 'black',face="bold"), 
  axis.text.y =element_text(size = 15,colour = 'black',face="bold"))+ annotate("text",x=50,y=1:2,label=PHF_tcs1,size=5.5,alpha=0.9,color="black")+coord_flip()
p_glmlf_tcs1 <- ggplot(LF_GLm, aes(lsmean,ora_health))+ geom_point(size=3.6)+
  geom_errorbarh(aes(xmax =LCL, xmin = UCL), height = 0.4)+labs(x='低频听力阈值最小二乘几何均值 (dB)',y='口腔健康')+xlim(5,40)+
  theme(
    axis.title.x = element_text(size = 20,color = "black",face="bold"),
    axis.title.y=element_text(size=20,color = "black",face="bold"),
    axis.text.x =element_text(size = 15,colour = 'black',face="bold"), 
    axis.text.y =element_text(size = 15,colour = 'black',face="bold"))+ annotate("text",x=25,y=1:2,label=PLF_tcs1,size=5.5,alpha=0.9,color="black")+coord_flip()
jpeg(file='三氯生检测限下.jpg',width=15000,height=8500,res=600)
ggdraw( xlim = c(0, 1.2)) +
  draw_plot(p_glmhf, x = 0, y = .2, width = .5, height = .6) +
  draw_plot(p_glmlf, x = 0.6, y = .2, width = .5, height = .6) +
  draw_plot_label(label = c("A","B"), size = 28,
                  x = c(0,0.6), y = c(0.9,0.9))

dev.off()

#########
#########
library(mediation)
all$OHQ845 <- as.numeric(all$OHQ845)
all <- mutate(all,OHQ845=ifelse(OHQ845==1,0,1))
hea<- filter(all,RIDAGEYR>=45)%>% na.omit(.)
heamodel <- mutate(all,inAnalysis= SEQN%in%hea$SEQN)
heaM <- svydesign(data=heamodel, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC6YR, nest=TRUE)
heamod <- subset(heaM, inAnalysis)



a <- svyglm(OHQ845~URXTRSgr2+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+URXUCR+stroke+HartF+BPQ020+DIQ010
       ,family = "binomial"
       ,design =heamod)
summary(a)
b <- svyglm(log(HF+12)~URXTRSgr2+OHQ845+Gender+RIDAGEYR+education+PIR+race+medicine+BMXBMI+smoke+URXUCR+stroke+HartF+BPQ020+DIQ010
            #,family = "quasibinomial"
            ,design  = heamod)
set.seed(123)
result = mediate(a,b,treat="URXTRSgr2",mediator = "OHQ845",robustSE = TRUE, sims = 1000)
summary(result)
