library(xlsx)
library(pROC)

f1 = read.csv('E:/work/ELISA/p36_0824統計/input/TableS3_Stage_forR.csv')



Names = unique(f1$分群名)
Groups = unique(f1$分群)
dict_name = vector(mode="list", length=length(Names))
names(dict_name) = Groups
for(i in 1:length(Groups)){
  dict_name[[i]] = Names [i]
}
dict_name
names(dict_name)

#20200206-OPMD+OSCCstage4
df = list()
gpAvsB = c()
gpA_v = c()
gpB_v = c()
AUC_v = c()
se_v = c()
p_v = c()
lower_v =c()
upper_v =c()

Neg_names = c(1)
conc_name = 'ELISA'
Pos_names = c(2,3,4,5)

for(n in Neg_names){
  for(pn in Pos_names){
    neg = f1[f1$分群==n,][,conc_name]
    neg <- neg[!is.na(neg)]
    pos = f1[f1$分群==pn,][,conc_name]
    pos <- pos[!is.na(pos)]
    y_true = c(rep(0,length(neg)), rep(1,length(pos)))
    y_pred = c(neg,pos)

    roc = roc(y_true, y_pred)
    ci = ci.auc(roc)
    lower = ci[1]
    auc = ci[2]
    upper = ci[3]
    se = (upper-auc)/qnorm(0.025,lower.tail = FALSE)
    z = (auc - 0.5)/se
    p = pnorm(z,  lower.tail = FALSE)
    gpA_number = length(neg)
    gpB_number = length(pos)
    gpA_median = median(neg)
    gpB_median = median(pos)
    gpA_range = range(neg)
    gpB_range = range(pos)

    gpAvsB = c(gpAvsB,paste(dict_name[toString(n)][[1]]," vs ",dict_name[toString(pn)][[1]]," (",gpA_number," vs ",gpB_number,")",sep=""))
    gpA_v = c(gpA_v,paste(gpA_median," (",gpA_range[1],"-",gpA_range[2],")",sep=""))
    gpB_v = c(gpB_v,paste(gpB_median," (",gpB_range[1],"-",gpB_range[2],")",sep=""))
    AUC_v = c(AUC_v,auc)
    se_v = c(se_v,se)
    p_v = c(p_v,p)
    lower_v =c(lower_v,lower)
    upper_v =c(upper_v,upper)
  }
}
df[['`group A vs B (case number)`']] =  gpAvsB
df[['`Group A Median(range)`']] =  gpA_v
df[['`Group B Median(range)`']] =  gpB_v
df[['AUC']] =  AUC_v
df[['標準誤差']] =  se_v
df[['漸近顯著性']] =  p_v
df[['Lower bound of AUC']] =  lower_v
df[['Upper bound of AUC']] =  upper_v

df = as.data.frame(df)
write.xlsx(df, file='E:/work/ELISA/p36_0824統計/output/Table2-Stage_AUC.xlsx', col.names = TRUE, row.names = FALSE, append = FALSE)
