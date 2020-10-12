
library(ggplot2)
library(latex2exp)
library(RColorBrewer)
library(ggpubr)
setwd("~/PosDoc/Coronavirus/LTCF/Code_september/cluster_data/")



iso = "total"

build  = "new"
beta = "06"

type_ind = "res"

################ testing AR
folder = paste("results_",build,"_",iso,"_S0/",sep = "")
data = read.table(paste(folder,"beta_0_",beta,"_sev_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
sum(data)/8000/120


type_data = "latent"
type_ind = "res"

n_boot = 1000

scenarios = c("S0","S1","S2a","S2b","S3a","S3b","S4a","S4b","S5a","S5b","S6a","S6b")

tempo = c()
result = c()
cenario = c()
for(j in 1:length(scenarios)){
  print(j)
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = rowMeans(data)
  t = seq(1,length(m))/24/5
  t = ceiling(t)
  
  media_d = rep(0,max(t))
  for(i in 1:length(t)){
    media_d[t[i]] = media_d[t[i]]+m[i]
  }
  tempo  = c(tempo,seq(1,max(t)))
  result = c(result,media_d)
  cenario = c(cenario,rep(scen,length(media_d)))
  
}

df = data.frame(tempo=tempo*5,media = result*1000/120,cenario)
write.table(df,paste("time_series_",build,"_",type_ind,".dat",sep = ""),row.names = F)
#rm(matrix_r)


type_data = "latent"
type_ind = "hcw"

n_boot = 1000

scenarios = c("S0","S1","S2a","S2b","S3a","S3b","S4a","S4b","S5a","S5b","S6a","S6b")

tempo = c()
result = c()
cenario = c()
for(j in 1:length(scenarios)){
  print(j)
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = rowMeans(data)
  t = seq(1,length(m))/24/5
  t = ceiling(t)
  
  media_d = rep(0,max(t))
  for(i in 1:length(t)){
    media_d[t[i]] = media_d[t[i]]+m[i]
  }
  tempo  = c(tempo,seq(1,max(t)))
  result = c(result,media_d)
  cenario = c(cenario,rep(scen,length(media_d)))
  
}

df = data.frame(tempo=tempo*5,media = result*1000/120,cenario)
write.table(df,paste("time_series_",build,"_",type_ind,".dat",sep = ""),row.names = F)
#rm(matrix_r)

#####################################################################################################################3
#######################################################################################################################
####################### Plots total number ###########################################################################
####################################################################################################################



library(ggplot2)
library(latex2exp)
library(RColorBrewer)
library(ggpubr)


scenarios = c("S0","S1","S2a","S2b","S3a","S3b","S4a","S4b","S5a","S5b","S6a","S6b")
n_boot = 1000


type_ind = "res"

type_data = "latent"
result = c()
cenario = c()

lab_ = "a"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/120,cenario)
write.table(df,paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)

type_data = "hosp"
result = c()
cenario = c()


lab_ = "b"
for(j in 1:length(scenarios)){
  print(j)
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/120,cenario)
write.table(df,paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)



type_data = "dead"
result = c()
cenario = c()

lab_ = "c"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/120,cenario)
write.table(df,paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)


type_ind = "hcw"

type_data = "latent"
result = c()
cenario = c()

lab_ = "d"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/68,cenario)
write.table(df,paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)


type_data = "hosp"
result = c()
cenario = c()


lab_ = "e"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/68,cenario)
write.table(df,paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)


type_data = "dead"
result = c()
cenario = c()

lab_ = "f"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/68,cenario)
write.table(df,paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)



###########################################################################################################################
#################################################### Plot percentage of silent infections
#######################################################################################################


library(ggplot2)
library(latex2exp)
library(RColorBrewer)
library(ggpubr)

setwd("/data/thomas-covid/LTCF/")


scenarios = c("S0","S1","S2a","S2b","S3a","S3b","S4a","S4b","S5a","S5b","S6a","S6b")
n_boot = 1000


type_ind = "res"

type_data = "latent"


n_boot = 1000

scenarios = c("S2a","S2b","S3a","S3b","S4a","S4b","S5a","S5b","S6a","S6b")

tempo = c()
result = c()
cenario = c()

for(j in 1:length(scenarios)){
  rs = c()
  print(j)
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_pre_res.dat",sep=""),h=F,stringsAsFactors = F)
  sr = colSums(data)+1
  found_r = read.table(paste(folder,"beta_0_",beta,"_iso_pre_res.dat",sep=""),h=F,stringsAsFactors = F)$V1
  
  data = read.table(paste(folder,"beta_0_",beta,"_pre_hcw.dat",sep=""),h=F,stringsAsFactors = F)
  sh = colSums(data)+1
  found_h = read.table(paste(folder,"beta_0_",beta,"_iso_pre.dat",sep=""),h=F,stringsAsFactors = F)$V1
  
  sh = sh+sr
  found_h = found_h+found_r
  
  data = read.table(paste(folder,"beta_0_",beta,"_asymp_res.dat",sep=""),h=F,stringsAsFactors = F)
  sr = colSums(data)+1
  found_r = read.table(paste(folder,"beta_0_",beta,"_iso_asymp_res.dat",sep=""),h=F,stringsAsFactors = F)$V1
  
  sh = sh+sr
  found_h = found_h+found_r
  
  data = read.table(paste(folder,"beta_0_",beta,"_asymp_hcw.dat",sep=""),h=F,stringsAsFactors = F)
  sr = colSums(data)+1
  found_r = read.table(paste(folder,"beta_0_",beta,"_iso_asymp.dat",sep=""),h=F,stringsAsFactors = F)$V1
  
  sh = sh+sr
  found_h = found_h+found_r
  
  for(i in 1:n_boot){
    pos = sample(seq(1,length(sr)),length(sr),replace = T)
    st = sum(sh[pos])
    if(st>0){
      rs[i] = sum(found_h[pos])/st
    }else
      rs[i] = NA
  }
  result = c(result,rs)
  cenario = c(cenario,rep(scen,length(rs)))
  
}


df = data.frame(proportion = result,scenario=cenario)
write.table(df,paste("proportion_",build,".dat",sep=""),row.names = F)
######################################################################################################33
######################33333 Calculating Reproductive Number



library(ggplot2)
library(latex2exp)
library(RColorBrewer)
library(ggpubr)

scenarios = c("S0","S1","S2a","S2b","S3a","S3b","S4a","S4b","S5a","S5b","S6a","S6b")
n_boot = 1000

R0_m = c()
R0_sd = c()

for(i in 1:length(scenarios)){
  scen = scenarios[i]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_R0.dat",sep=""),h=F,stringsAsFactors = F)$V1
  R0_m[i] = mean(data)
  R0_sd[i] = sd(data)
}

df = data.frame(scenarios,R0_m,R0_sd)
write.csv(df,paste("reproduction_number_",build,".csv",sep=""),row.names=F)
plot(df$R0_m)


######################################################################################################
######################################### Average reduction
#########################################################################################################


library(ggplot2)
library(latex2exp)
library(RColorBrewer)
library(ggpubr)

scenarios = c("S1","S2a","S2b","S3a","S3b","S4a","S4b","S5a","S5b","S6a","S6b")
n_boot = 1000

type_data = "latent"
type_ind = "res"
folder = paste("results_",build,"_",iso,"_S0/",sep = "")
data = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T,stringsAsFactors = F)

d0 = data$result[data$cenario=="S0"]

l_res = c()

for(i in 1:length(scenarios)){
  scen = scenarios[i]
  d = data$result[data$cenario==scen]
  d = -(d-d0)/d0
  m =median(d[!is.infinite(d)])
  s1 = quantile(d[!is.infinite(d)],0.975,name=F)
  s2 = quantile(d[!is.infinite(d)],0.025,name=F)
  l_res[i] = paste(round(m,digits=3)," ","(",round(s2,digits=3),",",round(s1,digits=3),")",sep="")
}

type_data = "hosp"
folder = paste("results_",build,"_",iso,"_S0/",sep = "")
data = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T,stringsAsFactors = F)

d0 = data$result[data$cenario=="S0"]

h_res = c()

for(i in 1:length(scenarios)){
  scen = scenarios[i]
  d = data$result[data$cenario==scen]
  d = -(d-d0)/d0
  m =median(d[!is.infinite(d)])
  s1 = quantile(d[!is.infinite(d)],0.975,name=F)
  s2 = quantile(d[!is.infinite(d)],0.025,name=F)
  h_res[i] = paste(round(m,digits=3)," ","(",round(s2,digits=3),",",round(s1,digits=3),")",sep="")
}

type_data = "dead"
folder = paste("results_",build,"_",iso,"_S0/",sep = "")
data = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T,stringsAsFactors = F)

d0 = data$result[data$cenario=="S0"]

d_res = c()

for(i in 1:length(scenarios)){
  scen = scenarios[i]
  d = data$result[data$cenario==scen]
  d = -(d-d0)/d0
  m =mean(d[!is.infinite(d)])
  m =median(d[!is.infinite(d)])
  s1 = quantile(d[!is.infinite(d)],0.975,name=F)
  s2 = quantile(d[!is.infinite(d)],0.025,name=F)
  d_res[i] = paste(round(m,digits=3)," ","(",round(s2,digits=3),",",round(s1,digits=3),")",sep="")
}

type_ind="hcw"

type_data = "latent"
folder = paste("results_",build,"_",iso,"_S0/",sep = "")
data = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T,stringsAsFactors = F)

d0 = data$result[data$cenario=="S0"]

l_hcw = c()

for(i in 1:length(scenarios)){
  scen = scenarios[i]
  d = data$result[data$cenario==scen]
  d = -(d-d0)/d0
  m =median(d[!is.infinite(d)])
  s1 = quantile(d[!is.infinite(d)],0.975,name=F)
  s2 = quantile(d[!is.infinite(d)],0.025,name=F)
  l_hcw[i] = paste(round(m,digits=3)," ","(",round(s2,digits=3),",",round(s1,digits=3),")",sep="")
}

type_data = "hosp"
folder = paste("results_",build,"_",iso,"_S0/",sep = "")
data = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T,stringsAsFactors = F)

d0 = data$result[data$cenario=="S0"]

h_hcw = c()

for(i in 1:length(scenarios)){
  scen = scenarios[i]
  d = data$result[data$cenario==scen]
  d = -(d-d0)/d0
  m =median(d[!is.infinite(d)])
  s1 = quantile(d[!is.infinite(d)],0.975,name=F)
  s2 = quantile(d[!is.infinite(d)],0.025,name=F)
  h_hcw[i] = paste(round(m,digits=3)," ","(",round(s2,digits=3),",",round(s1,digits=3),")",sep="")
}

type_data = "dead"
folder = paste("results_",build,"_",iso,"_S0/",sep = "")
data = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T,stringsAsFactors = F)

d0 = data$result[data$cenario=="S0"]

d_hcw = c()

for(i in 1:length(scenarios)){
  scen = scenarios[i]
  d = data$result[data$cenario==scen]
  d = -(d-d0)/d0
  m =median(d[!is.infinite(d)])
  s1 = quantile(d[!is.infinite(d)],0.975,name=F)
  s2 = quantile(d[!is.infinite(d)],0.025,name=F)
  d_hcw[i] = paste(round(m,digits=3)," ","(",round(s2,digits=3),",",round(s1,digits=3),")",sep="")
}




df = data.frame(scenarios,l_res,h_res,d_res,l_hcw,h_hcw,d_hcw)
write.csv(df,paste("reduction_",build,".csv",sep=""),row.names=F)



######################################################################################
#### SA


scenarios = c("SA60","SA70","SA80","SA90")
n_boot = 1000


type_ind = "res"

type_data = "latent"
result = c()
cenario = c()

lab_ = "a"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/120,cenario)
write.table(df,paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)

type_data = "hosp"
result = c()
cenario = c()


lab_ = "b"
for(j in 1:length(scenarios)){
  print(j)
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/120,cenario)
write.table(df,paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)



type_data = "dead"
result = c()
cenario = c()

lab_ = "c"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/120,cenario)
write.table(df,paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)


type_ind = "hcw"

type_data = "latent"
result = c()
cenario = c()

lab_ = "d"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/68,cenario)
write.table(df,paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)


type_data = "hosp"
result = c()
cenario = c()


lab_ = "e"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/68,cenario)
write.table(df,paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)


type_data = "dead"
result = c()
cenario = c()

lab_ = "f"
for(j in 1:length(scenarios)){
  scen = scenarios[j]
  folder = paste("results_",build,"_",iso,"_",scen,"/",sep = "")
  data = read.table(paste(folder,"beta_0_",beta,"_",type_data,"_",type_ind,".dat",sep=""),h=F,stringsAsFactors = F)
  m = colSums(data)
  media = c()
  for(i in 1:n_boot){
    media[i] = mean(sample(m,1000,replace = T))
  }
  result = c(result,media)
  cenario = c(cenario,rep(scen,length(media)))
}

df = data.frame(result=result*1000/68,cenario)
write.table(df,paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),row.names = F)





##################################################################
### Plot time series

build = "old"

df = read.table(paste("time_series_",build,"_res.dat",sep=""),h=T,stringsAsFactors = F)

baseline = "S0"
df1 = df[df$cenario == baseline,]

curves = c("S1","S2a","S2b","S3a","S3b")

cv = c(brewer.pal(5,"Set1"),brewer.pal(7,"Set1")[7])

df2 = df[df$cenario %in% curves,]
lab_ = "a"

mp1=ggplot()+geom_col(data=df1,aes(x=tempo,y=media),color="lightgrey",fill="lightgrey",alpha = 0.6,size = 0)+
  geom_line(data=df2,aes(x=tempo,y=media,color = as.factor(cenario)),alpha = 0.9,size = 2)+
  scale_x_continuous(name = TeX('Time (days)'),limits = c(0,150))+
  scale_y_continuous(name = TeX("Incidence of infections"))+
  scale_color_manual(values = cv[1:length(curves)],breaks = curves,labels = curves,name = "Scenarios")+
  guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          # panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "top",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp1

curves = c("S4a","S4b","S5a","S5b","S6a","S6b")

df2 = df[df$cenario %in% curves,]
lab_ = "b"

mp2=ggplot()+geom_col(data=df1,aes(x=tempo,y=media),color="lightgrey",fill="lightgrey",alpha = 0.6,size = 0)+
  geom_line(data=df2,aes(x=tempo,y=media,color = as.factor(cenario)),alpha = 0.9,size = 2)+
  scale_x_continuous(name = TeX('Time (days)'),limits = c(0,150))+
  scale_y_continuous(name = TeX("Incidence of infections"))+
  scale_color_manual(values = cv[1:length(curves)],breaks = curves,labels = curves,name = "Scenarios")+
  guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          # panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "top",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp2


paste("time_series_",build,"_hcw.dat",sep="")


baseline = "S0"
df1 = df[df$cenario == baseline,]

curves = c("S1","S2a","S2b","S3a","S3b")


df2 = df[df$cenario %in% curves,]
lab_ = "c"

mp3=ggplot()+geom_col(data=df1,aes(x=tempo,y=media),color="lightgrey",fill="lightgrey",alpha = 0.6,size = 0)+
  geom_line(data=df2,aes(x=tempo,y=media,color = as.factor(cenario)),alpha = 0.9,size = 2)+
  scale_x_continuous(name = TeX('Time (days)'),limits = c(0,150))+
  scale_y_continuous(name = TeX("Incidence of infections"))+
  scale_color_manual(values =cv[1:length(curves)],breaks = curves,labels = curves,name = "Scenarios")+
  guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          # panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "top",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp3

curves = c("S4a","S4b","S5a","S5b","S6a","S6b")

df2 = df[df$cenario %in% curves,]
lab_ = "d"

mp4=ggplot()+geom_col(data=df1,aes(x=tempo,y=media),color="lightgrey",fill="lightgrey",alpha = 0.6,size = 0)+
  geom_line(data=df2,aes(x=tempo,y=media,color = as.factor(cenario)),alpha = 0.9,size = 2)+
  scale_x_continuous(name = TeX('Time (days)'),limits = c(0,150))+
  scale_y_continuous(name = TeX("Incidence of infections"))+
  scale_color_manual(values = cv[1:length(curves)],breaks = curves,labels = curves,name = "Scenarios")+
  guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          # panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "top",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp4


legend <- cowplot::get_legend(mp1+guides(colour = guide_legend(override.aes = list(size=10),nrow=2,ncol = 3))+ theme(legend.position = "bottom",legend.text = element_text(size = 40),
                                                                                                                     legend.title = element_text(size = 40)))

legend2 <- cowplot::get_legend(mp2+guides(nrow=2,colour = guide_legend(override.aes = list(size=10)))+ theme(legend.position = "bottom",legend.text = element_text(size = 40),
                                                                                                             legend.title = element_text(size = 40)))

mp = ggarrange(mp1+theme(legend.position = "none"),NULL,mp2+theme(legend.position = "none")+rremove("y.title")+rremove("x.title"),mp3+theme(legend.position = "none"),NULL,mp4+rremove("y.title")+theme(legend.position = "none"),legend,NULL,legend2,ncol = 3,nrow = 3,widths = c(1,0.1,1,1,0.1,1),heights = c(1,1,0.5))

ggsave(paste("time_series_",build,".png",sep=""),mp,device="png",width = 50,height = 50,units = "cm",dpi=300)


########################################################3
########### here the plots for comparison #######


type_ind = "res"

type_data = "latent"
result = c()
cenario = c()

lab_ = "a"

df = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
mp1=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "red",fill = "red",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "red",fill="red",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of infections'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"))+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp1



type_data = "hosp"
result = c()
cenario = c()


lab_ = "b"

df = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
mp2=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "blue",fill = "blue",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "blue",fill="blue",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of hospitalizations'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"))+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp2

type_data = "dead"
result = c()
cenario = c()

lab_ = "c"

df = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
mp3=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "black",fill = "black",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "black",fill="black",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of deaths'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"))+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp3

type_ind = "hcw"

type_data = "latent"
result = c()
cenario = c()

lab_ = "d"

df = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
mp4=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "red",fill = "red",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "red",fill="red",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of infections'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"))+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp4



type_data = "hosp"
result = c()
cenario = c()


lab_ = "e"

df = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
mp5=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "blue",fill = "blue",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "blue",fill="blue",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of hospitalizations'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"))+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp5




type_data = "dead"
result = c()
cenario = c()

lab_ = "f"

df = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
mp6=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "black",fill = "black",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "black",fill="black",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of deaths'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"))+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp6


#mp=ggarrange(mp1,NULL,mp2+rremove("y.text")+rremove("y.title"),NULL,mp3+rremove("y.text")+rremove("y.title"), ncol = 5, nrow = 1,widths = c(1,0.2,1,0.2,1),vjust = 0.5)
#ggsave(paste("boxplot_",build,"_",type_ind,".png",sep=""),mp,device="png",width = 70,height = 25,units = "cm",dpi=300)
ts = theme(axis.text.x = element_text(colour="black", size = 35),
           axis.text.y = element_text(colour="black", size = 35),
           axis.title.y = element_text(colour="black", size = 40),
           axis.title.x = element_text(colour="black", size = 40))
mp=ggarrange(mp1+ts+rremove("x.title"),NULL,mp2+ts+rremove("y.text")+rremove("y.title")+rremove("x.title"),NULL,mp3+ts+rremove("y.text")+rremove("y.title")+rremove("x.title"),mp4+ts,NULL,mp5+ts+rremove("y.text")+rremove("y.title"),NULL,mp6+ts+rremove("y.text")+rremove("y.title"), ncol = 5, nrow = 2,widths = c(1,0.2,1,0.2,1,1,0.2,1,0.2,1),vjust = 0.5)
ggsave(paste("boxplot_",build,".png",sep=""),mp,device="png",width = 70,height = 60,units = "cm",dpi=300)


################################################ Sensitivity analysis




type_ind = "res"

type_data = "latent"
result = c()
cenario = c("S2b","SA60","SA70","SA80","SA90")
label_y = c("NP","60% sensitivity","70% sensitivity","80% sensitivity","90% sensitivity")
lab_ = "a"

%label_y = c(labels_y$NP,labels_y$`60\\% sensitivity`,labels_y$`70\\% sensitivity`,labels_y$`80\\% sensitivity`,labels_y$`90\\% sensitivity`)

df0 = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df0 = data.frame(df0[df0$cenario=="S2b",])

df = read.table(paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
df = rbind(df0,df)
mp1=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "red",fill = "red",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "red",fill="red",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of infections'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"),breaks = cenario,labels = label_y)+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp1



type_data = "hosp"


lab_ = "b"

df0 = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df0 = data.frame(df0[df0$cenario=="S2b",])

df = read.table(paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
df = rbind(df0,df)
mp2=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "blue",fill = "blue",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "blue",fill="blue",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of hospitalizations'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"),breaks = cenario,labels = label_y)+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp2

type_data = "dead"

lab_ = "c"

df0 = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df0 = data.frame(df0[df0$cenario=="S2b",])

df = read.table(paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
df = rbind(df0,df)
mp3=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "black",fill = "black",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "black",fill="black",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of deaths'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"),breaks = cenario,labels = label_y)+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp3

type_ind = "hcw"

type_data = "latent"

lab_ = "d"

df0 = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df0 = data.frame(df0[df0$cenario=="S2b",])

df = read.table(paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
df = rbind(df0,df)
mp4=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "red",fill = "red",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "red",fill="red",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of infections'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"),breaks = cenario,labels = label_y)+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp4



type_data = "hosp"


lab_ = "e"

df0 = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df0 = data.frame(df0[df0$cenario=="S2b",])

df = read.table(paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
df = rbind(df0,df)
mp5=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "blue",fill = "blue",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "blue",fill="blue",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of hospitalizations'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"),breaks = cenario,labels = label_y)+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp5




type_data = "dead"

lab_ = "f"

df0 = read.table(paste("average_bootstrap_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df0 = data.frame(df0[df0$cenario=="S2b",])

df = read.table(paste("average_bootstrap_SA_",type_ind,"_",build,"-",type_data,".dat",sep=""),h=T)
df = data.frame(df)
df = rbind(df0,df)
mp6=ggplot()+
  stat_summary(data=df,mapping = aes(x = result, y = cenario),fun = median,geom = "col",alpha = 0.4,color = "black",fill = "black",size = 0)+
  geom_boxplot(data=df,aes(x=result,y=cenario),color = "black",fill="black",outlier.shape = NA,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Number of deaths'),expand = c(0,0))+
  scale_y_discrete(name = TeX("Scenarios"),breaks = cenario,labels = label_y)+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )+labs(tag = paste(lab_,")",sep=""))
mp6


#mp=ggarrange(mp1,NULL,mp2+rremove("y.text")+rremove("y.title"),NULL,mp3+rremove("y.text")+rremove("y.title"), ncol = 5, nrow = 1,widths = c(1,0.2,1,0.2,1),vjust = 0.5)
#ggsave(paste("boxplot_",build,"_",type_ind,".png",sep=""),mp,device="png",width = 70,height = 25,units = "cm",dpi=300)
ts = theme(axis.text.x = element_text(colour="black", size = 35),
           axis.text.y = element_text(colour="black", size = 35),
           axis.title.y = element_text(colour="black", size = 40),
           axis.title.x = element_text(colour="black", size = 40))
mp=ggarrange(mp1+ts+rremove("x.title"),NULL,mp2+ts+rremove("y.text")+rremove("y.title")+rremove("x.title"),NULL,mp3+ts+rremove("y.text")+rremove("y.title")+rremove("x.title"),NULL,mp4+ts,NULL,mp5+ts+rremove("y.text")+rremove("y.title"),NULL,mp6+ts+rremove("y.text")+rremove("y.title"),NULL, ncol = 6, nrow = 2,widths = c(1,0.1,0.7,0.1,0.7,0.1,1,0.1,0.7,0.1,0.7,0.1),vjust = 0.5)
ggsave(paste("boxplot_SA_",build,".png",sep=""),mp,device="png",width = 80,height = 60,units = "cm",dpi=300)


###################################################33proportion of silent infections that are detected.


scenarios = c("S2a","S2b","S3a","S3b","S4a","S4b","S5a","S5b","S6a","S6b")


df = read.table(paste("proportion","_",build,".dat",sep=""),h=T)
df = data.frame(df)

mp1=ggplot()+
  stat_summary(data=df,mapping = aes(x = proportion, y = scenario),fun = median,geom = "col",alpha = 0.4,color = "red",fill = "red",size = 0)+
  geom_boxplot(data=df,aes(x=proportion,y=scenario),color = "red",fill="red",outlier.shape = 1,alpha = 0.7,size = 2)+
  scale_x_continuous(name = TeX('Proportion of detected silent infections'),expand = c(0,0),limits = c(0,0.19))+
  scale_y_discrete(name = TeX("Scenarios"),breaks = scenarios)+
  #scale_color_manual(values = brewer.pal(12,"Paired"),breaks = scenarios,labels = scenarios,name = "Scenarios")+
  #guides(color = guide_legend(reverse = FALSE, override.aes = list(size=5), byrow=TRUE,nrow = 2))+
  #scale_linetype_manual(values = c(1,4,10),breaks = c("not_testing","test_7","test_14"))+
  # scale_fill_manual(values = c("red","blue"),name="Individual",breaks = c("res","hcw"),labels = c("Residents","HCW"))+
  theme_bw()+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(colour = "black"),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25),
          axis.title.y = element_text(colour="black", size = 35),
          axis.title.x = element_text(colour="black", size = 35),
          plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.position = "none",
          plot.tag.position = "topleft",
          #axis.line = element_line(size=0.5, colour = "black")
  )#+labs(tag = paste(lab_,")",sep=""))
mp1
ggsave(paste("boxplot_proportion_",build,".png",sep=""),mp1,device="png",width = 30,height = 30,units = "cm",dpi=300)


