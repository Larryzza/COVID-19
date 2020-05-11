#setwd()
##### Tencent
data<-read.csv("cityflu_data_tencent.csv")
##### baidu
data<-read.csv("cityflu_data_baidu.csv")


colnames(data)[1]<-"place"
#data<-data[which(data$flu>10),]


#sum(data$flu)*(rate)
#data$flu<-round(data$flu*rate)

ff<-function(x)-log(dbinom(num_oversea,x,p))-d0-0.5*qchisq(0.95,1)
for(a in c(1:length(data$flu))){
  #a<-1
  num_oversea<-data[a,2]
  mean_time_detection<-10
  catchment_pop<-19e6
  daily_outbound<-data[a,3]
  daily_prob<-daily_outbound/catchment_pop
  p<-daily_prob*mean_time_detection 
  total<-round(num_oversea/p)
  d0<- -log(dbinom(num_oversea,total,p))
  temp<-NULL
  for(i in 1:150000){
    if(ff(i)*ff(i+1)<0){
      temp<-c(temp,i)
    }
  }
  data$forec[a]<-total
  data$P[a]<-p
}

####loglik
resul<-NULL
for(a in 1:20000){
  logli.sum<-NULL
  m<-mean_time_detection/catchment_pop
  for(i in c(1:length(data$flu))){
    ####estimate total cases here
    logli<-log(dbinom(data[i,2],a,(data[i,3]*m)))
    logli.sum<-c(logli.sum,logli)
    #print(logli)
    #print(i)
  }
  logli.sum<-cbind(sum(logli.sum),a)
  resul<-rbind(resul,logli.sum)
} 

##plot
'resul.baidu<-data.frame(resul,type="b) Baidu")
resul.tencent<-data.frame(resul,type="a) Tencent")
#dev.off()
#par(mfrow=c(1,2))
#par(las=1)
#plot(x=resul[,2],y=resul[,1],xlab="# of 2019-nCoV cases",ylab="log-likelihood", type="l")
#title(main = list("(a) Tencent data (θ = 1, λ = 1415)", cex = 1.5,col = "black", font = 1))
#title(main = list("(b) Baidu data (θ = 1, λ = 5574)", cex = 1.5,col = "black", font = 1))
resul.com<-rbind(resul.tencent,resul.baidu)
ggplot ()+
  geom_line(data = resul.com, aes(x=a, y=V1)) +
  facet_wrap(~type,scales="free")+
  labs(x= "# of COVID-19 cases",y = "log-likelihood")+ #,title = "effectiveness of different actions"
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),
        axis.text=element_text(size=21),axis.title.x=element_text(size=21),
        axis.title.y=element_text(size=21),legend.text=element_text(size=13),legend.title = element_text(size=13),
        strip.text.x = element_text(size = 23))
ggsave(filename = "est.pdf",width = 17,height = 8)'


resul[which.max(resul[,1]),]
k<-as.numeric(resul[which.max(resul[,1]),2])
print(k)
see<-resul[which.max(resul[,1]),1]-0.5*qchisq(0.95,1)

for(i in 1:(length(resul[,1])-1)){
  #i<-1
  x1<-as.numeric(resul[i,1])-see
  mark<-as.numeric(resul[i,2])
  x2<-as.numeric(resul[i+1,1])-see
  if(x1*x2<0){
    print(mark)
  }
}
data$est_total<-k
data$est_num<-data$P*data$est_total
data$est_log<-log(dbinom(ceiling(data$est_num),data$est_total,data$P))
#AIC=(2k-2L)/n 
#k<-5
#n<-24
data$AIC<-0
for(a in c(1:length(data$flu))){
  data$AIC[a]<-(2*1-2*data$est_log[a])
}
sum(data$AIC[-c(23,4)])
sum(data$AIC)



######scaling test
'
b.pop<-data[,c(1,3)]
t.pop<-data[,c(1,3)]
pop.com<-join(t.pop,b.pop,type="left",by="place")
colnames(pop.com)[c(2,3)]<-c("pop.tencent","pop.baidu")
cor.test(pop.com$pop.tencent,pop.com$pop.baidu,method = c("spearman"))
plot(pop.com[-25,2],pop.com[-25,3])
line<-lm(pop.com[-25,3]~pop.com[-25,2])
summary(line)
line$coefficients
pop.new<-pop.com
pop.new$pop.tencent<-pop.new$pop.tencent*line$coefficients[2]+line$coefficients[1]
pop.com$situation<-"a) Before processing"
pop.new$situation<-"b) After processing"
pop.add<-rbind(pop.com,pop.new)
ggplot ()+
  geom_point(data = pop.add[-c(25,50),], aes(x=pop.baidu, y=pop.tencent)) +
  facet_wrap(~situation,scales="free")+
  labs(x= "pop.baidu",y = "pop.tencent")+ #,title = "effectiveness of different actions"
  scale_x_log10(limits = c(50, 500000))+
  scale_y_log10(limits = c(50, 500000))+
  geom_abline(intercept = 0, slope = 1,linetype="dashed",color="black")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),
        axis.text=element_text(size=21),axis.title.x=element_text(size=21),
        axis.title.y=element_text(size=21),legend.text=element_text(size=13),legend.title = element_text(size=13),
        strip.text.x = element_text(size = 23))
ggsave(filename = "scaling.pdf",width = 15,height = 7)


#rate<-1
data<-data[-25,]
data<-data[-27,]
data$forec<-0
data$P<-0
'