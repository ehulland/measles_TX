#############################################################################
############ Read in relevant packages ######################################
#############################################################################
library(EpiEstim)
library(ggplot2)
library(data.table)
library(gt)
library(ggh4x)

#############################################################################
#####  TX measles by day of report      ####################################
#############################################################################

dat<-fread("~//Data/Measles TX 2025 from web.csv")
dat2<-copy(dat)
dat2[,Date:=as.Date(Date, format='%m/%d/%Y')]

#create a grid of ALL parameter options for supplement (SI and then also R0 for vax coverage)
si<-c(7,9,10, 11, 12,13, 14, 14.5, 21)
r0<-c(10,seq(11,27,1),14.5)
grid<-expand.grid(si=si,r0=r0)
r_samp_df<-data.frame(R_025=NA, R_median=NA, R_975=NA,Date=as.Date('2025-01-01',format='%Y-%m-%d'), si=NA, r0=NA, VC_median=NA, VC_025=NA, VC_975=NA)
medians<-data.frame(si=NA, r0=NA, median_rt=NA, VC_median=NA, R_025=NA, R_975=NA, VC_025=NA, VC_975=NA)
for(i in 1:nrow(grid)){
  si<-grid$si[i]
  r0<-grid$r0[i]
  
  #iterating 
  res_parametric_si<- estimate_R(dat2$Daily_cases,
                                 method="parametric_si",
                                 config = make_config(list(
                                   mean_si = si,
                                   std_si = 3.25,
                                   mean_prior = 14.5,
                                   std_prior = 1.5)))
  
  
  res_parametric_si$dates<-dat2$Date
  R_samp<-res_parametric_si$R
  
  #weekly window, so estimates start on 8th day of data
  R_samp$Date<-res_parametric_si$dates[8:length(res_parametric_si$dates)]
  #keep only 2.5th, 97.5th and median percentiles
  r_samp<-R_samp[,c(5,8,11,12)]
  names(r_samp)<-c("R_025","R_median","R_975", "Date")
  r_samp$si<-si
  r_samp$r0<-r0
  
  setDT(r_samp)
  #calculate vaccine coverage using R posterior estimates, 95% vaccine effectiveness, and iterating over values of R0
  r_samp[,VC_median:=(1-(R_median/r0))/0.95*100]
  r_samp[,VC_025:=(1-(R_025/r0))/0.95*100]
  r_samp[,VC_975:=(1-(R_975/r0))/0.95*100]
  
  
  #get all estimates for the most recent day of data: 
  median<-data.frame(si=si, r0=r0)
  median$median_rt<-r_samp[Date==max(Date), R_median]
  median$R_025<-r_samp[Date==max(Date), R_025]
  median$R_975<-r_samp[Date==max(Date), R_975]
  
  median$VC_median<-r_samp[Date==max(Date),VC_median]
  median$VC_025<-r_samp[Date==max(Date), VC_975]
  median$VC_975<-r_samp[Date==max(Date), VC_025]
  
  medians<-rbind(medians, median)
  r_samp_df<-rbind(r_samp_df,r_samp)
  
}

#focus on the period after which reporting became routine: 2025-02-25
setDT(r_samp_df)
r_samp_df2<-r_samp_df[!is.na(si) & Date>="2025-02-25"]

png('~//Measles_TX_bySI_all.png')
ggplot(data=r_samp_df2)+geom_line(aes(x=Date, y=R_median,col=factor(si)),lwd=3)+
  geom_ribbon(aes(ymin=R_025, ymax=R_975, x=Date,fill=factor(si)), alpha=0.5)+theme_bw()+
  scale_color_viridis_d("SI")+scale_fill_viridis_d("SI")+ylab(expression("R"[t]))+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  coord_cartesian(clip = "off", expand=F)
dev.off()

r_samp_df2[,r0_lab:=paste0("R(0):", r0)]
r_samp_df2[VC_025<0,VC_025:=0]
r_samp_df2[VC_975<0,VC_975:=0]
r_samp_df2[VC_median<0,VC_median:=0]

png("~//measles_TX_vaccine_coverage_cropped_all.png", width=1200, height=800)
ggplot(data=r_samp_df2)+geom_line(aes(x=Date, y=VC_median,col=factor(si)),lwd=3)+
  geom_ribbon(aes(ymax=VC_025, ymin=VC_975, x=Date,fill=factor(si)), alpha=0.5)+theme_bw()+
  scale_color_viridis_d("SI")+scale_fill_viridis_d("SI")+facet_grid(.~r0_lab)+ylab("Vaccination coverage (%)")+
  ylim(c(0,100))+geom_hline(aes(x=Date, yintercept=93), lty=2)+ theme(axis.text.x = element_text(angle = 45, vjust=0.5))+
  theme(axis.text.y=element_text(size=12), axis.title=element_text(size=14,face="bold"))+
  scale_x_date(date_breaks='1 week', date_labels="%b %d")+theme(strip.background =element_rect(fill="white"))+  coord_cartesian(clip = "off", expand=F)
dev.off()


setDT(medians)
medians2<-medians[!is.na(r0)]
medians2[,gap:=93-VC_median]
medians2[,gap_u:=93-VC_025]
medians2[,gap_l:=93-VC_975]
medians2[,si_lab:=paste("SI",si)]
medians2[,si:=NULL]
medians2[,blank_rowname := purrr::map(list(rep(" ", 1)), gt::html)]
medians2[,si_lab:=factor(si_lab, c("SI 7", "SI 9","SI 10", "SI 11","SI 12", "SI 13", "SI 14","SI 14.5", "SI 21"))]
medians2[,r0f:=factor(r0, c("10", "11","12","13","14","14.5","15","16","17","18", "19", "20","21","22","23","24","25","26","27"))]


#ggplot bar plot grouped

medians2[,size:=ifelse(r0 %in% c(11,14.5,18), 1, 0)]
png("~//Measles_TX_barchart_vc_all.png", width=800, height=500)
ggplot(medians2) +
  geom_col(aes(x=si_lab, y=VC_median, fill = factor(r0f)), position = position_dodge(0.8), 
           color = "black", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(x=si_lab,ymin = VC_025, ymax = VC_975, group = factor(r0f)), 
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_viridis_d(expression("R"[0])) +
  theme_bw()+xlab('Serial interval')+
  labs(x = NULL, y = "Vaccination coverage (%)") + geom_hline(aes(yintercept=93),lty=2)+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  theme(axis.text.x = element_text(angle = 45, vjust=0.5))+  coord_cartesian(clip = "off", expand=F)+
  ylim(c(0,100))
dev.off()


#line plot with all points from1 11 to 18 for main manuscript Figure 2
medians3<-medians2[r0>=11 & r0<=18,]
png("~//Measles_TX_linechart_vc_multi.png", width=650, height=600)
ggplot(data=medians3[si_lab %in% c("SI 9","SI 14.5")])+geom_line(aes(x=r0, y=VC_median,col=factor(si_lab)),lwd=3)+
  geom_ribbon(aes(ymax=VC_025, ymin=VC_975, x=r0,fill=factor(si_lab)), alpha=0.5)+theme_bw()+
  geom_point(aes(x=r0,y=VC_median, col=factor(si_lab), size=factor(size)))+
  scale_color_viridis_d("Serial interval", direction = 1)+scale_fill_viridis_d("Serial interval", direction=1)+ylab("Vaccination coverage (%)")+
  geom_hline(aes(yintercept=93), lty=2)+scale_size_manual(guide="none", values=c(NA,5))+
  theme_bw()+annotate(geom="text",x=17, y=94, label="Herd immunity threshold (93%)")+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  labs(x = expression("R"[0]), y = "Vaccination coverage (%)") + geom_hline(aes(yintercept=93),lty=2)+ theme(axis.line = element_line())+
  scale_y_continuous(limits =c(50, 100),
                     breaks = seq(50, 100, by = 10),
                     labels = c(0, seq(60, 100, by = 10)))+
  annotate("segment", x = c(10.73, 10.73), xend = c(10.87,10.87), y = c(53, 54), yend = c(56, 57))+
  annotate("segment", x = c(10.8,10.8), xend = c(10.8, 10.8), y = c(50, 52.9), yend = c(56.1, 100))+
  annotate("segment", x = c(18.1,18.1), xend = c(18.1, 18.1), y = c(50, 50), yend = c(100, 100))+
  annotate("segment", x = c(10.8,10.8), xend = c(18.1, 18.1), y = c(50, 50), yend = c(50,50))+
  annotate("segment", x = c(10.8,10.8), xend = c(18.1, 18.1), y = c(100, 100), yend = c(100,100))+
  annotate("segment", x = c(10.73, 10.73), xend = c(10.87,10.87), y = c(53, 54), yend = c(56, 57))+
  coord_cartesian(clip = "off", xlim=c(10.8,18.1), expand=F)+
  theme_minimal()
dev.off()

#now do one with ALL r0 values:

png("~//Measles_TX_linechart_vc_multi_ALL.png", width=600, height=400)
ggplot(data=medians2[si_lab %in% c("SI 9","SI 14.5")])+geom_line(aes(x=r0, y=VC_median,col=factor(si_lab)),lwd=3)+
  geom_ribbon(aes(ymax=VC_025, ymin=VC_975, x=r0,fill=factor(si_lab)), alpha=0.5)+theme_bw()+
  geom_point(aes(x=r0,y=VC_median, col=factor(si_lab), size=factor(size)))+
  scale_color_viridis_d("Serial Interval")+scale_fill_viridis_d("Serial Interval")+ylab("Vaccination coverage (%)")+
  geom_hline(aes(yintercept=93), lty=2)+scale_size_manual(guide="none", values=c(1,3))+
  theme_bw()+annotate(geom="text",x=13, y=94, label="Herd immunity threshold (93%)")+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  labs(x = expression("R"[0]), y = "Vaccination coverage (%)") + geom_hline(aes(yintercept=93),lty=2)+ theme(axis.line = element_line())+
  scale_y_continuous(limits =c(50, 100),
                     breaks = seq(50, 100, by = 10),
                     labels = c(0, seq(60, 100, by = 10)))+
  annotate("segment", x = c(9.73, 9.73), xend = c(10.07,10.07), y = c(53, 54), yend = c(56, 57))+
  annotate("segment", x = c(9.9,9.9), xend = c(9.9, 9.9), y = c(50, 52.9), yend = c(56.1, 100))+
  annotate("segment", x = c(27.1,27.1), xend = c(27.1, 27.1), y = c(50, 50), yend = c(100, 100))+
  annotate("segment", x = c(9.9,9.9), xend = c(27.1, 27.1), y = c(50, 50), yend = c(50,50))+
  annotate("segment", x = c(9.9,9.9), xend = c(27.1, 27.1), y = c(100, 100), yend = c(100,100))+
  annotate("segment", x = c(9.73, 9.73), xend = c(10.07,10.07), y = c(53, 54), yend = c(56, 57))+
  coord_cartesian(clip = "off", xlim=c(9.9,27.1), expand=F)+
  theme_minimal()
dev.off()
