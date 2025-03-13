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

dat<-fread("~/Measles TX 2025 - From web.csv")
dat2<-copy(dat)
dat2[,Date:=as.Date(Date, format='%m/%d/%Y')]

#make weekly incidence for plotting
dat2[,Week := as.Date(cut(dat2$Date, "week"))]
dat2[,Weekly_cases:=sum(Daily_cases), by=Week]
dat3<-dat2[,.(Weekly_cases, Week)]
dat3<-dat3[!duplicated(dat3)]
dat3[,cum_sum:=cumsum(Weekly_cases)]

#plot incident vs cumulative weekly cases (Figure 1 in text)

#scaling for two axes
ylim.prim <- c(0, 60)
ylim.sec <- c(0,200) 
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

png("~/Measles_TX_cases_by_week.png", width=700, height=600)
ggplot(data=dat3)+geom_col(aes(x=Week, y=Weekly_cases), fill='#440154', alpha=0.8)+theme_bw()+
  geom_line(aes(x=Week, y= a + cum_sum*b), col='#fde725', lwd=1.5)  + geom_vline(xintercept=as.Date("2025-02-25"), lty=2, col='#21918c', lwd=1.2)+
  scale_y_continuous("Weekly reported measles cases", sec.axis = sec_axis(~ (. - a)/b, name = "Cumulative reported measles cases")) +
  coord_cartesian(clip = "off", expand=F)+
  scale_x_date(date_breaks="1 week", date_labels = "%b %d")+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),axis.title=element_text(size=14,face="bold"))
dev.off()

#Priors derived from the literature
#SIs: 9, 14.5
#SD = 3.25
#Prior R0 = 14.5 (median value)
#STD = 1.5 (to cover range of 11 to 18)

#R0 values for vaccination coverage calculations: 11, 14.5 (median), 18


#create a grid of all parameter options (SI and then also R0 for vax coverage)
si<-c(9,14.5)
r0<-c(11,14.5,18)
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
#plot the R(t) estimates by the two SIs
png('~/Measles_TX_SI9_SI14.png')
ggplot(data=r_samp_df2)+geom_line(aes(x=Date, y=R_median,col=factor(si)),lwd=3)+
  geom_ribbon(aes(ymin=R_025, ymax=R_975, x=Date,fill=factor(si)), alpha=0.5)+theme_bw()+
  scale_color_viridis_d("Serial interval")+scale_fill_viridis_d("Serial interval")+ylab(expression("Effective reproductive number R"[t]))+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),axis.title=element_text(size=14,face="bold"))+
  coord_cartesian(clip = "off", expand=F)
dev.off()

#table of R(t), VC, and coverage gap (from a threshold of 93%) in a table for our six scenarios for most recent day of data
setDT(medians)
medians2<-medians[!is.na(r0)]
medians2[,gap:=93-VC_median]
medians2[,gap_u:=93-VC_025]
medians2[,gap_l:=93-VC_975]
medians2[,si_lab:=paste("SI",si)]
medians2[,si:=NULL]
medians2[,blank_rowname := purrr::map(list(rep(" ", 1)), gt::html)]

medians2|>
  gt(rowname_col = "blank_rowname", groupname_col = "si_lab") |>
  cols_merge_uncert(
    R_025,
    R_975,
    rows = everything(),
    sep = " - ",
    autohide = TRUE
  )|>
  cols_merge_uncert(
    VC_025,
    VC_975,
    rows = everything(),
    sep = " - ",
    autohide = TRUE
  )|>
  cols_merge_uncert(
    gap_l,
    gap_u,
    rows = everything(),
    sep = " - ",
    autohide = TRUE
  )|>
  cols_width(
    c("VC_025","R_025","gap_l") ~ px(150)
  )|>
  cols_label( r0="R(0)",
              median_rt="Estimated Median R(t)", 
              R_025='95% UI R(t)',
              VC_median="Estimated Median Vaccine Coverage",
              gap="Median gap in coverage to reach 93%",
              VC_025='95% UI Vaccine Coverage',
              gap_l="95% UI gap in coverage")|>
  tab_header(
    title = "Vaccine Coverage for most recent R(t) values\nby Serial Interval and R(0)",
    subtitle = "Using a prior R(0) value  14.5 (SD = 1.5) and vaccine efficacy estimated at 95%")|>
  cols_align(
    align = c("center"),
    columns = everything())|>
  cols_move(
    columns = VC_median,
    after = R_025
  ) |>
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_stub()
  ) |>
  fmt_number(columns = c(VC_median,median_rt,gap, VC_025, VC_975, R_025,R_975,gap_u,gap_l), decimal=1 , suffixing = TRUE)|>
  gtsave("~/clatest_VC_table_six_scenarios.png",expand=10, vwidth=1800, vheight=1000)



setDT(medians)
medians2<-medians[!is.na(r0)]
medians2[,gap:=93-VC_median]
medians2[,gap_u:=93-VC_025]
medians2[,gap_l:=93-VC_975]
medians2[,si_lab:=paste("SI",si)]
medians2[,si_lab0:=paste("",si)]
medians2[,si:=NULL]

#grouped bar plot
medians2[,si_lab0:=factor(si_lab0, c(" 9"," 14.5"))]
png("~/Measles_TX_barchart_vc.png")
ggplot(medians2) +
  geom_col(aes(x=si_lab0, y=VC_median, fill = factor(r0)), position = position_dodge(0.8), 
           color = "black", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(x=si_lab0,ymin = VC_025, ymax = VC_975, group = factor(r0)), 
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_viridis_d(expression("Basic reproduction number R"[0])) +
  theme_bw()+xlab('Serial Interval')+
  labs(y = "Vaccination Coverage (%)") + geom_hline(aes(yintercept=93),lty=2)+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  ylim(c(0,100))+coord_cartesian(clip = "off", expand=F)
dev.off()


#Line plot of R0 versus VC by SI
png("~/Measles_TX_linechart_vc.png", width=650, height=600)
ggplot(data=medians2)+geom_line(aes(x=r0, y=VC_median,col=factor(si_lab0)),lwd=3)+
  geom_ribbon(aes(ymax=VC_025, ymin=VC_975, x=r0,fill=factor(si_lab0)), alpha=0.5)+theme_bw()+
  geom_point(aes(x=r0,y=VC_median, col=factor(si_lab0), size=1.5))+
  scale_color_viridis_d("Serial interval")+scale_fill_viridis_d("Serial interval")+ylab("Vaccine coverage (%)")+
  geom_hline(aes(yintercept=93), lty=2)+
  theme_bw()+scale_size_continuous(guide="none")+annotate(geom="text",x=17, y=94, label="Herd immunity threshold (93%)")+
  labs(x = expression("Basic reproduction number R"[0]), y = "Vaccination coverage (%)") + geom_hline(aes(yintercept=93),lty=2)+ theme(axis.line = element_line())+
  
  #add y-axis discontinuity piece using the below:
  scale_y_continuous(limits =c(50, 100),
                     breaks = seq(50, 100, by = 10),
                     labels = c(0, seq(60, 100, by = 10)))+
  annotate("segment", x = c(10.73, 10.73), xend = c(10.87,10.87), y = c(53, 54), yend = c(56, 57))+
  annotate("segment", x = c(10.8,10.8), xend = c(10.8, 10.8), y = c(50, 52.9), yend = c(56.1, 100))+
  annotate("segment", x = c(18.1,18.1), xend = c(18.1, 18.1), y = c(50, 50), yend = c(100, 100))+
  annotate("rect", xmin = c(10.7999,10.7999), xmax = c(10.8, 10.8), ymin = c(54.55,54.55), ymax = c(55.4, 55.4), col='white', fill='white')+
  annotate("segment", x = c(10.8,10.8), xend = c(18.1, 18.1), y = c(50, 50), yend = c(50,50))+
  annotate("segment", x = c(10.8,10.8), xend = c(18.1, 18.1), y = c(100, 100), yend = c(100,100))+
  annotate("segment", x = c(10.73, 10.73), xend = c(10.87,10.87), y = c(53, 54), yend = c(56, 57))+
  coord_cartesian(clip = "off", xlim=c(10.8,18.1), expand=F)+
  theme_minimal()
dev.off()
#line plot of VC over time for various R0 values (SI of 14.5 ) 


png("~/Measles_TX_timeseries_VC2.png", width=650, height=600)
ggplot(data=r_samp_df2[si %in% c(14.5)])+geom_line(aes(x=Date, y=VC_median,col=factor(r0)),lwd=3)+
  geom_ribbon(aes(ymax=VC_025, ymin=VC_975, x=Date,fill=factor(r0)), alpha=0.5)+theme_bw()+
  scale_color_viridis_d(expression("Basic reproduction number R"[0]))+scale_fill_viridis_d(expression("Basic reproduction number R"[0]))+ylab("Vaccine coverage (%)")+
  geom_hline(aes(yintercept=93), lty=2)+
  theme_bw()+annotate(geom="text",x=as.Date('2025-02-28'), y=94, label="Herd immunity threshold (93%)")+
  scale_y_continuous(limits =c(00, 100))+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  coord_cartesian(clip = "off", xlim=c(as.Date("2025-02-25"),as.Date("2025-03-11")), expand=F)
dev.off()


#facet grid
r_samp_df2[,si_lab0:=paste("SI", si)]
r_samp_df2[,si_lab0:=factor(si_lab0, c("SI 9", "SI 14.5"))]

png("~/Measles_TX_timeseries_VC_facets.png", width=750, height=600)
ggplot(data=r_samp_df2)+geom_line(aes(x=Date, y=VC_median,col=factor(r0)),lwd=3)+
  geom_ribbon(aes(ymax=VC_025, ymin=VC_975, x=Date,fill=factor(r0)), alpha=0.5)+theme_bw()+
  scale_color_viridis_d(expression("Basic Reproduction Number R"[0]))+scale_fill_viridis_d(expression("Basic Reproduction Number R"[0]))+ylab("Vaccine coverage (%)")+
  geom_hline(aes(yintercept=93), lty=2)+facet_wrap(.~si_lab0)+
  theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=12),axis.title=element_text(size=16,face="bold"))+
  theme_bw()+annotate(geom="text",x=as.Date('2025-03-02'), y=94, label="Herd immunity threshold (93%)")+scale_x_date(breaks="1 week", date_labels = '%b %d')+
  scale_y_continuous(limits =c(0, 100))+theme(strip.background =element_rect(fill="white"), strip.text =element_text(size=12, face='bold'))+ coord_cartesian(clip = "off", expand=F)
dev.off()
