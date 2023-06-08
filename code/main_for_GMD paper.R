# SCRIPT TO REPRODUCE THE WORK IN MENGOLI ET. AL. 2023 GMD PAPER SUBMISSION

# RUN THE P MODEL AT SUB-DAILY TIMESCALE AT 67 SITES (EACH SITE ONE RUN)----
# USING THE SUB-DAILY VERSION OF THE P MODEL (Mengoli et al. 2022, JAMES Journal)
# model input from fluxnet 2015: air temperature (Ta), vapour pressure deficit (VPD), ambient CO2, light (PPFD)
# model input from MODIS: fraction of absorbed photosynthetically active radiation (fAPAR)
# model output: half-hourly gross primary production, GPP (well-watered predicted GPP: GPPpOpt1)
# updated script to use: main_running_mean_GMD paper,
# which refers to the updated functions in: .R_for_main_running_mean_GMD folder
# for the users: use the example settings files located in settings_main_script_rm folder
#                with the example fluxnet and MODIS data, in data/Data_for_sub-dailyPmodel folder

# RUN THE SPLASHv1 MODEL----
# USING THE PUBLISHED VERSION OF THE MODEL (Devis et. al. 2017)
# at this link: https://bitbucket.org/labprentice/splash/src/52d9454b566daf4f7e4f02b9d795aedd99fb53ff/releases/v1.0/r_version/?at=master
#               https://doi.org/10.5281/zenodo.376293
# model input from CRU TS4.06: lat, long, elevation, air temperature, precipitation, 
#                              and cloud cover to derive the fraction of sunlight hours (sf = 1- CLD) (Harris et al. 2014) 
# model output: daily soil moisture (swc), which then needs to be divided by the bucket size (150 mm) to obtain the relative swc
#               daily potential evapotranspiration (PET), that is needed to then compute the aridity index (AI = PET/Prec)

# AGGREGATE THE HALF-HOURLY GPP TO DAILY (EACH SITE SEPARATELY)----
# USING THE AGGREGATE FUNCTION
# script to use: dailyAggr.R
# the function generates a list: 'df_mean' for the daily averages;'df_nr' for the number of data used for averaging
# call the function, apply it to the dataset and then use the values in df_mean
# 
setwd('C:/pModelPlus/pModelPlus')
source ('./R/dailyAggr.R')    

dataOr1 = dailyAggr(dataOr)  
dataOr1 = dataOr1$df_mean    # data.frame with daily average

# MERGE THE DATASET----
# each site, merge the DAILY P model dataset (output) with the DAILY SPLASHv1 dataset (output)
# to have in each dataset (each site) GPPs from the P model, swc and AI from the SPLASHv1 model
# example name for the dataset: AU-GWW_test_3a2013-2014_CRUOPT_1

# APPLY THE GPP QC FILTER TO THE OBSERVED GPP (ALL SITES IN ONE RUN)----
# TO THEN COMPUTE THE BETA THETA RATIO
# the code below applies a TH of 0.8 for the NEE_CUT_REF_QC variable in fluxnet
# for the users: place all the dataset (of each site) in the list.files path 
#                example name for the dataset: AU-GWW_test_3a2013-2014_CRUOPT_1
#                the datasets must be daily and must have the observed and modelled GPP
#                place all the daily fluxnet dataset in the list.file2 path
#                two daily fluxnet datasets are provided in the FLUXNET_data_daily folder

rm(list = ls())
cat('\014')

list.files ('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/dataOr_splash_ratio_NewSites/')

for (lf in list.files('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/dataOr_splash_ratio_NewSites/')){
  
  sitecode = substr(lf,1,6)  
  
  df = read.csv(paste0('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/dataOr_splash_ratio_NewSites/', lf))
  
  list.files2 = list.files('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/FLUXNET_sites/DAILY', pattern = sitecode)
  
  dfc = read.csv(paste0('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/FLUXNET_sites/DAILY/', list.files2))
  dfc$YEAR = NA
  for (n in seq (1, nrow(dfc))){
    dfc$YEAR[n] = as.numeric(substr(as.character(dfc$TIMESTAMP[n]),1,4))
  }
  df$QC = NA
  for ( yy in sort(unique (df$YEAR))){
    pos1 = which(df$YEAR == yy) 
    pos2 = which(dfc$YEAR == yy)
    if ('NEE_CUT_REF_QC' %in% colnames(dfc)) {
      df$QC[pos1] = dfc$NEE_CUT_REF_QC[pos2]
    } else {
      df$QC[pos1] = dfc$NEE_VUT_REF_QC[pos2]
    }
  }
  
  pos_QC = which(df$QC < 0.8)
  if (length(pos_QC) > 0 ) df = df[-1*pos_QC,]
  rm(pos_QC)
  
  write.csv(df, paste0('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/dataOr_splash_ratio_QC_NewSites/', lf), row.names = F)
}

# PROCESS THE BETA THETA RATIO (ALL SITES IN ONE RUN)----
# TO REMOVE THOSE UNREALISTIC TOO HIGH VALUES 
# AND TO COMPUTE THE BETA THETA RATIO (OBSERVED GPP / WELL-WATERED PREDICTED GPP)
# the code below uses two filters: 
# removal of the positions where SWC = 1 and where predicted GPP is too low (i.e. 5th percentile)
# for the users: place all the datasets (of each site) in the working directory
#                where the example name for the dataset is: AU-GWW_test_3a2013-2014_CRUOPT_1
#                the code below will then merge all datasets in only one dataset, dataOr_t, with the beta theta ratio (BetaThetaObs)
#                in addition, the code needs the AI dataset as input 
#                and in data folder user can find an example file with AI values for two sites

# LOAD dataOr with daily predicted GPP (GPPpOpt1) and the observed GPP filtered for the QC
setwd('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/dataOr_splash_ratio_QC_NewSites')

# CREATE a df to see how many rows have been deleted
df_stat = data.frame("nomesito"= NA, "TH_05"= NA, "Min_GPP" = NA, "nrow"= NA, "ntot"= NA, "nrow2"= NA, "MAX_BETA" = NA, "AI" = NA,
                     'pos_na_swc' = NA,# % of rows with: which(is.na(dataOr$swc) == 1)
                     'pos_na_GPPpOpt1' = NA,	# % of rows with: which(is.na(dataOr$GPPpOpt1) == 1)
                     'pos_na_GPP_DT_CUT_REF' = NA,	# % of rows with: which(is.na(dataOr$GPP_DT_CUT_REF) == 1)
                     'pos_zero_GPPpOpt1' = NA,	# % of rows with: which(dataOr$GPPpOpt1 == 0)
                     'pos_neg_GPPpOpt1' = NA,	# % of rows with: which((dataOr$GPPpOpt1) <0)
                     'pos_over0p99_swc' = NA,	# % of rows with:which(dataOr$swc >= 0.99)
                     'pos_lt_quantile0p05_GPPpOpt1' = NA	# % of rows with:con which(dataOr$GPPpOpt1 < quantile(dataOr$GPPpOpt1, na.rm = T, probs = 0.05, names= F))
)

if(exists("dataOr_t")) rm(dataOr_t)
# for ( nomesito in c('AR-SLu_test_3a2009-2011_CRUOPT_1.csv',
#                     'AR-Vir_test_3a2009-2012_CRUOPT_1.csv',
#                     'AU-ASM_test_3a2010-2013_CRUOPT_1.csv')){
for ( nomesito in list.files ('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/dataOr_splash_ratio_QC_NewSites')){
  
  cat(sprintf('%s\n',substr(nomesito,1,6)))
  
  dataOr = read.csv(nomesito)
  nr_tot = nrow(dataOr)
  dataOr$TIME = ymd(as.character(dataOr$TIME))
  
  # dataOr$AI = sum(dataOr$AET)/sum(dataOr$PET)
  
  # LOAD 'AI' DATASET: one column with the site ID (e.g.AR-SLu)
  #                       another column with the aridity index values (AI = PET/Prec)
  #                       20 years records for computing the AI, where Prec are from CRU dataset while the PET is simulated by SPLASHv1
  AI_sites = read.csv('C:/Users/giuli/Desktop/Tiger team/AI_decades/AI_33+9+25Sites.csv', sep = ';')
  
  dataOr$AI = NA
  pos = which(AI_sites$sito == substr(nomesito,1,6))
  if (length(pos) > 0 ) AI = AI_sites$AI[pos]
  rm(pos)
  
  dataOr$AI = AI
  rm(AI)
  
  # A) NA
  pos_na_swc = which(is.na(dataOr$swc) == 1)
  if (length(pos_na_swc) > 0 ) dataOr = dataOr[-1*pos_na_swc,]
  
  pos_na_GPPpOpt1 = which(is.na(dataOr$GPPpOpt1) == 1)
  if (length(pos_na_GPPpOpt1) > 0 ) dataOr = dataOr[-1*pos_na_GPPpOpt1,]
  
  pos_na_GPP_DT_CUT_REF = which(is.na(dataOr$GPP_DT_CUT_REF) == 1)
  if (length(pos_na_GPP_DT_CUT_REF) > 0 ) dataOr = dataOr[-1*pos_na_GPP_DT_CUT_REF,]
  
  # B) INF
  pos_zero_GPPpOpt1 = which(dataOr$GPPpOpt1 == 0)
  if (length(pos_zero_GPPpOpt1) > 0 ) dataOr = dataOr[-1*pos_zero_GPPpOpt1,]
  
  # C) NEG
  pos_neg_GPPpOpt1 = which((dataOr$GPPpOpt1) <0)
  if (length(pos_neg_GPPpOpt1) > 0 ) dataOr = dataOr[-1*pos_neg_GPPpOpt1,]
  
  # swc == 1
  pos_over0p99_swc = which(dataOr$swc >= 0.99)
  if (length(pos_over0p99_swc) > 0 ) dataOr = dataOr[-1*pos_over0p99_swc,]
  
  # 5% bottom
  TH_05 = quantile(dataOr$GPPpOpt1, na.rm = T, probs = 0.05, names= F)
  MIN_GPPp = min(dataOr$GPPpOpt1, na.rm = T)
  
  pos_lt_quantile0p05_GPPpOpt1 = which(dataOr$GPPpOpt1 < quantile(dataOr$GPPpOpt1, na.rm = T, probs = 0.05, names= F))
  
  if (length(pos_lt_quantile0p05_GPPpOpt1) > 0 ) dataOr = dataOr[-1*pos_lt_quantile0p05_GPPpOpt1,]
  
  # # BETA THETA RATIO
  dataOr$BetaThetaObs = dataOr$GPP_DT_CUT_REF/ dataOr$GPPpOpt1
  
  df_stat = rbind(df_stat, data.frame("nomesito"= nomesito,
                                      "TH_05"= TH_05,
                                      "Min_GPP" = MIN_GPPp,
                                      "nrow" = nrow(dataOr),
                                      "ntot"= nr_tot,
                                      "nrow2"= nrow(dataOr)*0.05,
                                      "MAX_BETA" = max(dataOr$BetaThetaObs, na.rm = T), 
                                      "AI" = sum(dataOr$AET)/sum(dataOr$PET),
                                      'pos_na_swc' = length(pos_na_swc),#which(is.na(dataOr$swc) == 1)
                                      'pos_na_GPPpOpt1' = length(pos_na_GPPpOpt1),	# which(is.na(dataOr$GPPpOpt1) == 1)
                                      'pos_na_GPP_DT_CUT_REF' = length(pos_na_GPP_DT_CUT_REF),	#which(is.na(dataOr$GPP_DT_CUT_REF) == 1)
                                      
                                      'pos_zero_GPPpOpt1' = length(pos_zero_GPPpOpt1),	#which(dataOr$GPPpOpt1 == 0)
                                      
                                      'pos_neg_GPPpOpt1' = length(pos_neg_GPPpOpt1),	#which((dataOr$GPPpOpt1) <0)
                                      
                                      'pos_over0p99_swc' = length(pos_over0p99_swc),	#which(dataOr$swc >= 0.99)
                                      'pos_lt_quantile0p05_GPPpOpt1' = length(pos_lt_quantile0p05_GPPpOpt1)	#which(dataOr$GPPpOpt1 < quantile(dataOr$GPPpOpt1, na.rm = T, probs = 0.05, names= F))
  )) 
  
  if ( exists("dataOr_t")){
    dataOr_t = rbind(dataOr_t, dataOr)
  } else{
    dataOr_t = dataOr
  }
  rm(pos_na_swc,pos_na_GPPpOpt1,pos_na_GPP_DT_CUT_REF,pos_zero_GPPpOpt1)
  rm(pos_neg_GPPpOpt1,pos_over0p99_swc,pos_lt_quantile0p05_GPPpOpt1,nr_tot) 
}
rm(MIN_GPPp, nomesito, TH_05)
rm(dataOr)

df_stat = df_stat[-1,]

# %
df_stat$pos_na_swc = (df_stat$pos_na_swc / df_stat$ntot)*100
df_stat$pos_na_GPPpOpt1 = (df_stat$pos_na_GPPpOpt1 / df_stat$ntot)*100
df_stat$pos_na_GPP_DT_CUT_REF = (df_stat$pos_na_GPP_DT_CUT_REF / df_stat$ntot)*100
df_stat$pos_zero_GPPpOpt1 = (df_stat$pos_zero_GPPpOpt1 / df_stat$ntot)*100
df_stat$pos_neg_GPPpOpt1 = (df_stat$pos_neg_GPPpOpt1 / df_stat$ntot)*100
df_stat$pos_over0p99_swc = (df_stat$pos_over0p99_swc / df_stat$ntot)*100
df_stat$pos_lt_quantile0p05_GPPpOpt1 = (df_stat$pos_lt_quantile0p05_GPPpOpt1 / df_stat$ntot)*100

range(dataOr_t$BetaThetaObs)

# SEGMENTED CURVES COMPUTATION (ALL SITES IN ONE RUN)----
# USING THE BREAKPOINT REGRESSION ANALYSIS (https://rpubs.com/MarkusLoew/12164: segmented package)
# TO OBTAIN THE MAXIMUM LEVEL (y_psi) AND THE BREAKING POINT (psi) OF THE BETA THETA RATIO
# the code below uses two approaches: 
# the fixed, where the intercept is imposed equal to zero, and the free approach

cat('\014')
# loading of packages
library(lubridate)
library(ggplot2)
library(cowplot)
library(scales)
#Sys.setlocale("LC_TIME", "English")
library(magrittr) 
library(dplyr)    
library(segmented)

# LOAD dataOr_t generated by the previous section
dataOr = dataOr_t
# sites
unique_sito = unique(dataOr$sito)
#regression results with the 'free' approach
df_tot_free = data.frame('sito'=NA, 'veg'= NA, 'clim'= NA,
                         'psi' = NA,'psi_stErr' = NA, 
                         'y_psi' = NA,
                         'nr_psi'= NA,
                         'slope1' = NA,
                         'slope1_stErr' = NA,
                         'intercept1' = NA,
                         'slope2' = NA,
                         'slope2_stErr' = NA,
                         'intercept2' = NA,
                         'interc' = NA,
                         'interc_stErr' = NA)
#regression results with the 'fixed' approach (i.e. intercept equal to zero)
df_tot_i0 = data.frame('sito'=NA, 'veg'= NA, 'clim'= NA,
                       'psi' = NA,'psi_stErr' = NA, 
                       'y_psi' = NA,
                       'nr_psi'= NA,
                       'slope1' = NA,
                       'slope1_stErr' = NA,
                       'intercept1' = NA,
                       'slope2' = NA,
                       'slope2_stErr' = NA,
                       'intercept2' = NA,
                       'interc' = NA,
                       'interc_stErr' = NA)

#regression fitted values
dataOr$fitted_i0 = NA
dataOr$fitted_free = NA
# adding the ID row to make the merge
dataOr$id_riga = seq(1,nrow(dataOr))

for (nome_sito in unique_sito) {
  
  dataOr13 = dataOr[which(dataOr$sito == nome_sito),]
  # filters: NA, INF, neg
  # A) NA
  pos_na = which(is.na(dataOr13$swc) == 1)
  if (length(pos_na) > 0 ) dataOr13 = dataOr13[-1*pos_na,]
  pos_na = which(is.na(dataOr13$BetaThetaObs) == 1)
  if (length(pos_na) > 0 ) dataOr13 = dataOr13[-1*pos_na,]
  rm(pos_na)
  
  # B) INF
  pos_inf = which(is.infinite(dataOr13$BetaThetaObs) == 1)
  if (length(pos_inf) > 0 ) dataOr13 = dataOr13[-1*pos_inf,]
  rm(pos_inf)
  
  # B) NEG
  pos_inf = which((dataOr13$BetaThetaObs) <0)
  if (length(pos_inf) > 0 ) dataOr13 = dataOr13[-1*pos_inf,]
  rm(pos_inf)
  
  # B) 
  pos_inf = which((dataOr13$BetaThetaObs) >1)
  if (length(pos_inf) > 0 ) dataOr13 = dataOr13[-1*pos_inf,]
  rm(pos_inf)
  
  # regression
  for (opt_intercept in c('fixed_zero','free')) {
    if ( opt_intercept == 'fixed_zero')# fixed intercept regression
      lin_regr_model = lm(BetaThetaObs ~ 0 + swc, data= dataOr13)  
    if ( opt_intercept == 'free')# free intercept regression
      lin_regr_model = lm(BetaThetaObs ~ swc, data= dataOr13)
    
    ch_point = tryCatch({ segmented(lin_regr_model,seg.Z = ~swc)}, 
                        error = function(e) {
                          return (NA)},
                        warning = function(w) {
                          # return (segmented(lin_regr_model,seg.Z = ~swc))})
                          return (NA)})
    
    if (is.list(ch_point)) {
      ch_point_summary = summary(ch_point)
      
      df_tot = data.frame('sito'=nome_sito,'veg'=unique(dataOr13$veg),'clim'=unique(dataOr13$clim),
                          'psi' = ch_point$psi[2],
                          'psi_stErr' = ch_point$psi[3], 
                          'nr_psi'= length(ch_point$psi),
                          'slope1' = slope(ch_point)$swc[1,1],
                          'slope1_stErr' = slope(ch_point)$swc[1,2],
                          'intercept1' = intercept(ch_point)$swc[1],
                          'slope2' = slope(ch_point)$swc[2,1],
                          'slope2_stErr' = slope(ch_point)$swc[2,2],
                          'intercept2' = intercept(ch_point)$swc[2],
                          'interc' = ch_point_summary[["Ttable"]][1,1],
                          'interc_stErr' =ch_point_summary[["Ttable"]][1,2] )
      # when swc > change_point values are constant
      dataOr13$fitted_values  = NA
      dataOr13$fitted_values  = df_tot$slope1 * dataOr13$swc + df_tot$intercept1
      # maximum level computation (using psi)
      y_psi = df_tot$slope1 * df_tot$psi + df_tot$intercept1
      df_tot$y_psi = y_psi
      # replace y_psi values
      over_psi = which(dataOr13$swc > df_tot$psi)
      if ( length(over_psi) > 0) dataOr13$fitted_values[over_psi] = y_psi
      rm(y_psi,over_psi)
      if ( opt_intercept == 'fixed_zero') {
        df_tot_i0 = rbind(df_tot_i0,df_tot)
        for (nr in seq(1,nrow(dataOr13))) {
          id_r = dataOr13$id_riga[nr]
          dataOr$fitted_i0[which(dataOr$id_riga == id_r)] = dataOr13$fitted_values[nr]
        }
        rm(nr,id_r)
      }
      if ( opt_intercept == 'free') {
        df_tot_free = rbind(df_tot_free,df_tot)
        for (nr in seq(1,nrow(dataOr13))) {
          id_r = dataOr13$id_riga[nr]
          dataOr$fitted_free[which(dataOr$id_riga == id_r)] = dataOr13$fitted_values[nr]
        }
        rm(nr,id_r)
      }
      rm(df_tot,ch_point_summary)  
    }
    rm(lin_regr_model,ch_point)
  }
  rm(opt_intercept)
}
rm(dataOr13,nome_sito,unique_sito)

df_tot_i0 = df_tot_i0[-1,]
df_tot_free = df_tot_free[-1,]

# NON-LINEAR REGRESSION ANALYSIS----
# TO OBTAIN THE QUANTITATIVE FUNCTION TO USE TO MODIFY THE PREDICTED GPP
# NB: to make the non linear regression analysis include the AI values in 'df_tot_i0'

# packages
library(drc)
library(nlme)
library(aomisc)

# remove the outliers 
pos = which(df_tot_i0$sito == 'AU-Lox')
df_tot_i0 = df_tot_i0[-1*pos,] 
pos = which(df_tot_i0$sito == 'AU-RDF')
df_tot_i0 = df_tot_i0[-1*pos,] 
pos = which(df_tot_i0$sito == 'US-Wkg')
df_tot_i0 = df_tot_i0[-1*pos,]

# power curves of psi vs AI (threshold point)
#              and of y_psi vs AI (maximum level)
model1 <- drm(psi ~ AI, fct = DRC.powerCurve(),
              data = df_tot_i0)
summary(model1)
model2 <- drm(y_psi ~ AI, fct = DRC.powerCurve(),
              data = df_tot_i0)
summary(model2)

# APPLY THE NEW SOIL MOISTURE STRESS FUNCTION (ALL SITES IN ONE RUN)----
# AND THE STOCKER ET AL 2020 GMD function 
# TO OBTAIN THE REVISED DAILY PREDICTED GPP (GPPpOpt2)
# NB: to run the code below, use the P model daily GPP (GPPpOpt1) 
# and the relative daily swc from SPLASHv1, merged in one dataframe (each site one dataset)
# without applying the above filters for the beta theta ratio computation or the QC on GPPobs
# for the users: place all the dataset (of each site) in the list.files path 
#                example name for the dataset: AU-GWW_test_3a2013-2014_CRUOPT_1
#                the datasets must be daily and must have the observed, modelled GPP and the swc
#                place all the daily fluxnet dataset in the list.files path
#                in addition, the code needs the AI dataset as input 

rm(list = ls())
cat('\014')
# loading of packages
library(lubridate)
library(ggplot2)
library(cowplot)
library(scales)
Sys.setlocale("LC_TIME", "English")
library(magrittr) 
library(dplyr)    

# CREATE an empty df and load DAILY DATASET 
dataOr_t = data.frame('TIMESTAMP_START_ok' = NA,
                      'TIMESTAMP_START' = NA,
                      'TIME'= NA,
                      'YEAR'= NA,
                      'MONTH'= NA,
                      'DAY'= NA,'HOUR'= NA,'MINUTE'= NA,'WEEK'= NA,'GPP_DT_CUT_REF'= NA,
                      'GPPpOpt1'= NA,'BetaThetaObs'= NA,'PET'= NA,'AET'= NA,'swc'= NA,'swcObs'= NA,'sf'= NA,
                      'sito' = NA,'cluster' = NA,'clim' = NA,'veg' = NA,'AI' = NA,'psi' = NA,'y_psi' = NA,'GPPpOpt2' = NA,
                      'GPPpOpt3' = NA,'ft' = NA,
                      'ft3' = NA,'y_max' = NA,'y_min' = NA,'GPP_DT_VUT_REF' = NA)

# setwd('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/dataOr_splash_ratio_NewSites')
for (nomefile in list.files('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/dataOr_splash_ratio_NewSites')){

  cat(sprintf('%s\n',substr(nomefile,1,6)))
  dataOr = read.csv(nomefile)
  
  for (n in seq(1,nrow(dataOr))) 
    dataOr[n,'TIME'] = sprintf('%04d%02d%02d',dataOr[n,'YEAR'],dataOr[n,'MONTH'],dataOr[n,'DAY'] )
  dataOr$TIME = as.numeric(dataOr$TIME)
  
  # LOAD 'AI' DATASET: one column with the site ID (e.g.AR-SLu)
  #                       another column with the aridity index values (AI = PET/Prec)
  #                       20 years records for computing the AI, where Prec are from CRU dataset while the PET is simulated by SPLASHv1
  AI_sites = read.csv('C:/Users/giuli/Desktop/Tiger team/AI_decades/AI_33+9+25Sites.csv', sep = ';')
  
  dataOr$AI = NA
  pos = which(AI_sites$sito == substr(nomefile,1,6))
  if (length(pos) > 0 ) AI = AI_sites$AI[pos]
  rm(pos)
  
  dataOr$AI = AI
  
  # STOCKER'S function (Stocker et. al. 2020 GMD)
  # CONSTANTS
  a = 0.0
  b = 0.73300
  x0 = 0.0
  x1 = 0.6
  
  alpha = 1/AI
  
  y0 <- a + b * alpha
  q = (1 - y0 ) / (x0 - x1)^2
  outstress = 1.0 - q * (dataOr$swc - x1)^2
  
  outstress <- pmin(pmax(outstress, 0), 1)
  outstress[ dataOr$swc > x1] <- 1.0
  rm(a,b,x0,x1)
  
  # FUNCTION TO APPLY (ft3)
  dataOr$ft3 = outstress
  dataOr$GPPpOpt3 = NA
  
  dataOr$GPPpOpt3 = dataOr$GPPpOpt1 * dataOr$ft3 
  
  # NEW FUNCTION (Mengoli et. al.)
    # CRITICAL THRESHOLD POINT (psi)
    x1 = rep(1, nrow(dataOr))
    # x2 = 0.34*dataOr$AI^(-0.53)
    x2 = 0.34*dataOr$AI^(-0.60)
    dataOr$psi = pmin(x1,x2)
    
    # MAXIMUM LEVEL (y_psi)
    y1 = rep(1, nrow(dataOr))
    # y2 = 0.62*dataOr$AI^(-0.39)
    y2 = 0.62*dataOr$AI^(-0.45)
    dataOr$y_psi = pmin(y1,y2)
    
    # FUNCTION TO APPLY (ft)
    dataOr$ft = outstress
    dataOr$GPPpOpt2 = NA
  
      # TWO CONDITIONS TO APPLY THE NEW FX
      # 1) THETA >= PSI
      pos = which(dataOr$swc >= dataOr$psi)
      if (length(pos) > 0 ) dataOr$ft[pos] = dataOr$y_psi[pos]
      
      dataOr$GPPpOpt2[pos] = dataOr$GPPpOpt1[pos] * dataOr$y_psi[pos]
      rm(pos)
      
      # 2) THETA <= PSI
      pos = which(dataOr$swc <= dataOr$psi)
      if (length(pos) > 0 ) dataOr$ft[pos] = (dataOr$y_psi[pos]/dataOr$psi[pos]) *dataOr$swc[pos]
      
      dataOr$GPPpOpt2[pos] = dataOr$GPPpOpt1[pos] * dataOr$ft[pos]
      rm(pos)
  
  
  # DAILY FLUXNET
  fluxnet_file = list.files('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/FLUXNET_sites/DAILY', pattern =substr(nomefile,1,6) )
  cat(sprintf('%s\n', fluxnet_file))
  dfc = read.csv(paste0('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/FLUXNET_sites/DAILY/',fluxnet_file))
  
  dataOr$y_max = NA
  dataOr$y_min = NA
  dataOr$GPP_DT_VUT_REF = NA
  
  for ( yy in sort(unique (dataOr$TIME))){
    pos1 = which(dataOr$TIME == yy) 
    pos2 = which(dfc$TIMESTAMP == yy)
    
    if ('GPP_DT_CUT_05' %in% colnames(dfc)) {
      dataOr$y_max[pos1] = dfc$GPP_DT_CUT_05[pos2]
      dataOr$y_min[pos1] = dfc$GPP_DT_CUT_95[pos2]
      dataOr$GPP_DT_VUT_REF[pos1] = dfc$GPP_DT_CUT_REF[pos2]
    } else {
      dataOr$y_max[pos1] = dfc$GPP_DT_VUT_05[pos2]
      dataOr$y_min[pos1] = dfc$GPP_DT_VUT_95[pos2]
      dataOr$GPP_DT_VUT_REF[pos1] = dfc$GPP_DT_VUT_REF[pos2]
    }
  }
  
  dataOr = dataOr[, c('TIMESTAMP_START_ok','TIMESTAMP_START','TIME','YEAR',
                      'MONTH','DAY','HOUR','MINUTE','WEEK','GPP_DT_CUT_REF',
                      'GPPpOpt1','BetaThetaObs','PET','AET','swc','swcObs','sf',
                      'sito','cluster','clim','veg','AI','psi','y_psi','GPPpOpt2',
                      'GPPpOpt3','ft',
                      'ft3','y_max','y_min','GPP_DT_VUT_REF')]
  
  write.csv (dataOr, sprintf('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/apply_function/dataOrCRU_%s.csv',substr(nomefile,1,6) ) , row.names = F,quote = F)
  
  dataOr_t = rbind(dataOr_t, dataOr)
  
}

dataOr_t = dataOr_t[-1,]

# APPLY THE GPP QC FILTER TO THE OBSERVED GPP (ALL SITES IN ONE RUN)----
# TO THEN COMPUTE THE STATISTICS
# for the user: the example name for the datasts: dataOrCRU_AU-GWW

rm(list = ls())
cat('\014')

list.files ('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/apply_function/')
for (lf in list.files('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/apply_function/')){
  
  sitecode = substr(lf,11,16)  
  
  df = read.csv(paste0('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/apply_function/', lf))   
  
  list.files2 = list.files('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/FLUXNET_sites/DAILY', pattern = sitecode)
  
  dfc = read.csv(paste0('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/FLUXNET_sites/DAILY/', list.files2))
  dfc$YEAR = NA
  for (n in seq (1, nrow(dfc))){
    dfc$YEAR[n] = as.numeric(substr(as.character(dfc$TIMESTAMP[n]),1,4))
  }
  df$QC = NA
  for ( yy in sort(unique (df$YEAR))){
    pos1 = which(df$YEAR == yy) 
    pos2 = which(dfc$YEAR == yy)
    if ('NEE_CUT_REF_QC' %in% colnames(dfc)) {
      df$QC[pos1] = dfc$NEE_CUT_REF_QC[pos2]
    } else {
      df$QC[pos1] = dfc$NEE_VUT_REF_QC[pos2]
    }
  }
  
  pos_QC = which(df$QC < 0.8)
  if (length(pos_QC) > 0 ) df = df[-1*pos_QC,]
  rm(pos_QC)
  
  write.csv(df, paste0('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/after_fx_QC/', lf), row.names = F)
}

# COMPUTE THE STATISTICS----
# script to analyse the results after the function application
# datasets already filtered for the poor quality of GPP

rm(list = ls())

# loading of packages
library(lubridate)
library(ggplot2)
library(cowplot)
library(scales)
library("hydroGOF")

# LOAD dataOr with predicted GPP (GPPpOpt1,GPPpOpt2, GPPpOpt3)
setwd('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/after_fx_QC')

df_stat = data.frame("nomesito"= NA,  "rmse1"= NA, "rmse2"= NA, "rmse3"= NA,
                     "r1" = NA, "r2" = NA, "r3" = NA,
                     "r2ww" = NA, "r2Mf"= NA, "r2Bf"= NA)

if(exists("dataOr_t")) rm(dataOr_t)

for ( nomesito in list.files ('C:/pModelPlus/drylands/FLUXNET2015_Feb2020/Pmodel_rm/after_fx_QC')){
  
  dataOr = read.csv(nomesito)
  dataOr$TIME = ymd(as.character(dataOr$TIME))
  
  # A) NA
  pos_na = which(is.na(dataOr$GPPpOpt1) == 1)
  if (length(pos_na) > 0 ) dataOr = dataOr[-1*pos_na,]
  rm(pos_na)
  pos_na = which(is.na(dataOr$GPPpOpt2) == 1)
  if (length(pos_na) > 0 ) dataOr = dataOr[-1*pos_na,]
  rm(pos_na)
  pos_na = which(is.na(dataOr$GPPpOpt3) == 1)
  if (length(pos_na) > 0 ) dataOr = dataOr[-1*pos_na,]
  rm(pos_na)
  pos_na = which(is.na(dataOr$GPP_DT_VUT_REF) == 1)
  if (length(pos_na) > 0 ) dataOr = dataOr[-1*pos_na,]
  rm(pos_na)
  
  # # B) NEG
  # pos_inf = which((dataOr$GPPpOpt1) <0)
  # if (length(pos_inf) > 0 ) dataOr = dataOr[-1*pos_inf,]
  # rm(pos_inf)
  # pos_inf = which((dataOr$GPPpOpt2) <0)
  # if (length(pos_inf) > 0 ) dataOr = dataOr[-1*pos_inf,]
  # rm(pos_inf)
  # pos_inf = which((dataOr$GPPpOpt3) <0)
  # if (length(pos_inf) > 0 ) dataOr = dataOr[-1*pos_inf,]
  # rm(pos_inf)
  # pos_inf = which((dataOr$GPP_DT_VUT_REF) <0)
  # if (length(pos_inf) > 0 ) dataOr = dataOr[-1*pos_inf,]
  # rm(pos_inf)
  
  
  GPPpOpt1 = dataOr$GPPpOpt1
  GPPpOpt2 = dataOr$GPPpOpt2
  GPPpOpt3 = dataOr$GPPpOpt3
  GPP_DT_CUT_REF = dataOr$GPP_DT_VUT_REF
  
  rmse1 = rmse(GPPpOpt1, GPP_DT_CUT_REF)
  rmse2 = rmse(GPPpOpt2, GPP_DT_CUT_REF)
  rmse3 = rmse(GPPpOpt3, GPP_DT_CUT_REF)
  r1 = unname(cor.test(GPPpOpt1, GPP_DT_CUT_REF)$estimate) 
  r2 = unname(cor.test(GPPpOpt2, GPP_DT_CUT_REF)$estimate)
  r3 = unname(cor.test(GPPpOpt3, GPP_DT_CUT_REF)$estimate) 
  r2ww = r1*r1
  r2Mf = r2*r2
  r2Bf = r3*r3
  
  nomesito = substr(nomesito,11,16)
  
  df_stat = rbind(df_stat, data.frame("nomesito"= nomesito,
                                      "rmse1"= rmse1, "rmse2"= rmse2, "rmse3"= rmse3,
                                      "r1" = r1, "r2" = r2, "r3" = r3,
                                      "r2ww" = r2ww, "r2Mf"= r2Mf, "r2Bf"= r2Bf
  )) 
  
  
  if ( exists("dataOr_t")){
    dataOr_t = rbind(dataOr_t, dataOr)
  } else{
    dataOr_t = dataOr
  }
}

df_stat = df_stat[-1,]

rm(r1,r2,r3,r2Bf,r2Mf,r2ww,rmse1,rmse2,rmse3, nomesito, GPP_DT_CUT_REF, GPPpOpt1, GPPpOpt2, GPPpOpt3)
