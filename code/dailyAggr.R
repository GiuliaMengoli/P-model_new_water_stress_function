# FUNCTION for compute the daily values starting from half-hourly or hourly data
dailyAggr <- function(data_in = data_in,showMsg = F) {

  col_to_skip = c("TIMESTAMP_START_ok","TIMESTAMP_START","TIME","YEAR","MONTH","DAY","HOUR","MINUTE","WEEK")

  data_daily_mean = data_in[1,]
  data_daily_nr = data_in[1,]

  while (nrow(data_in) > 0) {
    pos_day = which(
      data_in$YEAR == data_in$YEAR[1] &
        data_in$MONTH == data_in$MONTH[1] &
        data_in$DAY == data_in$DAY[1]
    )
    if (showMsg)
      cat(sprintf('year: %d; month: %d; day: %d (row to process: %d)\n',
                  data_in$YEAR[1],data_in$MONTH[1],data_in$DAY[1],nrow(data_in)))

    tmp_mean = data_in[pos_day[1],]
    tmp_nr = data_in[pos_day[1],]

    for (col in colnames(data_in)) {
      if (col %in% col_to_skip) next
      if (!is.numeric(tmp_mean[,col])) next
      tmp_nr[,col] = length(pos_day) - sum(is.na(data_in[pos_day,col]))
      tmp_mean[,col] = NA
      if (tmp_nr[,col] > 0) tmp_mean[,col] = mean(data_in[pos_day,col],na.rm = T)
    }
    rm(col)
    data_daily_mean = rbind(data_daily_mean,tmp_mean)
    data_daily_nr = rbind(data_daily_nr,tmp_nr)
    data_in = data_in[-1*pos_day,]
    rm(tmp_nr,tmp_mean,pos_day)
  }
  data_daily_mean = data_daily_mean[-1,]
  data_daily_nr = data_daily_nr[-1,]

  # GPP CONVERSION: micromoli/m2/s to gC/m2/day
  col_gpp = grep('^GPP',colnames(data_daily_mean))
  if (length(col_gpp) > 0) {
    cat('\n')
    # ((peso molecolare C) / (peso molecolare CO2))*(44*10^-6) * (3600 s/ora) *
    # 24(ore/day)
    for (cy_col_gpp in col_gpp) {
      cat(sprintf('convert daily mean from umolCO2 m-2 s-1 to gC m-2 d-1 variable: %s\n',
                  colnames(data_daily_mean)[cy_col_gpp]))
      data_daily_mean[,cy_col_gpp] = data_daily_mean[,cy_col_gpp] *((12/44)*(44*10^-6) * 3600 * 24)
    }
  }

  return(list('df_mean' = data_daily_mean,'df_nr' = data_daily_nr))

}


