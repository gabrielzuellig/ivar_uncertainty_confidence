library(gdata)
library(data.table)

datain <- read.xls('Shadow rates all economies_DeRezende.xlsx', sheet = 3)
datain$yearmon <- substr(as.character(datain$datesd), 1, 6)
datain <- data.table(datain)
data <- datain[,j=list(OIS.1m.rate = mean(OIS.1m.rate, na.rm=T),
                       ShadowRatep2 = mean(Shadow.rate..p.2, na.rm=T),
                       ShadowRatep3 = mean(Shadow.rate..p.3, na.rm=T)),
               by=list(yearmon)]
write.fwf(data, file='Shadow rates all economies_DeRezende_monthly.csv')