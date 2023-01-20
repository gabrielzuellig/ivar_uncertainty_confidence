

indat <- read.table('datashares_selection.csv', header=T, sep=',')

# consumer outlook
cons <- indat[,grep('CONS', colnames(indat))]

Q <- cons[,grep('2', colnames(cons))]

mean <- Q[,1]+0.5*Q[,2]-0.5*Q[,4]+(-1)*Q[,5]
outdat <- sqrt(Q[,1]*(1 - mean)^2 + Q[,2]*(0.5 - mean)^2 + Q[,4]*(-0.5 - mean)^2 + Q[,5]*(-1 - mean)^2)

plot(1:431, outdat, 'l')


# industry
indu <- indat[,grep('INDU', colnames(indat))]
questions <- c(1:5)
disp <- matrix(NA, nrow=nrow(indu), ncol=length(questions))

for (q in 1:length(questions)){
  Q <- indu[,grep(as.character(questions[q]), colnames(indu))]
  fracplus <- Q[, 1]/100
  fracmin <- Q[,3]/100
  disp[,q] <- 100*(fracplus + fracmin - (fracplus - fracmin)^2)^(1/2)
  plot(1:431, disp[,q], 'l')
}

outdat <- cbind(outdat, rowMeans(disp))
plot(1:431, outdat[,2], 'l')

# export
colnames(outdat) <- c('ConsDisp','InduDisp')
outdat <- ts(outdat, freq=12, start=c(1985,1))
write.table(outdat, 'out.csv', sep=';')

