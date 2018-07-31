library(caret)
library(corrplot)
library(PRROC)

setwd('~/Desktop/EXP-16/')
### this T20_classifier was NOT preprocessed already in Rstudio in our server by scaling and centering
T20_classifier <- read.table('./classifier.txt', header=T)
head(T20_classifier)
dim(T20_classifier)
T20_classifier <- T20_classifier[!is.na(T20_classifier$PAPAprop), ]
dim(T20_classifier)
lapply(T20_classifier, class)



elastic_net_cat <- function(x, probs, lambdaType){
  # some clean up of the data
  x <- x[x$padj <= 0.05, ]
  x$padj <- NULL
  x$log2FoldChange <- NULL
  x$PROTlen <- NULL
  x$FI <- NULL
  x$FInumaa <- NULL
  x$c2 <- NULL
  x$ptb <- NULL
  
  x[,1:4] <- as.data.frame(lapply(x[,1:4], scale) )
  
  # make factors for relevant variables
  for (i in 5:ncol(x) ){
    x[,i] <- as.factor(x[,i])
  }
  
  # make dummy variables
  dmy <- dummyVars(" ~ .", data = x, fullRank = T)
  train_transformed <- data.frame(predict(dmy, newdata = x))
  
  
  # Converting the dependent variable back to categorical
  train_transformed$Residency.1 <- as.factor(train_transformed$Residency.1)
  
  # partition data
  index <- createDataPartition(train_transformed$Residency.1, p=probs, list=FALSE)
  trainSet <- train_transformed[ index,]
  testSet <- train_transformed[-index,]

  # parameter tuning using caret (for alpha)
  grid <- expand.grid(.alpha=seq(0.01,.5,by=.02),.lambda=seq(0,1,by=.02))
  custom <- trainControl(method = 'repeatedCV',
                         number = 10, repeats =5,
                         verboseIter = F)
  
  outcomeName <- 'Residency.1'
  predictors <- names(trainSet)[!names(trainSet) %in% outcomeName]
  en <- train(x = trainSet[,predictors], y = trainSet[,outcomeName], 
              method ="glmnet", trControl = custom, tuneGrid = grid)
  
  print(en$bestTune)
  
  ### find optimal lambda & 1SE lambda using cv from glmnet
  outcomeName<-'Residency.1'
  predictors <- names(testSet)[!names(testSet) %in% outcomeName]
  fitCV <- cv.glmnet(x = as.matrix(testSet[,predictors]), y = testSet[,outcomeName], 
                     family = 'binomial', type.measure = 'auc', nfolds = 10, 
                     alpha = en$bestTune$alpha )
  
  output <- data.frame(coefficients = coef(fitCV, s = lambdaType)[,1])  
  return(output) 
}


col_size = 100
results <- matrix(0, nrow = 15, ncol = col_size)
for (i in 1:col_size){
  dat <- elastic_net_cat(T20_classifier, 0.7, 'lambda.1se')
  results[,i] <- dat$coefficients
  row.names(results) <- row.names(dat)
}

dim(results)
head(results)

results.final <- data.frame(mean_coef = rowMeans(results), row.names = row.names(results))

sem <- matrix(0, nrow = 15, ncol = 1)
for (i in 1:nrow(results)){
  sem[i,]  <- sd(results[i,])/sqrt(length(results[i,]))
}

results.final$sem <- sem
head(results.final)
results.final <- results.final[-1, ]
results.final
## Makes a barplot of coefficients
par(pin=c(2,2), tcl=0.25, ps=9, family="Helvetica")
library(RColorBrewer)

barPlotter <- function(x){
  x$dir <- NA
  x[x$mean_coef > 0, 'dir'] = 'pos'
  x[x$mean_coef <= 0, 'dir'] = 'neg'
  cos_up <- x[x$dir == 'pos',1]
  cos_up <- data.frame(coef = cos_up, names = row.names(x)[x$dir=='pos'])
  
  cos_down <- x[x$dir == 'neg',1]
  cos_down <- data.frame(coef = cos_down, names = row.names(x)[x$dir=='neg'])
  
  cos_up <- cos_up[order(cos_up$coef, decreasing=F), ]
  cos_down <- cos_down[order(cos_down$coef, decreasing=F),]
  
  both <- rbind(cos_down, cos_up)
  both$sem <- x$sem[match(both$names, row.names(x) )]
  
  cols <- colorRampPalette(c("palevioletred",'lightgray' ,"slateblue"))
  barCenters <- barplot(both$coef, horiz=T, col =cols(nrow(both)),
                        names=both$names, xlim=c(min(both$coef)*1.5, max(both$coef)*1.5), xlab=c('Coefficients') )
  
  segments(both$coef - both$sem, barCenters, both$coef + both$sem, barCenters,
           lwd = 1.5)
  
  #arrows(both$coef - both$sem, barCenters, both$coef + both$sem, barCenters, 
  #       lwd = 1.5, angle = 90,
  #       code = 3, length = 0.05)
  return(both)
}

barPlotter(results.final)
