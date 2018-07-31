## ref: https://statcompute.wordpress.com/2017/09/03/variable-selection-with-elastic-net/

pkgs <- list("glmnet", "doParallel", "foreach", "PRROC")
lapply(pkgs, require, character.only = T)

setwd('~/Desktop/EXP-16/')
T20_classifier <- read.table('./classifier.txt', header=T)
head(T20_classifier)
dim(T20_classifier)
T20_classifier <- T20_classifier[!is.na(T20_classifier$PAPAprop), ]
dim(T20_classifier)
lapply(T20_classifier, class)


##### elastic net model 
##### elastic net model 
##### elastic net model 
T20_classifier_sig <- T20_classifier

T20_classifier_sig <- T20_classifier[T20_classifier$padj<=0.05, ]
T20_classifier_sig$padj <- NULL
T20_classifier_sig$log2FoldChange <- NULL
T20_classifier_sig$PROTlen <- NULL
T20_classifier_sig$FI <- NULL
T20_classifier_sig$FInumaa <- NULL
T20_classifier_sig$c2 <- NULL
T20_classifier_sig$ptb <- NULL
head(T20_classifier_sig)
dim(T20_classifier_sig)

# scale relevant variables
T20_classifier_sig[,1:4] <- as.data.frame(lapply(T20_classifier_sig[,1:4], scale) )

## make factors for relevant variables
for (i in 5:ncol(T20_classifier_sig) ){
  T20_classifier_sig[,i] <- as.factor(T20_classifier_sig[,i])
}

str(T20_classifier_sig)

# make dummy variables
set.seed(321)
dmy <- dummyVars(" ~ .", data = T20_classifier_sig,fullRank = T)
train_transformed <- data.frame(predict(dmy, newdata = T20_classifier_sig))
head(train_transformed)
str(train_transformed)

#Converting the dependent variable back to categorical
train_transformed$Residency.1 <- as.factor(train_transformed$Residency.1)
str(train_transformed)

# partition data
index <- createDataPartition(train_transformed$Residency.1, p=0.70, list=FALSE)
trainSet <- train_transformed[ index,]
testSet <- train_transformed[-index,]

#parameter tuning using caret (for alpha)
grid<-expand.grid(.alpha=seq(0,1,by=.05),.lambda=seq(0,1,by=.05))
custom <- trainControl(method = 'repeatedCV',
                       number = 10, repeats =5,
                       verboseIter = F)
outcomeName <- 'Residency.1'
predictors <- names(trainSet)[!names(trainSet) %in% outcomeName]

en <- train(x = trainSet[,predictors], y = trainSet[,outcomeName], 
            method ="glmnet", trControl = custom, tuneGrid = grid)
en

#plot important variables
plot(varImp(en))
attributes(en)
en$bestTune
plot(en)


### find optimal lambda & 1SE lambda using cv from glmnet
outcomeName<-'Residency.1'
predictors <- names(testSet)[!names(testSet) %in% outcomeName]
fitCV <- cv.glmnet(x = as.matrix(testSet[,predictors]), y = testSet[,outcomeName], 
                   family = 'binomial', type.measure = 'auc', nfolds = 10, 
                   alpha = en$bestTune$alpha )

plot(fitCV)
fitCV$lambda.1se
fitCV$lambda.min
coef(fitCV, s = 'lambda.1se')
coef(fitCV, s = 'lambda.min')

# use model on test sample
predCV <- predict(fitCV, newx = as.matrix(testSet[,predictors]),
                  s = "lambda.1se",
                  type = "response")

actuals <- test$Residency

fg <- predCV[actuals == 1]
bg <- predCV[actuals == 0]
fg <- fg[!is.na(fg)]
bg <- bg[!is.na(bg)]

# ROC Curve
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
par(pin=c(2,2), tcl=0.25, ps=9, family="Helvetica")
plot(roc, col='black', lwd=2)
lines(c(0, 1), c(0, 1))

# PR Curve
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(pr, col ='black', lwd=2)
