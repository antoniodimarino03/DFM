library(dfms)
library(xts)
QD <- read.csv("FREDQD.csv")
df <- QD[, colSums(is.na(QD)) == 0]
df <- df[-1, ]
n_obs <- nrow(df) - 1  
n_test <- 8  
n_train <- n_obs - n_test  
train_set <- df[1:(1 + n_train), ]
test_set <- df[(1 + n_train + 1):nrow(df), ]
write.csv(train_set, "train.csv", row.names = FALSE)
write.csv(test_set, "test.csv", row.names = FALSE)
x <- read.csv("TRAINSTAN.csv", header=F)
ic <- ICr(x)
screeplot(ic)
