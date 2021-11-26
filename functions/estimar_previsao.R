library(lubridate)
library(dplyr)

previsao <- read.table("DATA/previsao_doses_12_08.txt", sep = "\t", encoding = "UTF-8")

previsao$V1 <- factor(previsao$V1)
az = previsao %>% filter(V1 == "AZ") %>% select(-V1) %>% colSums()
cv = previsao %>% filter(V1 == "CV") %>% select(-V1) %>% colSums()
pf = previsao %>% filter(V1 == "PF") %>% select(-V1) %>% colSums()

df <- data.frame(coronavac = cv,
                 astrazeneca = az,
                 pfizer = pf)

seqdays = seq(as.Date("2021-08-01"),as.Date("2021-12-31"),by = "days")

perda = 0.95

df2 <- data.frame(matrix(0,length(seqdays),3))
colnames(df2) <- c("coronavac","astrazeneca","pfizer")
df2$date <- seqdays

this_month <- unique(month(df2$date))[1]
next_month <- unique(month(df2$date))[2]

df2[month(df2$date) == this_month,1:3] <- matrix( unlist(round(df[1,]/sum(month(df2$date) == this_month)*perda)), sum(month(df2$date) == this_month),3,byrow = T)
                                        
df2[month(df2$date) == next_month,1:3] <- matrix(unlist(round(df[2,]/sum(month(df2$date) == next_month)*perda)), sum(month(df2$date) == next_month),3,byrow = T)
                                        
df2[month(df2$date) > next_month,1:3] <- matrix(unlist(round(df[3,]/sum(month(df2$date) > next_month)*perda)), sum(month(df2$date) > next_month),3,byrow = T)
                                       
plot(df2$date, df2$coronavac, type = "l", ylim = c(0,max(df2[,1:3])), 
     xlim = as.Date(c("2021-07-01","2021-12-31")))
lines(df2$date, df2$astrazeneca, col = "red")
lines(df2$date, df2$pfizer, col = "blue")

write.csv(df2, file = "rollout_vacinas.csv", row.names = FALSE, quote = FALSE)
