
################################################################
## Prepare dataset for NG copula 

## Written by Himchan Jeong, Rosy Oh 

## Last updated : 28 Mar 2024

################################################################


#---------------------------------------------------------------
## train dataset 

load("data.rdata")


#### CO and CN ####
pdata = subset(data, CoverageCO*CoverageCN>0)
pdata$exposure = 1
summed = aggregate(cbind(FreqCO, FreqCN, ClaimCO,ClaimCN,
                         exposure, TypeCity, TypeCounty, 
                         TypeSchool, TypeTown, TypeVillage) ~ 
                     PolicyNum, data= pdata, sum)
# head(summed)

summedCO = summed[,-c(3,5)]
summedCO$Class <- "O"
colnames(summedCO)[2:3] <- c("N", "S")
# head(summedCO)

summedCN <- summed[,-c(2,4)]
summedCN$Class <- "N"
colnames(summedCN)[2:3] <- c("N", "S")

tmp_train = rbind(summedCO, summedCN)
tmp_train[,5:9] <- tmp_train[,5:9]/tmp_train[,4] 
train_dat = tmp_train%>%
  mutate(Class = factor(Class),
         M = S/N)%>%
  filter(S!=1)
 

train_dat$M[is.nan(train_dat$M)] <- 0


#----------------------------------------------------------
## test dataset 

load("dataout.rdata")
# head(dataout)

tmp = dataout %>% 
  select(PolicyNum, FreqCO, FreqCN, ClaimCO, ClaimCN, 
         TypeCity, TypeCounty, TypeSchool:TypeVillage)

tmp_o = tmp%>%
  select(-ends_with("CN"))%>%
  rename(N=FreqCO, S=ClaimCO)%>%
  mutate(Class="O")

tmp_n = tmp%>%
  select(-ends_with("CO"))%>%
  rename(N=FreqCN, S=ClaimCN)%>%
  mutate(Class="N")


test_dat = rbind(tmp_o, tmp_n)%>%
  mutate(M = S/N)%>%
  filter(PolicyNum %in% c(train_dat$PolicyNum))
test_dat$M[is.nan(test_dat$M)] <- 0

# dim(test_dat)   # 716  10

