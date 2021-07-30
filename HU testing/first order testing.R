# load in libraries 
library(readxl) # read in excel files
library(dplyr) # preform aggregation 
library(xlsx)  # export files into an excel workbook with sheets



# read in the data file by sheet name
df <- read_excel("20210728 combine_df2_first order kinetic.xlsx",sheet = "nmol")

# remove the extra info based on rows and columns we do not need 
df <- df[-c(1:13),-c(1)]

# Change colname of one column to x --> time in seconds 
colnames(df)[colnames(df) == "...2"] <- "x"

# remove all columns that have missing values NA
df<-df %>%   select(where(~!any(is.na(.))))
dim(df)
str(df)



# save as a .csv file to use later 
write.csv(df, "df.csv")
library(readr)
df <- read_csv("df.csv")
str(df)

# https://stackoverflow.com/questions/51181888/r-extracting-model-coefficients-from-a-nested-list-list-columns
# code outline

storage <- list()    # create an empty list for the coeff
storage2 <- list()   # create an empty list for the r^2 values
for(i in names(df)[-1]){
  storage[[i]] <- lm(get(i) ~ x, df)
  storage2[[i]] <- summary(storage[[i]])$r.squared
}
storage  # get the linear reg equation
storage2 # get the respected r^2 values of each column(y)

# put the r^2 all into a organized data frame to export into excel file
r_2 <-as.data.frame(do.call(rbind, storage2));colnames(r_2)[colnames(r_2) == "V1"] <- "r_2"   # Convert list to data frame rows
R_2 <- r_2
R_2 <- R_2 %>% filter(r_2 >= 0.8)

R_2_bottom <- r_2
R_2_bottom <- R_2_bottom  %>% filter(r_2 <= 0.2)



# put the linear regression into an organized data frame to export into excel file
coeff <- as.data.frame(do.call(rbind,storage))
cof<- as.data.frame(coeff$coefficients)
cof<- as.data.frame(t(cof))

# save the data into a data frame
write.xlsx(r_2, file = "calculations.xlsx", sheetName = "all squared")
write.xlsx(R_2, file = "calculations.xlsx", sheetName = "significant squared", append=TRUE)
write.xlsx(cof, file = "calculations.xlsx", sheetName = "all coeff", append=TRUE)




library('xlsx')
all_squared <- read_excel("calculations.xlsx", sheet = "all squared")
all_coeff <- read_excel("calculations.xlsx", sheet = "all coeff")[,c(1:3)]
sig_squared <- read_excel("calculations.xlsx", sheet = "significant squared")
info <- read_excel("20210728 combine_df2_first order kinetic.xlsx", sheet = "info")


# rename some columns
all_squared = all_squared %>% rename(UniProt = `...1`)
all_coeff = all_coeff %>% rename(UniProt = `...1`)
sig_squared  = sig_squared %>% rename(UniProt = `...1`)

r_cof <- merge(info, all_squared, by = "UniProt")
r_cof <- merge(r_cof, all_coeff, by = "UniProt")

r_sig <- merge(info, sig_squared,  by= "UniProt")
r_sig <- merge(r_sig, all_coeff, by = 'UniProt')


write.xlsx(r_cof, file = "cal_info.xlsx", sheetName = "all info")
write.xlsx(r_sig, file = "cal_info.xlsx", sheetName = "sig info", append=TRUE)


# adding df to the info 
df2 <-read_excel("df.xlsx")
head(df2)

# merge the file
df2 <- merge(r_cof,df2,by = "UniProt" )
# save it again
write.xlsx(df2 , file = "df.xlsx", sheetName = "df2_info")





hist(r_2$r_2, xlab = "R Squared value",
     main = "Histogram of all 521 Proteins")
lines(density(r_2$r_2), lwd = 2, col = 'red')

hist(abs(r_sig$`x (sploe)`), xlab = "K",
     main = "Histogram of K, 121 most sig")
lines(density(abs(r_sig$`x (sploe)`)), lwd = 2, col = 'red')




