# 04 - Compose results. for SnesSurv.
# Marco Andrello
# 12/01/2020

rm(list=ls())

library(tidyverse)

# names_var <- c("surv","growth","fruiting","fruits","recruit")
names_var <- c("surv")


# Produce analogous to Table 2
load("Coefficients_SensSurv.RData")
nrepl <- dim(sd_rand[[1]])[2]

## Fixed factors
table_res_coef <- array(NA,c(8,length(names_var)*3))
colnames(table_res_coef) <- rep("",length(names_var)*3)
colnames(table_res_coef)[(((c(1:length(names_var))-1)*3)+1)] <- paste0(names_var,"_Median")
colnames(table_res_coef)[(((c(1:length(names_var))-1)*3)+2)] <- name.low <-paste0(names_var,"_low")
colnames(table_res_coef)[(c(1:length(names_var))*3)] <- name.high <- paste0(names_var,"_high")

# Calculate median and CI
for(i.var in 1 : length(names_var)) {
    table_res_coef[,(((i.var-1)*3)+1)] <- apply(coef_fixed[[i.var]],1,median,na.rm=T)
    table_res_coef[,(((i.var-1)*3)+2)] <- apply(coef_fixed[[i.var]],1,quantile,c(0.025),na.rm=T)
    table_res_coef[,(i.var*3)] <- apply(coef_fixed[[i.var]],1,quantile,c(0.975),na.rm=T)
}
# Write csv
table_res_coef %>% as_tibble %>% round(2) %>% add_column(Predictor=row.names(coef_fixed[[1]])) %>% write.csv("Results_coef_fixed_SensSurv.csv")


## Random factors


# We will rearrange the 8 different rows in this order:
# Site on Intercept,
# Year on Intercept,
# Year on slope of SoilVeg2,
# Correlation between Year effect on the Intercept and on the slope of SoilVeg2
# Plot on Intercept,
# Plot on slope of Plant size,
# Correlation between Plot effect on the Intercept and on the slope of Plant size
# Residual (fruits) or theta (growth, recruit size)

# Set up matching table between the rows of table_res_sd and those of stat_sd_rand

matching_table <- c(4,1,2,3,NA,NA,NA,NA)


# Set up table with 8 rows and 3 columns (mean, lower and upper bound of 95% CI) for each of the 5 variables
table_res_sd <- array(NA,c(8,length(names_var)*3))
colnames(table_res_sd) <- rep("",length(names_var)*3)
colnames(table_res_sd)[(((c(1:length(names_var))-1)*3)+1)] <- paste0(names_var,"_Median")
colnames(table_res_sd)[(((c(1:length(names_var))-1)*3)+2)] <- name.low <-paste0(names_var,"_low")
colnames(table_res_sd)[(c(1:length(names_var))*3)] <- name.high <- paste0(names_var,"_high")
rownames(table_res_sd) <- c("Site (Intercept)",
                            "Year (Intercept)",
                            "Year (SoilVeg2)",
                            "Year (Intercept * SoilVeg2)",
                            "Plot (Intercept)",
                            "Plot (Plant size)",
                            "Plot effect (Intercept * Plant size)",
                            "Residual")


# Calculate median and CI
stat_sd_rand <- list()
for(i.var in 1 : length(names_var)) {
    stat_sd_rand[[i.var]] <- array(NA,c(dim(sd_rand[[i.var]])[1],3))
    rownames(stat_sd_rand[[i.var]]) <- rownames(sd_rand[[i.var]])
    stat_sd_rand[[i.var]][,1] <- apply(sd_rand[[i.var]],1,median,na.rm=T)
    stat_sd_rand[[i.var]][,2] <- apply(sd_rand[[i.var]],1,quantile,c(0.025),na.rm=T)
    stat_sd_rand[[i.var]][,3] <- apply(sd_rand[[i.var]],1,quantile,c(0.975),na.rm=T)
}

# Push values of stat_sd_rand in the correct rows of table_res_sd for each variable
for(i.var in 1 : length(names_var)) {
    for (i.row in 1 : 8) {
        id.row <- matching_table[i.row]
        if (!is.na(id.row)) {
            table_res_sd[i.row,(((i.var-1)*3)+1)] <- stat_sd_rand[[i.var]][id.row,1]
            table_res_sd[i.row,(((i.var-1)*3)+2)] <- stat_sd_rand[[i.var]][id.row,2]
            table_res_sd[i.row,(i.var*3)] <- stat_sd_rand[[i.var]][id.row,3]
        }
    }
}

table_res_sd %>% as_tibble %>% round(2) %>% add_column(Predictor=row.names(table_res_sd)) %>% write.csv("Results_sd_rand_SensSurv.csv")


## R-squared
table_res_r <- array(NA,c(2,length(names_var)*3))
colnames(table_res_r) <- rep("",length(names_var)*3)
colnames(table_res_r)[(((c(1:length(names_var))-1)*3)+1)] <- paste0(names_var,"_Mean")
colnames(table_res_r)[(((c(1:length(names_var))-1)*3)+2)] <- name.low <-paste0(names_var,"_low")
colnames(table_res_r)[(c(1:length(names_var))*3)] <- name.high <- paste0(names_var,"_high")
# Calculate median and CI
for(i.var in 1 : length(names_var)) {
    table_res_r[,(((i.var-1)*3)+1)] <- apply(r_squared[[i.var]],1,median,na.rm=T)
    table_res_r[,(((i.var-1)*3)+2)] <- apply(r_squared[[i.var]],1,quantile,c(0.025),na.rm=T)
    table_res_r[,(i.var*3)] <- apply(r_squared[[i.var]],1,quantile,c(0.975),na.rm=T)
}
# Write csv
table_res_r %>% as_tibble %>% round(2) %>% add_column(r.squared=row.names(r_squared[[1]])) %>% write.csv("Results_r_squared_SensSurv.csv")

# Then compare in Excel



