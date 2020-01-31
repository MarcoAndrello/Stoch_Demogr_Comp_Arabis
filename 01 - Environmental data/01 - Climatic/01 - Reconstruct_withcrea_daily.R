# Reconstruct temperature records of the iButtons using data from CREA weather stations
# Pierre de Villemereuil and Marco Andrello
# 23/10/2019


rm(list=ls())

library(tidyverse)
library(lubridate)
library(hms)
library(magrittr)
library(broom)

## ************************** Preliminary stuff

## Getting CREA information
load("crea_data.Rdata")

## Set some useful variables
levpop      <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
listquad    <- list()
listquad[["BRU"]] <- c(1, 2)
listquad[["CHA"]] <- c(1, 2, 3)
listquad[["VIL"]] <- c(1, 2, 3)
listquad[["LAU"]] <- c(1, 2, 3)
listquad[["GAL"]] <- c(1, 2, 3)
listquad[["PIC"]] <- c(1, 2, 3, 4)

## Compute some characteristics
# Some functions to summarise the weather
amplitude       <- . %>% {max(.) - min(.)} 
mean.dtemp      <- . %>% mean(.)
max.dtemp       <- . %>% max(.)
min.dtemp       <- . %>% min(.)
freeze          <- . %>% {any(. < 0)} %>% as.integer()
cum.dpostemp    <- . %>% {sum(.[. > 0])} 
cum.dprec       <- . %>% sum(.)

# Rounding the hours
rounded_hour <- function(t) {
    h <- hour(t)
    h <- ifelse(minute(t) >= 30,
                h + 1,
                h)
    return(hms(hour = h))
    #     return(hms(str_glue("{h}:00:00")))
    }

# Reset the year for all dates
reset_year <- function(vec_dates) {
    year(vec_dates) <- year(today())
    return(vec_dates)
}

# Checking for negative temperatures
is.frost <- function(vec) {
    ifelse(length(which(vec < 0)) > 0, 1, 0)
}

# Add predicted values to a df
add_predict <- function(df, model, prefix) {
    # Predict values to the df
    pred <- predict(model,
                    newdata  = df,
                    se.fit = TRUE)
    # Append predicted values to the df
    df[[str_glue("{prefix}_predict")]] <- pred[["fit"]]
    df[[str_glue("{prefix}_se")]] <- pred[["se.fit"]]
    # Enforce a tibble, 'cause I prefer it
    df <- as_tibble(df)
    # Return the df
    return(df)
}

# Number of missing values
n_missing <- . %>% is.na() %>% which() %>% length()

# Brute-force compute number of freezing days
compute_nfreeze <- function(probs, freeze, type) {
    # If no predicted values, just return nb days
    if (all(type != "Predicted")) {
        out <- tibble(NbDays    = sum(freeze, na.rm = TRUE),
                      Prob      = 1)
        return(out)
    } else {
        # Trimming out directly observed values
        probs <- probs[type == "Predicted"]
        
        # Known number of freezing days
        known_nbdays <- sum(freeze[type != "Predicted"], na.rm = TRUE)
        
        # Number of days
        k <- length(probs)
        # Number of replicates
        nrep <- 1000
        # Vector to store nb of freezing days
        nbfreeze <- numeric(nrep)
        
        for (i in 1:nrep) {
            # Simulating freezing or not
            Y <- rbinom(k, size = 1, prob = probs)
            # Compute nb of freezing days
            nbfreeze[i] <- sum(Y, na.rm = TRUE)
        }
        nbfreeze <- known_nbdays + nbfreeze
        
        # Format the output
        out <- 
            nbfreeze %>% 
            table() %>% 
            as_tibble() %>%
            set_colnames(c("NbDays", "Prob")) %>%
            mutate(Prob = Prob / nrep) %>%
            filter(Prob > 0.05)
        
        return(out)
    }
}

## ************************** Starting the loop

allmodels   <- tibble()
reconstruct <- tibble()
for (pop in levpop) {
    for (quad in listquad[[pop]]) {
        print(str_glue("Working on {pop}-{quad}"))
        
        ## Vercors or Lautaret?
        # Get the correct stations and month for each massif
        if (pop %in% c("BRU","CHA","VIL")) {
            stations    <- c("207", "209")
            monthkeep   <- 7 #c(6,7)
        } else {
            stations    <- c("217", "218")
            monthkeep   <- 7 #c(6,7)
        }
        
        ## ---------------------- Cleaning and formatting the data
        
        ## Get the iButtons data
        # "clean" means "imputation" from other quadrats was already done
        # Biased temperature measurements (PIC2-2009, PIC2-2012, PIC4-2009) were set to NA
        data <-
            read_csv(file = str_glue("iButtons clean/données.horaires/With 2015/{pop}-Q{quad}.csv"))
        # Some formatting
        data <- 
            data %>%
            mutate(Jour     = dmy(Jour),
                   Heure    = rounded_hour(Heure))
        # Creating a new df to work with (useless now, but kept because I'm lazy)
        df <-
            data %>%
            mutate(Time     = ymd_hms(str_c(Jour, " ", Heure)),
                   Year     = year(Jour)) %>%
            rename(Day  = Jour)
        # Finally, transfer time to UTC
        df[["Time_UTC"]] <- df[["Time"]]
        hour(df[["Time_UTC"]]) <- hour(df[["Time_UTC"]]) - 2
        
        ## Getting the month of focus
        df <-
            df %>%
            filter(month(Time) %in% monthkeep, month(Time_UTC) %in% monthkeep)
        
        ## Compute daily parameters
        daily <-
            df %>%
            group_by(Day) %>%
            summarise(mtemp     = mean.dtemp(Temp),
                      maxtemp   = max.dtemp(Temp),
                      mintemp   = min.dtemp(Temp),
                      amp       = amplitude(Temp),
                      freeze    = freeze(Temp)) %>%
            mutate(Year = year(Day),
                   Snow = as.integer(amp < 1.5 & mtemp < 3))
        
        ## ------------------------------ Matching with the weather information
        
        ## Formatting the weather into columns of temperatures
        # Selecting the data
        weather <- 
            crea %>% 
            filter(Station %in% stations) %>%
            select(Station, Time, contains("Temp")) %>%
            mutate(Day = str_glue("{year(Time)}-{month(Time)}-{day(Time)}") %>%
                         ymd()) %>%
            select(Day, Time, everything()) %>%
            filter(month(Day) %in% monthkeep)
        
        # Daily weather
        daily_weather <-
            weather %>%
            group_by(Day,Station) %>%
            # Getting the average, min and max temp for each day
            summarise(mtemp_m5cm    = mean.dtemp(Temp_m5cm),
                      mtemp_Sol     = mean.dtemp(Temp_Sol),
                      mtemp_30cm    = mean.dtemp(Temp_30cm),
                      mtemp_2m      = mean.dtemp(Temp_2m),
                      max_m5cm      = max.dtemp(Temp_m5cm),
                      max_Sol       = max.dtemp(Temp_Sol),
                      max_30cm      = max.dtemp(Temp_30cm),
                      max_2m        = max.dtemp(Temp_2m),
                      min_m5cm      = min.dtemp(Temp_m5cm),
                      min_Sol       = min.dtemp(Temp_Sol),
                      min_30cm      = min.dtemp(Temp_30cm),
                      min_2m        = min.dtemp(Temp_2m),
                      amp_m5cm      = amplitude(Temp_m5cm),
                      amp_Sol       = amplitude(Temp_Sol),
                      amp_30cm      = amplitude(Temp_30cm),
                      amp_2m        = amplitude(Temp_2m),
                      freeze_m5cm   = freeze(Temp_m5cm),
                      freeze_Sol    = freeze(Temp_Sol),
                      freeze_30cm   = freeze(Temp_30cm),
                      freeze_2m     = freeze(Temp_2m)) %>%
            # Adding year information
            mutate(Year = year(Day)) %>%
            ungroup()
    
        # Formatting the data in a spreaded df
        daily_weather <-
            left_join(daily_weather %>%
                          filter(Station == stations[1]) %>%
                          select(-Station),
                      daily_weather %>% 
                          filter(Station == stations[2]) %>%
                          select(-Station),
                      by = "Day",
                      suffix = c("_1", "_2")) %>%
            mutate(Year = year(Day))

        ## Matching with the data
        daily <-
            daily %>%
            full_join(daily_weather, by = c("Year", "Day")) %>%
            select(Year, Day, everything()) %>%
            arrange(Year, Day)
        
        ## Checking issues with snow
        issue_snow <-
            daily[["Snow"]] %>%
            replace_na("") %>%
            str_c(collapse = "") %>%
            str_detect("[1]{2,}")
        if (issue_snow) {
            warning(str_glue("Issue with snow at {pop}-{quad}."))
        }
        
        ## Modelling the temperatures
        # Creating a function to model mean, min or max temp (see "var")
        modelling <-
            . %>%
            # Removing days of snow when needed
            when(
                issue_snow  ~ filter(., Snow != 1),
                !issue_snow ~ .
            ) %>%
            # Selecting only some variables
            select(Year, 
                   contains("mtemp"),
                   contains("max"),
                   contains("min")) %>%
            # Removing underground measures
            select(-contains("m5cm")) %>%
            # Factorising year (used for the model, but removed as it impedes the reconstruction)
            mutate(Year = factor(Year)) %>%
            # Removing missing data
            drop_na()  %>%
            mutate(var = !!var) %>%
            # Perform the lm() and do a variable selection on-the-fly
            do(
                model = step(lm(var ~ mtemp_Sol_1 + mtemp_30cm_1 + mtemp_2m_1 +
                                      mtemp_Sol_2 + mtemp_30cm_2 + mtemp_2m_2 +
                                      max_Sol_1 + max_30cm_1 + max_2m_1 +
                                      max_Sol_2 + max_30cm_2 + max_2m_2 +
                                      min_Sol_1 + min_30cm_1 + min_2m_1 +
                                      min_Sol_2 + min_30cm_2 + min_2m_2,
                                data = .),
                             direction = "both", trace = 0)
            )
        
        # Running the models
        var <- quo(mtemp)
        models_mtemp <- modelling(daily)
        var <- quo(maxtemp)
        models_max <- modelling(daily)
        var <- quo(mintemp)
        models_min <- modelling(daily)
        
        # Getting summary statistics and combining the results
        add_glance <-
            . %>% 
            {bind_cols(., map_dfr(.[["model"]], glance))}
        models_mtemp <- add_glance(models_mtemp)
        models_max   <- add_glance(models_max)
        models_min   <- add_glance(models_min)
        allmodels <- 
            allmodels %>%
            bind_rows(add_column(models_mtemp,
                                 Var = "MeanTemp",
                                 Pop = pop, 
                                 Quad = as.integer(quad), 
                                 .before = 1),
                      add_column(models_max, 
                                 Var = "MaxTemp",
                                 Pop = pop, 
                                 Quad = as.integer(quad), 
                                 .before = 1),
                      add_column(models_min, 
                                 Var = "MinTemp",
                                 Pop = pop, 
                                 Quad = as.integer(quad), 
                                 .before = 1))
        
        ## Now reconstructing the data
        # A template with all days and year to work on
        # For reconstructing June in Lautaret and July in Galibier:
            # nbday <- ifelse(monthkeep == 6, 30, 31)
            # template <- tibble(Year = rep(2009:2014, each = nbday),
            #                    Day  = str_glue("{Year}-{monthkeep}-{rep(1:nbday, 6)}") %>%
            #                           ymd())
        # For reconstructing June and July in all sites:
            # nbday <- 30 + 31
            # template <- tibble(Year  = rep(2009:2014, each = nbday),
            #                    Month = rep(c(rep(6, 30), rep(7, 31)), 6),
            #                    Day  = str_glue("{Year}-{Month}-{rep(c(1:30,1:31), 6)}") %>%
            #                           ymd())
        # For reconstructing July in all sites:
        nbday <- 31
        template <- tibble(Year = rep(2009:2014, each = nbday),
                           Day  = str_glue("{Year}-{monthkeep}-{rep(1:nbday, 6)}") %>%
                                  ymd())
        template <- tibble(Year = rep(2008:2014, each = nbday),
                           Day  = str_glue("{Year}-{monthkeep}-{rep(1:nbday, 7)}") %>%
                               ymd())
                           
        # Reconstructing the data
        reconstruct <-
            daily %>%
            # Augment the data according to the template
            right_join(template, by = c("Year", "Day")) %>% 
            # Adding predictions for mean, max and min temp
            add_predict(models_mtemp[["model"]][[1]], "mtemp") %>%
            add_predict(models_max[["model"]][[1]], "maxtemp") %>%
            add_predict(models_min[["model"]][[1]], "mintemp") %>%
            # Selecting and renaming some variables
            select(Year, Day,
                   MeanTemp_Orig    = mtemp,
                   MeanTemp_Pred    = mtemp_predict,
                   MeanTemp_PredSE  = mtemp_se,
                   MinTemp_Orig     = mintemp,
                   MinTemp_Pred     = mintemp_predict,
                   MinTemp_PredSE   = mintemp_se,
                   MaxTemp_Orig     = maxtemp,
                   MaxTemp_Pred     = maxtemp_predict,
                   MaxTemp_PredSE   = maxtemp_se,
                   AmpTemp_Orig     = amp) %>%
            # Computing amplitude and prob. of freeze from there
            mutate(AmpTemp_Pred     = MaxTemp_Pred - MinTemp_Pred,
                   AmpTemp_PredSE   = sqrt(MaxTemp_PredSE^2 + MinTemp_PredSE^2),
                   Freeze_Orig      = as.integer(MinTemp_Orig < 0),
                   Freeze_Pred      = as.integer(MinTemp_Pred < 0),
                   Prob_Freeze      = pnorm(0, MinTemp_Pred, MinTemp_PredSE)) %>%
            # Adding some informations
            mutate(Pop  = pop,
                   Quad = as.integer(quad),
                   Year = as.integer(Year)) %>%
            # Sorting the variables
            select(Pop, Quad, Year, Day, everything()) %>%
            # Append to reconstruct
            {bind_rows(reconstruct, .)}
    }
}

## Saving the results
save(allmodels, file = "models_daily.Rdata")
save(reconstruct, file = "reconstruct.Rdata")

## ************************** Outputting the models

## Loading the models
load("models_daily.Rdata")

## Formatting the output
output <-
    allmodels %>%
    # Handle the formulae
    mutate(formula  = map(model, formula) %>%
                      map(deparse) %>%
                      map_chr(~str_c(., collapse = "")) %>%
                      str_replace_all("[ ]+", " ") %>%
                      str_replace_all("var [ ]*", "")) %>%
    select(Pop,
           Quad,
           Var,
           formula,
           R2 = r.squared,
           pval = p.value,
           df,
           df.res = df.residual) %>%
    arrange(Pop, Quad, Var)

## Checking that all models are significant
test <- any(output[["pval"]] > 0.05)
if (test) { stop("At least one model is not significant!") }

## Saving the output
write_csv(output, path = "Temperature_models_daily.csv")


## ************************** Formatting the reconstructed data

## Loading the reconstructed data
load("reconstruct.Rdata")

## Creating the "consensus" column
reconstruct <-
    reconstruct %>%
    # Consensus for temperature
    add_column(MeanTemp = ifelse(is.na(.[["MeanTemp_Orig"]]),
                                 .[["MeanTemp_Pred"]],
                                 .[["MeanTemp_Orig"]]),
               MeanTemp_Type = case_when(
                   is.na(.[["MeanTemp_Orig"]]) & is.na(.[["MeanTemp_Pred"]]) ~ "Missing",
                   is.na(.[["MeanTemp_Orig"]])                               ~ "Predicted",
                   TRUE                                                      ~ "Original"
               ),
               .before  = "MeanTemp_Orig") %>%
    add_column(MinTemp  = ifelse(is.na(.[["MinTemp_Orig"]]),
                                 .[["MinTemp_Pred"]],
                                 .[["MinTemp_Orig"]]),
               MinTemp_Type = case_when(
                   is.na(.[["MinTemp_Orig"]]) & is.na(.[["MinTemp_Pred"]]) ~ "Missing",
                   is.na(.[["MinTemp_Orig"]])                              ~ "Predicted",
                   TRUE                                                    ~ "Original"
               ),
               .before  = "MinTemp_Orig") %>%
    add_column(MaxTemp  = ifelse(is.na(.[["MaxTemp_Orig"]]),
                                 .[["MaxTemp_Pred"]],
                                 .[["MaxTemp_Orig"]]),
               MaxTemp_Type = case_when(
                   is.na(.[["MaxTemp_Orig"]]) & is.na(.[["MaxTemp_Pred"]]) ~ "Missing",
                   is.na(.[["MaxTemp_Orig"]])                              ~ "Predicted",
                   TRUE                                                    ~ "Original"
               ),
               .before  = "MaxTemp_Orig") %>%
    add_column(AmpTemp  = ifelse(is.na(.[["AmpTemp_Orig"]]),
                                 .[["AmpTemp_Pred"]],
                                 .[["AmpTemp_Orig"]]),
               AmpTemp_Type = case_when(
                   is.na(.[["AmpTemp_Orig"]]) & is.na(.[["AmpTemp_Pred"]]) ~ "Missing",
                   is.na(.[["AmpTemp_Orig"]])                              ~ "Predicted",
                   TRUE                                                    ~ "Original"
               ),
               .before  = "AmpTemp_Orig") %>%
    add_column(Freeze  = ifelse(is.na(.[["Freeze_Orig"]]),
                                .[["Freeze_Pred"]],
                                .[["Freeze_Orig"]]),
               Freeze_Type = case_when(
                   is.na(.[["Freeze_Orig"]]) & is.na(.[["Freeze_Pred"]]) ~ "Missing",
                   is.na(.[["Freeze_Orig"]])                             ~ "Predicted",
                   TRUE                                                  ~ "Original"
               ),
               .before  = "Freeze_Orig")

## Checking the number of missing values
missing <-
    reconstruct %>%
    group_by(Pop, Quad, Year) %>%
    summarise(n     = n_missing(MeanTemp),
              prop  = n / n())
test <- any(missing[["n"]] >= 5)
if (test) { stop("Many missing values for a year") }

## Saving the output
save(reconstruct,file="reconstruct_consensus.RData")
output <-
    reconstruct %>%
    mutate_if(is.double,
              funs(format(., digits = 2)))
write_csv(output, path = "Reconstruct_daily.csv")


## ************************** Generate the monthly estimates

## Compute the monthly estimates
monthly <-
    reconstruct %>%
    group_by(Pop, Quad, Year) %>%
    summarise(AvgDTemp  = mean(MeanTemp, na.rm = TRUE),
              AvgDTempSE= ifelse(MeanTemp_Type != "Predicted",
                                 0,
                                 MeanTemp_PredSE) %>%
                  raise_to_power(2) %>%
                  sum() %>%
                  divide_by(n()) %>%
                  sqrt(),
              AvgMinTemp  = mean(MinTemp, na.rm = TRUE),
              AvgMinTempSE= ifelse(MinTemp_Type != "Predicted",
                                 0,
                                 MinTemp_PredSE) %>%
                  raise_to_power(2) %>%
                  sum() %>%
                  divide_by(n()) %>%
                  sqrt(),
              AvgMaxTemp  = mean(MaxTemp, na.rm = TRUE),
              AvgMaxTempSE= ifelse(MaxTemp_Type != "Predicted",
                                   0,
                                   MaxTemp_PredSE) %>%
                  raise_to_power(2) %>%
                  sum() %>%
                  divide_by(n()) %>%
                  sqrt(),
              AvgAmp    = mean(AmpTemp, na.rm = TRUE),
              AvgAmpSE  = ifelse(AmpTemp_Type != "Predicted",
                                 0,
                                 AmpTemp_PredSE) %>%
                  raise_to_power(2) %>%
                  sum() %>%
                  divide_by(n()) %>%
                  sqrt(),
              NbFreeze  = sum(Freeze, na.rm = TRUE),
              NbFreeze2 = compute_nfreeze(Prob_Freeze, Freeze_Orig, Freeze_Type) %>%
                  list())


save(monthly, file = "reconstruct_monthly.Rdata")

out <- 
    monthly %>%
    mutate(NbFreeze = ifelse(NbFreeze2 %>% map_int(nrow) %>% {. > 1},
                             NA,
                             NbFreeze)) %>%
    select(-NbFreeze2)
write_csv(out, path = "Reconstruct_monthly.csv")

## Output freezing day probs for all sites relevant
which_uncertain_freeze <- monthly[["NbFreeze2"]] %>% map_int(nrow) %>% {which(. > 1)} %T>% print()
for (i in which_uncertain_freeze) {
    out <- monthly[["NbFreeze2"]][[i]]
    pop <- monthly[["Pop"]][[i]]
    quad <- monthly[["Quad"]][[i]]
    year <- monthly[["Year"]][[i]]
    write_csv(out, path = str_glue("NbFreeze_{pop}_{quad}_{year}.csv"))
}

