library(mosaic)

Epidemicmodel <- function(Population, Spreadrate, persymptoms, pertreated, initialdisease, numberofdays) {
  df <- data.frame(Day=integer(),
                   Healthy=integer(),
                   Symptoms=integer(),
                   Diseased=integer(),
                   Catch_Disease=integer(),
                   Treat = integer(),
                   cured = integer())
  df[1,] <- c(0, Population - initialdisease, persymptoms * (Population - initialdisease), round(initialdisease), 
             round(initialdisease * Spreadrate * (persymptoms * (Population - initialdisease))), round(pertreated * initialdisease), 0)
  for(i in 1:numberofdays) {
          df[i+1, ]$Day <- i
          df[i+1, ]$Healthy <- df[i, ]$Healthy - df[i, ]$Catch_Disease
          df[i+1, ]$Symptoms <- round(persymptoms * df[i+1, ]$Healthy)
          df[i+1, ]$Diseased <- df[i, ]$Diseased + df[i, ]$Catch_Disease - df[i, ]$Treat
          df[i+1, ]$Catch_Disease <- round(df[i+1,]$Diseased * Spreadrate * (persymptoms *(df[i + 1,]$Healthy)))
          df[i+1, ]$Treat <- round(pertreated * df[i+1,]$Diseased)
          df[i+1, ]$cured <- df[i,]$cured + df[i, ]$Treat
  
  }
return df }




Final <- Epidemicmodel(500, 0.01, 0.2, 0.7, 5, 50)

library(ggplot2)
ggplot(data  = Final, aes(x = Day)) +
  geom_line(aes(y = Healthy), color = "Green") + 
  geom_line(aes(y = Diseased), color = "Red") + 
  geom_line(aes(y = cured), color = "Blue")
         

