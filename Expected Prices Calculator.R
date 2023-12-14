# Install and load necessary packages

library(readxl)
library(openxlsx)
library(dplyr)


df <- read_excel("Datatest6.xlsx")

#Define the intervals and create a new column with group indices
intervals <- seq(0, max(df$St), by = 100)
df$Group <- cut(df$St, breaks = c(intervals, Inf), labels = seq_along(intervals), right = FALSE)

#Calculate the metrics
avg_zt <- df %>%
  group_by(Group) %>%
  summarize(Avg_Zt = mean(Zt, na.rm = TRUE),
            Probability = n() / nrow(df),
            Avg_St = mean(St, na.rm = TRUE),
            Adjusted_Probability = Probability * Avg_Zt,
            Adjusted_Price = Avg_St * Avg_Zt)

#Merge the average values back to the original data frame
df <- left_join(df, avg_zt, by = "Group")

expected_st <- df %>%
  group_by(Group) %>%
  summarize(Expected_St = sum(Probability * Avg_St))

df$Return <- (df$St - 487.16) / 487.16

df$Expected_Return <- (df$Adjusted_Price - 487.16) / 487.16

expected_st <- df %>%
  group_by(Group) %>%
  summarize(Expected_St = sum(Adjusted_Probability * Avg_St))

df <- left_join(df, expected_st, by = "Group")

write.xlsx(df, "R Codes", rowNames = FALSE)









