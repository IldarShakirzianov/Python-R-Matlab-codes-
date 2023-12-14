# Install and load the required package
library(readxl)

# Replace 'your_file.xlsx' with the actual file path
file_path <- "dataheston.xlsx"

# Read data from Excel file
data <- read_excel(file_path)

# Extract the first column (St)
prices <- data$St

# Define breaks with a step of 100
breaks <- seq(min(prices), max(prices) + 50, by = 50)

# Create the histogram
hist(prices, breaks = breaks, col = "skyblue", border = "black", 
     xlab = "S_t", ylab = "Frequency", main = "Stock Price Histogram")