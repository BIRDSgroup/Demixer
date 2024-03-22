library(dplyr)
library(tidyverse)
library(rstatix)
library(ggpubr)

# Step 1: Load the CSV file
data1 <- read.csv('dumpR1.csv') #Mix
data2 <- read.csv('dumpR2.csv') #NonMix

# Step 2: Extract the columns (assuming they are named 'Column1' and 'Column2')
#data_selected <- select(data, Mix, Non-Mix)

# Step 3: Run the Wilcoxon test
# Assuming these are independent samples
test_result <- wilcox.test(data1$Mix, data2$NonMix,alternative="greater")

# Display the test result
print(test_result$p.value)
