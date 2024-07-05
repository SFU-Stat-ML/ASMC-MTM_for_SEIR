# Load required packages
library(dplyr)
library(ggplot2)
library(MASS)
library(splines)

set.seed(444)

# Read the data
real_data <- read.csv("C://SFU//Research//Comparison//AdaptiveASMC//Real Data//Real Data Analysis//covid_19_clean_complete.csv")

# Extracting Sri Lankan Daily Covid-19 Data
sri_lanka_data <- real_data %>%
  filter(Country.Region == "Sri Lanka") %>%
  dplyr::select(Date, Confirmed) %>%
  arrange(Date) %>%  
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) 

# Create a new column for daily new confirmed cases
sri_lanka_data <- sri_lanka_data %>%
  arrange(Date) %>%  # Ensure data is in date order
  mutate(Previous_Confirmed = lag(Confirmed, 1, default = 0),  
         New_Confirmed = Confirmed - Previous_Confirmed) %>%
  slice(-1)  

sri_lanka_data <- sri_lanka_data[-(1:37), ]

# Compute the 'times' column
sri_lanka_data <- sri_lanka_data %>%
  arrange(Date) %>%  
  mutate(times = as.integer(Date - min(Date)))  

# Calculate summary statistics
summary_stats <- sri_lanka_data %>%
  summarise(
    Total_Cases = sum(New_Confirmed),  
    Average_Cases = mean(New_Confirmed),  
    Median_Cases = median(New_Confirmed),  
    SD_Cases = sd(New_Confirmed),  
    Min_Cases = min(New_Confirmed),  
    Max_Cases = max(New_Confirmed)   
  )

# Print the summary statistics
print(summary_stats)

# Calculate sample mean and variance
mean_cases <- mean(sri_lanka_data$New_Confirmed)
variance_cases <- var(sri_lanka_data$New_Confirmed)

# Estimate negative binomial parameters
(size <- mean_cases^2 / (variance_cases - mean_cases))

daily_cases <- sri_lanka_data$New_Confirmed

data <- data.frame(day = 1:length(daily_cases), daily_cases = daily_cases)

# Fit the negative binomial model with day as a predictor
fit_nb <- glm.nb(daily_cases ~ day, data = data)

# Extract the estimated dispersion parameter (phi)
phi <- fit_nb$theta

# Calculate the mean (mu)
mu <- mean(daily_cases)

# Print the estimated parameters
cat("Estimated mean (mu):", mu, "\n")
cat("Estimated dispersion parameter (phi):", phi, "\n")

# Generate negative binomial fit data points using the estimated parameters
set.seed(444)  # For reproducibility
nbfit <- rnbinom(n = nrow(sri_lanka_data), size = phi, mu = mu)

sri_lanka_data <- sri_lanka_data %>%
  mutate(nbfit = nbfit)

# Plotting number of cases over time with the negative binomial fit using months
slnb <- ggplot(sri_lanka_data, aes(x = Date)) +
  geom_point(aes(y = New_Confirmed, color = "Actual Cases"), size = 2) +  
  geom_line(aes(y = nbfit, color = "NB Fit"), size = 1.5, linetype = "twodash") +  
  #scale_x_date(date_labels = "%b", date_breaks = "1 month") +  
  scale_color_manual(values = c("Actual Cases" = "#FF7256", "NB Fit" = "#458B00")) +  
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +  
  labs(title = "", 
       x = "Time (months)", 
       y = "Daily Confirmed Cases",
       color = "Legend") +  
  theme_minimal() +
  theme(
    legend.position = "bottom",  
    axis.title = element_text(size = 18),  
    axis.text = element_text(size = 16),  
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 18)  
  )

print(slnb)
ggsave("SL_NB_Fit.jpeg", plot = slnb, width = 12, height = 8, dpi = 300)
