# Load required libraries
library(dplyr)
library(did)
library(ggplot2)

# Set seed for reproducibility
set.seed(12345)

# PARAMETERS: Number of individuals and time periods˙
N <- 500    # Number of individuals
T_periods <- 25  # Number of time periods

# Create individual-level data
ids <- 1:N
# Randomly assign treatment status (50% treated)
treated <- rbinom(N, 1, 0.5)
# For treated individuals, assign a treatment start time (between period 3 and 7)
treatment_time <- ifelse(treated == 1, sample(3:7, sum(treated == 1), replace = TRUE), Inf)

df_ind <- data.frame(id = ids,
                     treated = treated,
                     treatment_time = treatment_time)

# Expand to panel (long format)
df <- expand.grid(id = ids, time = 1:T_periods) %>%
  left_join(df_ind, by = "id") %>%
  arrange(id, time)

# Individual and time effects
alpha_i <- rnorm(N, mean = 0, sd = 1)      # individual-specific intercepts
gamma_t <- rnorm(T_periods, mean = 0, sd = 0.5)  # time-specific shocks

df <- df %>%
  mutate(alpha = rep(alpha_i, each = T_periods),
         gamma = rep(gamma_t, times = N))

# SIMULATING THE LATENT STOCK (S) AND FLOW (Y)
# Here we assume:
# - Every period, the stock S grows by a base rate.
# - After treatment onset, S grows by an additional small amount (beta) each period.
# - The flow outcome Y is defined as the period-to-period change in S.
base_growth <- 1
beta <- 0.1  # very small per-period treatment effect

# Initialize latent stock S matrix
S <- matrix(NA, nrow = N, ncol = T_periods)
for(i in 1:N){
  # Starting value incorporates individual and time effects + noise
  S[i,1] <- alpha_i[i] + gamma_t[1] + rnorm(1, 0, 1)
  for(t in 2:T_periods){
    treat_effect <- 0
    if(t >= df_ind$treatment_time[i]) {
      treat_effect <- beta
    }
    # Stock evolves by adding base growth, any treatment effect, and noise
    S[i,t] <- S[i, t-1] + base_growth + treat_effect + rnorm(1, 0, 1)
  }
}

# Convert S matrix to long format and attach to data frame
S_long <- as.vector(t(S))
df$S <- S_long

# Define the flow outcome Y as the first difference of S.
df <- df %>%
  group_by(id) %>%
  arrange(time) %>%
  mutate(Y = ifelse(time == 1, S, S - lag(S))) %>%
  ungroup()

# In this setting, Y (the flow) has only a tiny treatment effect (beta), likely lost in the noise.
# However, the latent stock S, which we treat as the cumulative outcome, accumulates this effect.
df <- df %>%
  group_by(id) %>%
  arrange(time) %>%
  mutate(Y_cum = S) %>%
  ungroup()

# For the did package, define 'G' as the treatment period (or 0 for never-treated).
df <- df %>% mutate(G = ifelse(is.finite(treatment_time), treatment_time, 0))

# ---------------------------------------------------------
# DID ESTIMATION USING CALLAWAY & SANT'ANNA APPROACH
# ---------------------------------------------------------

# Estimate treatment effects on the flow outcome (Y)
att_Y <- att_gt(yname = "Y",
                tname = "time",
                idname = "id",
                gname = "G",
                data = df,
                control_group = "nevertreated")
agg_Y <- aggte(att_Y, type = "dynamic")
cat("Summary for flow outcome (Y):\n")
print(summary(agg_Y))

# Estimate treatment effects on the cumulative outcome (Y_cum)
att_Ycum <- att_gt(yname = "Y_cum",
                   tname = "time",
                   idname = "id",
                   gname = "G",
                   data = df,
                   control_group = "nevertreated")
agg_Ycum <- aggte(att_Ycum, type = "dynamic")
cat("\nSummary for cumulative outcome (Y_cum):\n")
print(summary(agg_Ycum))

# Optional: Plot event study results for visual inspection
event_study_Y <- ggdid(agg_Y)
event_study_Ycum <- ggdid(agg_Ycum)
print(event_study_Y)
print(event_study_Ycum)

# retrieve Beta
df_plot <- data.frame(
  time = agg_Ycum$egt,         # Event time
  att = agg_Ycum$att.egt,      # Treatment effect estimates
  se = agg_Ycum$se.egt         # Standard errors
)

# Only select treatment periods
df_plot <- df_plot[df_plot$time>0,]
# For each time period, divide by time
df_plot$att_adj <- df_plot$att/df_plot$time

# Make plot - Beta is retrieved
ggplot(df_plot, aes(x = time, y = att_adj)) +
  geom_hline(yintercept = 0.1, col='red', linetype='dashed', size=1.1) +
  geom_hline(yintercept = 0, col='grey') +
  geom_point() +
  geom_errorbar(aes(ymin = att_adj - 1.96 * se, ymax = att_adj + 1.96 * se), width = 0.2, alpha=0.4) +
  theme_bw() +
  labs(title = "Beta retrieval from Cumulative Outcomes",
       x = "Time Since Treatment", y = "ATT (Treatment Effect)") 

# In this case, the slope of regressin att change over time also gives the treatment effect:
time <-  df_plot$time
atts <- df_plot$att
ols <- lm(atts ~ time)
# Extract the slope:
slope <- unname(ols$coefficients['time'])
# Ratio of slope over Beta is approximately 1:
slope/beta

# Create a data frame for the lines:
line_data <- data.frame(
  intercept = c(0.1, 0.1),
  slope = c(slope, beta),
  label = c("OLS Slope", "True Beta") # Legend labels
)

# Plot:
ggplot(df_plot, aes(x = time, y = att)) +
  geom_abline(data = line_data, aes(intercept = intercept, slope = slope, color = label), 
              linetype = 'dashed', size = 1.1, alpha = 0.5) +
  geom_hline(yintercept = 0, col = 'grey') +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("OLS Slope" = "red", "True Beta" = "blue")) +  # Define colors in legend
  labs(title = "Regression of ATT over time, OLS slope vs True Beta",
       x = "Time Since Treatment", y = "ATT (Treatment Effect)", color = "Effects:") +  # Legend title
  theme(legend.position = "bottom")  # Adjust legend position if needed












