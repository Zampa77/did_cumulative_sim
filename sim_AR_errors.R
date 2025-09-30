# Load required libraries
library(dplyr)
library(did)
library(ggplot2)

set.seed(12345)

# PARAMETERS
N <- 500
T_periods <- 30
rho <- 0.3       # AR(1) coefficient for error correlation within units
sigma_u <- 1      # SD of innovation in the AR(1) process

# Individual‐level setup
ids <- 1:N
treated <- rbinom(N, 1, 0.5)
treatment_time <- ifelse(treated == 1, sample(3:7, sum(treated), replace = TRUE), Inf)
df_ind <- data.frame(id = ids, treated, treatment_time)

# Expand to panel
df <- expand.grid(id = ids, time = 1:T_periods) %>%
  left_join(df_ind, by = "id") %>%
  arrange(id, time)

# Fixed effects
alpha_i <- rnorm(N, 0, 1)
gamma_t <- rnorm(T_periods, 0, 0.5)
df <- df %>%
  mutate(alpha = rep(alpha_i, each = T_periods),
         gamma = rep(gamma_t, times = N))

# Simulate AR(1) errors within each unit
Eps <- matrix(NA, nrow = N, ncol = T_periods)
for(i in 1:N) {
  # Initialize at stationary distribution:
  Eps[i,1] <- rnorm(1, 0, sigma_u / sqrt(1 - rho^2))
  for(t in 2:T_periods) {
    Eps[i,t] <- rho * Eps[i,t-1] + rnorm(1, 0, sigma_u)
  }
}

# Latent stock S and flow Y
base_growth <- 1
beta <- 0.1

S <- matrix(NA, N, T_periods)
for(i in 1:N){
  S[i,1] <- alpha_i[i] + gamma_t[1] + Eps[i,1]
  for(t in 2:T_periods){
    treat_eff <- if(t >= df_ind$treatment_time[i]) beta else 0
    S[i,t] <- S[i,t-1] + base_growth + treat_eff + Eps[i,t]
  }
}

# Attach to df
df$S <- as.vector(t(S))
df <- df %>%
  group_by(id) %>%
  arrange(time) %>%
  mutate(Y     = ifelse(time == 1, S, S - lag(S)),
         Y_cum = S,
         G     = ifelse(is.finite(treatment_time), treatment_time, 0)) %>%
  ungroup()

# DID estimation
att_Y <- att_gt(yname = "Y",    tname = "time", idname = "id",
                gname = "G", data = df, control_group = "nevertreated")
agg_Y <- aggte(att_Y, type = "dynamic")
cat("Flow outcome (Y) with AR(1) errors:\n")
print(summary(agg_Y))

att_Ycum <- att_gt(yname = "Y_cum", tname = "time", idname = "id",
                   gname = "G", data = df, control_group = "nevertreated")
agg_Ycum <- aggte(att_Ycum, type = "dynamic")
cat("\nCumulative outcome (Y_cum) with AR(1) errors:\n")
print(summary(agg_Ycum))


# Optional: Plot event study results for visual inspection
event_study_Y <- ggdid(agg_Y)
event_study_Ycum <- ggdid(agg_Ycum)
print(event_study_Y)
print(event_study_Ycum)