# Simple Mean-Field Ising Model 
library(ggplot2)
library(dplyr)
library(patchwork)

# M = tanh((J*M + H) / (k*T))
# Let J=1, k=1. Thus: M = tanh((M + H)/T)
# Critical Temperature is T_c = 1.0

solve_M <- function(T, H, tol=1e-6, max_iter=1000) {
  M <- 1.0 
  if(T == 0) return(1.0)
  for(i in 1:max_iter) {
    M_new <- tanh((M + H)/T)
    if(abs(M_new - M) < tol) return(M_new)
    M <- M_new
  }
  return(M)
}

T_c <- 1.0
T_seq <- seq(0.1, 2.0, by=0.01)
H_seq <- c(0, 0.05, 0.2)
df <- expand.grid(T = T_seq, H = H_seq)
df$M <- mapply(solve_M, df$T, df$H)

# Numerical derivative for Susceptibility Chi = dM/dH
delta_H <- 1e-4
df$M_plus <- mapply(solve_M, df$T, df$H + delta_H)
df$Chi <- (df$M_plus - df$M) / delta_H

df$Field <- as.factor(df$H)
df <- df %>% mutate(T_dist = T - T_c, abs_T_dist = abs(T_dist))

# Plot 1: Magnetization vs T - T_c
p1 <- ggplot(df, aes(x=T_dist, y=M, color=Field)) +
  geom_line(size=1) +
  geom_vline(xintercept=0, color="red", linetype="dashed") +
  theme_minimal(base_size=11) +
  labs(title="Order Parameter", x="T - T_c", y="Magnetization (M)")

# Plot 2: Susceptibility vs T - T_c (Linear)
p2_linear <- ggplot(df, aes(x=T_dist, y=Chi, color=Field)) +
  geom_line(size=1) +
  coord_cartesian(ylim=c(0, 25)) +
  geom_vline(xintercept=0, color="red", linetype="dashed") +
  theme_minimal(base_size=11) +
  labs(title="Susceptibility", x="T - T_c", y="Susceptibility (Chi)")

# Plot 3: Susceptibility vs |T - T_c| in log-log scale (only for H=0)
df_chi <- df %>% filter(H == 0, abs_T_dist > 0.005, T_dist > 0)
p3_log <- ggplot(df_chi, aes(x=abs_T_dist, y=Chi)) +
  geom_point(color="purple", size=1.5) +
  geom_smooth(method="lm", se=FALSE, color="blue", linewidth=0.8) +
  scale_x_log10() + scale_y_log10() +
  theme_minimal(base_size=11) +
  labs(title="Susceptibility Exp. (H=0)", x="log(|T - T_c|)", y="log(Chi)")

final_plot <- p1 + p2_linear + p3_log + plot_layout(ncol=3)
ggsave("ising_plots.pdf", final_plot, width=12, height=4)
