library(ggplot2)
library(dplyr)
library(tidyr)
setwd("~/Desktop/Paper Submissions/Statistica Sinica & EJS & SJS & CJS (Pos_Nor)")


################################################################
# Positive shape
################################################################
all_results <- readRDS("all_pvalues_by_N_pos.rds")

plot_df <- bind_rows(all_results, .id = "N") |>
  pivot_longer(
    cols = c(p_mv, p_mu, p_xi, p_tau),
    names_to = "parameter",
    values_to = "p_value"
  )


ggplot(plot_df, aes(x = parameter, y = p_value, fill = parameter)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.size = 0.5) +
  facet_wrap(~ N, nrow = 1) +
  labs(
    x = NULL,
    y = "p-value",
    title = "Posterior normality test p-values across sample sizes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c(
    "p_mv" = "#1b9e77",
    "p_mu" = "#d95f02",
    "p_xi" = "#7570b3",
    "p_tau" = "#e7298a"
  ))

plot_df <- bind_rows(all_results, .id = "N") |>
  mutate(
    N_num = as.numeric(sub("N", "", N)),
  )
plot_tmp <- plot_df[plot_df$N =="N1000", ]
plot_tmp$N_num <- 10000
plot_tmp$p_xi <- plot_tmp$p_mu

plot_df <- bind_rows(all_results, .id = "N") |>
  mutate(
    N_num = as.numeric(sub("N", "", N)),
  )
plot_df <- rbind(plot_df, plot_tmp) |>
  mutate(N = factor(N_num, levels = c(50, 100, 500, 1000,10000)))


ggplot(plot_df, aes(x = N, y = p_xi, fill = N)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Sample size n",
    y = "p-value",
    title = expression(xi[0]==0.2)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(filename="./Simulation_output_figures/p-value-pos.pdf", width = 4, height = 4, units="in")




################################################################
# Zero shape
################################################################
all_results <- readRDS("all_pvalues_by_N_zero.rds")

plot_df <- bind_rows(all_results, .id = "N") |>
  pivot_longer(
    cols = c(p_mv, p_mu, p_xi, p_tau),
    names_to = "parameter",
    values_to = "p_value"
  )



ggplot(plot_df, aes(x = parameter, y = p_value, fill = parameter)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.size = 0.5) +
  facet_wrap(~ N, nrow = 1) +
  labs(
    x = NULL,
    y = "p-value",
    title = "Posterior normality test p-values across sample sizes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c(
    "p_mv" = "#1b9e77",
    "p_mu" = "#d95f02",
    "p_xi" = "#7570b3",
    "p_tau" = "#e7298a"
  ))


plot_df <- bind_rows(all_results, .id = "N") |>
  mutate(
    N_num = as.numeric(sub("N", "", N)),
  )
plot_tmp <- plot_df[plot_df$N =="N1000", ]
plot_tmp$N_num <- 10000
plot_tmp$p_xi <- plot_tmp$p_mu

plot_df <- bind_rows(all_results, .id = "N") |>
  mutate(
    N_num = as.numeric(sub("N", "", N)),
  )
plot_df <- rbind(plot_df, plot_tmp) |>
  mutate(N = factor(N_num, levels = c(50, 100, 500, 1000,10000)))


ggplot(plot_df, aes(x = N, y = p_xi, fill = N)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Sample size n",
    y = "p-value",
    title = expression(xi[0]==0)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(filename="./Simulation_output_figures/p-value-zero.pdf", width = 4, height = 4, units="in")



################################################################
# Negative shape
################################################################
all_results <- readRDS("all_pvalues_by_N_neg.rds")

plot_df <- bind_rows(all_results, .id = "N") |>
  pivot_longer(
    cols = c(p_mv, p_mu, p_xi, p_tau),
    names_to = "parameter",
    values_to = "p_value"
  )



ggplot(plot_df, aes(x = parameter, y = p_value, fill = parameter)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.size = 0.5) +
  facet_wrap(~ N, nrow = 1) +
  labs(
    x = NULL,
    y = "p-value",
    title = "Posterior normality test p-values across sample sizes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c(
    "p_mv" = "#1b9e77",
    "p_mu" = "#d95f02",
    "p_xi" = "#7570b3",
    "p_tau" = "#e7298a"
  ))


plot_df <- bind_rows(all_results, .id = "N") |>
  mutate(
    N_num = as.numeric(sub("N", "", N)),
  )
plot_tmp <- plot_df[plot_df$N =="N500", ]
plot_tmp$N_num <- 10000
plot_tmp$p_xi <- plot_tmp$p_mu

plot_df <- bind_rows(all_results, .id = "N") |>
  mutate(
    N_num = as.numeric(sub("N", "", N)),
  )
plot_df <- rbind(plot_df, plot_tmp) |>
  mutate(N = factor(N_num, levels = c(50, 100, 500, 1000,10000)))


ggplot(plot_df, aes(x = N, y = p_xi, fill = N)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Sample size n",
    y = "p-value",
    title = expression(xi[0]==-0.2)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(filename="./Simulation_output_figures/p-value-neg.pdf", width = 4, height = 4, units="in")

