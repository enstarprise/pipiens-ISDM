
### --- POPULATION PROCESS: ---
#> Priors on coefficients:
#> These represent the additive effect on the log-scale for a 
#> one-unit change in a predictor. 
#> This is where regularization is most important.

#> beta_x ~ Normal(0, sigma = 0.5 to 1)
#> A Normal(0, 0.5) prior implies: 95% of the prior mass is between [-1, 1].
#> strongly regularized towards smaller, more plausible effects 
#> but can learn larger ones if the data is strong enough
hist(rnorm(100, 0, 0.5))
exp(-1);exp(1) # multiplicative effect of 0.3678 and 2.718

hist(rnorm(100, 0, 5))
exp(-3);exp(4) # multiplicative effect of 0.0498 and 54.598

hist(rnorm(100, 0, 5))
exp(-15);exp(15) # multiplicative effect of 3.059x10^-7 and 3.2M

hist(rnorm(100, 0, 10))


#> As descirbed by Manica et al: Multivariate Normal distribution
#> (MVNorm(0, 1000*I) was  chosen as minimally informative prior for the
#> parameters of the population process (β)

library(mnormt)
hist(rmnorm(100, 0, 1000))

library(ggplot2)
# prior samples
set.seed(123)
n_samples <- 10000
beta_samples <- rnorm(n_samples, 0, 1)
odds_ratio_samples <- exp(beta_samples)

plot_data <- data.frame(
  beta = beta_samples,
  odds_ratio = odds_ratio_samples
)

# plot 1: on coefficient scale (log-odds)
p1 <- ggplot(plot_data, aes(x = beta)) +
  geom_density(fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(-10, 10), linetype = "dashed", color = "red") +
  labs(title = "Prior: Coefficient (Log-Odds Scale)",
       subtitle = "Normal(0, 1)",
       x = "β coefficient (log-odds)",
       y = "Density") +
  theme_minimal()

# plot 2: on odds ratio scale (multiplicative)
p2 <- ggplot(plot_data, aes(x = odds_ratio)) +
  geom_density(fill = "forestgreen", alpha = 0.7) +
  geom_vline(xintercept = c(4.539993e-05, 22026.47), linetype = "dashed", color = "red") +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100, 500, 10000),
                labels = c("0.01x", "0.1x", "1x", "10x","100x", "500x", "10000x")) +
  labs(title = "Prior: Effect on Abundance (Fold Change)",
       subtitle = "exp(β) ~ Log-Normal(0, 1)",
       x = "Fold Change in Abundance",
       y = "Density") +
  theme_minimal()

library(patchwork)
p1 / p2 # combine plots


##### SURVEY DETECTION PROCESS
# What your prior expects vs what literature shows
library(ggplot2)

# thee prior expectation
prior_traps <- rnorm(10000, logit(0.10), 1)  # trap effects around mean
prior_probs <- plogis(prior_traps)    # must ! convert to probabilities

# literature values
literature <- data.frame(
  trap = c("Lower end prob.", "Higher end prob."),
  logit = c(-0.33, 1.52),
  prob = c(0.42, 0.82)
)

ggplot() +
  geom_density(aes(x = prior_probs), fill = "blue", alpha = 0.5) +
  # geom_vline(data = literature, aes(xintercept = prob, color = trap), 
  #            size = 2, linetype = "dashed") +
  # scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Prior Expectation vs Literature Reality",
       x = "Detection Probability",
       y = "Density") +
  theme_minimal()

############ SPATIAL RANDOM EFFECT PRIORS
# data exploration
library(invgamma)

gamma <- rgamma(1000,4, 1)
gamma <- as.data.frame(gamma)
mean(gamma$gamma)
median(gamma$gamma)

library(invgamma)
library(ggplot2)

inv_gamma(2,5) alpha_gp
invgamma <- rinvgamma(1000, 2, 5)
invgamma <- as.data.frame(invgamma)
mean(invgamma$invgamma) # 5
median(invgamma$invgamma) # 3



## LENGTH SCALE
data <- rgamma(1000, 4, 1)
df <- data.frame(values = data)

ggplot(df, aes(x = values)) +
  geom_histogram(fill = "skyblue",
                 color = "grey99",
                 bins = 30) +
  geom_vline(aes(xintercept = mean(values),
                 color = "Mean"),
             linetype = "solid",
             size = 1) +
  geom_vline(aes(xintercept = median(values),
                 color = "Median"),
             linetype = "dashed",
             size = 1) +
  scale_color_manual(values = c("Mean" = "red",
                                "Median" = "darkgreen")) +
  labs(title = "Gamma Distribution (shape=4, rate=1)",
       x = "Values",
       y = "Count") +
  theme_minimal()



## MARGINAL SD
data <- rinvgamma(1000, 2, 5)
df <- data.frame(values = data)

ggplot(df, aes(x = values)) +
  geom_histogram(fill = "olivedrab",
                 color = "grey99",
                 bins = 30) +
  geom_vline(aes(xintercept = mean(values),
                 color = "Mean"),
             linetype = "solid",
             size = 1) +
  geom_vline(aes(xintercept = median(values),
                 color = "Median"),
             linetype = "dashed",
             size = 1) +
  scale_color_manual(values = c("Mean" = "red",
                                "Median" = "black")) +
  labs(title = "Inverse Gamma Distribution (shape=2, rate=5)",
       x = "Values",
       y = "Count") +
  theme_minimal()
