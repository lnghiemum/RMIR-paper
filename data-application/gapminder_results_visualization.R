# Visualiation for the analysis gapminder data
rm(list = ls())
library(xtable)
library(ggpubr)
library(tidyverse)
library(grid)

load("data-application/gapminder_applications_withbinarycovariates.Rdata")
projection <- function(x, a = diag(nrow(x))){
  ldr::projection(x, a)
} 

fit_mcem_common <- data_application$RMIR_fit
countryinfo = data_application$country_info
df_merge <- data_application$data
estGamma = data_application$RMIR_pred_ranef
meanCSbasisRMIR = data_application$RMIR_fit$meanCSbasis
meanCS_df <- diag(projection(meanCSbasisRMIR))
separateMIRfit_unstr = data_application$SMIR

# Semi-orthogonal basis for the overall fixed effects central subspaces
print(xtable(t(as.matrix(qr.Q(qr(meanCSbasisRMIR)))), digits = 3), include.rownames = FALSE)

# Estimates for the covariance of random effects
print(xtable(data_application$RMIR_fit$estSigmastar, digits = 3), include.rownames = F)

# Figure 2 in the main paper
n = length(unique(data_application$data$geo))
country = unique(countryinfo$country)

# Create a matrix that stores the diagonal values for the projection matrix associated 
# with each cluster-specific central subspace
diagprojgamma <- 
  sapply(estGamma, function(j) diag(projection(j))) %>% t() %>%
  data.frame() %>%
  mutate(geo = country) %>%
  rename("Log Income" = X1, "Sex Ratio" = X2, "Infant mortality" = X3, "Emissions" = X4, 
         "Children per woman" = X5, "Inequality" = X6, "LDC" = X7, "West" = X8) %>%
  dplyr::select(geo, everything()) %>%
  pivot_longer(cols = -c("geo"), values_to = "coef", names_to = "Variable")

a <- diagprojgamma %>% group_by(Variable) %>% 
  filter(coef %in% sort(coef)[1:2] | coef %in% sort(coef, decreasing = TRUE)[1])


heatmap_RMIR <- ggplot(
  diagprojgamma %>% filter(geo %in% c(a$geo)) %>%
    mutate(geo = factor(geo, levels = c(unique(a$geo)))) %>%
    mutate(Variable = factor(Variable, 
                             levels = c("Log Income", "Sex Ratio", "Infant mortality", "Emissions", "Children per woman", "Inequality", "LDC", "West"))), 
  aes(x = Variable, y = geo, fill = coef) ) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = format(round(coef, 3), nsmall = 2)), color = "black", size = 8) +
  scale_fill_gradientn(colors = rev(terrain.colors(10)), guide = guide_colorbar(title = "")) +
  theme_classic(base_size = 20) +
  labs(y = "Central Subspaces")
ggsave(plot = heatmap_RMIR, filename = "Figure2.pdf", width = 20, height = 10)

# Figure S2 in the supplementary materials
df_plots1 <- data_application$data %>%
  mutate(meanDR1 = data_application$RMIR_meanDR[, 1]) %>%
  mutate(meanDR2 = data_application$RMIR_meanDR[, 2]) %>%
  dplyr::select(life_expectancy_female, west_and_rest, un_sdg_ldc, 
                meanDR1, meanDR2) %>%
  mutate(ldc = factor(un_sdg_ldc, levels = c("un_not_least_developed", "un_least_developed"), 
                      labels = c("No", "Yes")),
         west_and_rest = factor(west_and_rest, levels = c("rest", "west"), 
                                labels = c("No", "Yes"))) %>%
  pivot_longer(cols = meanDR1:meanDR2, names_to = "Direction", values_to = "value")

ggplot(df_plots1 %>% filter(life_expectancy_female > 30), 
       aes(x = value, y = life_expectancy_female)) +
  geom_point(aes(shape = west_and_rest, colour = ldc), size = 6) +
  theme_minimal(base_size = 20) +
  labs(x = "", y = "Female Life Expectancy") +
  facet_wrap(vars(Direction), scales = "free_x",
             labeller = as_labeller(c(`meanDR1` = "First Mean Direction",
                                      `meanDR2` = "Second Mean Direction"))) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 20)) +
  guides(colour = guide_legend(title = "Least Developed?"),
         shape = guide_legend(title = "Western?")) 
ggsave(filename = "FigureS2.pdf", width = 20, height = 10)

y = data_application$data[, 1:3]
y_with_country = left_join(y, countryinfo, by = "geo")

# Figure S3 in the supplementary materials
df_DR <- data.frame(data_application$data,
                    dr = Reduce("rbind", data_application$RMIR_DR), 
                    country = y_with_country$country) %>%
  mutate(ldc = factor(un_sdg_ldc, levels = c("un_not_least_developed", "un_least_developed"), 
                      labels = c("No", "Yes")),
         west_and_rest = factor(west_and_rest, levels = c("rest", "west"), 
                                labels = c("No", "Yes")))



country_list <- c("Laos", "Australia", "Canada", "Jamaica", "Namibia")

p1 <- ggplot(df_DR %>% filter(country %in% country_list), 
             aes(x = dr.1, y = life_expectancy_female)) +
  geom_point(aes(color = ldc, shape = west_and_rest), size = 6) + 
  geom_smooth(method = "loess", se = FALSE, size = 3) +
  theme_minimal(base_size = 24) +
  xlab("First direction") +
  ylab("Female life expectancy") +
  facet_wrap(vars(country), nrow = 1, scales = "free") +
  theme(strip.text = element_text(size = 20)) +
  guides(colour = guide_legend(title = "Least Developed?"),
         shape = guide_legend(title = "Western?")) 


p2 <- ggplot(df_DR %>% filter(country %in% country_list), aes(x = dr.2, y = life_expectancy_female)) +
  geom_point(aes(color = ldc, shape = west_and_rest), size = 6) + 
  geom_smooth(method = "loess", se = FALSE, size = 3) +
  theme_minimal(base_size =24) +
  xlab("Second direction") +
  ylab("Female life expectancy") +
  facet_wrap(vars(country), nrow = 1, scales = "free") + 
  theme(strip.text = element_text(size = 20)) +
  guides(colour = guide_legend(title = "Least Developed?"),
         shape = guide_legend(title = "Western?")) 

figure <- ggpubr::ggarrange(p1 + rremove("ylab"), 
                            p2 + rremove("ylab"), ncol = 1, nrow = 2, common.legend = TRUE, 
                            legend = "bottom")

annotate_figure(figure, left = textGrob("Female Life Expectancy", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))

ggsave(filename = "FigureS3.pdf", width = 20, height = 10)



