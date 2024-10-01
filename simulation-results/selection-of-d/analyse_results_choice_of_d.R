## Analysis of results
library(gtools)
library(tidyverse)
Settings <- expand.grid(n = c(100, 500, 750, 1000),
                         d = 1:2,
                         Sigmatype = 1:3)


dir.files <- gtools::mixedsort(list.files("./results_choice_d_selection/", full.names = T))
numCases <- 24

dresults_overall <- sapply(1:numCases, function(ii){
  bindex <- ii
  eindex <- ii
  dres <- lapply(dir.files[bindex:eindex], function(onefile){
    load(onefile)
    index_error <- sapply(Sim_onesigma, function(j) class(j) == "try-error")
    sapply(Sim_onesigma[!index_error], function(j) j$dselected)
  })
  mp <- Reduce("cbind", dres)
  trued <- Settings$d[ii]
  apply(mp, 1, function(k) sum(k == trued))/ncol(mp)
})


dresults_df <- data.frame(Settings,
                     t(dresults_overall))
dresults_df$Sigmatype <- factor(dresults_df$Sigmatype)

#### Final plots
plot_choice_d <- dresults_df %>%
  arrange(Sigmatype, d, n) %>%
  rename(s_aic = saic, s_bic = sbic, g_aic = gaic, g_bic = gbic) %>%
  pivot_longer(cols = c("s_aic", "s_bic", "g_aic", "g_bic"),
               names_to = c("Method", "Criterion"), 
               values_to = "prop", 
               names_sep = "_") %>%
  mutate(Sigmatype = factor(Sigmatype, 
                            levels = c(1, 2, 3), 
                            labels = c("Diagonal", "AR(1)", "Exchangeable")),
         Criterion = factor(Criterion, 
                         levels = c("aic", "bic"), 
                         labels = c("AIC", "BIC"))) %>%
  ggplot(aes(x = n, y = prop, col = Sigmatype)) +
  geom_line(aes(linetype = Criterion), linewidth = 2) + 
  theme_bw(base_size = 20) +
  xlim(c(100, 1000))+
  labs(x = bquote(n), y = "") +
  facet_grid(rows = vars(d), cols = vars(Method),
             scales = "free",
             labeller = labeller(d = c(`1` = "d = 1", `2` = "d = 2"),
                                 Method = c(`s` = "SPFC", `g` = "GPFC"))) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 20),
        legend.key.size=unit(3,"lines")) +
  guides(colour = guide_legend(title = bquote(tilde(Sigma))))

ggsave(plot_choice_d, 
       file = "FigureS1.pdf", width = 12, height = 8)




