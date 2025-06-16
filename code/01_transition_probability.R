library(readr)
library(stats)
library(tidyverse)
library(tidylog)
library(GGally)
library(patchwork)
library(gtsummary)
library(gt)
library(gtools)
library(kableExtra)
library(Matrix)
library(matrixcalc)
library(ggh4x)
library(svglite)
library(ggsvg)
library(gridExtra)
library(nflplotR)
library(furrr)
library(mgcv)
library(dtms)

stylable <- function(df){
  df %>%
    datatable(extensions = "Buttons",
              editable = FALSE,
              rownames = FALSE, 
              style = "auto",
              class = "cell-border stripe",
              options = list(dom = 'BSlftipr',
                             buttons = 'csv',
                             paging=FALSE,
                             ordering=FALSE,
                             scroller=FALSE,
                             info = FALSE,
                             scrollCollapse = T,
                             searching = FALSE))
}

round2 = function(x, d=0) {
  p = 10^d
  return((x * p * 2 + 1) %/% 2 / p)
}

theme_set(theme_bw() +
            theme(title =element_text(size=12),
                  axis.title.x = element_text(colour = "black", size = 12),
                  axis.title.y = element_text(colour = "black", size = 12),
                  axis.text.x = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 0.5, angle = 0),
                  axis.text.y = element_text(colour = "black", size = 12, hjust = 1),
                  axis.line = element_line(colour = "black", linewidth = rel(1)),
                  strip.text = element_text(size = 12, colour = "black"),
                  strip.background=element_blank(),
                  legend.position = "bottom",
                  legend.text = element_text(colour = "black", size = 12),
                  legend.background = element_blank(),
                  legend.box.background = element_blank(),
                  legend.axis.line = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_rect(fill='transparent'),
                  plot.background = element_rect(fill='transparent', color=NA),
                  plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")))


# Transition probability --------------------------------------------------

## smoothed probability ----------------------------------------------------

bootdata_aggregated <- read_rds("input/aggregted_bootdata.rds")

fit.dat_all <-
  bootdata_aggregated |>
  # filter(bootid==0) |> # for trial
  group_by(bootid, country, sex, status_start, status_end) %>%
  nest() %>% 
  mutate(data=map(.x=data, .f=~mutate(.x, pred_p = predict(
    gam(data=.x, Cases ~ s(age, bs="ps") + offset(log(Tpop)), 
        family = poisson(), weights = weight_c))))) %>%
  unnest(data)|> 
  mutate(pred_p=exp(pred_p)/Tpop) |> 
  ungroup() %>% 
  write_rds("intermediate/predicteddata_bootstrap.rds")

fit.dat_all <- read_rds("intermediate/predicteddata_bootstrap.rds")

fig_combine <-
  fit.dat_all %>%
  filter(bootid==0) %>% # bootid==0 means original data  
  left_join(fit.dat_all %>%
              filter(bootid!=0) %>% 
              group_by(country, sex, status_start, status_end, age) %>% 
              summarise(conf.low=quantile(pred_p, probs=0.025),
                        conf.high=quantile(pred_p, probs=0.975)),
            by=c("country", "sex", "status_start", "status_end", "age")) %>% 
  mutate(sex=factor(sex,
                    levels=c("men", "women"),
                    labels=c("Men", "Women")),
         country=factor(country, levels=c("sweden", "japan"), labels=c("Sweden", "Japan")),
         group=case_when(status_start==status_end ~ "Remaining",
                         status_start==1 ~ "From no care",
                         status_start==2 ~ "From home care",
                         status_start==3 ~ "From care home") %>%
           factor(levels=c("From no care", "From home care", "From care home", "Remaining")),
         status_end=factor(status_end),
         conf.high=if_else(conf.high > 0.5 & 
                             group %in% c("From no care", "From home care", "From care home"), 
                           0.5, conf.high),
         across(c(pred_p, conf.low, conf.high), ~.x*100)) %>%
  ggplot(data=., aes(y=pred_p, x=age, group=paste0(status_end, country), linetype=country)) +
  geom_line(aes(color=status_end), linewidth=0.5, alpha=0.8) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=status_end), alpha=0.3) +
  scale_x_continuous(breaks = seq(75, 105, 5), expand=c(0,0))+
  scale_fill_manual(values = c("1"="#add8e6", "2"="#ffb6c1", "3"="#eeee00", "4"="black"),
                    labels = c("No care", "Home care", "Care home", "Death")) +
  scale_color_manual(values = c("1"="#add8e6", "2"="#ffb6c1", "3"="#eeee00", "4"="black"),
                     labels = c("No care", "Home care", "Care home", "Death")) +
  labs(y="Probability, %", x="Age", fill="To", color="To", linetype="Country")+
  facet_grid2(sex~group, scales = "free", independent = "all") +
  scale_y_facet(
    group %in% "Remaining", limits = c(50, 100), breaks = seq(50, 100, 10), expand=c(0,0)
  ) +
  scale_y_facet(
    group %in% c("From no care", "From home care", "From care home"), limits = c(0, 50), 
    breaks = seq(0, 50, 10), expand=c(0,0)
  )

# fig_combine
ggsave(filename = "output/supl.fig2_smooth_probability.pdf", plot=fig_combine, width=10, height = 6)


## observed probability ----------------------------------------------------

bootdata_aggregated <- read_rds("input/aggregted_bootdata.rds")

fig_combine <-
  bootdata_aggregated |>
  filter(bootid==0) |> # bootid==0 means original data 
  mutate(sex=factor(sex,
                    levels=c("men", "women"),
                    labels=c("Men", "Women")),
         country=factor(country, levels=c("sweden", "japan"), labels=c("Sweden", "Japan")),
         group=case_when(status_start==status_end ~ "Remaining",
                         status_start==1 ~ "From no care",
                         status_start==2 ~ "From home care",
                         status_start==3 ~ "From care home") %>%
           factor(levels=c("From no care", "From home care", "From care home", "Remaining")),
         status_end=factor(status_end),
         p=p*100) %>%
  ggplot(data=., aes(y=p, x=age, group=paste0(status_end, country), linetype=country)) +
  geom_line(aes(color=status_end), linewidth=0.5) +
  scale_x_continuous(breaks = seq(75, 105, 5), expand=c(0,0))+
  scale_color_manual(values = c("1"="#add8e6", "2"="#ffb6c1", "3"="#eeee00", "4"="black"),
                     labels = c("No care", "Home care", "Care home", "Death")) +
  labs(y="Probability, %", x="Age", fill="To", color="To", linetype="Country")+
  facet_grid2(sex~group, scales = "free", independent = "all") +
  scale_y_facet(
    group %in% "Remaining", limits = c(50, 100), breaks = seq(50, 100, 10), expand=c(0,0)
  ) +
  scale_y_facet(
    group %in% c("From no care", "From home care", "From care home"), limits = c(0, 50), 
    breaks = seq(0, 50, 10), expand=c(0,0)
  )

# fig_combine
ggsave(filename = "output/supl.fig3_observed_probability.pdf", plot=fig_combine, width=10, height = 6)
