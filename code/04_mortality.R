
library(readr)
library(stats)
library(tidyverse)
library(tidylog)
library(gtools)
library(patchwork)
library(gtsummary)
library(ggsvg)
library(svglite)
library(MASS, exclude = "select")
library(ggrepel)
require(kableExtra)
library(arrow)
library(nflplotR) # ggpreview
library(lemon)
library(ggh4x)


stylable <- function(df){
  df %>%
    kableExtra::kable()  %>%
    kable_styling(bootstrap_options = "striped", 
                  position = "left",
                  full_width = FALSE)
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
                  legend.position = "none",
                  strip.text = element_blank(),
                  legend.text = element_text(colour = "black", size = 9),
                  panel.grid.major.x = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(colour = "black", linewidth = rel(1)),
                  plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")))



data <- read_parquet("input/mr_care_age_boot_jpn_2017_2020.parquet") %>% 
  mutate(country="Japan") %>% 
  bind_rows(read_parquet("input/mr_care_age_boot_swe_2017_2020.parquet") %>% 
              mutate(country="Sweden")) %>% 
  rename(care=basecarestatus, age=age_group) %>%
  mutate(care = as.character(case_when(care == "no care" ~ "N", 
                                       care == "home care" ~ "H", 
                                       care == "care home" ~ "C")) %>% 
           factor(levels=c("N", "H", "C")), 
         age=factor(age, levels=c("75-79", "80-84", "85-89", "90-94", "95-99", "100+")),
         sex = factor(sex, levels=c("men", "women"), labels=c("Men", "Women")),
         country=factor(country, levels=c("Sweden", "Japan")),
         mr=death/py*1000) %>% 
  select(bootid, country, age, sex, care, death, n, py, mr)

# Read SVG
japan <- paste(readLines("mark/japan.svg"), collapse = "\n")
sweden <- paste(readLines("mark/sweden.svg"), collapse = "\n")

icons_df <- data.frame(
  country = c('Sweden', 'Japan'),
  svg  = c(sweden ,  japan),
  stringsAsFactors = FALSE)

age_weight <- data %>% 
  filter(bootid==0) %>% 
  group_by(country, age) %>% 
  summarise(n=sum(n)) %>% 
  group_by(country) %>% 
  mutate(total=sum(n),
         proportion=n/total) %>% 
  pivot_wider(id_cols = age,
              names_from = country,
              values_from = proportion) %>% 
  mutate(weight=(Sweden+Japan)/2) %>% 
  select(age, weight)

# Follow up time ----------------------------------------------------------

data %>% 
  filter(bootid==0) %>% 
  summarise(across(c(py, n), ~sum(.x)), .by=c(country)) %>% 
  mutate(ave_follow=py/n)

data %>% 
  filter(bootid==0) %>% 
  summarise(across(c(py, n), ~sum(.x))) %>% 
  mutate(ave_follow=py/n)

# Age distribution   ------------------------------------------------------

table <- data %>% 
  filter(bootid==0) %>% 
  group_by(country, sex, care, age) %>% 
  summarise(n=sum(n)) %>% 
  group_by(country, sex, care) %>% 
  mutate(total=sum(n),
         proportion=n/total*100)

table %>% 
  group_by(country, sex, care) %>% 
  summarise(total=sum(n)) %>% 
  arrange(sex, care) %>% 
  stylable()

# Mortality ---------------------------------------------------------------

## Standardization ---------------------------------------------------------

stad_table <- data %>% 
  filter(bootid==0) %>% 
  mutate(mr_i=death/py,
         var_i=death/(py^2),
         conf.low=mr_i + qnorm(0.025)*sqrt(var_i),
         conf.high=mr_i + qnorm(0.975)*sqrt(var_i)) %>% 
  left_join(age_weight, by="age") %>% 
  group_by(sex, care, country) %>% 
  summarise(std.mr=sum(weight*mr_i),
            std.var=sum(weight^2*var_i),
            std.se=sqrt(std.var)) %>% 
  pivot_wider(id_cols=c(sex, care),
              names_from=country,
              values_from=starts_with("std")) %>% 
  mutate(conf.low_Japan=std.mr_Japan + qnorm(0.025)*sqrt(std.var_Japan),
         conf.high_Japan=std.mr_Japan + qnorm(0.975)*sqrt(std.var_Japan),
         conf.low_Sweden=std.mr_Sweden + qnorm(0.025)*sqrt(std.var_Sweden),
         conf.high_Sweden=std.mr_Sweden + qnorm(0.975)*sqrt(std.var_Sweden),
         smrd=std.mr_Sweden - std.mr_Japan,
         std.d_se=std.se_Japan + std.se_Sweden,
         conf.low.d=smrd + qnorm(0.025)*std.d_se,
         conf.high.d=smrd + qnorm(0.975)*std.d_se,
         smrr=std.mr_Sweden / std.mr_Japan,
         std.r_se=sqrt(std.var_Japan/std.mr_Japan^2 + std.var_Sweden/std.mr_Sweden^2),
         conf.low.r=exp(log(smrr) + qnorm(0.025)*std.r_se),
         conf.high.r=exp(log(smrr) + qnorm(0.975)*std.r_se),
         across(c(std.mr_Japan, conf.low_Japan, conf.high_Japan, std.mr_Sweden, conf.low_Sweden, conf.high_Sweden, smrd, std.d_se, conf.low.d, conf.high.d), ~.x*10^3),
         care=factor(care, 
                     levels=c("N", "H", "C"),
                     labels=c("No care", "Home care", "Care home"))) %>% 
  ungroup() %>% 
  select(sex, care, 
         std.mr_Japan, conf.low_Japan, conf.high_Japan, 
         std.mr_Sweden, conf.low_Sweden, conf.high_Sweden, 
         smrd, conf.low.d, conf.high.d, 
         smrr, conf.low.r, conf.high.r)

standardization <- stad_table %>% 
  select(sex, care, std.mr_Japan, conf.low_Japan, conf.high_Japan, 
         std.mr_Sweden, conf.low_Sweden, conf.high_Sweden) %>% 
  pivot_longer(cols=c(std.mr_Japan, conf.low_Japan, conf.high_Japan, 
                      std.mr_Sweden, conf.low_Sweden, conf.high_Sweden),
               names_to=c("type", "country"),
               names_sep = "_",
               values_to = "score") %>% 
  pivot_wider(names_from = type,
              names_vary = "fastest",
              values_from = score) %>% 
  left_join(icons_df, by="country") %>% 
  rename(age_standardized_MR=std.mr)


graph_men <- standardization |>   
  filter(sex=="Men") %>% 
  ggplot(data = ., aes(x=care, y=age_standardized_MR)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, group=country), position=position_dodge(width=0.3), linewidth=0.2, width=0.3) +
  geom_point_svg(aes(svg = country), position=position_dodge(width=0.3), size = 3) +
  labs(y="Age-standardized death rate \nper 1000 person-year", x="",
       title="Men",
       svg="Country") +
  scale_y_continuous(limits=c(0, 425),
                     breaks = seq(0, 400, 50),
                     expand = c(0,0)) +
  scale_svg_discrete_manual(
    aesthetics = 'svg', 
    values = c(Sweden = sweden, Japan = japan),
    guide = guide_legend(override.aes = list(size = 1))) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = -20, b = 0, l = 0)),
        legend.position = "right",
        axis.ticks.x = element_blank())

gtable_men <- stad_table %>% 
  filter(sex=="Men") %>% 
  mutate(SR=paste0(round2(smrr, d=2) %>% sprintf("%.2f", .),
                   " (",
                   round2(conf.low.r, d=2) %>% sprintf("%.2f", .),
                   ", ",
                   round2(conf.high.r, d=2) %>% sprintf("%.2f", .),
                   ")"),
         SD=paste0(round2(smrd, d=1) %>% sprintf("%.1f", .),
                   " (",
                   round2(conf.low.d, d=1) %>% sprintf("%.1f", .),
                   ", ",
                   round2(conf.high.d, d=1) %>% sprintf("%.1f", .),
                   ")")) %>% 
  select(care, `SR (95%CI)`=SR, `SD (95%CI)`=SD) %>% 
  pivot_longer(cols=-care,
               names_to = " ") %>% 
  mutate(` `=factor(` `, levels=c("SD (95%CI)", "SR (95%CI)"))) %>% 
  ggplot(data=., aes(x = care, y = ` `)) +
  geom_text(aes(label = value), size=3) +
  scale_y_discrete() +
  labs(y = NULL, x = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 10, vjust = 0.5, hjust = 0.5, angle = 0),
        panel.grid = element_blank(), 
        strip.text = element_blank(), 
        panel.spacing.x = unit(0, "mm"))

# women
graph_women <- standardization |>   
  filter(sex=="Women") %>% 
  ggplot(data = ., aes(x=care, y=age_standardized_MR)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, group=country), position=position_dodge(width=0.3), linewidth=0.2, width=0.3) +
  geom_point_svg(aes(svg = country), position=position_dodge(width=0.3), size = 3) +
  labs(y=NULL, x=NULL,
       title="Women",
       svg="Country") +
  scale_y_continuous(limits=c(0, 425),
                     breaks = seq(0, 400, 50),
                     expand = c(0,0)) +
  scale_svg_discrete_manual(
    aesthetics = 'svg', 
    values = c(Sweden = sweden, Japan = japan),
    guide = guide_legend(override.aes = list(size = 1))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        legend.position = "right",
        axis.ticks.x = element_blank())

gtable_women <- stad_table %>% 
  filter(sex=="Women") %>% 
  mutate(SR=paste0(round2(smrr, d=2) %>% sprintf("%.2f", .),
                   " (",
                   round2(conf.low.r, d=2) %>% sprintf("%.2f", .),
                   ", ",
                   round2(conf.high.r, d=2) %>% sprintf("%.2f", .),
                   ")"),
         SD=paste0(round2(smrd, d=1) %>% sprintf("%.1f", .),
                   " (",
                   round2(conf.low.d, d=1) %>% sprintf("%.1f", .),
                   ", ",
                   round2(conf.high.d, d=1) %>% sprintf("%.1f", .),
                   ")")) %>% 
  select(care, ` `=SR, `  `=SD) %>% 
  pivot_longer(cols=-care,
               names_to = " ") %>% 
  mutate(` `=factor(` `, levels=c("  ", " "))) %>% 
  ggplot(data=., aes(x = care, y = ` `)) +
  geom_text(aes(label = value), size=3) +
  scale_y_discrete() +
  labs(y = NULL, x = NULL) +
  theme_classic() +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 8, vjust = 0.5, hjust = 0.5, angle = 0),
        panel.grid = element_blank(), 
        strip.text = element_blank(), 
        panel.spacing.x = unit(0, "mm"))

gcombine <- graph_men + graph_women + 
  plot_spacer() + plot_spacer()+
  gtable_men + gtable_women + 
  plot_layout(guides = "collect", heights = c(5, -0.4, 1), nrow=3, ncol=2) & theme(legend.position = 'bottom')


ggsave(filename = "output/supl.fig4_standardized_mr.svg", plot=gcombine, width=8.5, height = 4)


## Age specific mortality --------------------------------------------------

sd <- stad_table %>% 
  mutate(SD=paste0("SMD (95%CI):\n", round2(smrd, d=1) %>% sprintf("%.1f", .),
                   " (",
                   round2(conf.low.d, d=1) %>% sprintf("%.1f", .),
                   ", ",
                   round2(conf.high.d, d=1) %>% sprintf("%.1f", .),
                   ")")) %>% 
  select(sex, care, SD)

mr_age <-
  data %>% 
  filter(bootid==0) %>% 
  group_by(age, sex, country, care) %>% 
  nest() %>% 
  mutate(conf.low=map(.x=data, .f=~poisson.test(.x$death)[["conf.int"]][[1]]),
         conf.high=map(.x=data, .f=~poisson.test(.x$death)[["conf.int"]][[2]])) %>% 
  unnest(data) %>% 
  unnest(conf.low) %>% 
  unnest(conf.high) %>% 
  mutate(across(c(conf.low, conf.high), .f=~.x/py*1000),
         care=factor(care, 
                     levels=c("N", "H", "C"),
                     labels=c("No care", "Home care", "Care home"))) %>% 
  left_join(icons_df, by="country") 

graph_men <- mr_age %>%    
  filter(sex=="Men") %>% 
  ggplot(data = ., aes(x=age, y=mr, group=country)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position=position_dodge(width=0.5), linewidth=0.3, width=0.5) +
  geom_point_svg(aes(svg = I(svg)), position=position_dodge(width=0.5), size = 3) +
  geom_text(data=sd %>% filter(sex=="Men"), aes(y=30, x=5, group=NULL, label = SD)) +
  labs(y="Death rate (log-scale) \nper 1000 person-year", x=NULL,
       subtitle="Men", color=NULL, svg=NULL) +
  coord_trans(y="log10", clip = 'off') +
  scale_y_continuous(limits=c(8, 1020.1),
                     breaks = c(8, 16, 32, 64, 128, 255, 510, 1020),
                     minor_breaks = NULL,
                     expand = c(0,0)) +
  scale_svg_discrete_manual(
    aesthetics = 'svg', 
    values = c(Sweden = sweden, Japan = japan),
    guide = guide_legend(override.aes = list(size = 1))) +
  theme(legend.position = "bottom")+
  facet_rep_grid(~care) + 
  theme_bw() +
  theme(title =element_text(size=12),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 0.5, angle = 0),
        axis.text.y = element_text(colour = "black", size = 12, hjust = 1),
        strip.text.x = element_text(size = 12, colour = "black", hjust = 0),
        strip.background=element_blank(),
        legend.position = "top",
        legend.text = element_text(colour = "black", size = 12),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = rel(1)),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

graph_women <- mr_age %>%    
  filter(sex=="Women") %>% 
  ggplot(data = ., aes(x=age, y=mr, group=country)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position=position_dodge(width=0.5), linewidth=0.3, width=0.5) +
  geom_point_svg(aes(svg = I(svg)), position=position_dodge(width=0.5), size = 3) +
  geom_text(data=sd %>% filter(sex=="Women"), aes(y=30, x=5, group=NULL, label = SD)) +
  labs(y="Death rate (log-scale) \nper 1000 person-year", x="Age",
       subtitle="Women", color=NULL, svg=NULL) +
  coord_trans(y="log10") +
  scale_y_continuous(limits=c(8, 1020.1),
                     breaks = c(8, 16, 32, 64, 128, 255, 510, 1020),
                     minor_breaks = NULL,
                     expand = c(0,0)) +
  scale_svg_discrete_manual(
    aesthetics = 'svg', 
    values = c(Sweden = sweden, Japan = japan),
    guide = guide_legend(override.aes = list(size = 1))) +
  theme_bw() +
  theme(title =element_text(size=12),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 0.5, angle = 0),
        axis.text.y = element_text(colour = "black", size = 12, hjust = 1),
        strip.text.x = element_text(size = 12, colour = "black", hjust = 0),
        strip.background=element_blank(),
        legend.position = "top",
        legend.text = element_text(colour = "black", size = 12),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = rel(1)),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))+
  facet_rep_grid(~care)

graph <- graph_men / graph_women +
  plot_layout(guides = "collect", 
              axis_titles = "collect") & 
  theme(legend.position = 'bottom')

# graph
ggsave(filename = "output/fig3_mr_age_sex_country_care.pdf", plot=graph, width=12, height = 8)
