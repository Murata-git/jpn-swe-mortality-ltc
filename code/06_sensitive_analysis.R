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

# Sensitive analysis: 2017-2022  


data_sensitive <- read_parquet("input/mr_care_age_boot_jpn_2017_2022.parquet") %>% 
  mutate(country="Japan") %>% 
  bind_rows(read_parquet("input/mr_care_age_boot_swe_2017_2022.parquet") %>% 
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

# Read multiple SVG somehow
japan <- paste(readLines("mark/japan.svg"), collapse = "\n")
sweden <- paste(readLines("mark/sweden.svg"), collapse = "\n")

icons_df <- data.frame(
  country = c('Sweden', 'Japan'),
  svg  = c(sweden ,  japan),
  stringsAsFactors = FALSE)

age_weight <- data_sensitive %>% 
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


# Mortality ---------------------------------------------------------------

## Standardized mortality   ------------------------------------------------

stad_table <- data_sensitive %>% 
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


##Age specific mortality  ------------------------------------------------

sd <- stad_table %>% 
  mutate(SD=paste0("SMD (95%CI):\n", round2(smrd, d=1) %>% sprintf("%.1f", .),
                   " (",
                   round2(conf.low.d, d=1) %>% sprintf("%.1f", .),
                   ", ",
                   round2(conf.high.d, d=1) %>% sprintf("%.1f", .),
                   ")")) %>% 
  select(sex, care, SD)

mr_age_sensitive <-
  data_sensitive %>% 
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

graph_men <- mr_age_sensitive %>%    
  filter(sex=="Men") %>% 
  ggplot(data = ., aes(x=age, y=mr, group=country)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position=position_dodge(width=0.5), linewidth=0.3, width=0.5) +
  geom_point_svg(aes(svg = I(svg)), position=position_dodge(width=0.5), size = 3) +
  geom_text(data=sd %>% filter(sex=="Men"), aes(y=30, x=5, group=NULL, label = SD)) +
  labs(y="Death rate (log-scale) \nper 1000 person-year", x="Age",
       subtitle="Men", color=NULL, svg=NULL) +
  coord_trans(y="log10") +
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

graph_women <- mr_age_sensitive %>%    
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

ggsave(filename = "output/supl.fig6_mr_age_sex_country_care_sensitive.pdf", plot=graph, width=12, height = 8)

# Horiuchi decomposition methods  -----------------------------------------

# Prepare Dataset
dat <- 
  data_sensitive %>% 
  rename(Dx = death, Nx = py) %>%
  mutate(Age = as.numeric(case_when(age == "75-79" ~ 75,
                                    age == "80-84" ~ 80,
                                    age == "85-89" ~ 85,
                                    age == "90-94" ~ 90,
                                    age == "95-99" ~ 95,
                                    age == "100+" ~ 100)))

# Calculate total rate
cdr <- 
  dat %>% 
  # filter(bootid==i) %>% 
  group_by(bootid, country, sex, Age) %>% 
  summarise(Dx = sum(Dx), 
            Nx = sum(Nx)) %>% 
  mutate(mx = Dx / Nx) %>% 
  pivot_wider(id_cols=c(bootid, sex, Age),
              names_from=country,
              values_from=mx) %>% 
  mutate(mr_dif=Sweden-Japan)

# cdr

# Decomposition of age specific mortality difference into care-specific mortality difference and LTC distribution difference
dec_c <- 
  dat %>%
  # filter(bootid==i) %>% 
  arrange(bootid, country, sex, Age) %>% 
  select(bootid, country, sex, Age, care, Dx, Nx, n) %>% 
  mutate(mx = Dx / Nx) %>% 
  group_by(bootid, sex, Age, country) %>% 
  mutate(p = n / sum(n)) %>% 
  pivot_wider(id_cols = c(bootid, sex, Age), 
              names_from = c(country, care), 
              values_from = c(mx, p, n)) %>% 
  # left_join(age_cx) %>% 
  mutate(ave_p_C = (p_Sweden_C+p_Japan_C)/2,
         ave_p_H = (p_Sweden_H+p_Japan_H)/2,
         ave_p_N = (p_Sweden_N+p_Japan_N)/2,
         ave_mx_C = (mx_Sweden_C+mx_Japan_C)/2,
         ave_mx_H = (mx_Sweden_H+mx_Japan_H)/2,
         ave_mx_N = (mx_Sweden_N+mx_Japan_N)/2)

full_cdr <- 
  dec_c %>% 
  mutate(mx_JPN = mx_Japan_N*p_Japan_N+mx_Japan_H*p_Japan_H+mx_Japan_C*p_Japan_C,
         mx_SWE = mx_Sweden_N*p_Sweden_N+mx_Sweden_H*p_Sweden_H+mx_Sweden_C*p_Sweden_C) %>% 
  mutate(diff = mx_SWE-mx_JPN)

######## Horiuchi
x <- dec_c %>% 
  ungroup() %>% 
  mutate(mort_N=(ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*mx_Sweden_N) - 
           (ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*mx_Japan_N),
         mort_H=(ave_p_C*ave_mx_C+ave_p_H*mx_Sweden_H+ave_p_N*ave_mx_N) - 
           (ave_p_C*ave_mx_C+ave_p_H*mx_Japan_H+ave_p_N*ave_mx_N),
         mort_C=(ave_p_C*mx_Sweden_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N) - 
           (ave_p_C*mx_Japan_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N),
         cx_N=(ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+p_Sweden_N*ave_mx_N) - 
           (ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+p_Japan_N*ave_mx_N),
         cx_H=(ave_p_C*ave_mx_C+p_Sweden_H*ave_mx_H+ave_p_N*ave_mx_N) - 
           (ave_p_C*ave_mx_C+p_Japan_H*ave_mx_H+ave_p_N*ave_mx_N),
         cx_C=(p_Sweden_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N) - 
           (p_Japan_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N),
         comp_effect=cx_N+cx_H+cx_C)

dec_result <- 
  x %>% 
  mutate(diff_dec = mort_N+mort_H+mort_C+comp_effect) %>%
  left_join(full_cdr %>% select(sex, Age, diff)) %>% 
  left_join(cdr %>% select(sex, Age, mr_dif)) %>% 
  mutate(residual=diff-diff_dec,
         residual2=mr_dif-diff_dec) %>% 
  write_rds("intermediate/boottable_horiuchi_dec_agespecific_sensitive.rds")


# diff_dec: total difference
# mort_N: mortality from no care
# mort_H: mortality from homeo care
# mort_C: mortality from care care
# comp_effect: distributional difference

horiuchi_table <- read_rds("intermediate/boottable_horiuchi_dec_agespecific_sensitive.rds") %>% 
  filter(bootid==0) %>% 
  select(-bootid) %>% 
  rename_with(~paste0(.x, ".estimate"), -c(Age,sex)) %>% 
  left_join(read_rds("intermediate/boottable_horiuchi_dec_agespecific_sensitive.rds") %>% 
              filter(bootid!=0) %>% 
              group_by(Age, sex) %>% 
              summarise(across(c(diff_dec, mort_N, mort_H, mort_C, comp_effect), 
                               ~quantile(.x, probs = c(0.025)),
                               .names="{.col}.conf_low"),
                        across(c(diff_dec, mort_N, mort_H, mort_C, comp_effect), 
                               ~quantile(.x, probs = c(0.975)),
                               .names="{.col}.conf_high")),
            by=c("Age", "sex")) %>% 
  select(Age, sex, contains("diff_dec"), contains("mort_N"), contains("mort_H"), contains("mort_C"), contains("comp_effect")) %>% 
  pivot_longer(cols=contains(c("diff_dec", "mort_N", "mort_H", "mort_C", "comp_effect")),
               names_to = c("type", ".value"),
               names_sep = "\\.") %>% 
  mutate(across(c("estimate", "conf_low", "conf_high"), ~.x*1000))

diff_dec <- read_rds("intermediate/boottable_horiuchi_dec_agespecific_sensitive.rds") %>% 
  filter(bootid==0) %>% 
  select(diff_dec, sex, Age) %>% 
  left_join(
    read_rds("intermediate/boottable_horiuchi_dec_agespecific_sensitive.rds") %>% 
      filter(bootid!=0) %>% 
      select(diff_dec, sex, Age) %>% 
      summarise(conf_low.diff_dec=quantile(diff_dec, probs = c(0.025)),
                conf_high.diff_dec=quantile(diff_dec, probs = c(0.975)),
                .by=c(sex, Age)))

horiuchi_table <- horiuchi_table %>% 
  left_join(horiuchi_table %>% 
              filter(type=="diff_dec") %>% 
              select(Age, sex, total=estimate),
            by=c("Age", "sex")) %>% 
  mutate(proportion=round2(estimate/total*100, d=1) %>% 
           sprintf("%.1f", .) %>% 
           paste0("%"),
         yposition=0,
         type=factor(type, levels=c("diff_dec", "mort_N", "mort_H", "mort_C", "comp_effect"))) %>% 
  left_join(diff_dec %>% 
              mutate(across(contains("diff_dec"), ~round2(.x*1000, 1))), by=c("sex", "Age")) %>% 
  {.->> horiuchi_table_restore} %>% 
  mutate(Age=as.numeric(factor(Age)))

fig_combine <- horiuchi_table %>% 
  mutate(conf_low=if_else(conf_low < -80, -80, conf_low),
         conf_high=if_else(conf_high > 160, 160, conf_high)) %>%
  ggplot(data=., aes(x=Age, y=estimate, fill=type)) +
  geom_hline(yintercept = 0, linewidth=0.3)+
  geom_bar(position = position_dodge2(width =0.5), stat="identity", width=0.7, color="black", linewidth=0.2) +
  geom_linerange(aes(ymax=conf_high, ymin=conf_low), position = position_dodge2(width =0.7), linewidth=0.2) +
  geom_segment(data=. %>% filter(sex=="Men"), aes(x = 5.72, xend = 5.72, y = -79, yend = -80), arrow = arrow(length = unit(1, "mm")), linewidth=0.01)+
  geom_segment(data=. %>% filter(sex=="Men"), aes(x = 5.86, xend = 5.86, y = -79, yend = -80), arrow = arrow(length = unit(1, "mm")), linewidth=0.01)+
  geom_segment(data=. %>% filter(sex=="Men"), aes(x = 6, xend = 6, y = -79, yend = -80), arrow = arrow(length = unit(1, "mm")), linewidth=0.01)+
  geom_segment(data=. %>% filter(sex=="Men"), aes(x = 5.72, xend = 5.72, y = 159, yend = 160), arrow = arrow(length = unit(1, "mm")), linewidth=0.01)+
  geom_segment(data=. %>% filter(sex=="Men"), aes(x = 6.14, xend = 6.14, y = 159, yend = 160), arrow = arrow(length = unit(1, "mm")), linewidth=0.01)+
  coord_cartesian(xlim = c(1, 6),    # set the x limits
                  ylim = c(-80, 160),  # set the y limits
                  clip = 'off')+
  scale_y_continuous(
    expand = c(0,0),
    breaks = seq(-80, 160, 20))+
  scale_x_continuous(breaks =1:6,
                     labels = c("75-79", "80-84", "85-89", "90-94", "95-99", "100+"),
                     expand = c(0.1,0.1))+
  scale_fill_manual(values = c("diff_dec"="grey", "mort_N"="#add8e6", "mort_H"="#ffb6c1", "mort_C"="#eeee00",  "comp_effect"="#2ecc71"),
                    breaks = c("diff_dec", "mort_N", "mort_H", "mort_C", "comp_effect"),
                    labels = c("Total", "No care mortality", "Home care mortality", "Care home mortality", "LTC distribution"))+
  labs(y="Difference per 1000 person-years\n(Sweden - Japan)", x="Age", fill=NULL, color=NULL, linetype=NULL, svg=NULL) +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks.x=element_blank(),
    axis.line.x=element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 0.5, angle = 0),
    axis.text.y = element_text(colour = "black", size = 10, vjust = 0.5, hjust = 1, angle = 0),
    strip.text = element_text(size = 12, colour = "black"),
    strip.background=element_blank(),
    legend.position = "bottom",
    plot.margin = unit(c(0.4,1,0.4,0.4), "cm"))+
  facet_grid2(.~sex,scales = "free", independent = "all")
# fig_combine

ggsave(filename = "output/supl.fig7_horiuchi_age_sensitive.svg", plot=fig_combine, width=8, height = 4)

horiuchi_table_restore %>% 
  mutate(difference=paste0(round2(estimate, d=2) %>% sprintf("%.2f", .),
                           " (",
                           round2(conf_low, d=2) %>% sprintf("%.2f", .),
                           ", ",
                           round2(conf_high, d=2) %>% sprintf("%.2f", .),
                           ")")) %>% 
  select(Age, sex, type, difference, proportion) %>% 
  stylable()




# Age-standardized mortality difference

weight <- dat |> 
  summarise(n=sum(n), .by=c(bootid, Age)) |> 
  mutate(total=sum(n), .by=bootid) |> 
  mutate(weight=n/total) |> 
  select(-c(n, total))

# Calculate total rate
cdr <- 
  dat %>% 
  summarise(across(c(Dx, Nx), ~sum(.x)), .by=c(bootid, country, Age, sex)) |> 
  left_join(weight, by=c("bootid", "Age")) %>%
  mutate(mx = Dx / Nx) %>%
  group_by(bootid, country, sex) %>% 
  summarise(std_mx=sum(mx*weight)) %>% 
  pivot_wider(id_cols=c(bootid, sex),
              names_from=country,
              values_from=std_mx) %>% 
  mutate(std_mr_dif=Sweden-Japan) %>%
  {.->>restore_cdr} |> 
  filter(bootid==0) |> 
  left_join(restore_cdr |> 
              filter(bootid!=0) %>% 
              group_by(sex) %>% 
              summarise(across(c(std_mr_dif), 
                               ~quantile(.x, probs = c(0.025)),
                               .names="{.col}.conf_low"),
                        across(c(std_mr_dif), 
                               ~quantile(.x, probs = c(0.975)),
                               .names="{.col}.conf_high")),
            by=c("sex"))
# cdr



dec_c <- 
  dat %>%
  # filter(bootid==i) %>% 
  arrange(bootid, country, sex, Age) %>% 
  select(bootid, country, sex, Age, care, Dx, Nx, n) %>% 
  mutate(mx = Dx / Nx) %>% 
  group_by(bootid, sex, Age, country) %>% 
  mutate(p = n / sum(n)) %>% 
  pivot_wider(id_cols = c(bootid, sex, Age), 
              names_from = c(country, care), 
              values_from = c(mx, p, n)) %>% 
  # left_join(age_cx) %>% 
  mutate(ave_p_C = (p_Sweden_C+p_Japan_C)/2,
         ave_p_H = (p_Sweden_H+p_Japan_H)/2,
         ave_p_N = (p_Sweden_N+p_Japan_N)/2,
         ave_mx_C = (mx_Sweden_C+mx_Japan_C)/2,
         ave_mx_H = (mx_Sweden_H+mx_Japan_H)/2,
         ave_mx_N = (mx_Sweden_N+mx_Japan_N)/2)

full_cdr <- 
  dec_c %>% 
  left_join(weight, by=c("bootid", "Age")) %>%
  group_by(bootid, sex) |> 
  summarise(mx_JPN = sum(mx_Japan_N*p_Japan_N*weight+mx_Japan_H*p_Japan_H*weight+mx_Japan_C*p_Japan_C*weight),
            mx_SWE = sum(mx_Sweden_N*p_Sweden_N*weight+mx_Sweden_H*p_Sweden_H*weight+mx_Sweden_C*p_Sweden_C*weight)) %>% 
  ungroup() |> 
  mutate(dif = mx_SWE-mx_JPN) %>% 
  {.->> full_cdr_restore} |> 
  filter(bootid==0) |> 
  select(sex, dif) |> 
  left_join(full_cdr_restore |> 
              filter(bootid!=0) |> 
              group_by(sex) %>% 
              summarise(across(c(dif), 
                               ~quantile(.x, probs = c(0.025)),
                               .names="{.col}.conf_low"),
                        across(c(dif), 
                               ~quantile(.x, probs = c(0.975)),
                               .names="{.col}.conf_high")),
            by=c("sex"))


######## Horiuchi

x <- dec_c %>% 
  left_join(weight, by=c("bootid", "Age")) %>%
  group_by(bootid, sex) %>% 
  summarise(mort_N=sum((ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*mx_Sweden_N)*weight) - 
              sum((ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*mx_Japan_N)*weight),
            mort_H=sum((ave_p_C*ave_mx_C+ave_p_H*mx_Sweden_H+ave_p_N*ave_mx_N)*weight) - 
              sum((ave_p_C*ave_mx_C+ave_p_H*mx_Japan_H+ave_p_N*ave_mx_N)*weight),
            mort_C=sum((ave_p_C*mx_Sweden_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N)*weight) - 
              sum((ave_p_C*mx_Japan_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N)*weight),
            cx_N=sum((ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+p_Sweden_N*ave_mx_N)*weight) - 
              sum((ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+p_Japan_N*ave_mx_N)*weight),
            cx_H=sum((ave_p_C*ave_mx_C+p_Sweden_H*ave_mx_H+ave_p_N*ave_mx_N)*weight) - 
              sum((ave_p_C*ave_mx_C+p_Japan_H*ave_mx_H+ave_p_N*ave_mx_N)*weight),
            cx_C=sum((p_Sweden_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N)*weight) - 
              sum((p_Japan_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N)*weight),
            comp_effect=cx_N+cx_H+cx_C)

look <- dec_c %>% 
  left_join(weight, by=c("bootid", "Age")) %>%
  filter(bootid==0) |> 
  mutate(mort_N=(ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*mx_Sweden_N)*weight - 
           (ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*mx_Japan_N)*weight,
         mort_H=(ave_p_C*ave_mx_C+ave_p_H*mx_Sweden_H+ave_p_N*ave_mx_N)*weight - 
           (ave_p_C*ave_mx_C+ave_p_H*mx_Japan_H+ave_p_N*ave_mx_N)*weight,
         mort_C=(ave_p_C*mx_Sweden_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N)*weight - 
           (ave_p_C*mx_Japan_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N)*weight,
         cx_N=(ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+p_Sweden_N*ave_mx_N)*weight - 
           (ave_p_C*ave_mx_C+ave_p_H*ave_mx_H+p_Japan_N*ave_mx_N)*weight,
         cx_H=(ave_p_C*ave_mx_C+p_Sweden_H*ave_mx_H+ave_p_N*ave_mx_N)*weight - 
           (ave_p_C*ave_mx_C+p_Japan_H*ave_mx_H+ave_p_N*ave_mx_N)*weight,
         cx_C=(p_Sweden_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N)*weight - 
           (p_Japan_C*ave_mx_C+ave_p_H*ave_mx_H+ave_p_N*ave_mx_N)*weight,
         comp_effect=cx_N+cx_H+cx_C)

dec_result <- 
  x %>% 
  mutate(diff_dec = mort_N+mort_H+mort_C+comp_effect) %>%
  filter(bootid==0) |> 
  left_join(x %>% 
              mutate(diff_dec = mort_N+mort_H+mort_C+comp_effect) %>%
              filter(bootid!=0) %>% 
              group_by(sex) %>% 
              summarise(across(c(diff_dec, mort_N, mort_H, mort_C, comp_effect), 
                               ~quantile(.x, probs = c(0.025)),
                               .names="{.col}.conf_low"),
                        across(c(diff_dec, mort_N, mort_H, mort_C, comp_effect), 
                               ~quantile(.x, probs = c(0.975)),
                               .names="{.col}.conf_high")),
            by=c("sex")) |> 
  left_join(full_cdr) %>% 
  rename_with(~ paste0(., ".estimate"), 
              .cols = c(diff_dec, mort_N, mort_H, mort_C, comp_effect, dif)) |> 
  mutate(residual=dif.estimate-diff_dec.estimate) |> 
  pivot_longer(cols=contains(c("diff_dec", "mort_N", "mort_H", "mort_C", "comp_effect")),
               names_to = c("type", ".value"),
               names_sep = "\\.") %>% 
  mutate(across(c("estimate", "conf_low", "conf_high"), ~.x*1000)) |> 
  ungroup() |> 
  select(sex, type, estimate, conf_low, conf_high, residual)

fig_combine <- dec_result %>% 
  # mutate(conf_low=if_else(conf_low < -80, -80, conf_low)) %>% 
  ggplot(data=., aes(x=type, y=estimate, fill=type)) +
  geom_hline(yintercept = 0, linewidth=0.3)+
  geom_bar(position = position_dodge2(width =2), stat="identity", width=0.9, color="black", linewidth=0.05) +
  geom_linerange(aes(ymax=conf_high, ymin=conf_low), position = position_dodge2(width =0.9), linewidth=0.2) +
  coord_cartesian(ylim = c(-5, 40),  # set the y limits
                  clip = 'off')+
  scale_y_continuous(
    expand = c(0,0),
    breaks = seq(-5, 40, 5))+
  scale_x_discrete(limits = c("diff_dec", "mort_N", "mort_H", "mort_C", "comp_effect"),
                   labels = c("Total", "No care mortality", "Home care mortality", "Care home mortality", "LTC distribution"))+
  scale_fill_manual(values = c("diff_dec"="grey", "mort_N"="#add8e6", "mort_H"="#ffb6c1", "mort_C"="#eeee00",  "comp_effect"="#2ecc71"),
                    breaks = c("diff_dec", "mort_N", "mort_H", "mort_C", "comp_effect"),
                    labels = c("Total", "No care mortality", "Home care mortality", "Care home mortality", "LTC distribution"))+
  labs(y="Age-standardized difference per 1000 person-years\n(Sweden - Japan)", x=NULL, fill=NULL, color=NULL, linetype=NULL, svg=NULL) +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks.x=element_blank(),
    axis.line.x=element_blank(),
    axis.title.y = element_text(colour = "black", size = 12),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 1, angle = 0),
    strip.text = element_text(size = 12, colour = "black"),
    strip.background=element_blank(),
    legend.position = "bottom",
    plot.margin = unit(c(0.4,1,0.4,0.4), "cm"))+
  facet_grid2(.~sex,scales = "free", independent = "all")
# ggpreview(fig_combine, width=7, height = 5)
ggsave(filename = "output/supl.fig8_horiuchi_agestandard_sensitive.svg", plot=fig_combine, width=7, height = 5)

