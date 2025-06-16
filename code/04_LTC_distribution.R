
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

# care proportion 


number_age_sex_care <- data %>% 
  filter(bootid==0) %>% 
  group_by(age, sex, country) %>% 
  mutate(total=sum(n)) %>% 
  ungroup() %>% 
  mutate(prop=n/total*100,
         prop_text=sprintf("%.0f", round2(prop, 0)),
         care=factor(care, 
                     levels=c("N", "H", "C"),
                     labels=c("No care", "Home care", "Care home"))) %>% 
  arrange(country, age, sex)

menplot <- number_age_sex_care %>% 
  filter(sex=="Men") %>% 
  group_by(country) %>% 
  nest() %>% 
  mutate(graph=map2(.x=data, .y=country, .f=~{ggplot(data=.x, aes(y=prop, x=age, fill=care)) +
      geom_bar(position = position_fill(reverse = TRUE), stat="identity", alpha=0.5) +
      geom_text(aes(label=prop_text), stat = "identity", position = position_fill(vjust=.5, reverse = TRUE), size = 3) +
      labs(y="%", x=NULL, fill=NULL, title = .y) + 
      scale_y_continuous(breaks = seq(0, 1, 0.2),
                         labels = seq(0, 100, 20),
                         expand = c(0,0)) +    
      scale_fill_manual(values=c("#add8e6", "#ee6363", "#eeee00")) +
      coord_flip()+
      theme_bw() +
      theme(plot.title =element_text(size=12, hjust=0.5),
            plot.subtitle = element_text(size=12, hjust=1.16),
            axis.title.x = element_text(colour = "black", size = 12),
            axis.title.y = element_text(colour = "black", size = 12),
            axis.text.x = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 0.5, angle = 0),
            axis.text.y = element_text(colour = "black", size = 12, hjust = 1),
            axis.ticks.y = element_blank(),
            legend.position = "none",
            strip.text.x = element_text(size = 12, colour = "black", hjust = 0),
            strip.background=element_blank(),
            panel.background = element_rect(fill='transparent'), 
            plot.background = element_rect(fill='transparent', color=NA),     
            legend.text = element_text(colour = "black", size = 12),
            panel.grid.major.y = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black", linewidth = rel(1)),
            plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"))}))


menplot$graph[[1]] <- menplot$graph[[1]] + 
  labs(subtitle="Age")+  
  annotate(geom = 'segment', y = -Inf, yend = -Inf, color = 'black', x = -Inf, xend = Inf, size = 0.8) +
  scale_y_reverse(breaks = seq(0, 1, 0.2),
                  labels = seq(0, 100, 20),
                  expand = c(0,0)) +
  theme(axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 0.5, angle = 0),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.5), "cm"),
        axis.line.x = element_line(colour = "black", linewidth = rel(1)))

womenplot <- number_age_sex_care %>% 
  filter(sex=="Women") %>% 
  group_by(country) %>% 
  nest() %>% 
  mutate(graph=map2(.x=data, .y=country, .f=~{ggplot(data=.x, aes(y=prop, x=age, fill=care)) +
      geom_bar(position = position_fill(reverse = TRUE), stat="identity", alpha=0.5) +
      geom_text(aes(label=prop_text), stat = "identity", position = position_fill(vjust=.5, reverse = TRUE), size = 3) +
      labs(y="%", x=NULL, fill=NULL, title = .y) + 
      scale_y_continuous(breaks = seq(0, 1, 0.2),
                         labels = seq(0, 100, 20),
                         expand = c(0,0)) +    
      scale_fill_manual(values=c("#add8e6", "#ee6363", "#eeee00")) +
      coord_flip()+
      theme_bw() +
      theme(plot.title =element_text(size=12, hjust=0.5),
            plot.subtitle = element_text(size=12, hjust=1.16),
            axis.title.x = element_text(colour = "black", size = 12),
            axis.title.y = element_text(colour = "black", size = 12),
            axis.text.x = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 0.5, angle = 0),
            axis.text.y = element_text(colour = "black", size = 12, hjust = 1),
            axis.ticks.y = element_blank(),
            legend.position = "none",
            strip.text.x = element_text(size = 12, colour = "black", hjust = 0),
            strip.background=element_blank(),
            legend.text = element_text(colour = "black", size = 12),
            panel.grid.major.y = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black", linewidth = rel(1)),
            plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"),
            panel.background = element_rect(fill='transparent'), 
            plot.background = element_rect(fill='transparent', color=NA) 
      )}))


womenplot$graph[[1]] <- womenplot$graph[[1]] + 
  labs(subtitle="Age")+
  annotate(geom = 'segment', y = -Inf, yend = -Inf, color = 'black', x = -Inf, xend = Inf, size = 0.8) +
  scale_y_reverse(breaks = seq(0, 1, 0.2),
                  labels = seq(0, 100, 20),
                  expand = c(0,0)) +
  theme(axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 0.5, angle = 0),
        plot.margin = unit(c(0.1,0.1,0.1,0.5), "cm"),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(colour = "black", linewidth = rel(1)))

carepoportionplot <- (menplot$graph[[1]] + menplot$graph[[2]] +
                        plot_layout(axes = "collect")) /
  (womenplot$graph[[1]] + womenplot$graph[[2]] +
     plot_layout(axes = "collect")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.1, 1)) &
  plot_annotation(tag_levels = list(c("Men", "", "Women", "")))

ggsave(filename = "output/fig3_proportion_ltc.pdf", plot=carepoportionplot, width=7, height = 6)


### Standardization

standardized_care <- data %>% 
  group_by(age, sex, country) %>% 
  mutate(total=sum(n)) %>% 
  ungroup() %>% 
  mutate(pr_i=n/total,
         var_i=pr_i*(1-pr_i)/total) %>% 
  arrange(country, age, sex) %>% 
  left_join(age_weight, by="age") %>% 
  group_by(sex, country, care) %>% 
  summarise(std.pr=sum(weight*pr_i),
            std.var=sum(weight^2*var_i),
            std.se=sqrt(std.var)) %>% 
  pivot_wider(id_cols=c(sex, care),
              names_from=country,
              values_from=starts_with("std")) %>% 
  mutate(conf.low_Japan=std.pr_Japan + qnorm(0.025)*sqrt(std.var_Japan),
         conf.high_Japan=std.pr_Japan + qnorm(0.975)*sqrt(std.var_Japan),
         conf.low_Sweden=std.pr_Sweden + qnorm(0.025)*sqrt(std.var_Sweden),
         conf.high_Sweden=std.pr_Sweden + qnorm(0.975)*sqrt(std.var_Sweden),
         std.difference=std.pr_Sweden - std.pr_Japan,
         std.d_se=std.se_Japan + std.se_Sweden,
         conf.low.d=std.difference + qnorm(0.025)*std.d_se,
         conf.high.d=std.difference + qnorm(0.975)*std.d_se,
         sprr=std.pr_Sweden / std.pr_Japan,
         std.r_se=sqrt(std.var_Japan/std.pr_Japan^2 + std.var_Sweden/std.pr_Sweden^2),
         conf.low.r=exp(log(sprr) + qnorm(0.025)*std.r_se),
         conf.high.r=exp(log(sprr) + qnorm(0.975)*std.r_se),
         across(c(std.pr_Japan, conf.low_Japan, conf.high_Japan, std.pr_Sweden, conf.low_Sweden, conf.high_Sweden, std.difference, std.d_se, conf.low.d, conf.high.d), ~.x*100),
         care=factor(care, 
                     levels=c("N", "H", "C"),
                     labels=c("No care", "Home care", "Care home"))) %>% 
  ungroup() %>% 
  select(sex, care,
         std.pr_Japan, conf.low_Japan, conf.high_Japan, 
         std.pr_Sweden, conf.low_Sweden, conf.high_Sweden, 
         std.difference, conf.low.d, conf.high.d, 
         sprr, conf.low.r, conf.high.r)

standardized_care %>% 
  mutate(std.difference=paste0(sprintf("%.1f", round2(std.difference, 1)),
                               " (",
                               sprintf("%.1f", round2(conf.low.d, 1)),
                               ", ",
                               sprintf("%.1f", round2(conf.high.d, 1)),
                               ")"),
         std.pr_Japan=paste0(sprintf("%.1f", round2(std.pr_Japan, 1)),
                             " (",
                             sprintf("%.1f", round2(conf.low_Japan, 1)),
                             ", ",
                             sprintf("%.1f", round2(conf.high_Japan, 1)),
                             ")"),
         std.pr_Sweden=paste0(sprintf("%.1f", round2(std.pr_Sweden, 1)),
                              " (",
                              sprintf("%.1f", round2(conf.low_Sweden, 1)),
                              ", ",
                              sprintf("%.1f", round2(conf.high_Sweden, 1)),
                              ")")) %>% 
  select(sex, care, 
         std.difference,
         std.pr_Japan, 
         std.pr_Sweden) %>% 
  stylable()
