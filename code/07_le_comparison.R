
library(readr)
library(stats)
library(tidyverse)
library(tidylog)
library(GGally)
library(patchwork)
library(scales)
library(DT)
library(readxl)
library(gtsummary)
library(gt)
library(gtools)
library(kableExtra)
library(arrow)
library(parallel)
library(foreach)
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
                             # autoWidth = TRUE,
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


# Life expectancy with dtms package  ---------------------------------------------

probs_all <- read_rds("input/aggregted_bootdata.rds")

# For loop
# i=1
age_list <- seq(75, 120 ,0.25)

combinedata <- tibble()

for (i in 101:1){
  
  trans <- dtms(transient=c("N", "H", "C"),
                absorbing="D",
                timestep=1/4,
                timescale=1:(120-i+1))
  
  probs <- probs_all %>%
    filter(bootid==0) %>% 
    filter(age >= age_list[[i]]) %>% 
    mutate(age=age*4-(298+i), # Multiply 4
           P=Cases/Tpop,
           across(c(status_start, status_end), ~case_when(.x==1 ~ "N",
                                                          .x==2 ~ "H",
                                                          .x==3 ~ "C",
                                                          .x==4 ~ "D") %>%
                    factor(levels=c("N", "H", "C", "D"))),
           from=paste0(status_start, "_", age),
           to=if_else(status_end=="D", status_end, paste0(status_end, "_", age+1))) %>% # used 1, not 0.25
    select(country, sex, from, time=age, to, P) %>% 
    group_by(country, sex) %>% 
    nest()
  
  class(probs$data[[1]]) <- class(probs$data[[2]]) <- class(probs$data[[3]]) <- class(probs$data[[4]]) <- 
    c("dtms_probs", "data.frame")
  
  probs <- probs %>% 
    mutate(t=map(.x=data, .f=~dtms_matrix(dtms=trans,
                                          probs=.x)))
  
  # probs$t[[1]][is.nan(probs$t[[1]])] = 0
  # probs$t[[3]][is.nan(probs$t[[3]])] = 0
  probs$t <- map(probs$t, function(.x) {
    .x[is.nan(.x)] = 0
    return(.x)  # Important: return the modified object
  })
  
  S <- probs_all %>% 
    filter(bootid==0) %>% 
    filter(age==age_list[[i]]) %>% 
    distinct(country, sex, status_start, .keep_all=TRUE) %>% 
    mutate(total=sum(Tpop), .by=c(country, sex)) %>% 
    mutate(status_start=case_when(status_start==1 ~ "N",
                                  status_start==2 ~ "H",
                                  status_start==3 ~ "C"),
           age=age*4-(298+i),
           prob=Tpop/total) %>% 
    pivot_wider(id_cols=c(bootid, country, sex),
                names_from=status_start,
                values_from=prob) %>% 
    rename_with(~paste0(.x, "_1"), c("N", "H", "C")) %>%
    select(-c(bootid)) %>% 
    group_by(country, sex) %>% 
    nest() %>% 
    mutate(data=map(.x=data, .f=~as.numeric(.x)))
  
  names(S$data[[1]]) <- names(S$data[[2]]) <- names(S$data[[3]]) <- names(S$data[[4]]) <- 
    c("N_1", "H_1", "C_1")
  
  combinedata <- combinedata %>% 
    bind_rows(probs %>% 
                select(-data) %>% 
                left_join(S) %>% 
                mutate(start_age=i*0.25+74.75))
  
}

table4fig_total75_100 <- combinedata %>% 
  mutate(expectancy=map2(.x=t, .y=data, 
                         .f=~dtms_expectancy(dtms=trans,
                                             matrix=.x,
                                             start_distr=.y) %>% 
                           as_tibble() %>% 
                           slice(4))) %>% 
  unnest(expectancy) %>% 
  select(-c(t, data)) %>% 
  pivot_longer(cols=c(N, H, C),
               names_to="statetime") %>% 
  mutate(statetime=factor(statetime,
                          levels=c("N", "H", "C")))

total_LE <- table4fig_total75_100 |>
  rename(age=start_age, total=TOTAL) %>% 
  distinct(country, sex, age, .keep_all=TRUE) %>% 
  filter(between(age, 75, 100))

lifetable <- read_table("hmd/lifetable_men_jpn_1x1.txt") |>
  mutate(sex="men") |>
  bind_rows(read_table("hmd/lifetable_women_jpn_1x1.txt") |>
              mutate(sex="women")) |>
  mutate(country="japan") |>
  bind_rows(read_table("hmd/lifetable_men_swe_1x1.txt") |>
              mutate(sex="men") |>
              bind_rows(read_table("hmd/lifetable_women_swe_1x1.txt") |>
                          mutate(sex="women")) |>
              mutate(country="sweden")) |>
  rename(year=Year, age=Age) |>
  filter(between(year, 2017, 2019)) |>
  mutate(age=if_else(age=="110+", 110, as.numeric(age))) |>
  filter(between(age, 75, 100)) |>
  summarise(ex=mean(ex), .by=c(age, sex, country))

fig_lecomparison <- total_LE |>
  select(age, sex, le=total, country) |>
  mutate(group=if_else(country=="japan", "LIFE", "AHC")) |>
  bind_rows(lifetable |>
              select(age, sex, le=ex, country) |>
              mutate(group="Human Mortality Database")) %>%
  mutate(group=factor(group, levels=c("AHC", "LIFE", "Human Mortality Database")),
         sex=factor(sex,
                    levels=c("men", "women"),
                    labels=c("Men", "Women")),
         country=factor(country, levels=c("sweden", "japan"),
                        labels=c("Sweden", "Japan"))) %>%
  ggplot(data=., aes(y=le, x=age, group=group, color=group)) +
  geom_line() +
  scale_y_continuous(limits=c(0, 16), breaks = seq(0, 16, 1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(75, 100, 5), expand = c(0,0)) +
  scale_color_manual(values = c("AHC"="#005293", "LIFE"="#BC002D", "Human Mortality Database"="black")) +
  labs(y="Total life expectancy", x="Age", group=NULL, color=NULL)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  facet_grid2(sex~country, scales = "free", independent = "all")

ggsave(filename = "output/supl.fig9_le_comparison.svg", plot=fig_lecomparison, width=8, height = 8, device = svglite)

diftable <- total_LE |>
  select(age, sex, le=total, country) |>
  distinct(country, sex, age, .keep_all=TRUE) %>% 
  mutate(group=if_else(country=="Japan", "LIFE", "AHC")) |>
  left_join(lifetable |>
              select(age, sex, le_hmd=ex, country),
            by=c("age", "sex", "country")) |>
  drop_na(le_hmd) %>%
  mutate(dif=le-le_hmd)

diftable %>% 
  mutate(abs_dif=abs(dif)) %>% 
  summarise(max_abs_dif=max(abs_dif), .by=c(country, sex))
