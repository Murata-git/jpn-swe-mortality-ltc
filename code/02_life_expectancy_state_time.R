
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

trans <- dtms(transient=c("N", "H", "C"),
              absorbing="D",
              timestep=1/4,
              timescale=1:120) 

probs <- probs_all %>%
  filter(age >= 75) %>% 
  {.->>restore} %>% 
  mutate(age=age*4-299, # Multiply 4 to make it integer  
         P=Cases/Tpop,
         se=sqrt(P*(1-P)/Tpop),
         across(c(status_start, status_end), ~case_when(.x==1 ~ "N",
                                                        .x==2 ~ "H",
                                                        .x==3 ~ "C",
                                                        .x==4 ~ "D") %>% 
                  factor(levels=c("N", "H", "C", "D"))),
         from=paste0(status_start, "_", age),
         to=if_else(status_end=="D", status_end, paste0(status_end, "_", age+1))) %>%
  select(bootid, country, sex, from, time=age, to, P) %>% 
  group_by(bootid, country, sex) %>% 
  nest()

probs$data <- map(probs$data, function(.x) {
  class(.x) <- c("dtms_probs", "data.frame")
  return(.x)
})

probs <- probs %>% 
  mutate(t=map(.x=data, .f=~dtms_matrix(dtms=trans,
                                        probs=.x)))

probs$t <- map(probs$t, function(.x) {
  .x[is.nan(.x)] = 0
  return(.x)
})

S <- restore %>% 
  filter(age==75) %>% 
  distinct(bootid, country, sex, status_start, .keep_all=TRUE) %>% 
  mutate(total=sum(Tpop), .by=c(bootid, country, sex)) %>% 
  mutate(status_start=case_when(status_start==1 ~ "N",
                                status_start==2 ~ "H",
                                status_start==3 ~ "C"),
         age=age*4-299,
         prob=Tpop/total) %>% 
  pivot_wider(id_cols=c(bootid, country, sex),
              names_from=status_start,
              values_from=prob) %>% 
  rename_with(~paste0(.x, "_1"), c("N", "H", "C")) %>% 
  group_by(bootid, country, sex) %>% 
  nest() %>% 
  mutate(data=map(.x=data, .f=~as.numeric(.x)))

S$data <- map(S$data, 
              function(.x){names(.x) <- c("N_1", "H_1", "C_1")
              return(.x)})

table4fig <- probs %>% 
  select(-data) %>% 
  left_join(S) %>% 
  mutate(start_age=75) %>% 
  mutate(expectancy=map2(.x=t, .y=data, 
                         .f=~dtms_expectancy(dtms=trans,
                                             matrix=.x,
                                             start_distr=.y) %>% 
                           as.data.frame() %>% 
                           rownames_to_column(var = "basestate"))) %>% 
  unnest(expectancy) %>% 
  select(-c(t, data)) %>% 
  pivot_longer(cols=c(N, H, C),
               names_to="statetime") %>% 
  mutate(statetime=factor(statetime,
                          levels=c("N", "H", "C"),
                          labels=c("No care", "Home care", "Care home")), 
         statetime= fct_rev(statetime),
         prop=value/TOTAL,
         basestate=factor(basestate, 
                          levels=c("start:N_1", "start:H_1", "start:C_1", "AVERAGE"), 
                          labels=c("N", "H", "C", "Total"))) %>% 
  write_rds("intermediate/le_dtms_boottable.rds")

fig_combine <- read_rds("intermediate/le_dtms_boottable.rds") %>% 
  filter(bootid==0) %>% 
  filter(basestate=="Total") %>% 
  mutate(text=round2(value,1)) %>% 
  ggplot(data=., aes(x=start_age, y=value, group=country, fill=statetime)) +
  geom_bar(data=. %>% filter(country=="japan") %>% filter(sex=="women"), aes(x=0.8), position="stack", stat="identity", width=0.2) +
  geom_text(data=. %>% filter(country=="japan") %>% filter(sex=="women"), aes(x=0.8, label = text), 
            size = 3, hjust = 0.5, position = position_stack(vjust = 0.5), size=2.5) +
  geom_bar(data=. %>% filter(country=="sweden") %>% filter(sex=="women"), aes(x=1.2), position="stack", stat="identity", width=0.2) +
  geom_text(data=. %>% filter(country=="sweden") %>% filter(sex=="women"), aes(x=1.2, label = text), 
            size = 3, hjust = 0.5, position = position_stack(vjust = 0.5), size=2.5) +
  geom_bar(data=. %>% filter(country=="japan") %>% filter(sex=="men"), aes(x=1.8), position="stack", stat="identity", width=0.2) +
  geom_text(data=. %>% filter(country=="japan") %>% filter(sex=="men"), aes(x=1.8, label = text), 
            size = 3, hjust = 0.5,position = position_stack(vjust = 0.5), size=2.5) +
  geom_bar(data=. %>% filter(country=="sweden") %>% filter(sex=="men"), aes(x=2.2), position="stack", stat="identity", width=0.2) +
  geom_text(data=. %>% filter(country=="sweden") %>% filter(sex=="men"), aes(x=2.2, label = text), 
            size = 3, hjust = 0.5, position = position_stack(vjust = 0.5), size=2.5) +
  coord_flip(xlim = c(0.7, 2.3),    # set the x limits
             ylim = c(0, 16),  # set the y limits
             clip = 'off')+
  scale_y_continuous(
    expand = c(0,0),
    breaks = seq(0, 16, 2)) +
  scale_x_continuous(breaks = seq(1,2),
                     labels = c("Women                   ", "Men                 "))+
  scale_fill_manual(values = c("No care"="#add8e6", "Home care"="#ffb6c1", "Care home"="#eeee00"),
                    breaks = c("No care", "Home care", "Care home")) +
  labs(y="Remaining life expectancy at age 75", x=NULL, fill=NULL, color=NULL, linetype=NULL, svg=NULL) +
  annotate("text", label="Japan", x = c(0.8, 1.8), y = -2.5, size=4.5, hjust = 0) +
  annotate("text", label="Sweden", x = c(1.2, 2.2), y = -2.5, size=4.5, hjust = 0) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y = element_text(colour = "black", size = 14, vjust = 0.5, hjust = 0, angle = 0),
    axis.text.x = element_text(colour = "black", size = 12, vjust = 0.5, hjust = 0.5, angle = 0),
    strip.text = element_text(size = 12, colour = "black"),
    plot.margin = unit(c(0.4,1,0.4,0.4), "cm"))

# ggpreview(plot=fig_combine, device = "png", width=7, height = 4)

ggsave(filename = "output/fig1_le_statetime_age75.pdf", plot=fig_combine, width=7, height = 4)
