# I think I'm just getting back the simulation + prior

pp <- x %>% 
  filter(term == "treatment", type == effect_type) %>% 
  arrange(postprob) %>% 
  {. ->> xint} %>%
  filter(postprob > 0.5) %>%
  slice(1) %>%
  pull(rep)

xint <- xint %>%
  mutate(row = row_number()) %>%
  filter(rep==pp) %>% 
  pull(row)

power <- x %>% 
  filter(term == "treatment", type == effect_type) %>% 
  mutate(check = ifelse(postprob > .5, 1, 0)) %>% 
  summarise(power = mean(check))

x %>% 
  filter(term == "treatment", type == effect_type) %>%
  ggplot(aes(x=reorder(rep, postprob), y=postprob)) + 
  geom_point(alpha=0.3) + 
  geom_vline(xintercept = xint) + 
  theme_bw() + 
  theme(plot.title = element_text(face="bold"),
        plot.title.position = "plot",
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Simulation index",
       y = paste0("P(effect is > ", min(x$b1_d), ")"),
       title = str_wrap(paste0(round(power*100, 0), "% of 250 simulations resulted in a posterior distribution with probability > 0.5 that the estimate is greater than d=", min(x$b1_d)), 80),
       subtitle = str_wrap(paste0("N = ", min(x$families)), 100)
  )
