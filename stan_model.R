my_data <- readr::read_rds("rds/my_data.rds")
m2 <- readr::read_rds("rds/jags_par.rds")

library(cmdstanr)

m2_stan <- cmdstan_model("STAN/integrated_model_simplified.stan",
                         cpp_options = list(stan_threads = TRUE))


stanm2_fit <- m2_stan$sample(data = my_data, 
                             chains = 4, threads_per_chain = 2,
                             parallel_chains = 4,
                             iter_warmup = 1500,
                             iter_sampling = 2000)

dras <- stanm2_fit$draws()
bayesplot::mcmc_trace(dras,pars = c("a"))

st_full <- stanm2_fit$summary()
stsum <- stanm2_fit$summary(variables = c("beta", "cc", "a", "zsum"))

stsum

# stsuma <- stanm2_fit$summary(variables = "z_obs")
dplyr::filter(stsum, grepl(c("beta|cc|^a|zsum"), variable))
dplyr::filter(st_full, grepl("psi", variable)) %>% pull(median) %>% rethinking::dens()
dplyr::filter(st_full, grepl("z\\[", variable)) %>% pull(median) %>% table

r <- agg_plane$x
values(r) <-    dplyr::filter(st_full, grepl("z\\[", variable)) %>% pull(median) 
plot(r)
calc(agg_plane$x, dplyr::filter(st_full, grepl("z\\[", variable)) %>% pull(median) )



readr::write_rds(stsum, "stan_par.rds")
theme_set(theme_dark())
diff1 <- 
  m2 %>% as_tibble(rownames = "variable") %>% janitor::clean_names() %>% 
  mutate(software="jags",q5 = lower95, q95=upper95,
         variable = stringr::str_replace(variable, "^a","a[1]" )) %>% 
  bind_rows(stsum %>% mutate(software='stan')) 

p1 <- diff1 %>% 
  filter(variable!="zsum")  %>% 
  ggplot(aes(variable, median, colour= software, shape = software)) + 
  geom_pointrange(aes(ymin=q5, ymax=q95), position = position_dodge(width = 0.2) )  +
  # facet_grid(variable~., scales='free') +
  scale_colour_viridis_d() +
  labs(x = "Variable", y = "Estimate")

p2 <- diff1 %>% 
  filter(variable=="zsum")  %>% 
  ggplot(aes(variable, median, colour= software, shape = software)) + 
  geom_pointrange(aes(ymin=q5, ymax=q95), position = position_dodge(width = 0.2) )  +
  scale_colour_viridis_d()+
  labs(x = "", y = "Estimate")


library(patchwork)
p1+p2 + patchwork::plot_layout(widths = c(0.9,0.1),guides = 'collect' ) +
  plot_annotation()
ggsave("jags_stan_compare.jpg")
