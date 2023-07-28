# Ajuste de modelos para supervivencia, altura del rebrote y producción de frutos,
# y elaboración de figuras.

# Estructura inicial de los modelos:
# y ~ year * trat + 
#     alt0 * trat + alt0 * year +
#     dom * trat + dom * year
 
# alt0 = altura inicial del arbusto,
# dom = situación de dominancia,
# trat = tratamiento de corte

# Tras ajustar los modelos, chequear residuos DHARMA en función de 
# año * trat * dom,
# alt0 * trat, 
# sp

# Complejizar los modelos si es necesario.

# --- Figuras ---
# Estructura general: facet_wrap donde los tratamientos van en las columnas.
# En el eje x, el año poscorte, 
# en las filas (3 subfiguras), mostrar datos/predicciones marginales, 
# en función de la altura inicial y en función de la dominancia.

# supervivencia: observados y predichos,
# altura del rebrote y frutos: histogramas de observados.


# Packages ----------------------------------------------------------------

library(tidyverse)  # imports several pakages: 
# ggplot for plots, 
# tidyr for pivot_longer, and
# magritr for pipe operator %>%
library(plyr)       # for revalue
library(car)        # for Anova() -type 2-
library(viridis)    # color-blind friendly colour palettes
library(mgcv)       # betar (para mirar residuos) y rmvn (para simular parametros)
library(VGAM)       # zero truncated negative binomial distribution (frutos)
library(countreg)   # para simular de la negbin truncada y calcular media y var.
# to install countreg:
# install.packages("countreg", repos="http://R-Forge.R-project.org") 
library(DHARMa)     # evaluación de residuos
library(grid)       # textGrob
library(gridExtra)  # grid.arrange
library(ggh4x)      # facet_nested

# Custom ggplot theme -----------------------------------------------------

# from https://rpubs.com/mclaire19/ggplot2-custom-themes

theme_mine <- function() { 
  font <- "Arial"   #assign font family up front
  marg <- 2 # figure margin in mm
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 12,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 1),               
      
      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 12),               
      
      # para separar el eje y de los nros
      axis.title.y = element_text(             
        margin = margin(t = 0, r = 2, b = 0, l = 0, "mm"),
        angle = 90),
      
      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 9),                #font size
      
      # legend.title = element_blank(),
      # legend.position = "bottom",
      legend.text = element_text(size = 9, family = font),
      
      strip.text = element_text(size = 10, family = font, color = "white"),
      strip.text.x = element_text(margin = margin(1.2,0,1.2,0, "mm")), # tamaño de la cajita
      strip.text.y = element_text(margin = margin(0,1.2,0,1.2, "mm")),
      strip.background = element_rect(fill = "gray10", color = "gray10"),
      
      plot.margin = unit(c(marg, marg, marg, marg), "mm")
    )
}

theme_set(theme_mine())



# Data --------------------------------------------------------------------

d <- read.csv("Pyr.csv") # no especifico la ruta porque está en el working directory

# categorize alt0
d$alt_cat <- cut(d$alt0, breaks = 3)
d$alt_cat <- plyr::revalue(as.character(d$alt_cat), replace = c(
  "(99.5,277]" = "[100; 277]", 
  "(277,453]" = "(277; 453]",
  "(453,631]" = "(453; 631]"
))
d$alt_cat <- factor(d$alt_cat, levels = unique(d$alt_cat))

# revalue trat
d$trat <- plyr::revalue(d$trat, replace = c(
  "a" = "Corte anual", "b" = "Corte bienal", "c" = "Corte al cuarto año"
))
d$trat <- factor(d$trat, levels = unique(d$trat))


# revalue dom
d$dom <- plyr::revalue(d$dom, replace = c(
  "dom" = "Dominantes", "sub" = "Subdominantes"
))
d$dom <- factor(d$dom, levels = unique(d$dom))


# longanize dataset
sup_cols <- grep("sup", names(d))
dlong <- pivot_longer(d, cols = sup_cols, 
                      values_to = "survival", names_to = "year")
dlong$year <- revalue(dlong$year,
                      replace = c(
                        "sup1" = "1",
                        "sup2" = "2",
                        "sup3" = "3",
                        "sup4" = "4",
                        "sup5" = "5"
                      )
) %>% as.numeric()

# relación entre alt0 y subdominancia?
plot(alt0 ~ dom, data = d)
summary(lm(alt0 ~ dom, data = d))
# Los subdominantes tienen una media 57 cm mayor. 
aggregate(alt0 ~ dom, data = d, FUN = summary)


# Survival model ------------------------------------------------------------

glm_surv <- glm(survival ~ year * trat + 
                           alt0 * year + alt0 * trat + 
                           dom * year + dom * trat, 
                family = binomial, data = dlong)

# Miramos residuos (deben tener distribución uniforme entre 0 y 1)
surv_res <- simulateResiduals(glm_surv, n = 2000, integerResponse = TRUE)
plot(surv_res) # OK
surv_res_raw <- resid(surv_res)
dlong$surv_res <- surv_res_raw

# en función de trat * año * dom
plotResiduals(surv_res_raw, 
              form = paste(dlong$trat, dlong$year, dlong$dom)) # OK

# en función de trat * año * alt0
ggplot(dlong, aes(x = alt0, y = surv_res)) + 
  geom_smooth(method = "gam", 
              method.args = list(family = betar(eps=.Machine$double.eps*500))) + 
  geom_point() + 
  facet_grid(rows = vars(year), cols = vars(trat)) # Bien

# en función de sp
plotResiduals(surv_res_raw, form = dlong$sp) # OK

# guardamos tabla de anova
aov_surv <- Anova(glm_surv, type = "II") # use this, not anova()
print(aov_surv)
# write.csv(aov_surv, "anova table - supervivencia.csv")


# Survival predictions ----------------------------------------------------


# get proportions to plot alongside model predictions

# marginal
sup_agg_marg <- aggregate(survival ~ year + trat, dlong, mean)
sup_agg_marg$label = "Promedio observado"

# alt0 (has to be binned)
sup_agg_alt <- aggregate(cbind(survival, alt0) ~ year + trat + alt_cat, 
                         dlong, mean)
# aggregate(survival ~ year + trat + alt_cat, dlong, length)
alt_means <- aggregate(alt0 ~ alt_cat, d, mean)

# dom
sup_agg_dom <- aggregate(survival ~ year + trat + dom, dlong, mean)

# Data to make predictions fixing height at its mean
p1_data <- expand.grid(year = seq(1, 5, length.out = 100),
                       trat = levels(dlong$trat),
                       alt0 = mean(d$alt0),
                       dom = levels(dlong$dom),
                       label = "Promedio estimado e\nIC del 95 %")

# Data to make predictions at three heights
p2_data <- expand.grid(year = seq(1, 5, length.out = 100),
                       trat = levels(dlong$trat),
                       alt0 = alt_means$alt0,
                       dom = levels(dlong$dom),
                       label = "Promedio estimado e\nIC del 95 %")

p2_data <- left_join(p2_data, alt_means, by = "alt0")


# predictions at linear predictor scale
p1_pred <- predict(glm_surv, p1_data, se.fit = TRUE)
p2_pred <- predict(glm_surv, p2_data, se.fit = TRUE)


# Predicitions for dominance plot (no simulation)
p1_data$p_mle <- plogis(p1_pred$fit)
p1_data$p_lower <- plogis(p1_pred$fit - qnorm(0.975) * p1_pred$se.fit)
p1_data$p_upper <- plogis(p1_pred$fit + qnorm(0.975) * p1_pred$se.fit)

# plot dominance
surv_dom <- 
  ggplot(p1_data, aes(x = year, y = p_mle, ymin = p_lower, ymax = p_upper,
                      colour = dom, fill = dom, group = dom)) +
  # predictions
  geom_ribbon(alpha = 0.3, color = NA) +
  geom_line() +
  # data 
  geom_point(data = sup_agg_dom, mapping = aes(x = year, y = survival, 
                                              colour = dom, shape = dom), 
             size = 3, inherit.aes = F) + 
  
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "D", end = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "D", end = 0.6) +
  
  facet_wrap(vars(trat), nrow = 1) +
  ylim(0, 1) + 
  xlab("Años poscorte") + 
  ylab("Supervivencia") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())
surv_dom
# Clear effect of position: dominant shrubs have higher survival. 


# Having dom as categorical predictor doesn't allow to fix it at its mean, 
# so we compute predictions for both dominance conditions and then average 
# predictions. Computing CIs is hard, so we achieve it by simulating at the 
# linear predictor scale.

# Simulate predictions
N_sim <- 3e4  

# get model matrices for predictions
p1_data$survival <- p2_data$survival <- 1
mm1 <- model.matrix(glm_surv, data = p1_data)
mm2 <- model.matrix(glm_surv, data = p2_data)

coef_sim <- mgcv::rmvn(N_sim, coef(glm_surv), vcov(glm_surv)) %>% t

p1_sim <- plogis(mm1 %*% coef_sim)
p2_sim <- plogis(mm2 %*% coef_sim)

# Aggregate by dominance

# First, get means only to create the space. Then, replace raw means with 
# weighted averages, based on the dominance distribution at each treatment.

# Get dominance distribution by year * trat, so predictions can match the 
# observed data
dom_dist <- (table(d$trat,d$dom) / rowSums(table(d$trat, d$dom))) #%>% as.matrix() %>% as.data.frame
dom_data <- data.frame(trat = factor(levels(p1_data$trat), levels = levels(p1_data$trat)), 
                       Dominantes = dom_dist[, "Dominantes"],
                       Subdominantes = dom_dist[, "Subdominantes"])

# Dom distribution by trat and alt_cat
d$Dominantes <- as.numeric(d$dom == "Dominantes")
d$Subdominantes <- as.numeric(d$dom == "Subdominantes")

dom_dist <- aggregate(cbind(Dominantes, Subdominantes) ~ trat, d, mean)
dom_dist_alt <- aggregate(cbind(Dominantes, Subdominantes) ~ trat + alt_cat, d, mean)
dom_dist2 <- left_join(dom_dist_alt, alt_means, by = "alt_cat")

p1_data <- left_join(p1_data, dom_dist, by = "trat")

p2_data <- left_join(p2_data, 
                     dom_dist2[, c("Dominantes", "Subdominantes", "alt0", "trat")], 
                     by = c("trat", "alt0"))

# simulations
p1_sim_agg <- p1_sim[p1_data$dom == "Dominantes", ] * 
                p1_data$Dominantes[p1_data$dom == "Dominantes"] +
              p1_sim[p1_data$dom == "Subdominantes", ] * 
                p1_data$Subdominantes[p1_data$dom == "Subdominantes"] 

p2_sim_agg <- p2_sim[p2_data$dom == "Dominantes", ] * 
                p2_data$Dominantes[p2_data$dom == "Dominantes"] +
              p2_sim[p2_data$dom == "Subdominantes", ] * 
                p2_data$Subdominantes[p2_data$dom == "Subdominantes"] 

# mle
p2_data$p_mle <- plogis(p2_pred$fit)

p1_data_agg <- aggregate(p_mle ~ year + trat + label, data = p1_data, mean)
p2_data_agg <- aggregate(p_mle ~ year + trat + alt0 + alt_cat, data = p2_data, mean)

p1_data_agg$p_mle <- p1_data[p1_data$dom == "Dominantes", "p_mle"] * 
                       p1_data$Dominantes[p1_data$dom == "Dominantes"] +
                     p1_data[p1_data$dom == "Subdominantes", "p_mle"] * 
                       p1_data$Subdominantes[p1_data$dom == "Subdominantes"] 

p2_data_agg$p_mle <- p2_data[p2_data$dom == "Dominantes", "p_mle"] * 
                       p2_data$Dominantes[p2_data$dom == "Dominantes"] +
                     p2_data[p2_data$dom == "Subdominantes", "p_mle"] * 
                       p2_data$Subdominantes[p2_data$dom == "Subdominantes"] 

# Get CIs
p1_ci <- apply(p1_sim_agg[, -(1:2)], 1, quantile, probs = c(0.025, 0.975)) %>% t %>% as.data.frame
colnames(p1_ci) <- c("p_lower", "p_upper")

p2_ci <- apply(p2_sim_agg[, -(1:3)], 1, quantile, probs = c(0.025, 0.975)) %>% t %>% as.data.frame
colnames(p2_ci) <- c("p_lower", "p_upper")

# merge data, mle and cis
p1_data_plot <- cbind(p1_data_agg, p1_ci)
p2_data_plot <- cbind(p2_data_agg, p2_ci)


# plot marginal predicions and data
surv_marg <-
  ggplot(p1_data_plot, aes(x = year, y = p_mle, ymin = p_lower, ymax = p_upper,
                      color = label, fill = label)) +
  # predictions
  geom_ribbon(alpha = 0.3, color = NA) +
  geom_line() +
  # data 
  geom_point(data = sup_agg_marg, mapping = aes(x = year, y = survival,
                                             shape = label), 
             size = 3, inherit.aes = F) + 
  # set graphic parameters
  scale_colour_manual(values = c("black")) +
  scale_fill_manual(values = c("black")) +
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  
  ylim(0, 1) + 
  xlab("Años poscorte") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"))
surv_marg


# plot altura 
surv_alt <- 
  ggplot(p2_data_plot, aes(x = year, y = p_mle, ymin = p_lower, ymax = p_upper,
                        colour = alt_cat, fill = alt_cat, group = alt_cat)) +
  # predictions
  geom_ribbon(alpha = 0.3, color = NA) +
  geom_line() +
  # data 
  geom_point(data = sup_agg_alt, mapping = aes(x = year, y = survival, 
                                              colour = alt_cat,
                                              shape = alt_cat), 
             size = 3, inherit.aes = F) + 
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "C", end = 0.6, 
                       name = "Clase de altura\ninicial (cm)") +
  scale_fill_viridis(discrete = TRUE, option = "C", end = 0.6,
                     name = "Clase de altura\ninicial (cm)") +
  scale_shape_manual(values = c(16, 17, 15), 
                     name = "Clase de altura\ninicial (cm)") +
  
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  ylim(0, 1) + 
  ylab("Supervivencia") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        # legend.key.size = unit(4.5, "mm"),
        legend.position = "right",
        legend.title = element_text())
surv_alt


# Merge survival plots 

surv_plot <- egg::ggarrange(
  surv_marg + ggtitle("a"), 
  surv_alt + ggtitle("b"), 
  surv_dom + ggtitle("c"),
  ncol = 1
)

ggsave("plots/survival2.png", surv_plot,
       width = 16, height = 13,
       units = "cm")


# Altura del rebrote modelo ------------------------------------------------

altcols <- which(names(d) %in% c("alt0", "alt1", "alt2", "alt3", "alt4", "alt5"))
altlong <- pivot_longer(d, cols = all_of(altcols), 
                        values_to = "alt", names_to = "year_alt")
altlong <- left_join(altlong, d[, c("ID", "alt0")], by = "ID")
altlong$year <- factor(altlong$year_alt, 
                           levels = c("alt0", "alt1", "alt2", "alt3", "alt4", "alt5"),
                           labels = 0:5)
# plot(alt ~ alt0, data = altlong) # OK

altlong_sub <- filter(altlong, year_alt != "alt0")

# fit model
glm_alt <- glm(alt ~ year * trat + 
                     alt0 * year + alt0 * trat + 
                     dom * year + dom * trat, 
                family = Gamma(link = "log"), data = altlong_sub)

# Miramos residuos de altura. 
alt_res <- simulateResiduals(glm_alt, n = 2000)
plot(alt_res) # bad overall fit

data_alt_res <- glm_alt$data[-glm_alt$na.action, ]
data_alt_res$alt_res <- resid(alt_res)

# en función de trat * año * dom
plotResiduals(alt_res, 
              form = paste(data_alt_res$trat, 
                           data_alt_res$year, 
                           data_alt_res$dom)) # OK

# en función de trat * año * alt0
ggplot(data_alt_res, aes(x = alt0, y = alt_res)) + 
  geom_smooth(method = "gam", 
              method.args = list(family = betar(eps=.Machine$double.eps*500))) + 
  geom_point() + 
  facet_grid(rows = vars(year), cols = vars(trat)) +
  ggtitle("modelo 1, sin cuadratico en altura")# Mal. 
# Casi siempre vemos una relación unimodal entre los residuos y altura

# en función de la especie
plotResiduals(alt_res, form = data_alt_res$sp) # subestima altura de atalantoides

# Modelo con altura cuadrática + interacciones
glm_alt_2 <- glm(alt ~ year * trat + 
                   alt0 * year + alt0 * trat + 
                   I(alt0 ^ 2) * year + I(alt0 ^ 2) * trat +
                   dom * year + dom * trat, 
                 family = Gamma(link = "log"), data = altlong_sub)

alt_res_2 <- simulateResiduals(glm_alt_2, n = 2000)
plot(alt_res_2) # bad overall fit, aunque mejoró un poco
data_alt_res$alt_res_2 <- resid(alt_res_2)

# en función de trat * año * dom
plotResiduals(alt_res_2,
              form = paste(data_alt_res$trat,
                           data_alt_res$year,
                           data_alt_res$dom)) # OK

# en función de altura * trat * año
ggplot(data_alt_res, aes(x = alt0, y = alt_res_2)) + 
  geom_smooth(method = "gam", 
              method.args = list(family = betar(eps=.Machine$double.eps*500))) +  # default is * 100
  geom_point() + 
  facet_grid(rows = vars(year), cols = vars(trat)) +
  ggtitle("modelo 6, cuadratico en alt con interacción")# Mal. 
# (tira warnings porque hay un par de ceros)
# mejoró bastante.

plotResiduals(alt_res_2, form = data_alt_res$sp) # sigue mal. 

# NOTA: modelos mixtos no mejoran el ajuste, ajustan muy mal (ver código
# extendido). La sobreestimación de la varianza para predichos altos mejora un
# poco utilizando un modelo gamma location-scale (con mgcv), en donde se modela 
# la scale en función de la altura y otras covariables. Sin embargo, las 
# predicciones de la media son muy similares, y para simplificar nos quedamos
# con este.

# Guardamos tabla de anova
aov_alt <- Anova(glm_alt_2, type = "II") # use this, not anova()
print(aov_alt)
# write.csv(aov_alt, "anova table - altura.csv")


# Modelo con especie
glm_alt_sp <- glm(alt ~ year * trat + 
                   sp +
                   alt0 * year + alt0 * trat +
                   I(alt0 ^ 2) * year + I(alt0 ^ 2) * trat +
                   dom * year + dom * trat, 
                 family = Gamma(link = "log"), data = altlong_sub)
alt_res_sp <- simulateResiduals(glm_alt_sp, n = 2000)
plot(alt_res_sp) # bad overall fit (sigue flojo)



# Altura del rebrote, predicciones en función de altura inicial -----------

nrep <- 120
altseq <- seq(min(altlong_sub$alt0), max(altlong_sub$alt0), length.out = nrep)

# datos para modelo sin sp
pdata_alt <- expand.grid(
  alt0 = altseq,
  year = unique(altlong_sub$year),
  trat = unique(altlong_sub$trat),
  dom = unique(altlong_sub$dom)
)

# datos para modelo con sp 
pdata_alt_sp <- expand.grid(
  alt0 = altseq,
  year = unique(altlong_sub$year),
  trat = unique(altlong_sub$trat),
  dom = unique(altlong_sub$dom),
  sp = unique(altlong_sub$sp)
)

predict_ic <- function(model, data = pdata_alt, name = "model_0") {
  p <- predict(model, newdata = data, se.fit = TRUE)
  preds <- data.frame(
    p_mle = model$family$linkinv(p$fit),
    p_lower = model$family$linkinv(p$fit - qnorm(0.975) * p$se.fit),
    p_upper = model$family$linkinv(p$fit + qnorm(0.975) * p$se.fit)
  )
  preds$name <- name
  return(cbind(data, preds))
}

# predicciones

alt_pred <- predict_ic(glm_alt_2, data = pdata_alt)
alt_pred_sp <- predict_ic(glm_alt_sp, data = pdata_alt_sp)

# plots

# (sin especie)
ggplot(alt_pred, aes(x = alt0, y = p_mle, ymin = p_lower,
                     ymax = p_upper, colour = dom, fill = dom)) +
  geom_ribbon(colour = NA, alpha = 0.25) +
  geom_line() +
  geom_point(data = altlong_sub, mapping = aes(x = alt0, y = alt, 
                                               colour = dom, fill = dom),
             inherit.aes = F) +
  scale_colour_viridis(discrete = TRUE, option = "D", end = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "D", end = 0.6) +
  facet_grid(rows = vars(year), cols = vars(trat)) +#, scales = "free_y") +
  xlab("Altura inicial (cm)") + 
  ylab("Altura del rebrote (cm)") +
  scale_y_continuous(
    sec.axis = sec_axis(trans= ~ . * 1, name = "Años poscorte")
  ) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank()) 

ggsave("plots/altura_alt0_yfixed.png",       ## OJO, quitar yfixed si se guarda con scales = free_y
       width = 16, height = 15, units = "cm")

# (con especie)
ggplot(alt_pred_sp, 
       aes(x = alt0, y = p_mle, ymin = p_lower, ymax = p_upper, 
           colour = sp, fill = sp)) +
  geom_ribbon(colour = NA, alpha = 0.2) +
  geom_line() +
  geom_point(data = altlong_sub, mapping = aes(x = alt0, y = alt, 
                                               colour = sp, fill = sp),
             inherit.aes = F, shape = 19, alpha = 0.7) +
  scale_colour_viridis(discrete = TRUE, option = "A", end = 0.6,
                       labels = c("*P. angustifolia*", "*P.* aff. *atalantioides*")) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.6,
                     labels = c("*P. angustifolia*", "*P.* aff. *atalantioides*")) +
  facet_nested(rows = vars(year), cols = vars(trat, dom)) +#, 
  #scales = "free_y") + 
  xlab("Altura inicial (cm)") + 
  ylab("Altura del rebrote (cm)") +
  scale_y_continuous(
    sec.axis = sec_axis(trans= ~ . * 1, name = "Años poscorte")
  ) + 
  theme(legend.position = "bottom",
        legend.text = ggtext::element_markdown(), # para usar codigo markdown en la legenda
        legend.title = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank()) #+ 

ggsave("plots/altura_sp_alt0_solo intercept de sp_yfixed.png",
       height = 16, width = 20, units = "cm")

# Altura del rebrote plots (histogramas) -----------------------------------

alt_marg <-
  ggplot(altlong, aes(x = year, y = alt)) +
  geom_boxplot(fill = "black", alpha = 0.3, width = 0.4) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"))
alt_marg

alt_alt <-
  ggplot(altlong, aes(x = year, y = alt, 
                      colour = alt_cat, fill = alt_cat)) +
  geom_boxplot(alpha = 0.3,
               position = position_dodge2()) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "C", end = 0.6, 
                       name = "Clase de altura\ninicial (cm)") +
  scale_fill_viridis(discrete = TRUE, option = "C", end = 0.6,
                     name = "Clase de altura\ninicial (cm)") +
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte") + 
  ylab("Altura (cm)") +
  theme(panel.grid.minor = element_blank(),
        # axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())
alt_alt


alt_dom <-
  ggplot(altlong, aes(x = year, y = alt, 
                      colour = dom, fill = dom)) +
  geom_boxplot(alpha = 0.3,
               position = position_dodge2()) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "D", end = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "D", end = 0.6) +
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())
alt_dom

alt_sp <-
  ggplot(altlong, aes(x = year, y = alt, 
                      colour = sp, fill = sp)) +
  geom_boxplot(alpha = 0.3,
               position = position_dodge2()) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "A", end = 0.6,
                       labels = c("P. angustifolia", "P. atalantoides")) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.6,
                     labels = c("P. angustifolia", "P. atalantoides")) +
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())
alt_sp

# Merge plots de altura

# sin especie

alt_plot <- egg::ggarrange(
  alt_marg + ggtitle("a"), 
  alt_alt + ggtitle("b"),
  alt_dom + ggtitle("c"),
  ncol = 1
)

ggsave("plots/altura.png", alt_plot,
       width = 16, height = 13,
       units = "cm")


# con especie

alt_plot_2 <- egg::ggarrange(
  alt_marg + ggtitle("a"), 
  alt_alt + ggtitle("b") + theme(axis.title.y = element_blank()),
  alt_dom + ggtitle("c") + theme(axis.title.x = element_blank()),
  alt_sp + ggtitle("d"),
  ncol = 1
)

eje_y <- textGrob("Altura (cm)", 
                  gp = gpar(angle = 90, fontsize = 12), 
                  rot = 90)
alt_plot_sp <-
grid.arrange(eje_y, alt_plot_2, ncol = 2, widths = c(1, 25))
plot(alt_plot_sp)

# ggsave("plots/altura_sp.png", alt_plot_sp,
#        width = 17, height = 16,
#        units = "cm")

# plot de altura solo

alt_sp_only <- alt_sp + ylab("Altura (cm)") + theme(axis.title.y = element_text())
# ggsave("plots/altura_sp_only.png", alt_sp_only, 
#        width = 17, height = 6, units = "cm")

# Frutos modelos -----------------------------------------------------------

# hurdle con neg bin truncada
# (preparar datos para modelos separados)

frut_names <- c("frutos0", "frutos1", "frutos2", "frutos3", "frutos4", "frutos5")
frutcols <- which(names(d) %in% frut_names)
frutlong <- pivot_longer(d, cols = all_of(frutcols), 
                         values_to = "frutos", names_to = "year_frut")
frutlong$year <- factor(frutlong$year_frut, 
                        levels = frut_names,
                        labels = 0:5)

frutlong_sub <- filter(frutlong, year != "0")
frutlong_sub <- frutlong_sub[!is.na(frutlong_sub$frutos), ]
frutlong_sub$frutos_bin <- as.numeric(frutlong_sub$frutos > 0)
with_fruits <- frutlong_sub$frutos_bin == 1

data_count <- frutlong_sub[with_fruits, ]
data_count$trat <- factor(as.character(data_count$trat), 
                          levels = unique(data_count$trat))
data_count$year <- factor(as.character(data_count$year), 
                          levels = unique(data_count$year))


# Fruit production probability (reg logística)
glm_frut_p <- glm(frutos_bin ~ year * trat + 
                    alt0 * trat + alt0 * year + 
                    dom * trat + dom * year, 
                  data = frutlong_sub, family = "binomial")

# Miramos residuos 
frut_p_res <- simulateResiduals(glm_frut_p, n = 2000, integerResponse = TRUE)
plot(frut_p_res) # perfect overall fit

plotResiduals(frut_p_res, form = glm_frut_p$data$alt0) # OK
plotResiduals(frut_p_res, 
              form = paste(glm_frut_p$data$dom, 
                           glm_frut_p$data$year)) # OK

# en función de trat * año * alt0
glm_frut_p$data$frutp_res <- frut_p_res$scaledResiduals

ggplot(glm_frut_p$data, aes(x = alt0, y = frutp_res)) + 
  geom_smooth(method = "gam", 
              method.args = list(family = betar(eps=.Machine$double.eps*500))) + 
  geom_point() + 
  facet_grid(rows = vars(year), cols = vars(trat)) # Bien

# en función de sp
plotResiduals(frut_p_res, form = glm_frut_p$data$sp) # OK

# Guardamos tabla de anova
aov_frut_p <- Anova(glm_frut_p, type = "II") # all significant
print(aov_frut_p)
# write.csv(aov_frut_p, "anova table - frutos (probabilidad).csv")



# Fruit count (positive neg binomial)

data_count_usecols <- c("frutos", "trat", "year", "alt0", "dom")
data_count_notna <- data_count[complete.cases(data_count[, data_count_usecols]), ]
glm_frut_n <- vglm(frutos ~ year + trat + alt0 + dom, 
                   data = data_count_notna, 
                   family = posnegbinomial)

# chequear residuos (hay que simular a mano)
# get linear predictors
fitted_lp <- predict(glm_frut_n, se.fit = TRUE) 
mu_nb <- exp(fitted_lp$fitted.values[, 1]) # mean for non-truncated NB distribution
k_param <- exp(fitted_lp$fitted.values[1, 2]) # dispersion parameter, = theta = size
# (it's a constant)

# Simulate data
set.seed(34234)
nsim <- 2000
nobs <- length(mu_nb)
y_sim <- matrix(NA, nobs, nsim)
for(i in 1:nobs) {
  y_sim[i, ] <- rztnbinom(n = nsim, mu = mu_nb[i], theta = k_param)   
}

# Get DHARMa residuals
frut_res_n <- createDHARMa(simulatedResponse = y_sim, 
                           observedResponse = data_count_notna$frutos, 
                           integerResponse = T)
plot(frut_res_n, quantreg = T) # OK

plotResiduals(frut_res_n, form = data_count_notna$alt0) # OK
plotResiduals(frut_res_n, form = data_count_notna$trat) # OK
plotResiduals(frut_res_n, form = data_count_notna$dom) # OK
plotResiduals(frut_res_n, form = data_count_notna$year) # maso pero zafa
plotResiduals(frut_res_n, form = data_count_notna$sp) # OK


# Miramos coeficientes
summary(glm_frut_n)
#                              Estimate Std. Error z value Pr(>|z|)    
#   (Intercept):1            2.483298   0.414940   5.985 2.17e-09 ***
#   (Intercept):2           -0.287587   0.128239  -2.243 0.024923 *  
#   year4                    1.786840   0.304257   5.873 4.29e-09 ***
#   year3                    1.517473   0.356920   4.252 2.12e-05 ***
#   tratCorte al cuarto año  0.993762   0.289605   3.431 0.000600 ***
#   alt0                     0.003860   0.001159   3.332 0.000863 ***
#   domSubdominantes        -0.922804   0.258437  -3.571 0.000356 ***


# Guardamos tabla de anova
aov_frut_n <- Anova(glm_frut_n, type = "II") # all significant
print(aov_frut_n) 
# write.csv(aov_frut_n, "anova table - frutos (conteo).csv")



# Frutos plots ------------------------------------------------------------


frut_marg <-
  ggplot(frutlong, aes(x = year, y = frutos + 1)) +
  geom_boxplot(fill = "black", alpha = 0.3, width = 0.4) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm")) #+ 
# scale_y_continuous(trans = "log10")
frut_marg

frut_alt <-
  ggplot(frutlong, aes(x = year, y = frutos, 
                       colour = alt_cat, fill = alt_cat)) +
  geom_boxplot(alpha = 0.3,
               position = position_dodge2()) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "C", end = 0.6, 
                       name = "Clase de altura\ninicial (cm)") +
  scale_fill_viridis(discrete = TRUE, option = "C", end = 0.6,
                     name = "Clase de altura\ninicial (cm)") +
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte") + 
  ylab("Número de frutos") +
  theme(panel.grid.minor = element_blank(),
        # axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())
frut_alt


frut_dom <-
  ggplot(frutlong, aes(x = year, y = frutos, 
                       colour = dom, fill = dom)) +
  geom_boxplot(alpha = 0.3,
               position = position_dodge2()) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "D", end = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "D", end = 0.6) +
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())
frut_dom

# Merge plots de altura

frut_plot <- egg::ggarrange(
  frut_marg + ggtitle("a"),
  frut_alt + ggtitle("b"),
  frut_dom + ggtitle("c"),
  ncol = 1
)

ggsave("plots/frutos.png", frut_plot,
       width = 16, height = 13,
       units = "cm")


# Modelos incluyendo especie - (para tabla de anovas) ----------------------

# supervivencia
glm_surv_sp <- glm(survival ~ year * trat + 
                  alt0 * year + alt0 * trat + 
                  dom * year + dom * trat + 
                  sp, 
                family = binomial, data = dlong)
aov_surv_sp <- Anova(glm_surv_sp, type = "II") 
print(aov_surv_sp)
# write.csv(aov_surv_sp, "anova table species - supervivencia.csv")

# altura
glm_alt_sp <- glm(alt ~ year * trat + 
                  alt0 * year + alt0 * trat + 
                  I(alt0 ^ 2) * year + I(alt0 ^ 2) * trat + 
                  dom * year + dom * trat +
                  sp, 
                 family = Gamma(link = "log"), data = altlong_sub)
aov_alt_sp <- Anova(glm_alt_sp, type = "II") # use this, not anova()
summary(glm_alt_sp)
print(aov_alt_sp)
# write.csv(aov_alt_sp, "anova table species - altura.csv")
# La especie da un efecto significativo en altura del rebrote


# producción de frutos

# Fruit production probability (reg logística)
glm_frut_p_sp <- glm(frutos_bin ~ year * trat + 
                    alt0 * trat + alt0 * year + 
                    dom * trat + dom * year + 
                    sp, 
                  data = frutlong_sub, family = "binomial")
aov_frut_p_sp <- Anova(glm_frut_p_sp, type = "II") # all significant
print(aov_frut_p_sp)
# write.csv(aov_frut_p_sp, "anova table species - frutos (probabilidad).csv")

# Fruit count (positive neg binomial)
glm_frut_n_sp <- vglm(frutos ~ year + trat + alt0 + dom +
                     sp, 
                   data = data_count, 
                   family = posnegbinomial)
glm_frut_table_n_sp <- summary(glm_frut_n)
aov_frut_n_sp <- Anova(glm_frut_n_sp, type = "II")
print(aov_frut_n_sp)
# write.csv(aov_frut_n_sp, "anova table species - frutos (conteo).csv")


