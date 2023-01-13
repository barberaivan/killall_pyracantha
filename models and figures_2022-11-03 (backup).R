# En caso de no usar proyectos de RStudio, setear directorio de trabajo,
# idealmente en la ruta donde están los datos. Por ej, correr algo así:
# setwd("home/daniel/papers/pyracantha")

# Packages ----------------------------------------------------------------

library(tidyverse)  # imports several pakages: 
                    # ggplot for plots, 
                    # tidyr for pivot_longer, and
                    # magritr for pipe operator %>%
library(plyr)       # for revalue
library(car)        # for Anova() -type 2-
library(viridis)    # color-blind friendly colour palettes
library(MASS)       # glm.nb (negative binomial glm)
library(VGAM)       # zero truncated negative binomial distribution (frutos)
library(countreg)   # para simular de la negbin truncada y calcular media y var.

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
        size = 16,                #set font size
        #face = 'bold',            #bold typeface
        hjust = -0.1,                #left align
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


# Survival as a function of year by treatment -----------------------------

glm1 <- glm(survival ~ year * trat, family = binomial, data = dlong)
summary(glm1)
Anova(glm1, type = "II") # use this, not anova()

# get proportions to plot alongside model predictions
dlong_agg <- aggregate(survival ~ year + trat, dlong, mean)

# make predictions
p1_data <- expand.grid(year = seq(1, 5, length.out = 100),
                       trat = levels(dlong$trat),
                       label = "Promedio estimado e\nIC del 95 %")

# predictions at linear predictor scale
p1_pred <- predict(glm1, p1_data, se.fit = TRUE)
p1_data$p_mle <- plogis(p1_pred$fit)
p1_data$p_lower <- plogis(p1_pred$fit - qnorm(0.975) * p1_pred$se.fit)
p1_data$p_upper <- plogis(p1_pred$fit + qnorm(0.975) * p1_pred$se.fit)

# plot

dlong_agg$label = "Promedio observado"

sup_marg <-
ggplot(p1_data, aes(x = year, y = p_mle, ymin = p_lower, ymax = p_upper,
                    color = label, fill = label)) +
  # predictions
  geom_ribbon(alpha = 0.3, color = NA) +
  geom_line() +
  # data 
  geom_point(data = dlong_agg, mapping = aes(x = year, y = survival,
                                             shape = label), 
             size = 3, inherit.aes = F) + 
  # set graphic parameters
  scale_colour_manual(values = c("black")) +
  scale_fill_manual(values = c("black")) +
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  
  ylim(0, 1) + 
  xlab("Años poscorte inicial") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"))
sup_marg


# Survival as a function of year and initial height by treatment ----------


glm2 <- glm(survival ~ year * trat + alt0 * trat + alt0 * year, 
            family = binomial, data = dlong)
summary(glm2)
Anova(glm2, type = "II") # use this, not anova()

# get proportions to plot alongside model predictions
# alt0 has to be binned
dlong_agg2 <- aggregate(cbind(survival, alt0) ~ year + trat + alt_cat, 
                        dlong, mean)
alt_means <- aggregate(alt0 ~ alt_cat, d, mean)

# get N
aggregate(alt0 ~ year + trat + alt_cat, dlong, length)

# make predictions
p2_data <- expand.grid(year = seq(1, 5, length.out = 100),
                       trat = unique(dlong$trat),
                       alt0 = alt_means$alt0,
                       agg_type = "Según altura inicial")

# predictions at linear predictor scale
p2_pred <- predict(glm2, p2_data, se.fit = TRUE)
p2_data$p_mle <- plogis(p2_pred$fit)
p2_data$p_lower <- plogis(p2_pred$fit - qnorm(0.975) * p2_pred$se.fit)
p2_data$p_upper <- plogis(p2_pred$fit + qnorm(0.975) * p2_pred$se.fit)

p2_data_2 <- left_join(p2_data, alt_means, by = "alt0")


# plot
sup_alt <- 
ggplot(p2_data_2, aes(x = year, y = p_mle, ymin = p_lower, ymax = p_upper,
                      colour = alt_cat, fill = alt_cat, group = alt_cat)) +
  # predictions
  geom_ribbon(alpha = 0.3, color = NA) +
  geom_line() +
  # data 
  geom_point(data = dlong_agg2, mapping = aes(x = year, y = survival, 
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
sup_alt

# Height shows a small positive effect if it's included without interactions,
# but its effect changes if an interaction with treatmeant or treatment and year
# are included.
# Does it make sense? With interactions, on treat a, shorter trees show higher 
# survival, and the opposite occurs in b and c.



# Survival as a function of year and position by treatment ----------------


glm3 <- glm(survival ~ year * trat + dom * trat, 
            family = binomial, data = dlong)
summary(glm3)
Anova(glm3, type = "II") # use this, not anova()

# get proportions to plot alongside model predictions
dlong_agg3 <- aggregate(survival ~ year + trat + dom, dlong, mean)


# make predictions
p3_data <- expand.grid(year = seq(1, 5, length.out = 100),
                       trat = unique(dlong$trat),
                       dom = unique(dlong$dom))

# predictions at linear predictor scale
p3_pred <- predict(glm3, p3_data, se.fit = TRUE)
p3_data$p_mle <- plogis(p3_pred$fit)
p3_data$p_lower <- plogis(p3_pred$fit - qnorm(0.975) * p3_pred$se.fit)
p3_data$p_upper <- plogis(p3_pred$fit + qnorm(0.975) * p3_pred$se.fit)

# plot
sup_dom <- 
ggplot(p3_data, aes(x = year, y = p_mle, ymin = p_lower, ymax = p_upper,
                    colour = dom, fill = dom, group = dom)) +
  # predictions
  geom_ribbon(alpha = 0.3, color = NA) +
  geom_line() +
  # data 
  geom_point(data = dlong_agg3, mapping = aes(x = year, y = survival, 
                                              colour = dom, shape = dom), 
             size = 3, inherit.aes = F) + 
  
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "D", end = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "D", end = 0.6) +
  
  facet_wrap(vars(trat), nrow = 1) +
  ylim(0, 1) + 
  xlab("Años poscorte inicial") + 
  ylab("Supervivencia") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())
sup_dom
# Clear effect of position: dominant shrubs have higher survival. 


# Merge survival plots ----------------------------------------------------

surv_plot <- egg::ggarrange(
  sup_marg, sup_alt, sup_dom,
  ncol = 1
  )

               # labels = c("A", "B", "C"),
               # label.args = list(gp = grid::gpar(font = 1, cex = 1),
               #                   hjust = 0.3))
ggsave("plots/survival.png", surv_plot, 
       width = 17, height = 13, 
       units = "cm")


# Modelos con alt0 y dom?
glm_sup_dom <- glm(survival ~ year * trat + dom, 
                   family = binomial, data = dlong)

glm_sup_dom_alt <- glm(survival ~ year * trat + dom + 
                                  alt0 * year + alt0 * trat, 
                       family = binomial, data = dlong)


summary(glm_sup_dom)
summary(glm_sup_dom_alt) 
# El coef de dominancia da más grande cuando se considera la alt inicial.
view(cov2cor(vcov(glm_sup_dom_alt)))
# Correlación muy baja entre subdominantes y los demás parámetros
















# Análisis de altura (marginal, por altura inicial y por dominancia -------


altcols <- which(names(d) %in% c("alt0", "alt1", "alt2", "alt3", "alt4", "alt5"))
altlong <- pivot_longer(d, cols = all_of(altcols), 
                        values_to = "alt", names_to = "year_alt")
altlong <- left_join(altlong, d[, c("ID", "alt0")], by = "ID")
altlong$year_alt <- factor(altlong$year_alt, 
                           levels = c("alt0", "alt1", "alt2", "alt3", "alt4", "alt5"),
                           labels = 0:5)
# plot(alt ~ alt0, data = altlong) # OK

altlong_sub <- filter(altlong, year_alt != "alt0")



# Plots de altura

alt_marg <-
  ggplot(altlong, aes(x = year_alt, y = alt)) +
  geom_boxplot(fill = "black", alpha = 0.3, width = 0.4) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte inicial") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"))
alt_marg

alt_alt <-
  ggplot(altlong, aes(x = year_alt, y = alt, 
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
  xlab("Años poscorte inicial") + 
  ylab("Altura (cm)") +
  theme(panel.grid.minor = element_blank(),
        # axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"))
alt_alt


alt_dom <-
  ggplot(altlong, aes(x = year_alt, y = alt, 
                      colour = dom, fill = dom)) +
  geom_boxplot(alpha = 0.3,
               position = position_dodge2()) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "D", end = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "D", end = 0.6) +
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte inicial") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"))
alt_dom

# Merge plots de altura

alt_plot <- egg::ggarrange(
  alt_marg, alt_alt, alt_dom,
  ncol = 1
)

ggsave("plots/altura.png", alt_plot, 
       width = 17, height = 13, 
       units = "cm")



# modelos de altura

# marginal
glm_alt_marg <- glm(alt ~ year_alt * trat, 
                    family = Gamma(link = "log"), 
                    data = altlong_sub)
Anova(glm_alt_marg, type = "II") # todo remil signif.

# con alt inicial
glm_alt_alt <- glm(alt ~ year_alt * trat + 
                         alt0 * year_alt + alt0 * trat, 
                    family = Gamma(link = "log"), 
                    data = altlong_sub)
Anova(glm_alt_alt, type = "II") # todo remil signif.

# con dominancia
glm_alt_dom <- glm(alt ~ year_alt * trat + 
                          dom * year_alt + dom * trat, 
                    family = Gamma(link = "log"), 
                    data = altlong_sub)
Anova(glm_alt_dom, type = "II") # todo remil signif, salvo año * dom




# Análisis de frutos ------------------------------------------------------

frut_names <- c("frutos0", "frutos1", "frutos2", "frutos3", "frutos4", "frutos5")
frutcols <- which(names(d) %in% frut_names)
frutlong <- pivot_longer(d, cols = all_of(frutcols), 
                         values_to = "frutos", names_to = "year_frut")
frutlong$year_frut <- factor(frutlong$year_frut, 
                           levels = frut_names,
                           labels = 0:5)


# Plots de frutos

frut_marg <-
  ggplot(frutlong, aes(x = year_frut, y = frutos + 1)) +
  geom_boxplot(fill = "black", alpha = 0.3, width = 0.4) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte inicial") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm")) #+ 
  # scale_y_continuous(trans = "log10")
frut_marg

frut_alt <-
  ggplot(frutlong, aes(x = year_frut, y = frutos, 
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
  xlab("Años poscorte inicial") + 
  ylab("Número de frutos") +
  theme(panel.grid.minor = element_blank(),
        # axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"))
frut_alt


frut_dom <-
  ggplot(frutlong, aes(x = year_frut, y = frutos, 
                      colour = dom, fill = dom)) +
  geom_boxplot(alpha = 0.3,
               position = position_dodge2()) +
  # geom_jitter(alpha = 0.5, width = 0.1) + 
  # set graphic parameters
  scale_colour_viridis(discrete = TRUE, option = "D", end = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "D", end = 0.6) +
  # facet_grid(cols = vars(trat), rows = vars(agg_type)) +
  facet_wrap(vars(trat), nrow = 1) +
  xlab("Años poscorte inicial") + 
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5, "mm"))
frut_dom

# Merge plots de altura

frut_plot <- egg::ggarrange(
  frut_marg, frut_alt, frut_dom,
  ncol = 1
)

ggsave("plots/frutos.png", frut_plot, 
       width = 17, height = 13, 
       units = "cm")



# modelos de frutos: hurdle con neg bin truncada

frutlong_sub <- filter(frutlong, year_frut != "0")
frutlong_sub <- frutlong_sub[!is.na(frutlong_sub$frutos), ]
frutlong_sub$frutos_bin <- as.numeric(frutlong_sub$frutos > 0)
with_fruits <- frutlong_sub$frutos_bin == 1

data_count <- frutlong_sub[with_fruits, ]
data_count$trat <- factor(as.character(data_count$trat), 
                          levels = unique(data_count$trat))
data_count$year_frut <- factor(as.character(data_count$year_frut), 
                          levels = unique(data_count$year_frut))

# marginal

# probability
glm_frut_marg_p <- glm(frutos_bin ~ year_frut * trat, 
                       data = frutlong_sub, family = "binomial")
Anova(glm_frut_marg_p, type = "II") # all significant
# number
glm_frut_marg_n <- vglm(frutos ~ year_frut + trat, 
                        data = data_count, 
                        family = posnegbinomial)
# poniendo interacción decía que no era full rank, porque el año 3 solo tiene 
# frutos en corte a los 4 años
aggregate(frutos ~ trat + year_frut, data = data_count,
          length)
# los N estaban bien
Anova(glm_frut_marg_n, type = "II") # all significant



# con altura inicial

# probability
glm_frut_alt_p <- glm(frutos_bin ~ year_frut * trat + 
                                    alt0 * year_frut + alt0 * trat, 
                       data = frutlong_sub, family = "binomial")
Anova(glm_frut_alt_p, type = "II") # la altura inicial no afecta la p de dar frutos
# number
glm_frut_alt_n <- vglm(frutos ~ year_frut + trat + alt0, 
                       data = data_count, 
                       family = posnegbinomial)
Anova(glm_frut_alt_n, type = "II") # la altura da marginalmente signif,
summary(glm_frut_alt_n)            # con efecto positivo. Tiene sentido.


# con altura dominancia

# probability
glm_frut_dom_p <- glm(frutos_bin ~ year_frut * trat + 
                        dom * year_frut + dom * trat, 
                      data = frutlong_sub, family = "binomial")
Anova(glm_frut_dom_p, type = "II") # sí, hay eff de dom
# number
glm_frut_dom_n <- vglm(frutos ~ year_frut + trat + dom, 
                       data = data_count, 
                       family = posnegbinomial)
Anova(glm_frut_dom_n, type = "II") # la dom da signif,
summary(glm_frut_dom_n)            # Los subdominantes producen menos!




# Me da fiaca explicar esto y mostrar tantos p-valores,
# quizas podemos decir que solo analizamos el nro de frutos porque lo otro 
# mostraba patrones clarísimos.







# COSAS VIEJAS, LA MAYORÍA ES BASURA PERO REVISAR -------------------------



# Initial height and position  --------------------------------------------

boxplot(alt0 ~ dom, data = d)
# are dominant shrubs shorter? weird


# Crecimiento del primer rebrote en funcion del tamaño inicial -------------

d$growth1 <- d$alt1 - d$tocon0
glm5 <- lm(growth1 ~ alt0, data = d)
Anova(glm5, type = "II")
summary(glm5)
# plot(glm5)
plot(growth1 ~ alt0, d)


# Producción de frutos total ~ trat ---------------------------------------

d$frut_sum <- rowSums(as.matrix(d[, c("frutos1", "frutos2", "frutos3", "frutos4", "frutos5")]), 
                      na.rm = TRUE)
names(d)
glm8 <- glm.nb(frut_sum ~ trat - 1, data = d)

unique(fitted(glm8))
summary(glm8)
Anova(glm8, type = "II")

# frutos en año 4 para trat c en funcion del tamaño inicial
glm6 <- MASS::glm.nb(frutos4 ~ alt0, data = d[d$trat == "c", ])
summary(glm6)
# plot(glm6)
plot(frutos4 ~ alt0, d[d$trat == "c", ])



# Tiempo de corte (suma) ~ trat -------------------------------------------


d$t_resprout <- rowSums(as.matrix(d[, c("t1", "t2", "t3", "t4", "t5")]), 
                        na.rm = TRUE)

# replace 0s because they are forbidden in a Gamma distribution
d$t_resprout[d$t_resprout == 0] <- 1e-5


glm7 <- glm(t_resprout ~ trat, data = d, family = Gamma(link = "log"))
summary(glm7)
Anova(glm7, type = "II")

# promedios por año:
# (a, b , c)
fitted(glm7) %>% unique %>% "/"(4)



# Tiempo de corte en función del tamaño -----------------------------------

timecols <- which(names(d) %in% c("t0", "t1", "t2", "t3", "t4", "t5"))
timelong <- pivot_longer(d, cols = timecols, 
                         values_to = "time", names_to = "year")

altcols <- which(names(d) %in% c("alt0", "alt1", "alt2", "alt3", "alt4", "alt5"))
altlong <- pivot_longer(d, cols = altcols, 
                        values_to = "alt", names_to = "year_alt")

timelong$alt <- altlong$alt
timelong$type <- NA
timelong$type[timelong$year == "t0"] <- "Adulto"
timelong$type[timelong$year != "t0"] <- "Rebrote"

ggplot(timelong, mapping= aes(x = alt, y = time, colour = type,
                              fill = type)) +
  geom_point(alpha = 0.4) + 
  geom_smooth(method = "glm", 
              method.args = list(family = Gamma(link = "log")),
              size = 0.5) + 
  scale_colour_viridis(discrete = TRUE, option = "A", end = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.6) +
  theme(panel.grid.minor = element_blank()) +
  xlab("Altura del arbusto (cm)") +
  ylab("Tiempo de corte (min)")
   
glm7.5 <- glm(time ~ alt * type, data = timelong, family = Gamma(link = "log"))
Anova(glm7.5, type = "II")


# altura del primer corte -------------------------------------------------


hist(d$tocon0)
summary(d$tocon0)
