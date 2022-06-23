# define plotting theme for ggplots
pt_theme <- function(){
  theme_classic() +
    theme(
      plot.margin = unit(c(1,.3,.3,.3), "cm"))
}



# compute 95% credible interval for binomial glmm
pt_get_effects.glmm_pa <- function(data, model, var_x){
  # general settings
  nsim <- 2000
  
  vars <- names(fixef(model))[-1]
  var_x_orig <- word(var_x, 1, sep = "\\.")

  newdata <- as.data.frame(
    setNames(replicate(length(vars), rep(0, 100), simplify = F), vars))
  newdata[, var_x] <- seq(
    min(data[[var_x]]), max(data[[var_x]]), length.out = 100)
  
  # back-transform z-transformed variables
  sd <- sd(data[[var_x_orig]])
  mu <- mean(data[[var_x_orig]])
  newdata[[var_x_orig]] <- newdata[[var_x]] * sd + mu
  
  # get 95% credible interval
  bsim <- sim(model, n.sim = nsim)
  fmla <- as.formula(paste("nest_pa ~", paste(vars, collapse = "+")))
  Xmat <- model.matrix(fmla[c(1,3)], data = newdata)
  newdata$fit <- plogis(Xmat %*% fixef(model))
  fitmat <- matrix(ncol = nsim, nrow = nrow(newdata))
  for(i in 1:nsim) fitmat[,i] <- plogis(Xmat %*% bsim@fixef[i,])
  newdata$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
  newdata$upr <- apply(fitmat, 1, quantile, prob = 0.975)
  
  return(newdata)
}




# plot vegetation effects from both scales in one plot
pt_plot_effects.pa.2scales <- function(dat_1m, dat_10cm, newdat_1m, newdat_10cm, var_x, lab_x, pval_10cm, pval_1m){
  # to show the plot back-transformed to the original scale
  var_x_orig <- word(var_x, 1, sep = "\\.")
  newdat_1m$x <- newdat_1m[[var_x_orig]]
  dat_1m$x <- dat_1m[[var_x_orig]]
  newdat_10cm$x <- newdat_10cm[[var_x_orig]]
  dat_10cm$x <- dat_10cm[[var_x_orig]]
  
  # add offset along y-axis to avoid overlapping of points from different scale
  dat_10cm <- dat_10cm %>% mutate(
    nest_pa = ifelse(nest_pa == 1, nest_pa + 0.02, nest_pa - 0.02)) %>%
    mutate(scale = "10cm")
  dat_1m <- dat_1m %>% mutate(nest_pa = ifelse(nest_pa == 1, nest_pa + 0.06, nest_pa - 0.06)) %>%
    mutate(scale = "1m")
  dat <- bind_rows(dat_10cm, dat_1m)
  
  # combine newdat
  newdat_10cm <- newdat_10cm %>% mutate(scale = "10cm")
  newdat_1m <- newdat_1m %>% mutate(scale = "1m")
  newdat <- bind_rows(newdat_10cm, newdat_1m)
  
  # general plot settings
  alpha_dots = .5; alpha_fill = .3; size_line = 1; size_point = 1.5
  position = position_jitter(w = 0, h = 0.015)
  col <- c("#0073C2FF", "#CD534CFF")
  
  # plot
  ggplot() +
    geom_ribbon(
      aes(x = x, y = fit, ymin = lwr, ymax = upr, fill = scale),
      newdat, alpha = 0.3) +
    geom_line(aes(x = x, y = fit, color = scale), data = newdat, size = 1) +
    geom_point(
      data = dat,
      mapping = aes(y = nest_pa, x = x, shape = scale, color = scale),
      position = position, size = size_point, alpha = alpha_dots) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    scale_shape_manual(values = c(21, 22), guide = FALSE) +
    scale_color_manual(
      name = "Scale",
      values = col,
      aesthetics = c("colour", "fill"),
      labels = c("1m" = expression(1~m^2), "10cm" = expression(10~cm^2))) +
    scale_linetype_manual(
      values = c(1,2),
      guide = "none") +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = "",
          linetype = c(1, 1), fill = col))) +
    annotate("text", x = min(dat_10cm$x, dat_1m$x), y = .5, hjust = "left", label = pval_10cm, parse = TRUE, color = col[1], size = 3) +
    annotate("text", x = min(dat_10cm$x, dat_1m$x), y = .4, hjust = "left", label = pval_1m, parse = TRUE, color = col[2], size = 3) +
    pt_theme() +
    theme(
      legend.position = c(0.14, 0.7)) +
    labs(y = "Nesting incidence", x = lab_x, linetype = element_blank())
}



# get model averaged summary nicely formated
pt_get_tab_mod_avg <- function(model, type){
  if(type == "cond.avg"){
    coefs <- summary(model)$coefmat.subset
    ci <- confint(model, full = FALSE)
  }
  if(type == "full.avg"){
    coefs <- summary(model)$coefmat.full
    ci <- confint(model, full = TRUE)
  }
  tab <- data.frame(row.names = NULL, rownames(coefs), coefs)
  colnames(tab) <- c("var", "est", "se", "adj.se", "z", "p")
  
  tab$ci.lwr <- unname(ci[, 1])
  tab$ci.upr <- unname(ci[, 2])
  
  imp <- data.frame(
    row.names = NULL,
    var = names(MuMIn::sw(model)),
    ri = MuMIn::sw(model))
  
  tab <- left_join(tab, imp, by = "var")
  
  tab_pretty <- tab %>%
    arrange(desc(abs(est))) %>%
    mutate(
      est = sprintf("%.2f", est),
      p = ifelse(p<0.001, "< 0.001", sprintf("%.3f", p)),
      lab_ci = paste0("[", sprintf("%.2f", ci.lwr), ", ", sprintf("%.2f", ci.upr), "]"),
      ri = sprintf("%.2f", ri)) %>%
    select(var, est, lab_ci, p, ri) %>%
    rename(
      Parameter = var,
      Estimate = est,
      `95% CI` = lab_ci,
      `p-value` = p,
      `Relative importance` = ri)
  
  return(tab_pretty)
}




# compute effects from 
pt_get_effects.MuMIn_pa <- function(data, models.top, var_x, type.avg){
  # use only models containing var_x
  if(type.avg == "cond.avg"){
    sel <- numeric()
    for(i in 1:length(models.top)){
      match <- var_x %in% str_match(models.top[[i]]@call$formula, var_x)
      if(match) sel <- c(sel, i)
    }
    allmods <- models.top[sel]
  }
  
  # use all models
  if(type.avg == "full.avg") allmods <- models.top
  
  # get model weights
  mod.avg <- model.avg(models.top)
  modweights <- Weights(mod.avg)[sel]
  
  vars <- colnames(mod.avg$coefficients)[-1]
  var_x_orig <- word(var_x, 1, sep = "\\.")
  var_trafo <- word(var_x, 2, sep = "\\.")
  
  sd <- sd(data[[var_x_orig]])
  mu <- mean(data[[var_x_orig]])
  
  vars <- c(vars, var_x_orig)
  newdat <- as.data.frame(setNames(replicate(length(vars), rep(0, 100), simplify = F), vars))
  newdat[[var_x]] <- seq(min(data[[var_x]]), max(data[[var_x]]), length.out = 100)
  
  if(var_trafo == "z")
    newdat[[var_x_orig]] <- newdat[[var_x]] * sd + mu
  if(var_trafo == "log1_z")
    newdat[[var_x_orig]] <- exp(newdat[[var_x]] * sd + mu) - 1
  
  # distribute simulations according to model weights
  nsim <- 2000
  nsimpm <- round(nsim*modweights)
  
  # prepare matrices to store random values from posterior distribution of model-averaged predicted value
  predmat <- matrix(nrow = nrow(newdat), ncol = sum(nsimpm))
  
  # fit each model in turn and simulate from the posterior distribution of the prediction 
  for(i in 1:length(nsimpm)){
    nsimm <- nsimpm[i]
    if(nsimm == 0) next
    termsinmodi <- colnames(coef(allmods[[i]])[[1]])[-c(1)]
    if(length(termsinmodi) < 1) termsinmodi <- "1"
    fmla <- as.formula(paste("nest_pa ~", paste(termsinmodi, collapse = "+"), " + (1|field)"))
    modi <- glmer(fmla, data = data, family = binomial(link = "logit"))
    bsim <- sim(modi, n.sim = nsimpm[i])
    fmla <- as.formula(paste("nest_pa ~", paste(termsinmodi, collapse = "+")))
    Xmat <- model.matrix(fmla[c(1,3)], data = newdat)
    for(r in 1:nsimpm[i]){
      predmat[, c(0, cumsum(nsimpm))[i] + r] <- Xmat %*% bsim@fixef[r,]
    }
  }
  newdat$fit <- plogis(apply(predmat, 1, mean))
  newdat$lwr <- plogis(apply(predmat, 1, quantile, prob = 0.025))
  newdat$upr <- plogis(apply(predmat, 1, quantile, prob = 0.975))
  
  return(newdat)
}




# plot effect sizes for nesting incidence data
pt_plot_effects.pa <- function(data, newdata, var_x, lab_x, lab_stats = "", sign,  breaks = waiver()){
  var_x_orig <- word(var_x, 1, sep = "\\.")
  newdata$x <- newdata[[var_x_orig]]
  data$x <- data[[var_x_orig]]
  
  alpha_dots = .5; alpha_fill = .3; size_line = 1; size_point = 1.5; shape = 21
  position = position_jitter(w = 0, h = 0.02)
  
  lt <- ifelse(sign == TRUE, 1, 2)
  
  ggplot() +
    geom_ribbon(aes(x = x, y = fit, ymin = lwr, ymax = upr), newdata, alpha = 0.3) +
    geom_line(aes(x = x, y = fit), data = newdata, size = 1, linetype = lt) +
    geom_point(
      data = data,
      mapping = aes(y = nest_pa, x = x),
      position = position, shape = shape, size = size_point, alpha = alpha_dots) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    annotate("text", x = min(data$x), y = 0.8, hjust = "left", label = lab_stats, parse = TRUE, size = 3) +    
    scale_x_continuous(breaks = breaks) +
    pt_theme() +
    labs(y = "Nesting incidence", x = lab_x)
}





# compute effect sizes from glm with frequentist confidence intervals
pt_get_effects.glm <- function(data, model, link, var_x){
  
  vars <- names(coef(model))[-1]
  var_x_orig <- word(var_x, 1, sep = "\\.")
  var_trafo <- word(var_x, 2, sep = "\\.")

  newdata <- as.data.frame(setNames(replicate(length(vars), rep(0, 100), simplify = F), vars))
  newdata[, var_x] <- seq(min(data[[var_x]]), max(data[[var_x]]), length.out = 100)
  
  # predict at scale of response variable with back transformation
  # based on variable name ending that identifies the transformation
  pred <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  
  if(link == "log"){
    newdata$fit <- exp(pred$fit)
    newdata$upr <- exp(pred$fit + 1.96 * pred$se.fit)
    newdata$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
  }
  if(link == "logit"){
    newdata$fit <- plogis(pred$fit)
    newdata$upr <- plogis(pred$fit + 1.96 * pred$se.fit)
    newdata$lwr <- plogis(pred$fit - 1.96 * pred$se.fit)
  }
  
  
  if(var_trafo == "z"){
    sd <- sd(data[[var_x_orig]])
    mu <- mean(data[[var_x_orig]])
    newdata[[var_x_orig]] <- newdata[[var_x]] * sd + mu
  }
  
  if(var_trafo == "log1_z"){
    sd <- sd(log(data[[var_x_orig]] + 1))
    mu <- mean(log(data[[var_x_orig]] +1))
    newdata[[var_x_orig]] <- exp(newdata[[var_x]] * sd + mu) - 1
  }
  return(newdata)
}





# plot effect sizes for nest abundance data
pt_plot_effects.n <- function(data, newdata, var_x, lab_x, lab_stats, sign, breaks = waiver()){
  var_x_orig <- word(var_x, 1, sep = "\\.")
  newdata$x <- newdata[[var_x_orig]]
  data$x <- data[[var_x_orig]]
  
  lt <- ifelse(sign == TRUE, 1, 2)
  
  ggplot() +
    geom_ribbon(aes(x = x, y = fit, ymin = lwr, ymax = upr), newdata, alpha = 0.3) +
    geom_line(aes(x = x, y = fit), data = newdata, size = 1, linetype = lt) +
    geom_point(
      data = data,
      mapping = aes(y = nest_n, x = x),
      shape = 21, size = 1.5, alpha = .8) +
    annotate("text", x = min(data$x), y = .95*max(newdata$upr, data$nest_n), hjust = "left", label = lab_stats, parse = TRUE, size = 3) +
    scale_x_continuous(breaks = breaks) +
    pt_theme() +
    labs(y = expression("Nests" ~ "per" ~ 400 ~ m^2 ~ "plot area"), x = lab_x)
}






# inext plotting function, based on ggiNEXT with some modification
pt_plot_inext <- function(data, lab_y){
  colors <- c("#CD534CFF", "#0073C2FF")
  shapes <- c(15, 16)
  
  pt_ggiNEXT(data, type = 1, se = TRUE, lwd = 1.1, size_obs = 3, lab_guides = "Soil management") +
    pt_theme() +
    scale_fill_manual(values = colors, guide = "none") +
    scale_color_manual(values = colors, guide = "none") +
    scale_shape_manual(values = shapes, guide = "none") +
    labs(y = lab_y, fill = "Soil management") + 
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 11))
}





# source code of ggiNEXT with slight modification for nicer visualization
pt_ggiNEXT <- function(x, type = 1, se = TRUE, facet.var = "none", color.var = "site", 
                       grey = FALSE, lwd = 1.5, size_obs = 5, lab_guides = "Guides") 
{
  TYPE <- c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if (is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1) 
    stop("invalid plot type")
  if (is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, 
                                               SPLIT) == -1) 
    stop("invalid facet variable")
  if (is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, 
                                               SPLIT) == -1) 
    stop("invalid color variable")
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if (facet.var == "order") 
    color.var <- "site"
  if (facet.var == "site") 
    color.var <- "order"
  options(warn = -1)
  z <- fortify(x, type = type)
  options(warn = 0)
  if (ncol(z) == 7) {
    se <- FALSE
  }
  datatype <- unique(z$datatype)
  if (color.var == "none") {
    if (levels(factor(z$order)) > 1 & "site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- z$shape <- paste(z$site, z$order, sep = "-")
    }
    else if ("site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$shape <- z$site
    }
    else if (levels(factor(z$order)) > 1) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$shape <- factor(z$order)
    }
    else {
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }
  else if (color.var == "order") {
    z$col <- z$shape <- factor(z$order)
  }
  else if (color.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- z$site
  }
  else if (color.var == "both") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- paste(z$site, z$order, sep = "-")
  }
  zz = z
  z$method[z$method == "observed"] = "interpolated"
  z$lty <- z$lty <- factor(z$method, levels = unique(c("interpolated", 
                                                       "extrapolated"), c("interpolation", "interpolation", 
                                                                          "extrapolation")), labels = c("Interpolated", "Extrapolated"))
  z$col <- factor(z$col)
  data.sub <- zz[which(zz$method == "observed"), ]
  g <- ggplot(z, aes_string(x = "x", y = "y", colour = "col")) + 
    geom_point(aes_string(shape = "shape"), size = size_obs, data = data.sub)
  g <- g + geom_line(aes_string(linetype = "lty"), lwd = lwd) + 
    guides(linetype = guide_legend(title = "Method"), colour = guide_legend(title = lab_guides), 
           fill = guide_legend(title = lab_guides), shape = guide_legend(title = lab_guides)) + 
    theme(legend.position = "bottom", legend.title = element_blank(), 
          text = element_text(size = 18), legend.key.width = unit(1.2, 
                                                                  "cm"))
  if (type == 2L) {
    g <- g + labs(x = "Number of sampling units", y = "Sample coverage")
    if (datatype == "abundance") 
      g <- g + labs(x = "Number of individuals", y = "Sample coverage")
  }
  else if (type == 3L) {
    g <- g + labs(x = "Sample coverage", y = "Species diversity")
  }
  else {
    g <- g + labs(x = "Number of sampling units", y = "Species diversity")
    if (datatype == "abundance") 
      g <- g + labs(x = "Number of individuals", y = "Species diversity")
  }
  if (se) 
    g <- g + geom_ribbon(aes_string(ymin = "y.lwr", ymax = "y.upr", 
                                    fill = "factor(col)", colour = "NULL"), alpha = 0.2)
  if (facet.var == "order") {
    if (length(levels(factor(z$order))) == 1 & type != 2) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")
    }
    else {
      g <- g + facet_wrap(~order, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = lab_guides, 
                                              ncol = length(levels(factor(z$order))), byrow = TRUE), 
                        fill = guide_legend(title = lab_guides))
      }
    }
  }
  if (facet.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites.")
    }
    else {
      g <- g + facet_wrap(~site, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = lab_guides, 
                                              nrow = length(levels(factor(z$order)))), fill = guide_legend(title = lab_guides))
      }
    }
  }
  if (facet.var == "both") {
    if (length(levels(factor(z$order))) == 1 | !"site" %in% 
        names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites or orders.")
    }
    else {
      g <- g + facet_wrap(site ~ order)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = lab_guides, 
                                              nrow = length(levels(factor(z$site))), byrow = TRUE), 
                        fill = guide_legend(title = lab_guides))
      }
    }
  }
  if (grey) {
    g <- g + theme_bw(base_size = 18) + scale_fill_grey(start = 0, 
                                                        end = 0.4) + scale_colour_grey(start = 0.2, end = 0.2) + 
      guides(linetype = guide_legend(title = "Method"), 
             colour = guide_legend(title = lab_guides), fill = guide_legend(title = lab_guides), 
             shape = guide_legend(title = lab_guides)) + theme(legend.position = "bottom", 
                                                               legend.title = element_blank())
  }
  g <- g + theme(legend.box = "vertical")
  return(g)
}
