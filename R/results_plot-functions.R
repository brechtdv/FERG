## FERG RESULTS :: FIGURES

## line_plot ........ median, 50%, 90%, 95% uncertainty intervals, per agent
## region_barplot ... median rate, per agent & region
## full_barplot ..... YLD/YLL, ages, scaled to 100%, per agent
## scatter_plot ..... DALY/case ~ DALY/100,000, per agent

## -------------------------------------------------------------------------#
## LINE PLOT ---------------------------------------------------------------#

line_plot <-
function(agents, samples, asc = TRUE,
         names = NULL, ylim = NULL, ylab = "DALY (global)", ylog = TRUE,
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
         col = NULL) {
  ## remove NA's
  is_na <- sapply(samples, function(x) any(is.na(x)))
  agents  <- agents[!is_na]
  samples <- samples[!is_na]

  ## calculate quantiles
  Q <-
    sapply(samples,
           function(x)
             quantile(x, c(.025, .05, .25, .50, .75, .95, .975)))

  order   <- order(Q[4, ])
  if (!asc) order <- rev(order)

  ## define hazard names
  if (is.null(names)) {
    names <- names(agents)
  }

  ## convert 'names' to factor
  names <- gsub("_", ". ", names)
  names <- factor(names, levels = names[order])

  df <-
    data.frame(agent = names,
               Q1 = Q[1, ], Q2 = Q[2, ], Q3 = Q[3, ],
               medians = Q[4, ],
               Q4 = Q[5, ], Q5 = Q[6, ], Q6 = Q[7, ])

  ## compute lower and upper whiskers
  ylim_lwr <- boxplot.stats(samples[[which.min(order)]])$stats[1]
  ylim_upr <- boxplot.stats(samples[[which.max(order)]])$stats[5]
  if (is.null(ylim)) ylim <- c(ylim_lwr, ylim_upr)

  p <-
    ggplot(df, aes_string(x = "agent", y = "medians")) +
      geom_linerange(aes_string(ymin = "Q1", ymax = "Q6"),
                     size = 1, col = "grey")

  if (!is.null(col)) {
    p <-
      p +
      geom_linerange(aes_string(ymin = "Q2", ymax = "Q5"),
                     size = 1, col = col) +
      geom_linerange(aes_string(ymin = "Q3", ymax = "Q4"),
                     size = 2, col = col)
  } else {
    p <-
      p +
      geom_linerange(aes_string(ymin = "Q2", ymax = "Q5"), size = 1) +
      geom_linerange(aes_string(ymin = "Q3", ymax = "Q4"), size = 2)
  }

  p <-
    p +
      geom_point(size = 1, colour = "white") +
      scale_x_discrete("") +
      coord_cartesian(ylim = ylim) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()) +
      theme(axis.text.x = axis.text.x)

  if (ylog) {
    p <-
      p +
      scale_y_log10(ylab,
                    breaks = 10^pretty(range(log10(ylim))),
                    labels = formatC(10^pretty(range(log10(ylim))),
                                     format = "fg", big.mark = ","))
    #p <- p + scale_y_log10(ylab); p

  } else {
    p <- p + scale_y_continuous(ylab, breaks = pretty(range(ylim)))
  }

  print(p)
}


## -------------------------------------------------------------------------#
## REGION BARPLOT ----------------------------------------------------------#

region_barplot <-
function(x, names = NULL, xlab = NULL, ylab = NULL, col = NULL, nrow = 1) {
  ## define hazard names
  if (is.null(names)) {
    names <- colnames(x)
  }

  ## define 'xlab'
  if (is.null(xlab)) {
    xlab <- ifelse(nrow(x) == 6, "WHO Region", "WHO Subregion")
  }

  ## define 'col'
  if (is.null(col)) {
    col <- hue_pal()(ncol(x))
  }

  ## create df needed by ggplot2
  df <-
    data.frame(x = rep(rownames(x), ncol(x)),
               y = c(x),
               hazard = factor(rep(names, each = nrow(x)),
                               levels = names))

  ## generate barplot
  ggplot(data = df, aes_string(x = "x", y = "y", fill = "hazard")) +
    geom_bar(stat = "identity") +
    scale_x_discrete(xlab) +
    scale_y_continuous(ylab) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = col) +
    guides(fill = guide_legend(title = NULL, nrow = nrow))
}


## -------------------------------------------------------------------------#
## 100% BARPLOT ------------------------------------------------------------#

full_barplot <-
function(..., labels = NULL, names = NULL, sort = c("none", "asc", "desc"),
         xlab = "Agent", ylab = "Proportion", col = NULL, nrow = 1,
         axis.text.x = element_text()) {
  ## check 'sort'
  sort <- match.arg(sort)

  ## calculate means
  z <- list(...)
  n_cat <- length(z)
  n_haz <- length(z[[1]])
  n_sim <- length(z[[1]][[1]])

  x <- matrix(nrow = n_haz, ncol = n_cat)
  for (i in seq(n_haz)) {
    total <-
      rowSums(
        matrix(
          unlist(
            lapply(z,
                   function(x) {
                     if (length(x[[i]]) == 1) {
                       rep(x[[i]], n_sim)
                     } else {
                       x[[i]]
                     }
                   })),
        ncol = n_cat))

    for (j in seq(n_cat)) {
      x[i, j] <- mean(z[[j]][[i]] / total)
    }
  }

  ## define 'labels'
  if (is.null(labels)) {
    labels <- seq(n_cat)
  }

  ## define hazard names
  if (is.null(names)) {
    names <- names(z[[1]])
  }

  ## define order of hazards
  haz_order <-
    switch(sort,
           none = seq(n_haz),
           asc = order(x[, 1]),
           desc = rev(order(x[, 1])))
  names <- factor(names, levels = names[haz_order])

  ## define 'col'
  if (is.null(col)) {
    col <- hue_pal()(ncol(x))
  }

  ## create df needed by ggplot2
  df <-
    data.frame(x = rep(names, n_cat),
               y = c(x),
               label = factor(rep(labels, each = n_haz),
                              levels = labels))

  ## generate barplot
  ggplot(data = df, aes_string(x = "x", y = "y", fill = "label")) +
    geom_bar(stat = "identity", width = 1, colour = "black") +
    scale_x_discrete(xlab) +
    scale_y_continuous(ylab) +
    theme_bw() +
    scale_fill_manual(values = col) +
    guides(fill = guide_legend(title = NULL)) +
    theme(#legend.position = "right",
          #legend.title = element_text(NULL),
          #legend.box = "vertical",
          #legend.direction = "vertical",
          axis.text.x = axis.text.x) # + coord_flip()
}


## -------------------------------------------------------------------------#
## SCATTER PLOT ------------------------------------------------------------#

scatter_plot <-
function(x_val, y_val, fb = FALSE, names = NULL, col = NULL,
         xlog = TRUE, ylog = TRUE, labels = FALSE) {
  ## define hazard names
  if (is.null(names)) {
    names <- colnames(x_val)
  }
  names <- factor(names, levels = names)
  haz <- substr(colnames(x_val), 3, 6)

  ## define 'col'
  if (is.null(col)) {
    col <- hue_pal()(ncol(x_val))
  }

  ## create df needed by ggplot2
  df <-
    data.frame(x = x_val[1, ], y = y_val[1, ],
               x_lwr = x_val[2, ], x_upr = x_val[3, ],
               y_lwr = y_val[2, ], y_upr = y_val[3, ],
               Hazard = names, haz = haz)

  gplot <-
    ggplot(df, aes_string(x = "x", y = "y")) +
      geom_segment(aes_string(x = "x", xend = "x", y = "y_lwr", yend = "y_upr")) +
      geom_segment(aes_string(y = "y", yend = "y", x = "x_lwr", xend = "x_upr")) +
      #geom_point(aes_string(col = "Hazard", shape = "Hazard")) +
      geom_point(aes_string(col = "Hazard")) +
      scale_colour_manual(values = col) +
      #theme_bw() +
      theme(legend.position = "bottom") +
      guides(col = guide_legend(title = NULL, nrow = 4))

  if (xlog) {
   gplot <- gplot + scale_x_log10("DALY per 100,000")

  } else {
    gplot <- gplot + scale_x_continuous("DALY per 100,000")
  }

  if (ylog) {
   gplot <- gplot + scale_y_log10("DALY per case")

  } else {
    gplot <- gplot + scale_y_continuous("DALY per case")
  }

  if(labels) 
    gplot <-
      gplot + 
      geom_text(aes_string(label = "haz"), size = 3, hjust = -.2, vjust = -.2)

  print(gplot)
}
