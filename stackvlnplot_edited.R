library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
base_idv_plot <- function(
  # A Seurat object
  obj,
  # The base plot should only made for one gene and we manually stitch later
  feature,
  # Points on the violins
  pt.size = 0.1,
  # If set as an identity, the violins will be split side-by-side like
  # vanilla VlnPlot(). split.plot will be automatically set as TRUE
  split.by = NULL,
  # If the plot is splitted, Wilcoxon test will be performed per cluster
  # between the two identities. If p.adj = TRUE, we adjust the p-value using
  # FDR (Benjamini-Hochberg). Otherwise, raw p-value is reported / used to
  # annotate the plot.
  p.adj = TRUE,
  deg.tbl = NULL,
  # Significance threshold
  p.cut.off = 0.05,
  # If TRUE, an asterisk will be added. Otherwise, p-value (adjusted or not)
  # will be shown two digits after the decimal point
  star = TRUE,
  # If you don't want an asterisk, change it
  star_symbol = "*",
  # How opaque the violins should be
  alpha = 0.5,
  # If not set, the viridis palette will be used. If a two element
  # color vector is provided here, it will be used instead.
  cols = NULL,
  # Whether to show y-axis labels
  x.label = FALSE,
  ...
    ) {

  # Automatically set split.plot if split.by is provided
  if (!is.null(split.by)) {
    split.plot <- TRUE
  } else {
    split.plot <- FALSE
  }

  # Generate a base plot from VlnPlot() and modify
  p <- VlnPlot(
    obj,
    # Make sure we only use one feature
    feature[1],
    # Return a raw ggplot object so it's easier to manipulate
    combine = FALSE,
    # Split the plot if an identity for splitting is provided
    split.by = split.by,
    # and set the split accordingly
    split.plot = split.plot,
    # Set point sizes
    pt.size = pt.size, ...
  )[[1]] +
    # VlnPlot() returns a list of ggplots when there are more than one feature
    # requested, but this makes it unpredictable, so we are forcing it to take
    # one feature only (see above) and retrieve only the first and only ggplot
    # object.
    # We remove the legend
    guides(fill = "none") +
    # And make the y-axis label the feature name
    labs(y = feature)

  # Add median indicators
  if (!is.null(split.by)) {
    # If the violins are split...
    p <- p +
      stat_summary(
        fun = median,
        geom = "crossbar",
        # Make bars per split and color them differently
        aes(group = split, color = split),
        # Make sure the bars "dodge" to the side they belong
        position = position_dodge(width = 1), width = 0.7
      )

    # Manually make violins (first ggplot layer) to be less opaque
    # since VlnPlot() does not allow setting this natively
    p$layers[[1]]$aes_params$alpha <- alpha
    p$layers[[2]]$aes_params$alpha <- 0.3
  } else {
    # If there's no split, just plot median bars per cluster
    p <- p +
      stat_summary(fun = median, geom = "crossbar")
  }

  # VlnPlot() is colored, so when we recolor the plot later, ggplot will
  # warn you about existing color scales being overwritten.
  # Here we manually purge existing scales to suppress that.
  p$scales$scales <- list()

  # Re-color the plot
  if (!is.null(cols) & length(cols) >= 2) {
    # If at least two colors are provided (= not NULL), we color it per user
    # setting
    p <- p +
      scale_fill_manual(
        # We only use the first two colors even if the user provides more since
        # split violins can only account for two levels
        values = cols[1:2]
      ) +
      scale_color_manual(
        values = cols[1:2]
      )
  } else {
    # If not specified, we use the viridis color palette
    p <- p +
      scale_fill_viridis_d() +
      scale_color_viridis_d()
  }


  # base_idv_plot() is not supposed to be used alone. Instead, we expect it
  # to generate multiple plots to be manually stacked together.
  # In such a case, we don't want x-axis labels (i.e., cluster labels) to
  # be shown for every plot but only for the bottom-most one.
  # To allow the user-facing function to do that, we allow a switch here.
  # That is, the x-axis label is only shown when x.label is TRUE.
  if (x.label) {
    p <- p +
      theme(
        axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1),
        legend.title = element_blank(),
        plot.title = element_blank()
      )
  } else {
    p <- p +
      theme(
        axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1),
        legend.title = element_blank(),
        plot.title = element_blank()
      )
  }

  if (!is.null(split.by)) {
    # The following is to perform Wilcoxon test between split within each cluster
    # so it should only be executed when split.by is set.

    # Extract raw data for Wilcoxon's test
    # ggplot saves raw data used to generate the plot in $data
    rawdata <- p$data

    if (is.null(deg.tbl)) {

      # Get all cluster identity involved so we can do Wilcox on each cluster
      clust_id <- unique(as.character(rawdata$ident))
      # Set the name so the output will be a named vector
      names(clust_id) <- clust_id

      # Retrieve p-value by running Wilcoxon test per cluster
      # vapply() is lapply() except that it checks the type of each element
      # by comparing to FUN.VALUE
      # I am specifying that the p values must be numeric otherwise vapply()
      # will fail -- failing early due to a check makes debugging easier.
      raw_p <- vapply(
        clust_id, function(cluster) {
          # For each clust_id, we only keep the rows belonging to
          # that cluster in the subset data
          cluster_data <- rawdata[rawdata$ident == cluster, ]

          # If a cluster comes from purely from one condition, the test cannot be
          # performed, so we check that here and if this happens, we return
          # a p-value of 1 because the comparison does not make sense.
          if (length(unique(cluster_data$split)) == 1) {
            return(1)
          }

          # We rename the gene expression column to always be named "exp".
          # It is originally the name of the feature, but this makes the later
          # wilcox.test() formula harder to write...
          colnames(cluster_data)[colnames(cluster_data) == feature] <- "exp"

          raw_p <- wilcox.test(data = cluster_data, exp ~ split)$p.value
          return(raw_p)
        },
        # Every element of the returned list must be one (1) numeric value
        # or this line will throw an error
        FUN.VALUE = numeric(1)
      )

      # Calculate FDR from raw p-values
      p_adj <- p.adjust(raw_p, method = "fdr")


      if (p.adj) {
        # If the user decided to use adjusted p-values (default)
        p_use <- p_adj
      } else {
        # Otherwise we use the raw Wilcoxon test results
        p_use <- raw_p
      }

      # We generate a p-value containing data.frame because ggplot2 likes
      # data.frames more than anything, and we will be using this
      # data.frame to add stars to the plot
      p_plot <- data.frame(
        p = p_use,
        ident = names(p_use)
      )

    } else {
      p_plot <- deg.tbl[deg.tbl$gene == feature , c("p_val_adj", "celltype")]
      colnames(p_plot) <- c("p", "ident")
    }

    # We add stars above the maximum value
    y_max <-  max(rawdata[[feature]])

    if (star) {
      # If we decide to plot a star
      p <- p +
        geom_text(
          # We take the p-value data.frame and only keep rows that
          # has a p-value (adjusted or not based on user choices)
          # below the cutoff
          data = subset(p_plot, p < p.cut.off),
          # Plot an asterisk for the ones passing the cut off by default
          label = star_symbol,
          aes(x = ident),
          y = y_max + 0.25,
          size = 8,
          # This layer is using a different data, so we don't want it
          # to inherit the aes() info from other layers of the plot.
          inherit.aes = FALSE
        ) +
        # We set the maximum of y-axis to be 0.5 above the highest measurement
        # so the stars will never be cropped.
        scale_y_continuous(limits = c(NA, y_max + 0.5))
    } else {
      p <- p +
        geom_text(
          data = subset(p_plot, p < p.cut.off),
          # If we don't want stars, we print out p-values in the place where
          # you'd otherwise plot stars
          aes(
            x = ident,
            label = round(p, digits = 2)
          ),
          y = y_max + 0.25,
          hjust = 0.5,
          size = rel(4),
          inherit.aes = FALSE
        ) +
        scale_y_continuous(limits = c(NA, y_max + 0.5))
    }
  }



  # Return the completed plot
  return(p)
  }

StackedVlnPlot <- function(
  # Receive the same arguments as the base plot except that now we can have
  # multiple features
  obj,
  features,
  pt.size = 0,
  split.by = NULL,
  p.adj = TRUE,
  star = TRUE,
  cols = NULL,
  p.cut.off = 0.05,
  alpha = 0.5,
  deg.tbl = NULL,
  title = NULL,
  star_symbol = "*",
  ...
) {
  # We count how many features are there so we know which plot is the last one
  # and should include x-axis labels
  feature_count <- length(features)

  # For each feature...
  gg_list <- lapply(
    # seq_along() counts along a vector (e.g., seq_along(c("A" ,"B", "C")) will
    # return c(1, 2, 3))
    seq_along(features), function(x) {
      # We check for each feature whether it is the last (i.e., index equals
      # the length)
      x.label = x == feature_count
      return(
        # Draw a base plot using provided setting
        base_idv_plot(
          obj = obj,
          feature = features[x],
          pt.size = pt.size,
          split.by = split.by,
          p.adj = p.adj,
          star = star,
          p.cut.off = p.cut.off,
          alpha = alpha,
          cols = cols,
          star_symbol = star_symbol,
          deg.tbl = deg.tbl,
          # and plot the x-axis labels for the last plot
          x.label = x.label, ...
        )
      )
    }
  )

  # Wrap the ggplot list into a column and show only one integrated legend
  composite_plot <- wrap_plots(gg_list, ncol = 1, guides = "collect")


  if (!is.null(title)) {
    # If a title is provided, we show it
    composite_plot <- composite_plot +
      plot_annotation(
        title = title,
        theme = theme(plot.title = element_text(size = 18, face = "bold")))
  }
  return(composite_plot)
}

FindSplitMarkers <- function(obj, split.by, split.ident = NULL) {
  mk <- lapply(
    # Get all active clusters and run FindMarkers() one by one
    levels(Idents(obj)),
    function(x) {
      # Get each subset
      subobj <- subset(obj, idents = x)
      # Use split.by as active identity
      Idents(subobj) <- split.by

      # Set ident 1 and 2 for marker finding
      if (is.null(split.ident)) {
        ident_lvls <- levels(Idents(subobj))
        if (length(ident_lvls) != 2) {
          # Error out if there are more than two levels for split.by
          stop("split.by must contain two levels.")
        }
      } else {
        ident_lvls <- split.ident
        if (length(ident_lvls) != 2) {
          # Error out if there are more than two levels for split.by
          stop("split.ident must contain two levels.")
        }
      }
      ident_1 <- ident_lvls[1]
      ident_2 <- ident_lvls[2]

      mk_tbl <- FindMarkers(subobj, ident_1, ident_2)
      mk_tbl$gene <- row.names(mk_tbl)
      mk_tbl$ident <- x
      return(mk_tbl)
    }
  )
  return(do.call(rbind, mk))
}

# This should work
StackedVlnPlot(
  mtl, features = goi, split.by = "genotype", pt.size = 0.0, cols =c("#F6416C", "#00B8A9"),
  title = "Title is here"
)
