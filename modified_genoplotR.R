library(genoPlotR)

annotation_grob <-
  function (annotation, ...) 
{
  if (!is.annotation(annotation)) 
    stop("An annotation object is required")
  grob_list <- gList()
  if (nrow(annotation) < 1) 
    return(grob_list)
  for (i in 1:nrow(annotation)) {
    label <- as.annotation(annotation[i, ])
    grob_list[[i]] <- label_grob(label, ...)
  }
  grob_list
}

label_grob <-
  function (label, cex = 0.8, y_val = 0.0) 
{
  y <- y_val
  w <- 0.1
  if (!is.annotation(label)) 
    stop("An annotation object is required")
  if (nrow(label) > 1) 
    stop("A single-line annotation is required")
  grob_list <- gList()
  if (!is.na(label$x2)) {
    bracket_coord <- genoPlotR:::bracket_coord(label$x1, label$x2, y = y, 
                                   w = w)
    grob_list[[2]] <- linesGrob(x = bracket_coord$x, y = bracket_coord$y, 
                                name = paste("annot", "line", gsub(" ", "_", label$text), 
                                             sep = "."), default.units = "native", gp = gpar(col = label$col))
    x <- mean(c(label$x1, label$x2))
    w <- w * 2
  }
  else {
    x <- label$x1
  }
  if (label$rot == 0) {
    just <- c(0.5, 0)
  }
  else {
    just <- c(-0.1, 0.5)
  }
  grob_list[[1]] <- textGrob(parse(text = label$text), x = x, y = y + w, 
                             just = just, name = paste("annot", "label", gsub(" ", 
                                                                              "_", label$text), sep = "."), rot = label$rot, default.units = "native", 
                             gp = gpar(col = label$col, cex = cex))
  grob_list
}

plot_gene_map <-
  function (dna_segs, comparisons = NULL, tree = NULL, tree_width = NULL, 
          tree_branch_labels_cex = NULL, tree_scale = FALSE, legend = NULL, 
          annotations = NULL, top_annotation_banner = 0, top_annotation_padding = 0, annotation_height = 1, 
          annotation_y_val = 0.0, annotation_cex = 0.8, seg_height = 1, seg_plots = NULL, seg_plot_height = 3, 
          seg_plot_height_unit = "lines", seg_plot_yaxis = 3, seg_plot_yaxis_cex = scale_cex, 
          xlims = NULL, offsets = NULL, minimum_gap_size = 0.05, fixed_gap_length = FALSE, 
          limit_to_longest_dna_seg = TRUE, main = NULL, main_pos = "centre", 
          dna_seg_labels = NULL, dna_seg_label_cex = 1, dna_seg_label_col = "black", 
          gene_type = NULL, arrow_head_len = 200, dna_seg_line = TRUE, 
          scale = TRUE, dna_seg_scale = !scale, n_scale_ticks = 7, 
          scale_cex = 0.6, global_color_scheme = c("auto", "auto", 
                                                   "blue_red", 0.5), override_color_schemes = FALSE, plot_new = TRUE, 
          debug = 0, ...) 
{
  if (missing(dna_segs)) 
    stop("Argument dna_segs must be provided")
  if (!is.list(dna_segs) || !all(sapply(dna_segs, is.dna_seg))) 
    stop("Argument dna_segs must be a list of dna_seg objects")
  n_dna_segs <- length(dna_segs)
  n_rows <- 3 * n_dna_segs - 1
  n_comparisons <- length(comparisons)
  if (n_comparisons > 1) {
    if (!is.list(comparisons) || !all(sapply(comparisons, 
                                             is.comparison))) 
      stop("Argument comparisons must be a list of comparison objects")
  }
  if (n_comparisons > 0 && !(n_dna_segs - n_comparisons == 
                             1)) 
    stop("Number of comparisons not correct")
  if (is.null(dna_seg_labels) && !is.null(dna_segs)) {
    dna_seg_labels <- names(dna_segs)
  }
  if (!is.null(dna_seg_labels) && !(length(dna_seg_labels) == 
                                    n_dna_segs)) 
    stop("Argument dna_seg_labels doesn't have the same length as dna_segs")
  if (length(dna_seg_label_col) == 1) {
    dna_seg_label_col <- rep(dna_seg_label_col, n_dna_segs)
  }
  else if (!length(dna_seg_label_col) == n_dna_segs) {
    stop("Length of argument dna_seg_label_col must be 1 or as dna_segs")
  }
  if (!is.null(tree)) {
    if (!inherits(tree, "phylog")) 
      stop("Argument tree should be of class phylog (ade4)")
    if (is.null(dna_seg_labels)) 
      stop("If tree is given, label names should be provided via named list dna_segs or dna_seg_labels")
    if (length(tree$leaves) != n_dna_segs) 
      stop("Number of leaves in the tree not equal to number of dna segs")
    if (!all(dna_seg_labels %in% names(tree$leaves))) 
      stop("Tree leaves not corresponding to dna_seg labels")
    if (is.null(tree_branch_labels_cex)) {
      if (length(grep("^[^I]", names(tree$nodes)))) {
        tree_branch_labels_cex <- 0.8
      }
      else {
        tree_branch_labels_cex <- 0
      }
    }
  }
  if (!is.null(seg_plots)) {
    if (is.seg_plot(seg_plots)) {
      s_plot <- seg_plots
      seg_plots <- c(list(s_plot), rep(list(NULL), n_dna_segs - 
                                         1))
    }
    else if (length(seg_plots) == n_dna_segs) {
      if (!all(sapply(seg_plots, function(x) is.seg_plot(x) || 
                      is.null(x)))) 
        stop("All elements of seg_plots should be NULL or seg_plots objects")
    }
    else stop("seg_plots must be of same length as dna_segs")
    seg_plot_h <- ifelse(sapply(seg_plots, is.null), 0, 
                         seg_plot_height)
  }
  else {
    seg_plot_h <- rep(0, n_dna_segs)
  }
  if (is.null(seg_plot_yaxis) || !is.numeric(seg_plot_yaxis) || 
      is.null(seg_plots)) 
    seg_plot_yaxis <- 0
  if (!is.null(annotations)) {
    if (is.annotation(annotations)) {
      annot <- annotations
      annotations <- c(list(annot), rep(list(NULL), n_dna_segs - 
                                          1))
    }
    else if (length(annotations) == n_dna_segs) {
      if (!all(sapply(annotations, function(x) is.annotation(x) || 
                      is.null(x)))) 
        stop("All elements of annotations should be NULL or annotation objects")
    }
    else stop("annotation must be of same length as dna_segs")
    annot_h <- ifelse(sapply(annotations, is.null), 0, annotation_height)
  }
  else {
    annot_h <- rep(0, n_dna_segs)
  }
  if (!is.null(xlims)) {
    if (!is.list(xlims)) 
      stop("xlims must be a list")
    if (length(xlims) != n_dna_segs) 
      stop("xlims must be of the same length as dna_segs")
    if (!all(sapply(xlims, function(x) (length(x)%%2) == 
                    0))) 
      stop("All elements of xlims should have an even number of elements")
    for (i in 1:length(xlims)) {
      xlim <- xlims[[i]]
      rng <- range(dna_segs[[i]])
      min <- rng[1] - 0.02 * diff(rng)
      max <- rng[2] + 0.02 * diff(rng)
      if (is.null(xlim)) 
        xlim <- c(min, max)
      xlim[xlim == -Inf] <- min
      xlim[xlim == Inf] <- max
      if (!is.numeric(xlim)) 
        stop("All elements of xlims should be numeric")
      xlim <- data.frame(matrix(xlim, ncol = 2, byrow = TRUE))
      names(xlim) <- c("x0", "x1")
      xlim$strand <- ifelse(xlim$x0 < xlim$x1, 1, -1)
      for (j in 1:nrow(xlim)) xlim[j, 1:2] <- sort(as.numeric(xlim[j, 
                                                                   1:2]))
      xlims[[i]] <- xlim
    }
  }
  else {
    xlims <- list()
    for (i in 1:n_dna_segs) {
      rng <- range(dna_segs[[i]])
      xlims[[i]] <- data.frame(x0 = rng[1] - 0.02 * diff(rng), 
                               x1 = rng[2] + 0.02 * diff(rng), strand = 1)
    }
  }
  if (!is.null(offsets) && length(offsets) != n_dna_segs) {
    stop("Length of offsets not equal to number of dna_segs")
  }
  if (main_pos == "centre") {
    main_x <- 0.5
    main_just <- "centre"
  }
  else if (main_pos == "left") {
    main_x <- 0
    main_just <- "left"
  }
  else if (main_pos == "right") {
    main_x <- 1
    main_just <- "right"
  }
  else {
    stop("main_pos should be one of centre, left, right")
  }
  if (is.logical(dna_seg_line)) {
    dna_seg_line <- as.character(dna_seg_line)
    dna_seg_line[dna_seg_line == "TRUE"] <- "black"
  }
  if (!is.character(dna_seg_line)) 
    stop("dna_seg_line should be eiher a logical or character giving color")
  if (length(dna_seg_line) == 1) {
    dna_seg_line <- rep(dna_seg_line, n_dna_segs)
  }
  else if (length(dna_seg_line) != n_dna_segs) {
    stop("dna_seg_line should be of length 1 or same length as dna_segs")
  }
  if (!is.null(gene_type) && !(gene_type %in% gene_types())) 
    stop(paste("gene_type muste be one of:", paste(gene_types(), 
                                                   collapse = ", ")))
  if (is.logical(dna_seg_scale)) {
    if (length(dna_seg_scale) == 1) {
      dna_seg_scale <- rep(dna_seg_scale, n_dna_segs)
    }
    else if (length(dna_seg_scale) != n_dna_segs) {
      stop("dna_seg_scale must be the same length dna_segs")
    }
  }
  else {
    stop("dna_seg_scale must be logical")
  }
  if (length(global_color_scheme) != 4) 
    stop("global_color_scheme should be length 4")
  glob_col_sch_2_vals <- c("increasing", "decreasing", "auto")
  if (length(grep(global_color_scheme[2], glob_col_sch_2_vals)) != 
      1) {
    stop(paste("Second argument to global_color_scheme should be one of", 
               paste(glob_col_sch_2_vals, collapse = ", ")))
  }
  else {
    global_color_scheme[2] <- grep(global_color_scheme[2], 
                                   glob_col_sch_2_vals, value = TRUE)
  }
  h <- rep(1, n_rows)
  h[seq(2, n_rows, by = 3)] <- 1 + scale_cex * dna_seg_scale + 
    annot_h
  h[seq(1, n_rows, by = 3)] <- seg_plot_h
  h[1] <- h[1] + top_annotation_banner
  dna_seg_heights <- unit(h, c(rep(c(seg_plot_height_unit, 
                                     "lines", "null"), n_dna_segs), "lines"))
  if (!is.null(gene_type)) {
    if (gene_type == "auto") {
      n_genes <- sapply(dna_segs, nrow)
      gene_type <- auto_gene_type(n_genes)
    }
    for (i in 1:n_dna_segs) {
      dna_segs[[i]]$gene_type <- gene_type
    }
  }
  if ((!any("col" %in% unlist(lapply(comparisons, names))) || 
       override_color_schemes) && !is.null(comparisons)) {
    num_cols <- lapply(comparisons, function(x) names(x)[sapply(x, 
                                                                is.numeric)])
    shared_num_cols <- names(which(table(unlist(num_cols)) == 
                                     length(num_cols)))
    shared_num_cols <- shared_num_cols[!shared_num_cols %in% 
                                         c("start1", "start2", "end1", "end2")]
    if (global_color_scheme[1] == "auto") {
      names_comp_1 <- names(comparisons[[1]])
      global_color_scheme[1] <- if ("per_id" %in% shared_num_cols) {
        "per_id"
      }
      else if ("e_value" %in% shared_num_cols) {
        "e_value"
      }
      else {
        names_comp_1[names_comp_1 %in% shared_num_cols][1]
      }
    }
    else if (!global_color_scheme[1] %in% shared_num_cols) {
      stop("One or all columns don't have the indicated column for global color scheme")
    }
    if (global_color_scheme[2] == "auto") {
      global_color_scheme[2] <- if (global_color_scheme[1] %in% 
                                    c("mism", "gaps", "e_value")) 
        TRUE
      else FALSE
    }
    else if (global_color_scheme[2] == "decreasing") {
      global_color_scheme[2] <- TRUE
    }
    else if (global_color_scheme[2] == "increasing") {
      global_color_scheme[2] <- FALSE
    }
    else {
      stop("Invalid value for global_color_scheme[2]")
    }
    range_col_from <- c(Inf, -Inf)
    for (i in 1:n_comparisons) {
      if (nrow(comparisons[[i]]) > 0) {
        range_col_from[1] <- min(c(range_col_from[1], 
                                   comparisons[[i]][[global_color_scheme[1]]]))
        range_col_from[2] <- max(c(range_col_from[2], 
                                   comparisons[[i]][[global_color_scheme[1]]]))
      }
    }
    for (i in 1:n_comparisons) {
      comparisons[[i]]$col <- apply_color_scheme(x = comparisons[[i]][[global_color_scheme[1]]], 
                                                 direction = comparisons[[i]]$direction, color_scheme = global_color_scheme[3], 
                                                 decreasing = global_color_scheme[2], rng = range_col_from, 
                                                 transparency = as.numeric(global_color_scheme[4]))
    }
  }
  for (i in 1:n_dna_segs) {
    xlims[[i]]$length <- xlims[[i]]$x1 - xlims[[i]]$x0
  }
  def_gap_length <- max(sapply(xlims, function(x) sum(x$length))) * 
    minimum_gap_size
  unpadded_lengths <- sapply(xlims, function(x) sum(x$length) + 
                               (nrow(x) - 1) * def_gap_length)
  longest_seg <- which.max(unpadded_lengths)
  max_length <- unpadded_lengths[longest_seg]
  scale_unit <- diff(pretty(c(0, max_length), n = n_scale_ticks + 
                              2)[1:2])
  seg_subplots <- list()
  dna_subsegs <- list()
  for (i in 1:n_dna_segs) {
    n_subsegs <- nrow(xlims[[i]])
    dna_subsegs[[i]] <- list()
    seg_subplots[[i]] <- list()
    for (j in 1:n_subsegs) {
      dna_subsegs[[i]][[j]] <- genoPlotR:::trim.dna_seg(dna_segs[[i]], 
                                            c(xlims[[i]]$x0[j], xlims[[i]]$x1[j]))
      if (seg_plot_h[[i]] > 0) {
        seg_subplots[[i]][[j]] <- genoPlotR::trim.seg_plot(seg_plots[[i]], 
                                                c(xlims[[i]]$x0[j], xlims[[i]]$x1[j]))
      }
    }
  }
  if (n_comparisons > 0) {
    for (i in 1:n_comparisons) {
      comp1 <- comparisons[[i]][0, ]
      for (j in 1:nrow(xlims[[i]])) {
        for (k in 1:nrow(xlims[[i + 1]])) {
          comp1 <- rbind(comp1, genoPlotR:::trim.comparison(comparisons[[i]], 
                                                xlim1 = as.numeric(xlims[[i]][j, c("x0", 
                                                                                   "x1")]), xlim2 = as.numeric(xlims[[i + 
                                                                                                                        1]][k, c("x0", "x1")])))
        }
      }
      comparisons[[i]] <- comp1
    }
  }
  if (is.null(offsets)) {
    prel_offsets <- lapply(xlims, function(x) c(0, rep(def_gap_length, 
                                                       nrow(x) - 1)))
    offsets <- minimize_comps(comparisons, xlims, unpadded_lengths, 
                              prel_offsets, fixed_gap_length)
  }
  else {
    offsets <- as.list(offsets)
    for (i in 1:n_dna_segs) {
      if (length(offsets[[i]]) == 1) {
        offsets[[i]] <- c(offsets[[i]], rep(def_gap_length, 
                                            nrow(xlims[[i]]) - 1))
      }
      else if (length(offsets[[i]]) != nrow(xlims[[i]])) {
        stop("The length of each element of offsets should be either one or equal to the number of subsegments in the corresponding segment.")
      }
    }
  }
  if (limit_to_longest_dna_seg) {
    for (i in 1:n_dna_segs) {
      tot_length <- sum(c(xlims[[i]]$length, offsets[[i]]))
      if (tot_length > max_length) {
        excess <- tot_length - max_length
        for (j in length(offsets[[i]]):1) {
          if ((offsets[[i]][j] - excess) < def_gap_length) {
            excess <- excess - offsets[[i]][j] + def_gap_length
            offsets[[i]][j] <- def_gap_length
          }
          else {
            offsets[[i]][j] <- offsets[[i]][j] - excess
            break
          }
        }
      }
    }
  }
  else {
    padded_lengths <- sapply(xlims, function(x) sum(x$length)) + 
      sapply(offsets, sum)
    max_length <- max(padded_lengths)
  }
  if (n_comparisons > 0) {
    for (i in 1:n_comparisons) {
      comparisons[[i]] <- genoPlotR:::calc_comp_coor(offsets[[i]], 
                                         xlims[[i]], comparisons[[i]], side = 1)
      comparisons[[i]] <- genoPlotR:::calc_comp_coor(offsets[[i + 
                                                    1]], xlims[[i + 1]], comparisons[[i]], side = 2)
    }
  }
  padded_lengths <- sapply(xlims, function(x) sum(x$length)) + 
    sapply(offsets, sum)
  max_length <- max(padded_lengths)
  longest_segment <- which.max(padded_lengths)
  dna_seg_grobs <- list()
  dna_seg_scale_grobs <- list()
  for (i in 1:n_dna_segs) {
    dna_seg_grobs[[i]] <- list()
    dna_seg_scale_grobs[[i]] <- list()
    for (j in 1:length(dna_subsegs[[i]])) {
      if (debug > 0 && debug < nrow(dna_subsegs[[i]][[j]])) 
        dna_subsegs[[i]][[j]] <- dna_subsegs[[i]][[j]][1:debug, 
        ]
      dna_seg_grobs[[i]][[j]] <- genoPlotR:::dna_seg_grob(dna_subsegs[[i]][[j]], 
                                              arrow_head_len, i, ...)
      dna_seg_scale_grobs[[i]][[j]] <- if (dna_seg_scale[[i]]) 
        dna_seg_scale_grob(range = xlims[[i]][j, c("x0", 
                                                   "x1")], cex = scale_cex, unit = scale_unit, 
                           i = i, j = j)
      else NULL
    }
  }
  seg_plot_grobs <- list()
  seg_plot_ylims <- list()
  for (i in 1:n_dna_segs) {
    seg_plot_grobs[[i]] <- list()
    xl_sg <- c(Inf, -Inf)
    for (j in 1:length(dna_subsegs[[i]])) {
      if (length(seg_plots[[i]]) > 0) {
        grb <- do.call(seg_subplots[[i]][[j]]$func, 
                       seg_subplots[[i]][[j]]$args)
        rng <- nice_ylim.seg_plot(seg_subplots[[i]][[j]])
        xl_sg[1] <- min(xl_sg[1], rng[1])
        xl_sg[2] <- max(xl_sg[2], rng[2])
        seg_plot_grobs[[i]][[j]] <- grb
      }
      else {
        seg_plot_grobs[[i]] <- NULL
      }
    }
    seg_plot_ylims[[i]] <- if (is.null(seg_plots[[i]]$ylim)) 
      xl_sg
    else seg_plots[[i]]$ylim
  }
  seg_plot_yaxis_grobs <- list()
  for (i in 1:n_dna_segs) {
    seg_plot_yaxis_grobs[[i]] <- if (length(seg_plots[[i]]) > 
                                     0 && seg_plot_yaxis > 0) 
      yaxis_grob(seg_plot_ylims[[i]], cex = seg_plot_yaxis_cex, 
                 n = seg_plot_yaxis, i)
    else NULL
  }
  comparison_grobs <- list()
  if (n_comparisons > 0) {
    for (i in 1:n_comparisons) {
      if (debug > 0 && debug < nrow(comparisons[[i]])) 
        comparisons[[i]] <- comparisons[[i]][1:debug, 
        ]
      comparison_grobs[[i]] <- genoPlotR:::comparison_grob(comparisons[[i]], 
                                               i)
    }
  }
  if (!is.null(annotations)) {
    annotation_grobs <- list()
    for (i in 1:n_dna_segs) {
      if (is.null(annotations[[i]])) {
        annotation_grobs[[i]] <- list(NULL)
      }
      else {
        annotation_grobs[[i]] <- list()
        for (j in 1:length(dna_subsegs[[i]])) {
          annot <- genoPlotR:::trim.annotation(annotations[[i]], 
                                   xlims[[i]][j, c("x0", "x1")])
          annotation_grobs[[i]][[j]] <- annotation_grob(annot, 
                                                        annotation_cex,
                                                        annotation_y_val)
        }
      }
    }
  }
  if (scale) {
    scale_grob <- genoPlotR:::scale_grob(max_length)
    scale_h <- 1
  }
  else {
    scale_h <- 0
  }
  if (!is.null(main)) {
    main_grob <- textGrob(x = main_x, label = main, gp = gpar(cex = 1.2), 
                          just = main_just)
    main_h <- 1.8
  }
  else {
    main_h <- 0
  }
  if (!is.null(tree)) {
    y <- permute_tree(tree, dna_seg_labels)
    tree_grob <- phylog_grob(tree, 1 - ((y - 1)/(n_dna_segs - 
                                                   1)), clabel.leaves = dna_seg_label_cex, clabel.nodes = tree_branch_labels_cex, 
                             tree.scale = tree_scale)
    tree_w <- unit(0.1, "npc") + tree_grob$width
  }
  else if (!is.null(dna_seg_labels)) {
    tree_grob <- genoPlotR:::dna_seg_label_grob(dna_seg_labels, cex = dna_seg_label_cex, 
                                    col = dna_seg_label_col)
    tree_w <- tree_grob$width
  }
  else {
    tree_grob <- NULL
    tree_w <- unit(0, "npc")
  }
  if (!is.null(tree_width)) 
    tree_w <- unit(tree_width, "inches")
  if (tree_scale) 
    scale_h <- 1
  if (plot_new) 
    grid.newpage()
  pushViewport(viewport(width = unit(1, "npc") - unit(1, "lines"), 
                        height = unit(1, "npc") - unit(1, "lines"), name = "oma"), 
               viewport(layout = grid.layout(2, 1, heights = unit(c(main_h, 
                                                                    1), c("lines", "null"))), name = "oma_layout"))
  if (!is.null(main)) {
    pushViewport(viewport(layout.pos.row = 1, name = "main"))
    grid.draw(main_grob)
    upViewport()
  }
  seg_plot_yaxis_w <- if (seg_plot_yaxis > 0 && !is.null(seg_plots[[longest_segment]])) 
    unit(3, "grobwidth", data = seg_plot_yaxis_grobs[[longest_segment]])
  else unit(0, "null")
  pushViewport(viewport(layout.pos.row = 2, layout = grid.layout(2, 
                                                                 3, heights = unit(c(1, scale_h), c("null", "lines")), 
                                                                 widths = unit.c(tree_w, unit(1, "null"), seg_plot_yaxis_w)), 
                        name = "frame"))
  if (scale) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2, 
                          xscale = c(0, max_length), name = "scale"))
    grid.draw(scale_grob)
    upViewport()
  }
  if (!is.null(tree_grob)) {
    annot_margin <- unit(if (is.null(annotations[[1]])) 
      0 + top_annotation_banner
      else annotation_height + top_annotation_banner, "lines")
    seg_plot_margin <- unit(if (is.null(seg_plot_grobs[[1]][[1]])) 
      0
      else seg_plot_height, seg_plot_height_unit)
    hli <- unit(0.5, "lines")
    dna_scale_margin <- unit(scale_cex * dna_seg_scale[[n_dna_segs]], 
                             "lines")
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, 
                          layout = grid.layout(6, 1, heights = unit.c(seg_plot_margin, 
                                                                      annot_margin, hli, unit(n_dna_segs * (1 + seg_plot_height), 
                                                                                              "null"), hli, dna_scale_margin)), name = "tree_outer"))
    pushViewport(viewport(layout.pos.row = 4, width = unit(1, 
                                                           "npc") - unit(1, "lines"), just = c("centre", "bottom"), 
                          name = "tree"))
    grid.draw(tree_grob$grob)
    upViewport(2)
  }
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, 
                        name = "plotarea_outer", clip = "on"), viewport(width = unit(1, 
                                                                                     "npc") - unit(1, "lines"), height = unit(1, "npc") - 
                                                                          unit(0, "lines"), name = "plotarea", clip = "off"))
  pushViewport(viewport(layout = grid.layout(n_rows, 1, heights = dna_seg_heights), 
                        name = "map"))
  if (n_comparisons > 0) {
    for (i in 1:n_comparisons) {
      pushViewport(viewport(layout.pos.row = 3 * i, yscale = c(0, 
                                                               1), xscale = c(0, max_length), name = paste("comparison", 
                                                                                                           i, sep = ".")))
      grid.draw(comparison_grobs[[i]])
      upViewport()
    }
  }
  for (i in 1:n_dna_segs) {
    n_dna_subsegs <- length(dna_subsegs[[i]])
    n_cols <- n_dna_subsegs * 2 + 1
    widths <- numeric(n_cols)
    widths[1:n_dna_subsegs * 2] <- xlims[[i]]$length
    widths[1:n_dna_subsegs * 2 - 1] <- offsets[[i]]
    widths_units <- unit(widths, rep("native", n_cols))
    heights <- unit(c(annot_h[i], seg_height, scale_cex * dna_seg_scale[i]), 
                    c("lines", "lines", "lines"))
    if (i == 1) {
      heights[1] <- heights[1] + unit(top_annotation_padding,"lines")
    }
    pushViewport(viewport(layout.pos.row = 3 * i - 2, layout = grid.layout(1, 
                                                                           n_cols, widths = widths_units, just = c("left", 
                                                                                                                   "centre")), xscale = c(0, max_length), name = paste("seg_plot", 
                                                                                                                                                                       i, sep = ".")))
    for (j in 1:n_dna_subsegs) {
      idx <- if (xlims[[i]]$strand[j] == 1) 
        c("x0", "x1")
      else c("x1", "x0")
      xscale <- as.numeric(xlims[[i]][j, idx])
      if (!is.null(seg_plots[[i]])) {
        pushViewport(viewport(layout.pos.col = j * 2, 
                              yscale = seg_plot_ylims[[i]], xscale = xscale, 
                              just = c("left", "centre"), name = paste("seg_subplot", 
                                                                       i, j, sep = ".")))
        grid.draw(seg_plot_grobs[[i]][[j]])
        upViewport()
      }
    }
    if (!is.null(seg_plots[[i]]) && seg_plot_yaxis > 0) {
      pushViewport(viewport(layout.pos.col = n_cols, yscale = seg_plot_ylims[[i]], 
                            width = unit(1, "grobwidth", data = seg_plot_yaxis_grobs[[i]]), 
                            just = c("left", "centre"), name = paste("seg_plot_yaxis", 
                                                                     i, sep = ".")))
      grid.draw(seg_plot_yaxis_grobs[[i]])
      upViewport()
    }
    upViewport()
    pushViewport(viewport(layout.pos.row = 3 * i - 1, layout = grid.layout(3, 
                                                                           n_cols, heights = heights, widths = widths_units, 
                                                                           just = c("left", "centre")), xscale = c(0, max_length), 
                          name = paste("scale_and_dna_seg", i, sep = ".")))
    for (j in 1:n_dna_subsegs) {
      idx <- if (xlims[[i]]$strand[j] == 1) 
        c("x0", "x1")
      else c("x1", "x0")
      xscale <- as.numeric(xlims[[i]][j, idx])
      if (!is.null(annotations[[i]])) {
        pushViewport(viewport(layout.pos.row = 1, layout.pos.col = j * 
                                2, yscale = c(0, 1), xscale = xscale, just = c("left", 
                                                                               "centre"), name = paste("annotation", i, j, 
                                                                                                       sep = ".")))
        grid.draw(annotation_grobs[[i]][[j]])
        upViewport()
      }
      if (dna_seg_scale[i]) {
        pushViewport(viewport(layout.pos.row = 3, layout.pos.col = j * 
                                2, yscale = c(0, 1), xscale = xscale, just = c("left", 
                                                                               "centre"), name = paste("dna_seg_scale", i, 
                                                                                                       j, sep = ".")))
        grid.draw(dna_seg_scale_grobs[[i]][[j]])
        upViewport()
      }
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j * 
                              2, yscale = c(0, 1), xscale = xscale, just = c("left", 
                                                                             "centre"), name = paste("dna_seg", i, j, sep = ".")))
      if (!dna_seg_line[i] == "FALSE") {
        grid.segments(x0 = unit(xlims[[i]]$x0[j], "native"), 
                      y0 = unit(0.5, "native"), x1 = unit(xlims[[i]]$x1[j], 
                                                          "native"), y1 = unit(0.5, "native"), name = paste("dna_seg_line", 
                                                                                                            i, j, sep = "."), gp = gpar(col = dna_seg_line[i]))
      }
      grid.draw(dna_seg_grobs[[i]][[j]])
      upViewport()
      if (j > 1) {
        pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j * 
                                2 - 1, yscale = c(0, 1), xscale = c(0, widths[2 * 
                                                                                j - 1]), just = c("centre", "centre"), name = paste("gap", 
                                                                                                                                    i, j, sep = ".")))
        grid.draw(gap_grob(w = def_gap_length, m = widths[2 * 
                                                            j - 1]/2, i, j))
        upViewport()
      }
    }
    upViewport()
  }
  upViewport(2)
  upViewport(2)
  upViewport(2)
}


# 
# # Output file
# pdf(file=output_fn,width=6,height=5)
# 
# # Plot
# plot_gene_map(dna_sequences,
#               comparisons = comparison_files,
#               annotations = annotations,
#               seg_height = 0.75,
#               annotation_height = 0.0,
#               annotation_cex = 1.0,
#               dna_seg_label_cex = 1.0,
#               dna_seg_labels = seq_labels,
#               offsets = offset_lengths,
#               top_annotation_padding = 0.25,
#               annotation_y_val = 0.5,
#               top_annotation_banner = 3.0)
# 
# dev.off()