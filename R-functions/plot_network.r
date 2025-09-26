plot_network <- function(
  mrIML_network,
  strength_threash = 0.05,
  network_title = "",
  nmds_layout = FALSE,
  degree = FALSE,
  group_colours = FALSE
) {
  # Filter to only edges above specified strength
  assoc_net_filtered <- mrIML_network %>%
    dplyr::filter(mean_strength > strength_threash)
  if (nrow(assoc_net_filtered) == 0) {
    warning(
      paste0("No effects above the specified threashold of ", strength_threash)
    )
    return()
  }
  # Convert to igraph for visualization
  g <- igraph::graph_from_data_frame(
    assoc_net_filtered,
    directed = TRUE
  )
  # Set edge and node attributes
  E(g)$Value <- assoc_net_filtered$mean_strength
  E(g)$Color <- ifelse(
    assoc_net_filtered$direction == "negative",
    yes = "red",
    no = "blue"
  )

  if (degree) {
    V(g)$sizeD <- (degree(g, mode = "out") * 2) + 3
  } else {
    V(g)$sizeD <- 1
  }

  if (group_colours) {
    clusters <- cluster_optimal(g)
    V(g)$color <- rainbow(length(clusters))[membership(clusters)]
  } else {
    V(g)$color <- "grey"
  }

  # Create network plot
  if (nmds_layout) {
    nmds_layout <- layout_with_mds(g, dim = 2)
    circle_layout <- layout_in_circle(g)
    gg <- ggnetwork(g, layout = circle_layout) %>%
      arrange(name)
  } else {
    gg <- ggnetwork::ggnetwork(g)
  }

  network_plot <- gg %>%
    ggplot2::ggplot(
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    ggnetwork::geom_edges(
      aes(color = Color, linewidth = Value),
      curvature = 0.1,
      arrow = arrow(length = unit(5, "pt"), type = "closed")
    ) +
    ggnetwork::geom_nodes(
      aes(color = color, size = sizeD)
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    ggnetwork::geom_nodelabel_repel(
      aes(label = name),
      data = gg,
      size = 3,
      segment.colour = "black",
      colour = "white",
      fill = "grey36"
    ) +
    ggplot2::labs(
      title = network_title,
      subtitle = "Blue = positive association, Red = negative association"
    )

  network_plot
}
