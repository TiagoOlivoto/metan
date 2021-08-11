#' Draw Venn diagrams
#' @description
#' `r badge('stable')`
#'
#' Produces ggplot2-based Venn plots for 2, 3 or 4 sets. A Venn diagram shows
#' all possible logical relationships between several sets of data.
#' @param ... A list or a comma-separated list of vectors in the same class. If
#'   vector contains duplicates they will be discarded. If the list doesn't have
#'   names the sets will be named as `"set_1"`, "`Set_2"`, `"Set_3"` and so on.
#'   If vectors are given in `...`, the set names will be named with the names
#'   of the objects provided.
#' @param names By default, the names of the sets are set as the names of the
#'   objects in `...` (`names = NULL`). Use `names` to override this default.
#' @param show_elements Show set elements instead of count. Defaults to `FALSE`.
#' @param show_sets Show set names instead of count. Defaults to `FALSE`.
#' @param fill Filling colors in circles. Defaults to the default ggplot2 color
#'   palette. A vector of length 1 will be recycled.
#' @param alpha Transparency for filling circles. Defaults to `0.5`.
#' @param stroke_color Stroke color for drawing circles.
#' @param stroke_alpha Transparency for drawing circles.
#' @param stroke_size Stroke size for drawing circles.
#' @param stroke_linetype Line type for drawing circles. Defaults to `"solid"`.
#' @param name_color Text color for set names. Defaults to `"black"`.
#' @param name_size Text size for set names.
#' @param text_color Text color for intersect contents.
#' @param text_size Text size for intersect contents.
#' @param label_sep The separator for labs when `show_elements = TRUE`. Defaults
#'   to `","`.
#' @importFrom dplyr rowwise
#' @md
#' @return A ggplot object.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' (A <- letters[1:4])
#' (B <- letters[2:5])
#' (C <- letters[3:7])
#' (D <- letters[4:12])
#'
#' # create a Venn plot
#' venn_plot(A, B)
#'
#' # Three sets
#' venn_plot(A, B, C)
#'
#' # Four sets
#' venn_plot(A, B, C, D)
#'
#'
#' # Use a list
#' dfs <- list(A = A, B = B, C = C, D = D)
#' venn_plot(dfs,
#'           show_elements = TRUE,
#'           fill = c("red", "blue", "green", "gray"),
#'           stroke_color = "black",
#'           alpha = 0.8,
#'           text_size = 8,
#'           label_sep = ".")
#' }
venn_plot <- function(...,
                      names = NULL,
                      show_elements = FALSE,
                      show_sets = FALSE,
                      fill = ggplot_color(4),
                      alpha = 0.5,
                      stroke_color = "white",
                      stroke_alpha = 1,
                      stroke_size = 1,
                      stroke_linetype = "solid",
                      name_color = "black",
                      name_size = 6,
                      text_color = "black",
                      text_size = 4,
                      label_sep = ",") {

  set_name <-
    as.character(
      sapply(quos(...), function(x){
        rlang::quo_get_expr(x)
      })
    )
  if(is.list(c(...))){
    sets <- as.list(...)
    if(is.null(names)){
      if(is.null(names(sets))){
        names(sets) <- paste("set", seq_len(length(sets)), sep = "_")
      }
    } else{
      if(length(sets) != length(names)){
        stop("Incompatible length of names", call. = FALSE)
      }
      names(sets) <- names
    }
  }else{
    sets <- list(...)
    if(is.null(names)){
      if(is.null(names(sets))){
        names(sets) <- set_name
      }
    } else{
      if(length(sets) != length(names)){
        stop("Incompatible length of names", call. = FALSE)
      }
      names(sets) <- names
    }
  }
  if (!(length(sets) %in% 2:4)) {
    stop("Number of sets should be 2 or 4.", call. = FALSE)
  }
  if (length(unique(sapply(sets, class))) != 1) {
    stop("Vectors must be in the same class.", call. = FALSE)
  }
  if (!(sapply(sets, class)[1] %in% c("integer", "numeric", "character"))) {
    stop("The list must contain only integers, numerics or characters.", call. = FALSE)
  }
  if (!(alpha >= 0 & alpha <= 1)) {
    stop("'alpha' should be between 0 and 1.", call. = FALSE)
  }
  if(length(fill) == 1){
    fill_color <- rep(fill, length(sets))
  } else{
    fill_color <- fill
  }
  prepare_data <- function(data,
                           show_elements = show_elements,
                           label_sep = label_sep) {
    # Adapted from ggvenn' ggvenn()
    # https://github.com/yanlinlin82/ggvenn/blob/master/R/ggvenn.R
    # Copyright (c) 2019-2020 Linlin Yan
    # Thanks to @yanlinlin82 (https://twitter.com/yanlinlin82)
    get_shape <- function(group,
                           x_offset = 0,
                           y_offset = 0,
                           radius = 1,
                           radius_b = radius,
                           theta_offset = 0,
                           length.out = 200) {
      tibble(group = group,
             theta = seq(0, 2 * pi, length.out = length.out)) %>%
        mutate(x_raw = radius * cos(theta),
               y_raw = radius_b * sin(theta),
               x = x_offset + x_raw * cos(theta_offset) - y_raw * sin(theta_offset),
               y = y_offset + x_raw * sin(theta_offset) + y_raw * cos(theta_offset))
    }
    columns <- names(data) %>% head(4)

    # Two sets
    if (length(columns) == 2) {
      d <-   rbind(get_shape(1L, -.6, 0, 1),
                   get_shape(2L,  .6, 0, 1))
      set_a <- data[[1]]
      set_b <- data[[2]]
      texts <-
        tribble(~name, ~x, ~y, ~label,
                "A",   -1,  0,  set_difference(set_a, set_b),
                "B",    1,  0,  set_difference(set_b, set_a),
                "AB",   0,  0,  set_intersect(set_a, set_b)) %>%
        rowwise() %>%
        mutate(n = length(label),
               label = ifelse(show_elements == TRUE, paste(label, collapse = label_sep), n),
               label = ifelse(show_sets == TRUE, name, label))

      labels <-
        tribble(~name, ~x,   ~y,  ~hjust, ~vjust,
                "A",   -.6, 1.2, 0.5,    0,
                "B",    .6, 1.2, 0.5,    0) %>%
        mutate(label = columns)
      return(list(shapes = d, texts = texts, labels = labels))
    }
    # Three sets
    if (length(columns) == 3) {
      d <-
        rbind(get_shape(1L, -.6, (sqrt(1) + 2) / 6, 1),
              get_shape(2L, .6,(sqrt(1) + 2) / 6, 1),
              get_shape(3L, 0, -(sqrt(1) + 1.5) / 6, 1))
      set_a <- data[[1]]
      set_b <- data[[2]]
      set_c <- data[[3]]
      texts <-
        tribble(~name, ~x,    ~y, ~label,
                "A",   -1,    0.75, set_a %>% set_difference(set_union(set_b, set_c)),
                "B",    1,    0.75, set_b %>% set_difference(set_union(set_a, set_c)),
                "C",    0,      -1, set_c %>% set_difference(set_union(set_a, set_b)),
                "AB",   0,     0.9, set_a %>% set_intersect(set_b) %>% set_difference(set_c),
                "AC",  -0.5, -0.15, set_a %>% set_intersect(set_c) %>% set_difference(set_b),
                "BC",   0.5, -0.15, set_b %>% set_intersect(set_c) %>% set_difference(set_a),
                "ABC",  0,     0.2, set_a %>% set_intersect(set_b, set_c)
        ) %>%
        rowwise() %>%
        mutate(n = length(label),
               label = ifelse(show_elements == TRUE, paste(label, collapse = label_sep), n),
               label = ifelse(show_sets == TRUE, name, label))
      labels <-
        tribble(~name, ~x,    ~y,  ~hjust, ~vjust,
                "A",   -0.8,  1.6, 0.5,    0,
                "B",    0.8,  1.6, 0.5,    0,
                "C",    0,   -1.5, 0.5,    1) %>%
        mutate(label = columns)
      return(list(shapes = d, texts = texts, labels = labels))
    }
    # Four sets
    if (length(columns) == 4) {
      d <-   rbind(get_shape(1L, -.7, -1/2, .75, 1.5, pi/4),
                   get_shape(2L, -.72+2/3, -1/6, .75, 1.5, pi/4),
                   get_shape(3L, .72-2/3, -1/6, .75, 1.5, -pi/4),
                   get_shape(4L, .7, -1/2, .75, 1.5, -pi/4))
      set_a <- data[[1]]
      set_b <- data[[2]]
      set_c <- data[[3]]
      set_d <- data[[4]]
      texts <-   tribble(
        ~name, ~x,    ~y,  ~label,
        "A",   -1.5,  0,   set_a %>% set_difference(set_union(set_b, set_c, set_d)),
        "B",   -0.6,  0.7, set_b %>% set_difference(set_union(set_a, set_c, set_d)),
        "C",    0.6,  0.7, set_c %>% set_difference(set_union(set_a, set_b, set_d)),
        "D",    1.5,  0,   set_d %>% set_difference(set_union(set_a, set_b, set_c)),
        "AB",  -0.9,  0.3, set_a %>% set_intersect(set_b) %>% set_difference(set_c, set_d),
        "BC",   0,    0.4, set_b %>% set_intersect(set_c) %>% set_difference(set_a, set_d),
        "CD",   0.9,  0.3, set_c %>% set_intersect(set_d) %>% set_difference(set_a, set_b),
        "AC",  -0.8, -0.9, set_a %>% set_intersect(set_c) %>% set_difference(set_b, set_d),
        "BD",   0.8, -0.9, set_b %>% set_intersect(set_d) %>% set_difference(set_a, set_c),
        "AD",   0,   -1.4, set_a %>% set_intersect(set_d) %>% set_difference(set_b, set_c),
        "ABC", -0.5, -0.2, set_a %>% set_intersect(set_b, set_c) %>% set_difference(set_d),
        "BCD",  0.5, -0.2, set_b %>% set_intersect(set_c, set_d) %>% set_difference(set_a),
        "ACD", -0.3, -1.1, set_a %>% set_intersect(set_c, set_d) %>% set_difference(set_b),
        "ABD",  0.3, -1.1, set_a %>% set_intersect(set_b, set_d) %>% set_difference(set_c),
        "ABCD", 0,   -0.7, set_intersect(set_a, set_b, set_c, set_d)
      ) %>%
        rowwise() %>%
        mutate(n = length(label),
               label = ifelse(show_elements == TRUE, paste(label, collapse = label_sep), n),
               label = ifelse(show_sets == TRUE, name, label))
      labels <-
        tribble(~name, ~x,   ~y,    ~hjust, ~vjust,
                "A",   -1.6,  0.85,      1,     1,
                "B",   -0.8,  1.1,     0.5,     0,
                "C",    0.8,  1.1,     0.5,     0,
                "D",    1.6,  0.85,      0,     1) %>%
        mutate(label = columns)
      return(list(shapes = d, texts = texts, labels = labels))
    }
  }
  venn <- prepare_data(sets, show_elements, label_sep)
  venn$shapes %>%
    mutate(group = LETTERS[group]) %>%
    ggplot() +
    geom_polygon(aes(x = x, y = y, group = group, fill = group),
                 alpha = alpha) +
    geom_polygon(aes(x = x, y = y, group = group),
                 fill = NA,
                 color = stroke_color,
                 size = stroke_size,
                 alpha = stroke_alpha,
                 linetype = stroke_linetype) +
    geom_text(data = venn$labels,
              aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
              color = name_color,
              size = name_size) +
    geom_text(data = venn$texts,
              aes(x = x, y = y, label = label, hjust = 0.5, vjust = 0.5),
              color = text_color,
              size = text_size) +
    scale_x_continuous(expand = expansion(0.1)) +
    scale_y_continuous(expand = expansion(0.1)) +
    scale_fill_manual(values = fill_color) +
    guides(fill = "none") +
    coord_fixed() +
    theme_void()
}
