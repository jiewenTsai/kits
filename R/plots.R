#' @include style.R
NULL

# 圖的顏色與符號（PAL、FONT 在 style.R）
# DIF 的樣式照方向：up = 對 focal 組有利（橘、朝上的三角），down = 對 ref 組有利（藍、朝下）；
# 男女的標籤固定：favors girls 是 up、favors boys 是 down（和以前相同）
COL_DIF <- c(none = PAL$gray, up = PAL$orange, down = PAL$blue, inconsistent = PAL$black)
SYM_DIF <- c(none = "circle", up = "triangle-up", down = "triangle-down", inconsistent = "square")
COL_TWO <- c(PAL$blue, PAL$orange)
COL_MANY <- c(PAL$blue, PAL$orange, "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
SYM_TWO <- c("circle", "square")

#' Age ticks for plots
#'
#' The default `age_ticks` of the plots against age: `plot()` of a [kits_item_table()], of a
#' [kits_norms()] by group and of a [kits_vertical()].
#'
#' @format A numeric vector: 3, 6, 12, 18, 24, 36, 48, 60, 72, 84, 96.
#' @examples
#' AGE_TICKS
#' @export
AGE_TICKS <- c(3, 6, 12, 18, 24, 36, 48, 60, 72, 84, 96)


# 題幹的短版（內部用，給成長階梯）：拿掉括號裡的例子；中文留前 14 個字；
# 英文留前 36 個字元，在字和字之間切，切掉的加刪節號。
short_stem <- function(stem) {
  stem <- ifelse(is.na(stem), "", stem)
  stem <- trimws(gsub("[\uff08(].*$", "", stem))
  english_short <- function(x) {
    if (nchar(x) <= 36) return(x)
    paste0(sub("\\s+\\S*$", "", substr(x, 1, 36)), "\u2026")
  }
  is_zh <- grepl("[\u4e00-\u9fa5]", stem)
  ifelse(is_zh, substr(stem, 1, 14), map_chr(stem, english_short))
}


# 滑鼠停在題目上看到的文字（rp：題目位置的反應機率，預設 RP50）
item_hover <- function(d, rp = 0.5) {
  domain <- lab_name("domain", d$domain)
  age <- d[[age_col(rp)]]
  paste0("<b>", d$item_id, "</b>  ", domain, "<br>", d$stem,
         "<br>", LAB$hover_step, if (rp != 0.5) sprintf(" (RP%g)", 100 * rp) else "", " = ", sprintf("%.1f", d$D),
         "<br>", sprintf(LAB$hover_age, 100 * rp), " = ", ifelse(is.na(age), LAB$hover_na, sprintf("%.1f", age)),
         "<br>", LAB$hover_reach, " = ", ifelse(is.na(d$age_reach), "\u2014", sprintf("%.1f", d$age_reach)),
         "<br>", LAB$hover_dif, ": ", LAB$dif[d$DIF])
}

# 每一格的標題：置中在格子上方（像 matplotlib 的 ax.set_title）
panel_title <- function(text) {
  list(list(x = 0.5, y = 1, xref = "paper", yref = "paper", text = text,
            showarrow = FALSE, xanchor = "center", yanchor = "bottom", font = list(size = 12)))
}

# 格子右下角的註記：這一格沒畫出來的題（圖上沒有位置）
not_drawn_note <- function(ids, label) {
  if (length(ids) == 0) return(list())
  list(list(x = 0.98, y = 0.02, xref = "paper", yref = "paper", text = paste0(label, paste(ids, collapse = ", ")),
            showarrow = FALSE, xanchor = "right", yanchor = "bottom", font = list(size = 10, color = PAL$gray)))
}

# 檢查題目表的 DIF 欄。NA 當成 "none"。
clean_dif <- function(dif) {
  dif[is.na(dif)] <- "none"
  bad <- unique(dif[!dif %in% c("none", "inconsistent") & !startsWith(dif, "favors ")])
  abort_if(length(bad) > 0, "Unknown DIF value(s): ", paste(bad, collapse = ", "),
           ". Use: none, inconsistent, or \"favors \" and a group (e.g. favors girls).")
  dif
}

# 每個 DIF 的值是哪一種樣式（none、up、down、inconsistent），照 dif_direction() 的方向：
# sides 是題目表記下的兩個標籤（up、down；kits_item_table() 的 attr "dif_labels"）。
# favors girls、favors boys 固定（up、down）；沒有記錄的其他組照名字排，第一個是 up
dif_style_keys <- function(values, sides = NULL) {
  values <- unique(values)
  keys <- stats::setNames(rep(NA_character_, length(values)), values)
  keys[values %in% c("none", "inconsistent")] <- values[values %in% c("none", "inconsistent")]
  keys[values == "favors girls"] <- "up"
  keys[values == "favors boys"] <- "down"
  if (!is.null(sides)) {
    for (k in c("up", "down")) keys[is.na(keys) & values == sides[[k]]] <- k
  }
  rest <- sort(values[is.na(keys)])
  keys[rest] <- c("up", "down", rep("up", max(0, length(rest) - 2)))[seq_along(rest)]
  keys
}

# 圖例的文字：男女和 none、inconsistent 用 LAB$dif；其他組例如 "Favors F"
dif_legend <- function(value) {
  if (value %in% names(LAB$dif)) return(LAB$dif[[value]])
  paste0(toupper(substr(value, 1, 1)), substring(value, 2))
}

# 畫的順序：none、up、down、inconsistent
dif_draw_order <- function(keys) {
  keys[order(match(keys, names(COL_DIF)), names(keys))]
}

# D 軸的範圍：有給就用；沒給就用所有格的資料算一個共同範圍（取到 10 的倍數）
shared_range <- function(values, y_range) {
  if (!is.null(y_range)) return(y_range)
  values <- values[is.finite(values)]
  abort_if(length(values) == 0, "No finite values to plot.")
  r <- range(values)
  c(floor(r[1] / 10) * 10, ceiling(r[2] / 10) * 10)
}

# 檢查領域的代碼在資料裡（domain 欄）
need_domains <- function(domains, data = NULL, arg = "data") {
  abort_if(length(domains) == 0, "`scores` is empty.")
  if (!is.null(data)) {
    none <- setdiff(domains, data$domain)
    abort_if(length(none) > 0, "No rows in `", arg, "` for domain(s): ", paste(none, collapse = ", "), ".")
  }
}

# 成長曲線和題目（內部用，給 plot(kits_item_table(), type = "growth")）：
# 一個領域一格，灰色的線是常模的 P10、P50、P90，點是題目（顏色和形狀是 DIF）
growth_plot <- function(norm, items, domains, age_ticks = AGE_TICKS, y_range = NULL) {
  rp <- rp_of(items)
  need_cols(norm, c("domain", "age_m", "P10", "P50", "P90"), "norm")
  need_cols(items, c("item_id", "domain", "x", "D", "empirical", "DIF", "stem", age_col(rp),
                     "age_reach"), "items")
  need_domains(domains, norm, "norm")
  items$DIF <- clean_dif(items$DIF)
  keys <- dif_style_keys(items$DIF, attr(items, "dif_labels"))
  norm <- filter(as.data.frame(norm), domain %in% domains)
  items <- filter(as.data.frame(items), domain %in% domains)
  y_range <- shared_range(c(norm$P10, norm$P90, items$D), y_range)
  no_x <- filter(items, is.na(x))$item_id
  if (length(no_x) > 0) {
    kits_inform("Not drawn (no age on the chart; harder than the highest norm median and no age at 50% passing): ",
                paste(no_x, collapse = ", "), class = "kits_message_data")
  }
  # 圖例每一類只出現一次：百分位線在第一格，DIF 的類別在第一個有這一類的格
  first_panel <- items |>
    filter(!is.na(x)) |>
    summarise(panel = min(match(domain, domains)), .by = DIF)
  legend_dif <- stats::setNames(first_panel$panel, first_panel$DIF)

  panels <- imap(domains, \(sc, panel) {
    growth_panel(filter(norm, domain %in% sc), filter(items, domain %in% sc), sc, rp, y_range,
                 legend_lines = panel == 1, legend_dif = names(legend_dif)[legend_dif == panel],
                 keys = keys)
  })

  rows <- ceiling(length(domains) / 2)
  plotly::subplot(panels, nrows = rows, shareX = TRUE, titleY = TRUE, margin = 0.05) |>
    plotly::layout(xaxis = list(title = LAB$age, tickvals = age_ticks),
                   xaxis2 = list(title = LAB$age, tickvals = age_ticks)) |>
    kits_style(sprintf(LAB$title_growth, 100 * rp, if (rp != 0.5) sprintf(" (items at RP%g)", 100 * rp) else ""))
}


# 成長曲線的一格：常模的 P10、P50、P90，和這個領域的題目（顏色、形狀是 DIF）
growth_panel <- function(norm, items, sc, rp, y_range, legend_lines, legend_dif, keys) {
  segments <- norm_segments(norm)
  drawn <- filter(items, !is.na(x))
  p <- plotly::plot_ly()
  for (s in names(segments)) {
    d <- segments[[s]]
    for (q in c("P10", "P50", "P90")) {
      p <- plotly::add_lines(p, x = d$age_m, y = d[[q]],
                             line = list(color = PAL$gray, width = if (q == "P50") LINE_W$bold else LINE_W$thin),
                             name = q, legendgroup = q, showlegend = legend_lines && s == names(segments)[1], hoverinfo = "skip")
    }
  }
  keys <- dif_draw_order(keys)
  for (lv in names(keys)) {
    d <- filter(drawn, DIF == lv)
    if (nrow(d) == 0) next
    key <- keys[[lv]]
    symbol <- if_else(d$empirical, SYM_DIF[[key]], paste0(SYM_DIF[[key]], "-open"))
    p <- plotly::add_markers(p, x = d$x, y = d$D,
                             name = dif_legend(lv), legendgroup = lv, showlegend = lv %in% legend_dif,
                             marker = list(color = COL_DIF[[key]], symbol = symbol, size = MARKER_SIZE),
                             hovertext = item_hover(d, rp), hoverinfo = "text")
  }
  note <- not_drawn_note(filter(items, is.na(x))$item_id, LAB$not_drawn)
  y_axis <- list(title = LAB$D, range = y_range)
  plotly::layout(p, annotations = c(domain_panel_title(sc), note), yaxis = y_axis)
}

# 成長階梯（內部用，給 plot(kits_item_table(), type = "ladder")）：
# 一欄一個領域，每題一條橫線在它的 D；右邊的刻度是常模中位數到達這個 D 的月齡
ladder_plot <- function(items, columns, marks) {
  rp <- rp_of(items)
  need_cols(items, c("item_id", "domain", "column", "D", "DIF", "stem", age_col(rp), "age_reach"),
            "items")
  need_cols(marks, c("age_m", "D"), "marks")
  it <- filter(as.data.frame(items), column %in% columns)
  abort_if(nrow(it) == 0, "No items in `items` have a `column` in `columns`.")
  keys <- dif_draw_order(dif_style_keys(clean_dif(it$DIF), attr(items, "dif_labels")))
  it <- it |>
    mutate(DIF = clean_dif(DIF),
           x = match(column, columns),
           label = paste0(sprintf("%.1f", D), "  ", short_stem(stem)))

  p <- plotly::plot_ly()
  for (lv in names(keys)) {
    d <- filter(it, DIF == lv)
    if (nrow(d) == 0) next
    key <- keys[[lv]]
    p <- p |>
      plotly::add_trace(x = as.vector(rbind(d$x - 0.2, d$x + 0.02, NA)), y = as.vector(rbind(d$D, d$D, NA)),
                        type = "scatter", mode = "lines",
                        line = list(color = COL_DIF[[key]], width = if (lv == "none") LINE_W$main else LINE_W$bold + 1),
                        name = dif_legend(lv), legendgroup = lv, hoverinfo = "skip") |>
      plotly::add_trace(x = d$x + 0.05, y = d$D, type = "scatter", mode = "text",
                        text = d$label, textposition = "middle right",
                        textfont = list(size = 10, color = COL_DIF[[key]]),
                        hovertext = item_hover(d, rp), hoverinfo = "text",
                        legendgroup = lv, showlegend = FALSE)
  }
  for (i in seq_len(nrow(marks))) {
    p <- plotly::add_lines(p, x = c(0.5, length(columns) + 0.6), y = rep(marks$D[i], 2),
                           line = list(color = PAL$light, width = LINE_W$thin,
                                       dash = "dot"), hoverinfo = "skip", showlegend = FALSE)
  }
  mark_text <- map2(marks$D, marks$age_m, \(D, age) {
    list(x = length(columns) + 0.62, y = D, text = paste0(age, LAB$month),
         showarrow = FALSE, xanchor = "left", font = list(size = 10, color = PAL$gray))
  })

  plotly::layout(p,
                 xaxis = list(tickvals = seq_along(columns), ticktext = lab_name("column", columns),
                              range = c(0.4, length(columns) + 0.9), showgrid = FALSE, zeroline = FALSE),
                 yaxis = list(title = paste0(LAB$D_step, if (rp != 0.5) sprintf(", RP%g", 100 * rp) else ""), dtick = 10),
                 annotations = mark_text) |>
    kits_style(LAB$title_ladder)
}


# 一格的標題：領域的名稱；沒有領域（一個領域）時沒有標題
domain_panel_title <- function(sc) {
  if (is.na(sc)) {
    return(list())
  }
  panel_title(lab_name("domain", sc))
}


# 常模的格子（內部用，給 plot() 照組的常模）：一列一個領域，一欄一組（性別）。
# norms 有 domain、sex 兩欄（沒有時 NA）和百分位的欄；最外面的兩條百分位畫虛線，P50 畫粗線
norm_grid_plot <- function(norms, scores, sexes, pct, age_ticks, y_range, title, y_title,
                           overlay = NULL) {
  if (is.null(scores)) scores <- unique(norms$domain)
  need_domains(scores, norms, "norms")
  norms <- as.data.frame(norms)
  empty <- expand_grid(sex = sexes, domain = scores) |>   # 領域變得最快，和錯誤訊息的順序相同
    anti_join(norms, by = c("domain", "sex"))
  abort_if(nrow(empty) > 0, "No rows in `norms` for: ", paste(empty$domain, empty$sex, collapse = ", "), ".")
  cols <- paste0("P", pct)
  shown <- filter(norms, domain %in% scores, sex %in% sexes)
  y_range <- shared_range(c(shown[[cols[1]]], shown[[cols[length(cols)]]]), y_range)
  grid <- expand_grid(sc = scores, i = seq_along(sexes))
  panels <- map2(grid$sc, grid$i, \(sc, i) {
    norm_grid_panel(norms, sc, sexes, i, cols, sc %in% scores[1] && i == 1, y_range, y_title,
                    overlay)
  })
  margin <- c(0.04, 0.04, 0.03, 0.03)
  plotly::subplot(panels, nrows = length(scores), shareX = TRUE, titleY = TRUE, margin = margin) |>
    plotly::layout(xaxis = list(title = LAB$age, tickvals = age_ticks),
                   xaxis2 = list(title = LAB$age, tickvals = age_ticks)) |>
    kits_style(title)
}

# 一格：一個領域、一組的百分位曲線（每一段各畫一次）
norm_grid_panel <- function(norms, sc, sexes, i, cols, legend, y_range, y_title, overlay = NULL) {
  sx <- sexes[i]
  colour <- COL_MANY[(i - 1) %% length(COL_MANY) + 1]
  segments <- norm_segments(filter(norms, domain %in% sc, sex %in% sx))
  outer <- cols[c(1, length(cols))]
  p <- plotly::plot_ly()
  if (!is.null(overlay$children)) {   # 這一組的孩子，畫在曲線後面
    own <- filter(overlay$children, domain %in% sc, sex %in% sx)
    p <- add_children(p, own, overlay$persons, overlay$min_cell)
  }
  for (s in names(segments)) {
    d <- segments[[s]]
    for (q in cols) {
      hover <- sprintf("%s, %s, %.1f%s: D = %.1f", norm_panel_name(sc, sx), q, d$age_m, LAB$month,
                       d[[q]])
      line <- list(color = colour, width = if (q == "P50") LINE_W$bold else LINE_W$thin,
                   dash = if (q %in% outer && length(cols) > 1) "dot" else "solid")
      p <- plotly::add_lines(p, x = d$age_m, y = d[[q]], name = q, legendgroup = q,
                             showlegend = legend && s == names(segments)[1], line = line,
                             hovertext = hover, hoverinfo = "text")
    }
  }
  y_axis <- list(title = y_title, range = y_range)
  plotly::layout(p, annotations = norm_panel_title(sc, sx), yaxis = y_axis)
}

# 一格的名稱：「領域, 性別」（滑鼠停的文字）；沒有的部分不寫
norm_panel_name <- function(sc, sx) {
  parts <- c(if (!is.na(sc)) lab_name("domain", sc), if (!is.na(sx)) lab_name("sex", sx))
  paste(parts, collapse = ", ")
}

# 一格的標題：「領域 (性別)」；只有一個常模（沒有領域、性別）時沒有標題
norm_panel_title <- function(sc, sx) {
  if (is.na(sc) && is.na(sx)) {
    return(list())
  }
  text <- if (is.na(sc)) lab_name("sex", sx) else lab_name("domain", sc)
  if (!is.na(sc) && !is.na(sx)) text <- paste0(text, " (", lab_name("sex", sx), ")")
  panel_title(text)
}

# 存 html 要 pandoc（htmlwidgets 的 selfcontained = TRUE）：找不到時先停下來，說怎麼設定
need_pandoc <- function() {
  found <- if (requireNamespace("rmarkdown", quietly = TRUE)) rmarkdown::pandoc_available() else nzchar(Sys.which("pandoc"))
  abort_if(!found, "Saving html needs pandoc, which was not found. RStudio and Quarto come with pandoc: ",
           "set the environment variable RSTUDIO_PANDOC to the folder that holds it, for example ",
           paste0("Sys.setenv(RSTUDIO_PANDOC = \"/usr/lib/rstudio/resources/app/bin/quarto/bin/tools\") or the tools folder ",
                  "of Quarto "),
           paste0("(on Windows \"C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools\"), or install pandoc from ",
                  "https://pandoc.org. "), "Then call kits_save_plot() again.")
}

#' Save a plot: interactive html and static png
#'
#' The html is an interactive plot (zoom, hover); the png is a browser screenshot of the same plot.
#' Saving html needs pandoc (RStudio and Quarto come with it); when it is not found, the function
#' stops (class `kits_error_input`) and says how to set `RSTUDIO_PANDOC`.
#' Saving png needs the webshot2 package and a Chromium-based browser (its path in
#' `CHROMOTE_CHROME`).
#'
#' @param p A plotly plot.
#' @param html The html file name; it must end in `.html`.
#' @param png The png file name. Omit it when no png is wanted.
#' @param width,height Width and height of the png (pixels).
#' @return The saved file names (invisibly).
#' @examples
#' set.seed(1)
#' s <- data.frame(cohort = rep(c("M03", "M36"), each = 200), wave_m = 36,
#'                 theta = rnorm(400, 68, 4))
#' p <- plot(kits_vertical(s), type = "median")
#' \donttest{
#' # Save html only (needs pandoc)
#' if (rmarkdown::pandoc_available()) kits_save_plot(p, file.path(tempdir(), "cohort.html"))
#' }
#' \dontrun{
#' # Also save png (needs webshot2 and a browser)
#' kits_save_plot(p, file.path(tempdir(), "cohort.html"), file.path(tempdir(), "cohort.png"),
#'           width = 800, height = 420)
#' }
#' @export
kits_save_plot <- function(p, html, png = NULL, width = 1200, height = 800) {
  abort_if(!inherits(p, "plotly"), "`p` must be a plotly plot (for example from plot()).")
  abort_if(!grepl("\\.html$", html), "`html` must end with .html")
  need_pandoc()
  name <- sub("\\.html$", "", basename(html))
  htmlwidgets::saveWidget(p, html, selfcontained = TRUE, title = name)
  unlink(sub("\\.html$", "_files", html), recursive = TRUE)
  if (is.null(png)) {
    return(invisible(html))
  }

  abort_if(!requireNamespace("webshot2", quietly = TRUE),
           "Saving png needs the webshot2 package: install.packages(\"webshot2\")")
  still <- p |>
    plotly::layout(margin = list(t = 90)) |>
    plotly::config(displayModeBar = FALSE)
  still$width <- width     # 圖的大小設在 htmlwidget 上（plotly 不再建議用 layout() 設）
  still$height <- height
  tmp <- file.path(tempdir(), paste0(name, ".html"))
  htmlwidgets::saveWidget(still, tmp, selfcontained = TRUE)
  webshot2::webshot(tmp, png, vwidth = width, vheight = height, zoom = 2, quiet = TRUE)
  invisible(c(html, png))
}
