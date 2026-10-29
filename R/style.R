# 所有圖共用的樣式：matplotlib.pyplot 的樣子。
# 白底、四邊黑框、刻度朝外、淺灰格線；主色是藍和橘（Okabe-Ito 的色盲友善色）。
# 圖例沒有框；圖例的字太長時放在圖的下方。
# 每個畫圖函數最後都交給 kits_style()，所以字型、框線、標題、圖例都一樣。

PAL <- list(blue = "#0072B2", orange = "#D55E00", black = "#000000", gray = "#7F7F7F", light = "#C7C7C7")
FONT <- list(family = "DejaVu Sans, Noto Sans CJK TC, Arial, sans-serif", size = 12, color = "#000000")
LINE_W <- list(thin = 1, main = 1.5, bold = 2.5)
MARKER_SIZE <- 7

# 一條軸的樣式；已經設定的（標題、範圍、刻度）優先
mpl_axis <- function(axis) {
  style <- list(showline = TRUE, mirror = TRUE, linecolor = PAL$black, linewidth = 1,
                ticks = "outside", ticklen = 5, tickwidth = 1, tickcolor = PAL$black,
                showgrid = TRUE, gridcolor = "#E5E5E5", gridwidth = 1, zeroline = FALSE)
  utils::modifyList(style, if (is.null(axis)) list() else axis)
}

# 套用樣式：每一條軸、背景、字型、標題（置中）、圖例（淺灰框）
kits_style <- function(p, title = NULL) {
  built <- plotly::plotly_build(p)
  lay <- built$x$layout
  for (ax in grep("^[xy]axis[0-9]*$", names(lay), value = TRUE)) {
    lay[[ax]] <- mpl_axis(lay[[ax]])
  }
  # 每一格畫一個黑框（共用軸的格子沒有自己的軸線，所以用框補上）
  axis_or <- function(axis, default) if (is.null(axis)) default else axis
  panels <- unique(data.frame(x = map_chr(built$x$data, \(tr) axis_or(tr$xaxis, "x")),
                              y = map_chr(built$x$data, \(tr) axis_or(tr$yaxis, "y"))))
  boxes <- map2(panels$x, panels$y, \(xa, ya) {
    list(type = "rect", x0 = 0, x1 = 1, y0 = 0, y1 = 1, xref = paste(xa, "domain"), yref = paste(ya, "domain"),
         line = list(color = PAL$black, width = 1), fillcolor = "rgba(0,0,0,0)", layer = "above")
  })
  lay$shapes <- c(if (is.null(lay$shapes)) list() else lay$shapes, boxes)
  lay$paper_bgcolor <- "#FFFFFF"
  lay$plot_bgcolor <- "#FFFFFF"
  lay$font <- FONT
  shown <- map_lgl(built$x$data, \(tr) !identical(tr$showlegend, FALSE) && !is.null(tr$name))
  legend_names <- map_chr(built$x$data[shown], \(tr) as.character(tr$name)[1])
  longest <- max(c(0, nchar(legend_names)))
  legend <- list(borderwidth = 0, bgcolor = "rgba(255,255,255,0)")
  if (longest > 28) {   # 字太長：圖例放在下方，橫排
    legend <- c(legend, list(orientation = "h", x = 0, xanchor = "left", y = -0.22, yanchor = "top"))
    lay$margin <- utils::modifyList(if (is.null(lay$margin)) list() else lay$margin, list(b = 130))
  }
  lay$legend <- utils::modifyList(legend, if (is.null(lay$legend)) list() else lay$legend)
  if (!is.null(title)) {
    lines <- lengths(regmatches(title, gregexpr("<br>", title))) + 1
    # 標題放在工具列（右上角，約 30 px）下面
    lay$title <- list(text = title, x = 0.5, xanchor = "center", xref = "paper", y = 1, yanchor = "top", yref = "container",
                      pad = list(t = 34), font = list(size = 14))
    lay$margin <- utils::modifyList(if (is.null(lay$margin)) list() else lay$margin, list(t = 56 + 20 * lines))
  }
  built$x$layout <- lay
  built
}
