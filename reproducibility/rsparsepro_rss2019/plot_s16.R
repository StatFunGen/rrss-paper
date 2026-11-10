plot_s16 <- function(comparison, out, preview = "") {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  rrss_lambda <- jsonlite::fromJSON(comparison)
  stopifnot(rrss_lambda$summary$cases == 36L,
            rrss_lambda$summary$maximum_alpha_difference < 1e-10,
            rrss_lambda$summary$identical_PIP_cases == 36L,
            nrow(rrss_lambda$single_effect) == 1296L,
            nrow(rrss_lambda$variant_pip) == 432L)
  rrss_draw_lambda <- function() {
    colors <- c(identity = "#0072B2", low_rank = "#009E73", adsp = "#D55E00")
    labels <- c(adsp = "ADSP sketch", low_rank = "Singular LD", identity = "Identity LD")
    commit_label <- "(susieR, GitHub commit ID: 9837e100, March 21, 2019)"
    layout(matrix(c(1, 2), nrow = 1))
    par(family = "Arial", las = 1, mgp = c(2.1, 0.5, 0), tcl = -0.22,
        oma = c(0, 0.3, 0.5, 0), cex.axis = 0.88, xaxs = "i", yaxs = "i")
    draw_pair <- function(d, letter, title, xlab, ylab) {
      par(mar = c(4.9, 4.15, 2.05, 0.65), pty = "s")
      plot(NA, xlim = c(-0.025, 1.025), ylim = c(-0.025, 1.025), asp = 1,
           xlab = "", ylab = "", axes = FALSE)
      axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
      axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
      box(bty = "l", lwd = 0.7)
      abline(0, 1, lty = 2, col = "#828282", lwd = 0.8)
      for (kind in names(colors)) {
        x <- d[d$kind == kind, ]
        is_adsp <- kind == "adsp"
        points(x$rss2019, x$rsparsepro,
               pch = if (is_adsp) 21 else if (kind == "low_rank") 17 else 16,
               cex = if (is_adsp) 0.88 else 0.65, lwd = 0.65,
               col = adjustcolor(colors[kind], alpha.f = if (is_adsp) 0.8 else 0.38),
               bg = adjustcolor(colors[kind], alpha.f = 0.28))
      }
      mtext(xlab, side = 1, line = 2.0, cex = 0.92)
      mtext(commit_label, side = 1, line = 3.35, cex = 0.76)
      mtext(ylab, side = 2, line = 2.3, las = 0, cex = 0.92)
      mtext(title, side = 3, line = 0.5, cex = 1.04, font = 2)
      mtext(letter, side = 3, line = 1.3, adj = 0, cex = 1.17, font = 2)
    }
    draw_pair(rrss_lambda$single_effect, "a", "Single effect PIP",
              expression("RSS-" * lambda * " (" * alpha[jl] * ")"),
              expression("RSparsePro (" * gamma[kg]^"*" * ")"))
    legend("topleft", inset = c(0.01, 0.01),
           legend = labels[c("adsp", "low_rank", "identity")],
           col = colors[c("adsp", "low_rank", "identity")],
           pch = c(21, 17, 16), pt.bg = adjustcolor(colors["adsp"], 0.28),
           pt.cex = c(0.88, 0.7, 0.7), cex = 0.83, bty = "n",
           y.intersp = 1.0, x.intersp = 0.55)
    draw_pair(rrss_lambda$variant_pip, "b", "Variant PIP", expression("RSS-" * lambda), "RSparsePro")

  }
  cairo_pdf(file.path(out, "Figure_S16.pdf"), width = 7.3, height = 3.1, family = "Arial", pointsize = 10)
  rrss_draw_lambda()
  dev.off()
  if (nzchar(preview)) {
    png(preview, width = 2190, height = 930, res = 300, type = "cairo", family = "Arial", pointsize = 10)
    rrss_draw_lambda()
    dev.off()
  }
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  plot_s16(args[1], args[2])
}
