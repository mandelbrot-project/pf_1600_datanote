library(circlize)

set.seed(999)
n = 1000
df = data.frame(sectors = sample(letters[1:8], n, replace = TRUE),
    x = rnorm(n), y = runif(n))

circos.par("track.height" = 0.1)
circos.initialize(df$sectors, x = df$x)

circos.track(df$sectors, y = df$y,
    panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, 
            CELL_META$cell.ylim[2] + mm_y(5), 
            CELL_META$sector.index)
        circos.axis(labels.cex = 0.6)
})
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(df$sectors, df$x, df$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "what", sector.index = "a", track.index = 1)

bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(df$sectors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)



circos.track(df$sectors, x = df$x, y = df$y,
    panel.fun = function(x, y) {
        ind = sample(length(x), 10)
        x2 = x[ind]
        y2 = y[ind]
        od = order(x2)
        circos.lines(x2[od], y2[od])
})

circos.update(sector.index = "d", track.index = 2, 
    bg.col = "#FF8080", bg.border = "black")
circos.points(x = -2:2, y = rep(0.5, 5), col = "white")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "updated", col = "white")


circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 0.1)
    n_breaks = length(breaks)
    circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
                breaks[-1], rep(ylim[2], n_breaks - 1),
                col = rand_color(n_breaks), border = NA)
})

circos.link("a", 0, "b", 0, h = 0.4)
circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red",
    border = "blue", h = 0.2)
circos.link("e", 0, "g", c(-1,1), col = "green", border = "black", lwd = 2, lty = 2)

circos.clear()


