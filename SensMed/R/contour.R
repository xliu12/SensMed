

# from dml.sensemakr https://github.com/carloscinelli/dml.sensemakr/blob/main/R/plot.R
contour_plot <- function(grid_values.x,
                         grid_values.y,
                         z_axis,
                         levels = NULL,
                         labels = NULL,
                         threshold = 0,
                         nlevels = 10,
                         grid.number = 70,
                         col.contour = "grey40",
                         col.thr.line = "red",
                         xlab = NULL,
                         ylab = NULL,
                         round = 0,
                         cex.lab = 0.8,
                         cex.axis = 0.8,
                         cex.main = 1,
                         asp = 1,
                         list.par = NULL){
    
    if(is.null(xlab)) xlab <- expression(paste("1-",R[alpha%~%alpha[s]]^2))
    if(is.null(ylab)) ylab <- expression(paste(R[Y-theta[s]%~%theta-theta[s]]^2))
    # if(is.null(ylab)) ylab <- expression(paste(eta[Y%~%A~"|"~DX]^2))
    
    default_levels <- pretty(range(z_axis), nlevels)
    too_close      <- abs(default_levels - threshold) < min(diff(default_levels)) * 0.25
    line_color     <- ifelse(too_close, "transparent", col.contour)
    line_type      <- ifelse(too_close, 1, 1)
    line_width     <- ifelse(too_close, 1, 1)
    
    # Plot contour plot:
    if (is.null(list.par)) {
        oldpar <- par(mar = c(5, 5, 4, 1) + .1, pty = "s")
        on.exit(par(oldpar))
    } else {
        if (!is.list(list.par)) stop("list.par needs to be a named list")
        oldpar <- do.call("par", list.par)
        on.exit(par(oldpar))
    }
    
    # contour ----
    contour(grid_values.x, grid_values.y, z_axis,
            nlevels = nlevels,
            levels = levels,
            labels = labels,
            xlab = xlab,
            ylab = ylab,
            cex.lab = cex.lab,
            cex.axis = cex.axis,
            cex.main = cex.main,
            asp = asp,
            col = line_color,
            lty = line_type,
            lwd = line_width)
    
    contour(grid_values.x, grid_values.y, z_axis,
            level = threshold,
            label = round(threshold, digits = round),
            add = TRUE,
            col = col.thr.line,
            lwd = 2,
            lty = 2,
            cex.lab = cex.lab,
            cex.axis = cex.axis,
            cex.main = cex.main,
            asp = asp)
}
