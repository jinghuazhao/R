# 10-9-1996 Terry Therneau to S-News (StatLib)

hsrcolor <- function(printer = "har7cps1", color)
{
        if(missing(color))
                color <- c("black", "red", "green", "blue", "magenta", "red4", 
                        "orange", "DarkGreen", "cyan2", "DarkViolet")
        temp <- match(color, dimnames(ps.colors.rgb)[[1]])
        if(any(is.na(temp)))
                stop("Invalid color name")
        if(is.null(printer) || printer == "") {
                ps.options(colors = ps.colors.rgb[color,  ])
                if(dev.cur() > 1)
                        ps.options.send(colors = ps.colors.rgb[color,  ])
        }
        else {
                cmd <- paste("/usr/local/bin/splusprint -d", printer)
                ps.options(colors = ps.colors.rgb[color,  ], command = cmd)
                if(dev.cur() > 1)
                        ps.options.send(colors = ps.colors.rgb[color,  ], 
                                command = cmd)
        }
        invisible(color)
}
