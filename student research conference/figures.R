library(grid)
x <- rep(seq(0, 0.8, 0.2), 5)
y <- rep(seq(0, 0.8, 0.2), each=5)
pdf("blank_grid.pdf")
grid.newpage()
vp1 <- viewport(x=0.1, y=0.1, w=0.8, h=0.8, 
                just=c("left", "bottom"), name="vp1")
grid.rect(x=x, y=y, height=1/5, width=1/5, hjust=0, vjust=0, vp=vp1, 
          gp=gpar(col=1))
dev.off()

pdf("grid_spot.pdf")
grid.newpage()
vp1 <- viewport(x=0.1, y=0.1, w=0.8, h=0.8, 
                just=c("left", "bottom"), name="vp1")
grid.rect(x=x, y=y, height=1/5, width=1/5, hjust=0, vjust=0, vp=vp1, 
          gp=gpar(col=1))
grid.points(x=0.74, y=0.629, pch=19, gp=gpar(col="red"), vp=vp1)
dev.off()
