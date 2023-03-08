library(ggplot2)

dat <- read.table(text="Program School1  School2 School3
 Program1        1        1       1
 Program2        1        0       1 
 Program3        1        1       0", header=TRUE, stringsAsFactors=FALSE)

dat_long <- reshape2::melt(dat)

# discrete vs continuous
dat_long$value <- factor(dat_long$value)

gg <- ggplot(dat_long)

# fill + legend, gray border
gg <- gg + geom_tile(aes(x=Program, y=variable, fill=value), color="#7f7f7f")

# custom fill colors
gg <- gg + scale_fill_manual(values=c("white", "black"))

# squares
gg <- gg + coord_equal()

# no labels
gg <- gg + labs(x=NULL, y=NULL)

# remove some chart junk
gg <- gg + theme_bw()
gg <- gg + theme(panel.grid=element_blank())
gg <- gg + theme(panel.border=element_blank())
gg
####################################################################

dat <- matrix(c(0,1,1,0,1,0,0,1), nrow=2)
dat <- dat[nrow(dat):1, ]
colnames(dat) <- colnames(dat, do.NULL = FALSE, prefix = "f")
rownames(dat) <- rownames(dat, do.NULL = FALSE, prefix = "ind")
#df <- as.data.frame(dat)
df <- reshape2::melt(dat)

df$value <- factor(df$value)

gg <- ggplot(df)

# fill + legend, gray border
gg <- gg + geom_tile(aes(x=Var2, y=Var1, fill=value), color="#7f7f7f")

# custom fill colors
gg <- gg + scale_fill_manual(values=c("white", "black"))

# squares
gg <- gg + coord_equal()

# no labels
gg <- gg + labs(x=NULL, y=NULL)

# remove some chart junk
gg <- gg + theme_bw()
gg <- gg + theme(panel.grid=element_blank())
gg <- gg + theme(panel.border=element_blank())
gg
