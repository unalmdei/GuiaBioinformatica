library("xkcd")
library("showtext")
library("sysfonts")
library("tibble")

head(DNase)
plot(DNase$conc, DNase$density)


plot(DNase$conc, DNase$density,
     ylab = attr(DNase, "labels")$y,
     xlab = paste(attr(DNase, "labels")$x, attr(DNase, "units")$x),
     pch = 3,
     col = "blue")


hist(DNase$density, breaks=25, main = "")
boxplot(density ~ Run, data = DNase)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Hiiragi2013")
library("Hiiragi2013")
data("x")
dim(Biobase::exprs(x))


head(pData(x), n = 2)


library("dplyr")
groups = group_by(pData(x), sampleGroup) %>%
  summarise(n = n(), color = unique(sampleColour))
groups






library("ggplot2")
ggplot(DNase, aes(x = conc, y = density)) + geom_point()


ggplot(groups, aes(x = sampleGroup, y = n)) +
  geom_bar(stat = "identity")


stopifnot(formals(ggplot2::geom_bar)$stat=="count")


groupColor = setNames(groups$color, groups$sampleGroup)


ggplot(groups, aes(x = sampleGroup, y = n, fill = sampleGroup)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = groupColor, name = "Groups") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




gg = ggplot(DNase, aes(x = conc, y = density)) + geom_point()


gg
print(gg)


ggplot2::ggsave("DNAse-histogram-demo.pdf", plot = gg)


file.remove("DNAse-histogram-demo.pdf")


library("mouse4302.db")




dftx = data.frame(t(Biobase::exprs(x)), pData(x))
ggplot( dftx, aes( x = X1426642_at, y = X1418765_at)) +
  geom_point( shape = 1 ) +
  geom_smooth( method = "loess" )


stopifnot(is(dftx, "data.frame"))


ggplot( dftx, aes( x = X1426642_at, y = X1418765_at ))  +
  geom_point( aes( color = sampleColour), shape = 19 ) +
  geom_smooth( method = "loess" ) +
  scale_color_discrete( guide = "none" )




library("mouse4302.db")


AnnotationDbi::select(mouse4302.db,
                      keys = c("1426642_at", "1418765_at"), keytype = "PROBEID",
                      columns = c("SYMBOL", "GENENAME"))


dfx = as.data.frame(Biobase::exprs(x))
ggplot(dfx, aes(x = `20 E3.25`)) + geom_histogram(binwidth = 0.2)


pb = ggplot(groups, aes(x = sampleGroup, y = n))


class(pb)
pb


pb = pb + geom_bar(stat = "identity")
pb = pb + aes(fill = sampleGroup)
pb = pb + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pb = pb + scale_fill_manual(values = groupColor, name = "Groups")
pb


pb.polar = pb + coord_polar() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") + ylab("")
pb.polar


selectedProbes = c( Fgf4 = "1420085_at", Gata4 = "1418863_at",
                    Gata6 = "1425463_at",  Sox2 = "1416967_at")




selectedProbes2 = AnnotationDbi::select(mouse4302.db,
                                        keys = selectedProbes, keytype = "PROBEID", columns = c("SYMBOL"))
stopifnot(identical(sort(selectedProbes2$SYMBOL), sort(names(selectedProbes))),
          all(selectedProbes[selectedProbes2$SYMBOL] == selectedProbes2$PROBEID))


library("reshape2")
genes = melt(Biobase::exprs(x)[selectedProbes, ],
             varnames = c("probe", "sample"))


genes$gene =
  names(selectedProbes)[match(genes$probe, selectedProbes)]
head(genes)


ggplot(genes, aes(x = gene, y = value)) +
  stat_summary(fun = mean, geom = "bar")


library("Hmisc")
ggplot(genes, aes( x = gene, y = value, fill = gene)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               width = 0.25)


p = ggplot(genes, aes( x = gene, y = value, fill = gene))
p + geom_boxplot()


p + geom_dotplot(binaxis = "y", binwidth = 1/6,
                 stackdir = "center", stackratio = 0.75,
                 aes(color = gene))
library("ggbeeswarm")
p + geom_beeswarm(aes(color = gene))


ggplot(genes, aes( x = value, color = gene)) + geom_density()


p + geom_violin()


library("ggridges")
ggplot(genes, aes(x = value, y = gene, fill = gene)) + 
  geom_density_ridges()


top42 = order(rowMeans(Biobase::exprs(x)), decreasing = TRUE)[1:42]
g42 = melt(Biobase::exprs(x)[rev(top42), ], varnames = c("probe", "sample"))
ggplot(g42, aes(x = value, y = probe)) 


ggplot(g42, aes(x = value, y = probe)) + 
  geom_density_ridges() + theme(legend.position = "none",
                                axis.title.y = element_blank(), axis.text.y = element_blank(),
                                axis.ticks.y = element_blank()) + xlim(13, 15) 


simdata = rnorm(70)
tibble(index = seq(along = simdata),
       sx = sort(simdata)) %>%
  ggplot(aes(x = sx, y = index)) + geom_step()


ggplot(genes, aes( x = value, color = gene)) + stat_ecdf()






ggplot(dfx, aes(x = `64 E4.5 (EPI)`)) + geom_histogram(bins = 100)
ggplot(dfx, aes(x = 2 ^ `64 E4.5 (EPI)`)) + 
  geom_histogram(binwidth = 20) + xlim(0, 1500)


scp = ggplot(dfx, aes(x = `59 E4.5 (PE)` ,
                      y = `92 E4.5 (FGF4-KO)`))
scp + geom_point()


scp  + geom_point(alpha = 0.1)


scp + geom_density2d()


scp + geom_density2d(h = 0.5, bins = 60)
library("RColorBrewer")
colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))

scp + stat_density2d(h = 0.5, bins = 60,
                     aes( fill = ..level..), geom = "polygon") +
  colorscale + coord_fixed()


scp + geom_hex() + coord_fixed()
scp + geom_hex(binwidth = c(0.2, 0.2)) + colorscale +
  coord_fixed()


library("ggthemes")
sunsp = tibble(year   = time(sunspot.year),
               number = as.numeric(sunspot.year))
sp = ggplot(sunsp, aes(x = year, y = number)) + geom_line()
sp
ratio = with(sunsp, bank_slopes(year, number))
sp + coord_fixed(ratio = ratio)






library("magrittr")
dftx$lineage %<>% sub("^$", "no", .)
dftx$lineage %<>% factor(levels = c("no", "EPI", "PE", "FGF4-KO"))

ggplot(dftx, aes( x = X1426642_at, y = X1418765_at)) +
  geom_point() + facet_grid( . ~ lineage )


ggplot( dftx,
        aes( x = X1426642_at, y = X1418765_at)) + geom_point() +
  facet_grid( Embryonic.day ~ lineage )


ggplot(mutate(dftx, Tdgf1 = cut(X1450989_at, breaks = 4)),
       aes( x = X1426642_at, y = X1418765_at)) + geom_point() +
  facet_wrap( ~ Tdgf1, ncol = 2 )






data("volcano")
volcanoData = list(
  x = 10 * seq_len(nrow(volcano)),
  y = 10 * seq_len(ncol(volcano)),
  z = volcano,
  col = terrain.colors(500)[cut(volcano, breaks = 500)]
)
library("rgl")
with(volcanoData, persp3d(x, y, z, color = col))


rgl.snapshot(filename = file.path("figure", "chap3-r_rglvolcano.png"))


.volcanocut = cut(volcano, breaks = 500)
stopifnot(!any(is.na(.volcanocut)), all(as.integer(.volcanocut) %in% 1:500))


par(mai = rep(0,4))
pie(rep(1, 8), col=1:8)




display.brewer.all()


head(brewer.pal.info)
table(brewer.pal.info$category)


brewer.pal(4, "RdYlGn")


mypalette  = colorRampPalette(
  c("darkorange3", "white","darkblue")
)(100)
head(mypalette)
par(mai = rep(0.1, 4))
image(matrix(1:100, nrow = 100, ncol = 10), col = mypalette,
      xaxt = "n", yaxt = "n", useRaster = TRUE)


library("pheatmap")
topGenes = order(rowVars(Biobase::exprs(x)), decreasing = TRUE)[1:500]
rowCenter = function(x) { x - rowMeans(x) }
pheatmap( rowCenter(Biobase::exprs(x)[ topGenes, ] ),
          show_rownames = FALSE, show_colnames = FALSE,
          breaks = seq(-5, +5, length = 101),
          annotation_col =
            pData(x)[, c("sampleGroup", "Embryonic.day", "ScanDate") ],
          annotation_colors = list(
            sampleGroup = groupColor,
            genotype = c(`FGF4-KO` = "chocolate1", `WT` = "azure2"),
            Embryonic.day = setNames(brewer.pal(9, "Blues")[c(3, 6, 9)],
                                     c("E3.25", "E3.5", "E4.5")),
            ScanDate = setNames(brewer.pal(nlevels(x$ScanDate), "YlGn"),
                                levels(x$ScanDate))
          ),
          cutree_rows = 4
)


groupColor[1]


hexvals = sapply(1:3, function(i) substr(groupColor[1], i*2, i*2+1))
decvals = strtoi(paste0("0x", hexvals))


library("colorspace")
library("grid")

plothcl = function(h, c, l, what, x0 = 0.5, y0 = 0.5, default.units = "npc", ...) {
  switch(what,
         "c" = {
           stopifnot(length(l)==1)
           n = length(c)
         },
         "l" = {
           stopifnot(length(c)==1)
           n = length(l)
         },
         stop("Sapperlot"))
  
  cr = seq(0.1, 0.5, length = n+1)
  dr = 0.05 / n
  
  for (j in seq_len(n)) {
    r = c(cr[j]+dr, cr[j+1]-dr)
    for(i in 1:(length(h)-1)){
      phi = seq(h[i], h[i+1], by=1)/180*pi
      px = x0 + c(r[1]*cos(phi), r[2]*rev(cos(phi)))
      py = y0 + c(r[1]*sin(phi), r[2]*rev(sin(phi)))
      mycol = switch(what,
                     "c" = hcl(h=mean(h[i+(0:1)]), c=c[j], l=l),
                     "l" = hcl(h=mean(h[i+(0:1)]), c=c, l=l[j]))
      grid.polygon(px, py, gp=gpar(col=mycol, fill=mycol),
                   default.units=default.units,...)
    }
  }
}


plothcl( h = seq(0, 360, by=3), c = seq(5, 75, by=10), l = 75,   what="c")
grid.newpage()
plothcl( h = seq(0, 360, by=3), c = 55, l = seq(20, 100, by=10), what="l")


gg = ggplot(tibble(A = Biobase::exprs(x)[, 1], M = rnorm(length(A))),
            aes(y = M))
gg + geom_point(aes(x = A), size = 0.2)
gg + geom_point(aes(x = rank(A)), size = 0.2)




volume = function(rho, nu)
  pi^(nu/2) * rho^nu / gamma(nu/2+1)

ggplot(tibble(nu    = 1:15,
              Omega = volume(1, nu)), aes(x = nu, y = Omega)) +
  geom_line() +
  xlab(expression(nu)) + ylab(expression(Omega)) +
  geom_text(label =
              "Omega(rho,nu)==frac(pi^frac(nu,2)~rho^nu, Gamma(frac(nu,2)+1))",
            parse = TRUE, x = 6, y = 1.5)


ggplot(genes, aes( x = value, color = gene)) + stat_ecdf() +
  theme(text = element_text(family = "Times"))






library("ggbio")
data("hg19IdeogramCyto", package = "biovizBase")
plotIdeogram(hg19IdeogramCyto, subchr = "chr1")


library("GenomicRanges")
data("darned_hg19_subset500", package = "biovizBase")
autoplot(darned_hg19_subset500, layout = "karyogram",
         aes(color = exReg, fill = exReg))


data("ideoCyto", package = "biovizBase")
dn = darned_hg19_subset500
seqlengths(dn) = seqlengths(ideoCyto$hg19)[names(seqlengths(dn))]
dn = keepSeqlevels(dn, paste0("chr", c(1:22, "X")))
autoplot(dn, layout = "karyogram", aes(color = exReg, fill = exReg))


darned_hg19_subset500[1:2,]


stopifnot(is(darned_hg19_subset500, "GRanges"), identical(start(darned_hg19_subset500),end(darned_hg19_subset500)))


ggcars = ggplot(mtcars, aes(x = hp, y = mpg)) + geom_point()
ggcars
ggcars + theme_bw()
ggcars + theme_minimal()
