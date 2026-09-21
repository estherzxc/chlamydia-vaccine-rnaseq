#Cell Ranger Results QC

setwd("/work/users/x/c/xczhou/Chlamydia/CellRanger_QC_New")

df <- read.csv("mapping_metrics.csv", stringsAsFactors = FALSE)

df$mapped_conf <- as.numeric(sub("%","", df$Mouse.Reads.Mapped.Confidently.to.Genome))
df$mapped_genome <- as.numeric(sub("%","", df$Mouse.Reads.Mapped.to.Genome))
df$mouse <- as.numeric(sub("^M([0-9]{2})_.*", "\\1", df$sample_folder))
df$day <- sub("^M[0-9]{2}_D([0-9]+)_.*", "\\1", df$sample_folder)
df$day[df$day == "0"] <- "-2"
df$day_label <- paste0(df$day, " dpi")

df$group <- ifelse(df$mouse <71, "Vaccinated", "Naive")

#QC Graphs, percent reads mapped to mouse genome
df$xgroup <- paste(df$group, df$day_label)
df$xgroup <- factor(df$xgroup,
                    levels = c(
                      "Naive -2 dpi",
                      "Naive 4 dpi",
                      "Naive 6 dpi",
                      "Naive 11 dpi",
                      "Vaccinated -2 dpi",
                      "Vaccinated 4 dpi",
                      "Vaccinated 6 dpi",
                      "Vaccinated 11 dpi"
                    ))
df$x <- as.numeric(df$xgroup)

#Make a shape variable
df$pch <- NA
df$pch[df$day_label == "-2 dpi"] <- 16   # filled circle
df$pch[df$day_label == "4 dpi"]  <- 15   # filled square
df$pch[df$day_label == "6 dpi"]  <- 17   # filled triangle
df$pch[df$day_label == "11 dpi"] <- 18   # filled diamond

par(mar = c(7, 4, 1, 1))

plot(df$x, df$mapped_genome,
     type = "n",
     xaxt = "n",
     xlab = "",
     ylab = "% Reads Mapped to Mouse Genome",
     main="Percent Reads Mapped to the Mouse Genome Across Samples",
     ylim = c(0, 100),
     xlim = c(0.5, 8.5),
     bty = "l")

points(df$x, df$mapped_genome,
       pch = df$pch,
       col = ifelse(df$group == "Vaccinated", "orange", "purple"),
       cex = 1.5)

axis(1, at = 1:8,
     labels = c("-2 dpi", "4 dpi", "6 dpi", "11 dpi",
                "-2 dpi", "4 dpi", "6 dpi", "11 dpi"),
     tick = FALSE)

# group labels underneath
mtext("Naive", side = 1, at = 2.5, line = 2.2, cex = 1.3)
mtext("Vaccinated", side = 1, at = 6.5, line = 2.2, cex = 1.3)

# underline group labels
segments(x0 = 0.7, y0 = -11, x1 = 4.3, y1 = -11, xpd = TRUE)
segments(x0 = 4.7, y0 = -11, x1 = 8.3, y1 = -11, xpd = TRUE)


#QC Graphs, percent reads confidently mapped to mouse genome
plot(df$x, df$mapped_conf,
     type="n",
     xaxt="n",
     xlab="",
     ylab="% Reads Mapped Confidently to Mouse Genome",
     main="Percent Reads Mapped Confidently to the Mouse Genome Across Samples",
     ylim=c(0,100),
     xlim=c(0.5,8.5),
     bty="l")
points(df$x, df$mapped_conf,
       pch=df$pch,
       col=ifelse(df$group=="Vaccinated","orange","purple"),
       cex=1.5)
axis(1, at=1:8,
     labels=c("-2 dpi","4 dpi","6 dpi","11 dpi",
              "-2 dpi","4 dpi","6 dpi","11 dpi"),
     tick=FALSE)

mtext("Naive", side=1, at=2.5, line=2.2, cex=1.3)
mtext("Vaccinated", side=1, at=6.5, line=2.2,cex=1.3)

segments(x0 = 0.7, y0 = -11, x1 = 4.3, y1 = -11, xpd = TRUE)
segments(x0 = 4.7, y0 = -11, x1 = 8.3, y1 = -11, xpd = TRUE)


#QC Graphs, Number of Reads
df$Number.of.Reads <- as.numeric(df$Number.of.Reads)
summary(df$Number.of.Reads)

df$reads_millions <- df$Number.of.Reads / 1e6

plot(df$x, df$reads_millions,
     type="n",
     xaxt="n",
     xlab="",
     ylab="Number of Reads (Millions)",
     main="Sequencing Depth Across Samples",
     ylim=c(0, max(df$reads_millions)*1.1),
     xlim=c(0.5,8.5),
     bty="l")
points(df$x, df$reads_millions,
       pch=df$pch,
       col=ifelse(df$group=="Vaccinated","orange","purple"),
       cex=1.5)
axis(1, at=1:8,
     labels=c("-2 dpi","4 dpi","6 dpi","11 dpi",
              "-2 dpi","4 dpi","6 dpi","11 dpi"),
     tick=FALSE)

mtext("Naive", side=1, at=2.5, line=2.2, cex=1.3)
mtext("Vaccinated", side=1, at=6.5, line=2.2, cex=1.3)

segments(0.7,-3.5,4.3,-3.5, xpd=TRUE)
segments(4.7,-3.5,8.3,-3.5, xpd=TRUE)

