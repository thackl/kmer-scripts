### kmer.R ##
## A collection of functions for kmer statistics analysis and plotting
## Author: Thomas Hackl
## Version: 2.1.1

##-- make_plots.R ------------------------------------------------------------------##

run_task <- function(task, ...){
    switch(task,
           scr=scr(...),
           kfr=kfr(...),
           kmerPlot=kmerPlot(...),
           unknown_task(task)
           )
}


unknown_task <- function(task){
    write(paste("Unknown task:", task), stderr());
    quit(status = 1);
}




##-- task kmerPlot --##

kmerPlot <- function(..., coverage=NULL, out="kmerPlot.pdf", useggplot=FALSE){

    coverage <- as.integer(coverage)
    sets <- c(...);

    ## kmer filter
    data <- list();
    lims <- list();
    for (set in sets){
        write(paste("Reading kmers", set), stderr());
        ## read data and determine peaks
        if(grepl(".jf$", set, perl=TRUE)){
            data[[set]]$raw  <- read.table(pipe(paste('jellyfish histo', set, sep=" ")), header=F);
            data[[set]] <- c(data[[set]],  kmerPeaks(data[[set]]$raw, breaks=300, trim.tail=F));
        }else{
            data[[set]]$raw  <- read.table(set, header=F);
            data[[set]] <- c(data[[set]], kmerPeaks(data[[set]]$raw, breaks=300, trim.tail=F));
        }

        ## determine limits
        lims$covs <- c(lims$covs, max(data[[set]]$peaks$covs[1]))
        lims$freqs <- c(lims$freqs, max(data[[set]]$peaks$freqs[1]))
        lims$frels <- c(lims$frels, max(data[[set]]$peaks$frels[1]))

        write(paste(set, max(data[[set]]$peaks$covs[1]), sep="\t"), stdout());
    }

    if(useggplot){
        library("ggplot2")

        df <- data.frame(
            set=rep(sets[1], length(data[[1]]$raw$V1)),
            x=data[[1]]$raw$V1,
            y=data[[1]]$raw$V2
        )

       str(df);

        for (set in sets[2:length(sets)]){
            df <- rbind(df, data.frame(
                set=rep(set, length(data[[set]]$raw$V1)),
                x=data[[set]]$raw$V1,
                y=data[[set]]$raw$V2
            ))
        }

       str(df);

        gg <- ggplot(df, aes(x=x, y=y, group=set, colour=set));
        gg <- gg + xlim(0, max(lims$covs)*3)
        gg <- gg + ylim(1,max(lims$frels)*1.5)
        gg <- gg + geom_line();
        gg <- gg + labs(x="coverage", y="frequency"
                        #,title="kmer-coverage of data sets"
                        )

        ggsave(gg, file=out, width=10, height=6);

    }else{
        pdf(out, width=10, height=6);

        write("Plotting filtered kmers", stderr());
        plot(c(1,1),
             type="n",
             main="kmer-coverage of data sets",
             ##log="xy",
             xlim=c(0,max(lims$covs)*3),
             ylim=c(1,max(lims$frels)*1.5),
             xlab="coverage",
             ylab="frequency"
             );


        i=1;
        for (set in sets){
            i=i+1;
            lines(data[[set]]$raw, col=cl[i], lwd=3);
                                        #lines(data[[set]]$data, col=cl[i], lwd=2, lty=5);
                                        #abline(v=data[[set]]$peaks$covs[1], lty=5, lwd=2, col=cl[5]);
            rel=T;
            peakSizes.add(data[[set]]$peaks, rel=rel);
            y <- ifelse(rel, max(data[[set]]$peaks$frels), max(data[[set]]$peaks$freqs));
            total.size=sum(apply(data[[set]]$raw[-1:-15,], MARGIN=1, FUN=prod)) ## ignore first 10 kmers
            gen.cov <- data[[set]]$peaks$covs[1];
            text(x=gen.cov, pos=4, y=y*1.4, labels=paste("total:", round( total.size/10^6/gen.cov),"Mbp",sep=" "), vfont=c("sans serif", "bold"))
        }


        if(!is.null(coverage)){
            abline(v=coverage, lwd=2, col=cl[5])
        }


        legend(
            "topright",
            basename(sets),
            lwd=3,
            lty=1,
            seg.len=2,
            col=cl[(1:length(sets))+1]
        );

        msg <- dev.off()
    }
}








##-- task scr --##

scr <-function(reads.hash, seeds.hash, do.plot=TRUE){

    ## from seed reads
    write("Reading seed kmers", stderr());
    seeds <- kmerPeaks(read.table(pipe(paste('jellyfish histo', seeds.hash, sep=" "))), freq.min=0);

    ## from entire input
    write("Reading total kmers", stderr());
    raw <- kmerPeaks(read.table(pipe(paste('jellyfish histo', reads.hash, sep=" "))), freq.min=0);

    ## TODO: warn/die if no peaks

    ## get raw peak closest to seed peak
    delta <- abs(raw$peaks$covs - seeds$peaks$covs[1]);
    ## get closest peak(s) - minimum distance, two in case of 2 peaks in of equal dist
    p1.i <- which(delta == sort(unique(delta))[1]);
    if(length(p1.i) < 2){ ## get second-closest peak(s)
        p2.i <- which(delta == sort(unique(delta))[2])
    }else{
        p2.i <- p1.i[2];
        p1.i <- p1.i[1];
    }

    double.peak <- FALSE;
    ## check if we are dealing with double peak
    for (i in p2.i){
        ratio = raw$peaks$covs[p1.i]/raw$peaks$covs[i];
        if(ratio > 1.8 && ratio < 2.2){ ## double peak
            double.peak = TRUE;
            p2.i <- p1.i;
            p1.i <- i;
        }else if(ratio > 0.4 && ratio < 0.6){
            double.peak = TRUE;
            p2.i <- i;
        }
    }

    p1 <- raw$peaks[p1.i,];
    p2 <- NULL;
    if(double.peak){
        p2 <- raw$peaks[p2.i,];
        pt.size <- round(p1$sizes/p1$covs + p2$sizes/p2$covs)
    }else{
        pt.size <- round(p1$sizes/p1$covs)
    }

    ## TODO: warn/die if pt.size > 400kb

    write(paste("coverage", round(p1$cov), sep="\t"), stdout());
    write(paste("size", pt.size, sep="\t"), stdout());

    if(do.plot){

        pdf("scr-seeds.pdf", width=10, height=6);

        write("Plotting seed and total kmers", stderr());

        ## get some plotting limits

        xmax <- max(seeds$peaks$covs[1])^1.5
        xmin <- max(seeds$peaks$covs[1])^(1/2)
        ymax <- max(seeds$peaks$freqs[1])*2.5


        barplot(
            seeds$data$freqs,
            diff(seeds$data$covs),
            space=0,
            main="kmer-coverage of plastid seed reads",
            xlab="coverage",
            ylab="frequency",
            xlim=c(xmin,xmax),
            ylim=c(0,ymax),
            log="x",
            col=cl[2],
            axes=FALSE
        );

        ## scale raw kmer frequency of seed peak surrounding area to fit the plot area
        raw.i <- which.min(abs(raw$data$covs - seeds$peaks$covs[1]))
        raw.sf <- seeds$peaks$freqs[1] / log(max(raw$data$frels[(raw.i-5):(raw.i+5)]))
        raw$peaks.scaled <- raw$peaks
        raw$peaks.scaled$frels <- log(raw$peaks$frels)*raw.sf

        axis(2)
        axis(1, pos=-3);
        lines(raw$data$covs,log(raw$data$frels)*raw.sf, col=cl[4], lwd=3);
        abline(v=seeds$peaks$covs[1], lwd=3, col=cl[5]);
        abline(v=p1$covs[1], lwd=3, col=cl[1]);
        if(! is.null(p2)){
            abline(v=p2$covs[1], lwd=3, lty=2, col=cl[1]);
        }
        peakSizes.add(raw$peaks.scaled, rel=T);

        legend(
            "topright",
            c("plastid seed kmers", "plastid seed coverage", "total data kmers (30%)", "total plastid coverage"),
            lwd=3,
            lty=c(1,1,1,1),
            seg.len=2,
            col=cl[c(2,5,4,1)]
        );

        msg <- dev.off()

    }
}


##-- task kfr --##

kfr <- function(..., coverage=NULL){

    coverage <- as.integer(coverage)
    sets <- c(...);

    ## kmer filter
    data <- list();
    lims <- list();
    for (set in sets){
        write(paste("Reading kmers", set), stderr());
        data[[set]]  <- kmerPeaks(read.table(pipe(paste('jellyfish histo -l 10', set, sep=" ")), header=F));
        lims$covs <- c(lims$covs, max(data[[set]]$peaks$covs[1]))
        lims$freqs <- c(lims$freqs, max(data[[set]]$peaks$freqs[1]))
        write(paste(set, max(data[[set]]$peaks$covs[1]), sep="\t"), stdout());
    }

    pdf("kfr.pdf", width=10, height=6);

    write("Plotting filtered kmers", stderr());
    plot(c(1,1),
         type="n",
         main="kmer-coverage of subsetted and filtered data sets",
         ##log="xy",
         xlim=c(0,max(lims$covs)*3),
         ylim=c(1,max(lims$freqs)*1.5),
         xlab="coverage",
         ylab="frequency"
         );


    i=1;
    for (set in sets){
        i=i+1;
        lines(data[[set]]$data, col=cl[i], lwd=3);
        abline(v=data[[set]]$peaks$covs[1], lty=5, lwd=2, col=cl[i]);
    }

    if(!is.null(coverage)){
        abline(v=coverage, lwd=2, col=cl[5])
    }
    ##add_psizes(kfr2.ex);

    legend(
	"topright",
    	sets,
    	lwd=3,
    	lty=1,
    	seg.len=2,
    	col=cl[(1:length(sets))+1]
    );

    msg <- dev.off()

}


##-- task rrm --##

rrm <- function(){

    ## ref map coverage
    rrm<-read.table("rrm-cov.tsv", header=T);
    for (df in split(rrm, rrm$id)){
    	id=df[1,1];
	if(id == "genome") next;
    	dt=df[,2:3]
  	dt.ex = get_extrema(dt, peaks=c(50,100,150,200));
    	plot(dt,
    	     type="n",
    	     main=paste("per-base coverage of reference (", id, ")", sep=""),
    	     xlab="coverage",
    	     ylab="frequency",
	     xlim=c(1,dt.ex$cov[length(dt.ex$cov)]*3),
	     ylim=c(0,dt[,2][dt.ex$cov]*2)
    	    );
    	lines(dt, col=cl[2], lwd=3);
	add_psizes(dt.ex);
    };

}

##-- gctile --##

gcmx <- function(..., coverage.max=300, out="kmerPlot.pdf"){
    library(reshape2);
    library(ggplot2);
    library(scales);

    mt.file <- c(...);
    mt <- read.table(mt.file, header=F);
    mt <- mt[,2:(dim(mt)[2]-1)]; # remove first and last col - aggregate of missing values
    #z.max <- max(mt[,(dim(mt)[2]/5):dim(mt)[2]]);
    z.max <- max(mt[, -1:-5]); # ignore first 5 kmers for max

    mt.df <- expand.grid(y=1:dim(mt)[1], x=1:dim(mt)[2]);
    mt.df$z <- unlist(mt);

    #cls <- rev(rainbow(6))
    #cls <- cls[-1]
    #cls <- c("darkblue", "cyan", "green", "yellow", "red");
    cls <- rev(c("#6a0000", "#d40000","#fb8b00", "#fddf01", "#90ff36", "#90fcfc", "#020061"));
    gg <- ggplot(mt.df, environment=environment())
    gg <- gg + geom_tile(aes(x=x, y=y, fill=z))
          #  scale_x_continuous(limits=c(0,500))
    gg <- gg + scale_fill_gradientn("kmer frequency", colours=cls, limits=c(0,z.max), na.value=cls[length(cls)]);
    gg <- gg + labs(x="coverage", y="GC coer kmer")
    if(coverage.max) gg <- gg + xlim(c(0,coverage.max))

    ggsave(gg, file=out, width=10, height=6);
}

##-- gccov --#
gccov <- function(..., out="gccov.pdf", length.min=1000, coverage.max=500,
                  tax.occ.min=1, bin.num=100, tax.ignore=FALSE, theme="gg",
                  palette="gg", length.min.scatter=0, jitter=FALSE,
                  sample.scatter=0, width=10, height=6
                  ){

    library(reshape2);
    library(ggplot2);
    library(scales);
    library(grid);
    library(gridExtra);
    library(RColorBrewer);
    library(colorspace);

    if (tax.occ.min < 1) tax.occ.min <- 1;

    ## read data
    write("reading table", stderr());
    df.file <- c(...);
    df.fh <- OpenRead(df.file) # prevent R peek bug on <() constructs
    df <- read.table(df.fh, header=F, fill=T, sep="\t");
    close(df.fh)

    with.tax <- FALSE
    if ( !tax.ignore &&
            length(df[1,]) == 5 ){
        colnames(df) <- c("contig","length","GC","coverage","tax");
        with.tax <- TRUE
    }else{
        colnames(df) <- c("contig","length","GC","coverage");
    }

    ## prepare data
    write("filtering data", stderr());
    df <- subset(df, GC > 0 & GC < 1 & length >= length.min); # ignore poly AAAA,GGGG, ..
    df <- df[order(df$length),];
    get_length.bin <- function(x){ as.integer(log10(x)) }

    df$length.bin <- sapply(df$length, get_length.bin)
    df$length.bin <- factor(df$length.bin, levels=sort(unique(df$length.bin), decreasing=T)) # order in which hist is stacked

    ## get some specs
    x.min <- min(df$GC)
    x.max <- max(df$GC)
    y.max <- ifelse(coverage.max, coverage.max, max(df$coverage)) # quantile(df$coverage, c(.5))*3

    write("setting up plots", stderr());
    aes.points <- c();
    scale.colour <- c();
    scale.size <- c();
    scale.shape <- c();
    scale.fill <- c();

    ## aestetics
    labs <- c(">10",">100",">1k", ">10k", ">100k", ">1M", ">10M", ">100M");
    breaks <- 1:8;
    cls <- c()
    if(palette == "gg"){
        cls <- gg_color_hue(length(breaks)-1);
        cls <- c("#525252", cls[c(5,3,1,6,7,2,4)]); # mixed gg + grey base
    }else{
        cls <- brewer.pal(length(breaks), palette)
    }
    sizes  <- c(0.2,0.3,.5,.9,3,6,10,15);
    shapes <- rep(1, 8)   # only one shape

    names(labs) <- breaks
    names(cls) <- breaks
    names(shapes) <- breaks
    names(sizes) <- breaks


    range <- as.integer(log10(range(df$length)))
    breaks <- breaks[range[1]:range[2]];
    labs <- labs[range[1]:range[2]];

    aes.points <- aes(x=GC, y=coverage,
                      colour=length.bin,
                      size=length.bin,
                      shape=length.bin)

    scale.shape <- scale_shape_manual("Contigs (bp)", labels=labs, values=shapes, limits=breaks)
    scale.size <- scale_size_manual("Contigs (bp)", labels=labs, breaks=breaks, values=sizes, limits=breaks)
    scale.colour <- scale_colour_manual("Contigs (bp)", labels=labs, breaks=breaks, values=cls, limits=breaks)

                                        # gh
    cls.fill <- rev(cls)
    cls.guide <- c()
    scale.fill <- scale_fill_manual(values=cls.fill);
    tax.length.level <- c();

    if(with.tax){ # tax column
        ## sort tax
        tax.df <- (as.data.frame(table(df$tax)))
        tax.df <- tax.df[order(tax.df$Freq, decreasing=T),]
        tax.df <- tax.df[tax.df$Freq >= tax.occ.min,]
        tax.levels <- as.vector(tax.df$Var1);
        tax.n <- length(tax.levels)
        tax.labels <- as.vector(apply(tax.df, 1, function(x){paste(x[1]," (",as.numeric(x[2]),")", sep="")}))
        df$tax <- factor(df$tax, levels=tax.levels)
        df <- df[complete.cases(df),] # remove low freq taxs with NA tax factor

        ## compute colours
        if(palette == "gg"){
            cls <- rev(gg_color_hue(tax.n));
        }else{
            cls <- brewer.pal(ifelse(tax.n>3, tax.n, 3), palette)[1:tax.n]
        }

        tax.cls <- c()
        tax.cls.fill <- c()

        for (i in 1:tax.n){
            cls.seq <- hex_col_seq(cls[i], n=8)
            cls.seq <- cls.seq[range[1]:range[2]]

                                        #x <- tax.levels[i]
                                        #x.id <-unique(as.integer(log10(df$length[df$tax==x])))
                                        #cls.seq <- cls.seq[x.id]

            tax.cls <- c(tax.cls, rev(cls.seq))
            tax.cls.fill <- c(tax.cls.fill, rev(cls.seq))
        }

        tax.length.level <- levels(factor(df$tax):df$length.bin)
        names(tax.cls) <- tax.length.level

                                        # guide colour
        cls.guide <- hex_col_seq(cls[1], n=8)
        cls.guide <- cls.guide[range[1]:range[2]]

        tax.breaks <- factor(sapply(tax.levels, function(x){ paste(x, as.integer(log10(max(df$length[df$tax==x]))), sep=":") }), level=tax.length.level)

        aes.points <- aes(x=GC, y=coverage,
                          colour=factor(tax):factor(length.bin),
                          size=length.bin,
                          shape=length.bin
                          )
                                        #scale.shape <- scale_shape_manual("Taxonomy", labels=tax.labels, values=shapes, limits=breaks)
                                        #scale.shape <- scale_shape("Taxonomy", solid=FALSE, labels=tax.labels)
        scale.colour <- scale_colour_manual("Taxonomy", breaks=tax.breaks, labels=tax.labels, values=tax.cls)
        scale.fill <- scale_fill_manual(values=tax.cls, breaks=tax.breaks);

    }


    ##- gg ----------------------------------------------------------------------##
    ## themes
    gg.theme <- theme(
        text=element_text(size=10),
        legend.position="none",
        plot.margin=unit(c(.5,.5,.5,0), "cm"))

    ##- themes ----------------------------------------------------------------------##
    if (theme == "bw"){
        gg.theme <- theme_bw() + gg.theme;
    } else if (theme == "classic"){
        gg.theme <- theme_classic() + gg.theme;
    }

    ## plots
    df.gg <- subset(df, length > length.min.scatter)
    if(sample.scatter){
        df.gg <- df.gg[sample(nrow(df.gg), nrow(df.gg)*sample.scatter), ]
        df.gg <- df.gg[order(df.gg$length),];
    }

    gg <- ggplot(df.gg) +
        geom_point(aes.points, position=ifelse(jitter, "jitter", "identity")) +
            scale.shape +
                scale.colour +
                    scale.size +
                        gg.theme +
                            xlim(x.min, x.max) +
                                ylim(0,y.max)


    guides.legend <- guides(
        colour = guide_legend(order=1, reverse=TRUE),
        size = guide_legend(order=1, reverse=TRUE),
        shape = guide_legend(order=1, reverse=TRUE)
    )

    if(with.tax){
        guides.legend <- guides(
            colour = guide_legend(order=1, override.aes=list(colour=cls)),
            size = guide_legend(order=2, reverse=TRUE),
            shape = guide_legend(order=2, reverse=TRUE, override.aes=list(colour=rev(cls.guide)))
        )
    }

    gg.legend <- get_legend(
        gg + theme(
            legend.justification=c(1,1),
            legend.position=c(0,1),
            legend.box.just = "right",
            legend.text=element_text(size=8),
            text = element_text(size=10)
        ) + guides.legend
    )


    ##- gh ----------------------------------------------------------------------##
    gh.theme <- theme(
        text=element_text(size=10),
        legend.position="none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=unit(c(.5,.5,.5,0), "cm"))

    if (theme == "bw"){
        gh.theme <- theme_bw() + gh.theme;
    }else if (theme == "classic"){
        gh.theme <- theme_classic() + gh.theme;
    }


    gh <- ggplot(df) +
        ylab("sum of length") +
            coord_flip() +
                labs(x=NULL) +
                    gh.theme +
                        scale_y_continuous(labels = scientific_format(digits=0)) +
                             xlim(0,y.max) + scale.fill

    if(with.tax){ # tax column
        gh <- gh + geom_bar(aes(x=coverage, weight=length, fill=factor(tax):factor(length.bin)), binwidth=y.max/bin.num)
    }else{
        gh <- gh + geom_bar(aes(x=coverage, weight=length, fill=factor(length.bin)), binwidth=y.max/bin.num)
    }


    ##- plotting ----------------------------------------------------------------------##
    write("plotting", stderr());

    if(grepl(".pdf$", out)){
        pdf(out, width=width, height=height);
    }else if(grepl(".png$", out)){
        png(out, width=width*100, height=height*100);
    }else{
        stop("only .pdf and .png output supported");
    }
    grid.arrange(gg, gh, gg.legend, nrow = 1, widths = c(0.65, .35, .0))
    dev.off();
}


hist_count <- function (x){
  dt <- as.data.frame(table(x))
  dt$x <- as.numeric(levels(dt$x))[dt$x]
  return(dt)
}


seqcov <- function(..., out="kmer-plot-seqcov.pdf"){
    library(reshape2);
    library(ggplot2);
    library(grid);
    library(gridExtra);

    tsv.file <- c(...)[1];
    tsv.id <- c(...)[2];
    cmd <- paste("grep -w '^",tsv.id,"' ",tsv.file," | reshape", sep="");
    print(cmd);
    d1 <- read.table(pipe(cmd), header=F);
    head(d1)

    d1$pos <- 1:length(d1[,1])
    colnames(d1) <- c("id", "rm", "rr", "am", "ar", "rc", "ac", "pos")

    gg1 <- ggplot(d1, aes(x=pos)) +
        geom_line(aes(y=rc, colour="raw counts")) +
            geom_line(aes(y=ac, colour="adjusted counts")) +
                geom_abline(intercept=41, slope=0) +
                    scale_y_continuous(limits=c(0, 600))


                                        # extract legend
    gg1.legend <- get_legend(gg1)
    gg1 <- gg1+theme(legend.position="none",
                     plot.margin=unit(c(.5,.5,.5,0), "cm"))

    h1 <- cbind(set="rc", hist_count(d1$rc))
    h1 <- rbind(h1, cbind(set="ac", hist_count(d1$ac)))

    gh1 <- ggplot(h1, aes(x=x, y=Freq, colour=set)) +
        geom_line() +
            coord_flip() +
                labs(x=NULL) +
                    theme(axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          axis.title.y=element_blank(),
                          plot.margin=unit(c(.5,.5,.5,0), "cm")) +
                    scale_x_continuous(limits=c(0, 600))

    pdf(out, width=10, height=6);
    grid.arrange(gg1, gh1, gg1.legend, nrow = 2, widths = c(0.65, 0.35), heights=c(.8,.2))
    dev.off();
}


asmcov <- function(..., out="kmerPlot.pdf", length.min=1000, coverage.max=500,
                   bin.num=100, anscombe=FALSE, theme="gg", palette="gg",
                   width=10, height=6, kmer.histo=""
                   ){

    library(reshape2);
    library(ggplot2);
    library(scales);
    library(grid);
    library(gridExtra);
    library(RColorBrewer);
    library(colorspace);

    ## read data
    files <- c(...);
    df <- data.frame(contig=character(0), length=numeric(0), GC=numeric(0), coverage=numeric(0), assembly=character(0));

    for (df.file in files){
        write(paste("reading table: ", df.file), stderr());

        df.fh <- OpenRead(df.file) # prevent R peek bug on <() constructs
        df.tmp <- read.table(df.fh, header=F, fill=T, sep="\t");
        close(df.fh)

        df.tmp[5] <- df.file
        colnames(df.tmp) <- c("contig","length","GC","coverage", "assembly");
        df <- rbind(df, df.tmp)
    }

    dk <- c();
    if(length(kmer.histo)){
        dk <- read.table(kmer.histo, header=F)
        colnames(dk) <- c("coverage", "count")
                                        #dk$frequency
    }

    ## prepare data
    write("filtering data", stderr());
    df <- subset(df, GC > 0 & GC < 1 & length >= length.min); # ignore poly AAAA,GGGG, ..

    get_length.bin <- function(x){ as.integer(log10(x)) }
    df$length.bin <- sapply(df$length, get_length.bin)
    df$length.bin <- factor(df$length.bin, levels=sort(unique(df$length.bin), decreasing=T)) # order in which hist is stacked

    write("setting up plots", stderr());

    ## aestetics
    y.max <- max(df$length)
    x.max <- ifelse(coverage.max, coverage.max, max(df$coverage))
    x.breaks <- c();

    breaks <- 1:8
    labs <- c(">10",">100",">1k", ">10k", ">100k", ">1M", ">10M", ">100M");
    names(labs) <- breaks

    cls <- c()
    if(palette == "gg"){
        cls <- gg_color_hue(length(breaks)-1);
        cls <- c("#525252", cls[c(5,3,1,6,7,2,4)]); # mixed gg + grey base
    }else{
        cls <- brewer.pal(length(breaks), palette)
    }
    names(cls) <- breaks

    range <- as.integer(log10(range(df$length)))
    breaks <- breaks[range[1]:range[2]];
    labs <- labs[range[1]:range[2]];

    scale.fill <- scale_fill_manual("Contigs (bp)", labels=labs, breaks=breaks, values=cls, limits=breaks)

    ##- themes ----------------------------------------------------------------------##
    gh.guides <- guides(fill = guide_legend(reverse=FALSE))

    gh.theme <- theme(text=element_text(size=10), panel.grid.minor = element_blank())

    if (theme == "bw"){
        gh.theme <- theme_bw() + gh.theme;
    }else if (theme == "classic"){
        gh.theme <- theme_classic() + gh.theme;
    }

    ##- plot ------------------------------------------------------------------------##
    gh.aes <- c();
    scale.x <- c();
    bin.width <- c();
    if(! anscombe){
        bin.width <- x.max/bin.num
        gh.aes <- aes(x=coverage, weight=length, fill=length.bin)
                                        #gh.aes <- aes(x=coverage, fill=length.bin)
        scale.x <- scale_x_continuous("coverage", limits=c(NA, x.max))
    }else{
        x.breaks <- anscombe_breaks(x.max);
        bin.width <- 1;
        gh.aes <- aes(x=anscombe(coverage), weight=length, fill=length.bin)
                                        #gh.aes <- aes(x=anscombe(coverage), fill=length.bin)
        scale.x <- scale_x_continuous(breaks=x.breaks, labels=anscombe_inv, "coverage", limits=c(NA, anscombe(x.max)))
    }

    gh <- ggplot(df) +
        geom_bar(gh.aes, stat="bin", binwidth=bin.width) +
            scale_y_continuous("sum of length [bp]", labels=scientific_format(digits=0)) +
                gh.theme +
                    gh.guides +
                        scale.x +
                            scale.fill +
                                facet_wrap(~assembly, ncol=1)

    if(length(kmer.histo)){
        gh <- gh + geom_line(
            data=dk[dk$coverage >10,], aes(x=anscombe(coverage), y=count*50, linetype="count * 50")) +
                scale_linetype("k-mers")

    }

                                        #gg.c <- ggplot(df, aes(x=anscombe(x))) +
                                        #  geom_bar(aes(weight=len, fill=bin), stat="bin", binwidth=1) +
                                        #  scale_x_continuous(breaks=breaks.anscombe, labels=anscombe_inv)


    ##- plotting ----------------------------------------------------------------------##
    write("plotting", stderr());

    if(grepl(".pdf$", out)){
        pdf(out, width=width, height=height);
    }else if(grepl(".png$", out)){
        png(out, width=width*100, height=height*100);
    }else{
        stop("only .pdf and .png output supported");
    }

    print(gh)

    dev.off();
}

anscombe <- function(x){2*sqrt(x + 3/8)}
anscombe_int <- function(x){round(2*sqrt(x + 3/8), digits=0)}
anscombe_inv <- function(x){(x/2)^2 - 3/8}

anscombe_breaks <- function(max, list=FALSE){
    ex <- as.integer(log10(max))+2
    b.n <- as.vector(sapply(1:ex, function(x){c(.1,.2,.5) * 10^x} ))
    b.n <- b.n[1:which(b.n > max)[1]]
    b.a <- sapply(b.n, anscombe)
    d.a <- max(b.a)/15

    r <- list();
    r$b.a <- b.a[1]
    r$b.n <- b.n[1]
    for(i in 2:length(b.a)){
        if(r$b.a[length(r$b.a)] + d.a < b.a[i]){
            r$b.a <- c(r$b.a, b.a[i])
            r$b.n <- c(r$b.n, b.n[i])
        }
    }

    if(list){
        return(r)
    }else{
        return(r$b.a)
    }
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##-- shared --##

kmerPeaks <- function(d, cov.dev=0.25, cov.min=15, k=3, freq.min=15, smooth="histGeom", breaks=100, trim.tail=TRUE){

    colnames(d) <- c("covs", "freqs");

    peaks.covs <- c()
    peaks.freqs <- c()
    peaks.frels <- c()
    peaks.sizes <- c()
    peaks.err.covs <- c()
    peaks.err.freqs <- c()
    peaks.err.frels <- c()
    peaks.err.freqs <- c()
    peaks.minor.covs <- c()
    peaks.minor.freqs <- c()
    peaks.minor.frels <- c()
    peaks.minor.sizes <- c()

    total.size <- NA

    ##pits.covs <- c()
    ##pits.freqs <- c()
    ##peak.max.cov <- NA
    ##peak.max.freq <- NA

    ## total.size
    total.size=sum(apply(d, MARGIN=1, FUN=prod))

    ## remove last cov - might be a max bin that contains any higher cov freqs as well
    d <- d[-dim(d)[1],]

    ## trim tail
    if(trim.tail){
        ## trim tail
        tail.start <- tail(which(runmed(d$freqs, k=7, endrule="median") > freq.min), n=1);
        if(length(tail.start)){
            d <- d[1:tail.start,]
        }
    }

    ## smooth da biatch
    d <- switch(smooth,
                medAdj = smooth.medAdj(d, k=k),
                med = smooth.med(d, k=k),
                hist = smooth.hist(d, breaks=breaks),
                histGeom = smooth.hist(d, breaks=breaks, geom=T),
                none = d,
                smooth.unknown(smooth)
                );

    ##print(str(d));
    ## get peaks
    maxs.i <- localMaxima(d$freqs)

    for (i in maxs.i){
        freq <- d$freqs[i]
        cov <- d$covs[i]
        d.env <- d[which(d$covs >= cov-cov*cov.dev & d$covs < cov+cov*cov.dev),]

        freq.env=d$freqs[which(d$covs >= cov & d$covs < cov+cov*cov.dev)]

        ## eval peaks in da hood
        if( length(freq.env)==1 || freq==max(freq.env) ){
            size <- peakSize(d, cov)
            if(length(peaks.covs) && sum(peaks.covs > cov-cov*cov.dev)){ ## the fight is on
                if(peaks.freqs[length(peaks.freqs)] <= freq*1.1){ ## demote the looser
                    peaks.minor.covs <- c(peaks.minor.covs, peaks.covs[length(peaks.covs)])
                    peaks.covs[length(peaks.covs)] <- cov

                    peaks.minor.freqs <- c(peaks.minor.freqs, peaks.freqs[length(peaks.freqs)])
                    peaks.freqs[length(peaks.freqs)] <- freq

                    peaks.minor.sizes <- c(peaks.minor.sizes, peaks.sizes[length(peaks.sizes)])
                    peaks.sizes[length(peaks.sizes)] <- size
                }
                ## wannabe
                peaks.minor.covs <- c(peaks.minor.covs, cov)
                peaks.minor.freqs <- c(peaks.minor.freqs, freq)
                peaks.minor.sizes <- c(peaks.minor.sizes, size)
            }else{
                peaks.covs <- c(peaks.covs, cov);
                peaks.freqs <- c(peaks.freqs, freq);
                peaks.sizes <- c(peaks.sizes, size);
            }
        }
    }

    ## trim min cov
    d <- d[which(d$covs >= cov.min),]

    ## handle error peak (cov 1 or 2)
    if(length(peaks.covs)){
        peaks.err.i <- peaks.covs < cov.min

        peaks.err.covs <- peaks.covs[peaks.err.i]
        peaks.covs <- peaks.covs[!peaks.err.i]

        peaks.err.freqs <- peaks.freqs[peaks.err.i]
        peaks.freqs <- peaks.freqs[!peaks.err.i]

        peaks.err.sizes <- peaks.sizes[peaks.err.i]
        peaks.sizes <- peaks.sizes[!peaks.err.i]
    }

    ## king of da hill
    if(! is.null(peaks.freqs)){
        peaks.frels <- d$frels[d$covs %in% peaks.covs]

        peaks.o <- order(peaks.freqs, decreasing=T)
        peaks.covs <- peaks.covs[peaks.o]
        peaks.freqs <- peaks.freqs[peaks.o]
        peaks.sizes <- peaks.sizes[peaks.o]
        peaks.frels <- peaks.frels[peaks.o]
    }

    if(! is.null(peaks.minor.freqs)){
        peaks.minor.frels <- d$frels[d$covs %in% peaks.minor.covs]

        peaks.minor.o <- order(peaks.minor.freqs, decreasing=T)
        peaks.minor.covs <- peaks.minor.covs[peaks.minor.o]
        peaks.minor.freqs <- peaks.minor.freqs[peaks.minor.o]
        peaks.minor.sizes <- peaks.minor.sizes[peaks.minor.o]
        peaks.minor.frels <- peaks.minor.frels[peaks.minor.o]

    }

    if(! is.null(peaks.err.freqs)){
        peaks.err.frels <- d$frels[d$covs %in% peaks.err.covs]

        peaks.err.o <- order(peaks.err.freqs, decreasing=T)
        peaks.err.covs <- peaks.err.covs[peaks.err.o]
        peaks.err.freqs <- peaks.err.freqs[peaks.err.o]
        peaks.err.sizes <- peaks.err.sizes[peaks.err.o]
        peaks.err.frels <- peaks.err.frels[peaks.err.o]
    }

    return(
        list(
            data = d,
            peaks = data.frame(
                covs=peaks.covs,
                freqs=peaks.freqs,
                sizes=peaks.sizes,
                frels=peaks.frels
            ),
            peaks.minor = data.frame(
                covs=peaks.minor.covs,
                freqs=peaks.minor.freqs,
                sizes=peaks.minor.sizes,
                frels=peaks.minor.frels
            ),
            peaks.error = data.frame(
                covs=peaks.err.covs,
                freqs=peaks.err.freqs,
                sizes=peaks.err.sizes,
                frels=peaks.err.frels
            ),
            total.size=total.size
        )
    );
}



## stolen from http://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
localMaxima <- function(x) {
  ## Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-Inf, x)) > 0L
  ## print(rle(y)$lengths)
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y
  }
  y
}


geomSeries <- function(base, max) {
    base^(0:floor(log(max, base)))
}

plotKmerPeaks <- function(...){

}

peakSize <- function(d, peak.cov){
    peak.i <- which(d$covs > peak.cov*0.5 & d$covs <= peak.cov*1.5)
    size <- sum(apply(d[peak.i,c("covs", "freqs")], MARGIN=1, FUN=prod));
    return(size);
}


smooth.unknown <- function(d, smooth){
    write(paste("unknown smooth algorithm", smooth), stderr());
    quit(status=1);
    return(d);
}

smooth.med <- function(d, k=3){
    d$freqs <- runmed(d$freqs, k=k, endrule="median");
    return(d)
}

smooth.medAdj <- function(d, k=3){
    ## smooth frequencies, gradually increase k size
    d$freqs[1:50] <- runmed(d$freqs[1:50], k=3, endrule="median");
    d$freqs[45:100] <- runmed(d$freqs[45:100], k=5, endrule="median");
    d$freqs[90:160] <- runmed(d$freqs[90:160], k=9, endrule="median");
    d$freqs[150:length(d$freqs)] <- runmed(d$freqs[150:length(d$freqs)], k=13, endrule="median");
    return(d)
}

smooth.hist <- function(d, breaks=200, geom=FALSE){

    cov.max <- tail(d$covs, n=1)

    ## geometric hist bin size
    if(geom){
        #### ## fixed number of bases
        base <- cov.max^(1/breaks);
        breaks <- cumsum(round(c(1, diff(base^(1:breaks)))));

        ## fixed base
        ##breaks <- geomSeries(2^(1/6), cov.max)

        ## handle too small breaks
        ## remove too small bins at the start
        breaks <- breaks[which(c(diff(breaks)>3, TRUE))];
        breaks <- c(0:floor(breaks[1]/3 -1)*3+1, breaks); ## add 3-bin

        breaks[length(breaks)] <- cov.max + 0.0001 ## make sure, max is contained
    }

    cov.hist <- hist(d$covs, breaks=breaks, plot=F)

    ## hist freqs for cov.hist :)
    cumsums <- cumsum(cov.hist$counts)
    freqs <- c();

    for(i in 1:length(cumsums)){
        f <- ifelse(i>1, cumsums[i-1]+1, 0)
        t <- cumsums[i]
        freqs <- c(freqs, sum(as.numeric(d$freqs[f:t])))
    }

    ## smooth freqs
    ##freqs <- runmed(freqs, k=5, endrule="median")

    ## scale
    if(cov.hist$equidist){
        sf <- cov.max/length(cov.hist$breaks)
        frels <- freqs/sf;
    }else{ ## scale geom
        break.intervals <- diff(cov.hist$breaks)
        frels <- freqs/break.intervals;
    }

    return(data.frame(covs=cov.hist$mids, freqs=freqs, frels=frels))
}

####---- DEPRECATED ----####

get_extrema <- function(data,peaks){

    cmins=c();
    fmins=c();
    for ( i in 1:length(peaks)){
        if(i>1){
            p1=peaks[i-1]
        }else{
            p1=0
        }
        p2=peaks[i];
        p2p=data[,1] >= p1 & data[,1] < p2;
        cmins=c(cmins,data[,1][p2p][which.min(data[,2][p2p])])
        fmins=c(fmins,data[,2][p2p][which.min(data[,2][p2p])])
    }

                                        ## peak max
    cmaxs=c();
    fmaxs=c();
    for ( i in 1:length(cmins) ){
        p1=cmins[i];
        if(i==length(cmins)){
            p2=data[,1][length(data[,1])]
        }else{
            p2=cmins[i+1]
        }

        p2p=data[,1] >= p1 & data[,1] < p2;
        cmaxs=c(cmaxs,data[,1][p2p][which.max(data[,2][p2p])])
        fmaxs=c(fmaxs,data[,2][p2p][which.max(data[,2][p2p])])
    }

    psizes=c();
    for ( i in 1:length(cmins) ){
        p1= cmins[i];
        if(i==length( cmins)){
            p2=data[,1][length(data[,1])]
        }else{
            p2= cmins[i+1]
        }

        p2p=data[,1] >= p1 & data[,1] < p2;
        psizes=c(psizes,sum(apply(data[p2p,], MARGIN=1, FUN=prod)))
    }

    cov=cmaxs[which.max(fmaxs)]
    freq=max(fmaxs)

    tbp=sum(apply(data, MARGIN=1, FUN=prod))

    return(list(cov.maxima=cmaxs, freq.maxima=fmaxs, cov.minima=cmins, freq.minima=fmins, peaksizes=psizes, cov=cov, freq=freq, total=tbp));
}

add_psizes <- function(data=extrema){
    cov=data$cov;
    y=data$freq*1.2
    for(i in 1:length(data$cov.maxima)){
        text(x= data$cov.maxima[i], y=y, labels=paste(c(round( data$peaksizes[i] /10^3/cov)),"kbp",sep=" "), vfont=c("sans serif", "bold"))
    }
    text(x=data$cov/10, pos=4, y=y*1.2, labels=paste("total:", round( data$total/10^3/data$cov),"kbp",sep=" "), vfont=c("sans serif", "bold"))
}

peakSizes.add <- function(p, rel=FALSE){

    for(i in 1:length(p$covs)){
        y <- ifelse(rel, p$frels[i], p$freqs[i]);
        points(p$covs[i], y, cex=1.3, pch=17);
        y <- y * 1.2;
        text(x= p$covs[i], y=y, labels=paste(c(round( p$sizes[i] /10^6/p$covs[i])),"Mbp",sep=" "), vfont=c("sans serif", "bold"))
    }
    ## only gets peaks, doesn't know total !!
    ##   y <- ifelse(rel, max(p$frels), max(p$freqs));
    ##    print(p);
    ##    text(x=p$covs[1], pos=4, y=y*1.2, labels=paste("total:", round( data$total/10^6/data$cov),"Mbp",sep=" "), vfont=c("sans serif", "bold"))
}


gg_color_hue <- function(n, l=65, c=100) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=l, c=c)[1:n]
}

hex_col_shade <- function(colour, n=5){
  colour.rgb <- col2rgb(colour)/255
  colour.shades <- sapply(seq(0.4,0.9,length=n), function(x){
    do.call(rgb, as.list(colour.rgb*x))
  })
  return(colour.shades)
}

hex_col_seq <- function(colour, n=5, c.=c(60,85), l=c(90,50)){
  c.hsv <- rgb2hsv(col2rgb(colour))
  sequential_hcl(n, h=c.hsv[1]*360, c.=c., l=l)
}

OpenRead <- function(arg) {
    ## http://stackoverflow.com/questions/15784373/process-substitution
    if (arg %in% c("-", "/dev/stdin")) {
        file("stdin", open = "r")
    } else if (grepl("^(/dev/fd/|/proc/self/fd/)", arg)) {
        fifo(arg, open = "r")
    } else {
        file(arg, open = "r")
    }
}

##-- main ---------------------------------------------------------------------##

## globals
cl=rainbow(5);

#params <- commandArgs(trailingOnly=T);

#if(length(params) > 0) do.call(run_task, as.list(params));

#if(! is.null(warnings())) warnings();
