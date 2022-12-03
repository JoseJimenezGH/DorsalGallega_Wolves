
# Modified from scrbook and oSCR packages
spiderplotJJ<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = "grey", lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)                            # Tamaño  # fondo  # grosor
        points(mu.x[[s]], mu.y[[s]], pch = 21, cex = 1.75, bg = clr, lwd=2.25)
    }
    par(op)
}


make.grid<-function (ll = NA, minx = NA, maxx = NA, miny = NA, maxy = NA,
nx = 40, ny = NULL, buffer = 0)
  {
    if (is.null(ny))
    ny <- nx
    if (!is.na(ll)) {
    minx <- min(ll[, 1])
    maxx <- max(ll[, 1])
    miny <- min(ll[, 2])
    maxy <- max(ll[, 2])
    bx <- (maxx - minx) * buffer
    by <- (maxy - miny) * buffer
    minx <- minx - bx
    maxx <- maxx + bx
    miny <- miny - by
    maxy <- maxy + by
  }
  x <- sort(rep(seq(minx, maxx, , nx), ny))
  y <- rep(seq(maxy, miny, , ny), nx)
  cbind(x, y)
}


rot<-function (m) 
{
    nr <- nrow(m)
    nc <- ncol(m)
    v <- matrix(NA, nrow = nc, ncol = nr)
    for (i in 1:nr) {
        v[, nr - (i - 1)] <- m[i, ]
    }
    v
}

make.scrFrame<-function (caphist, traps, indCovs = NULL, trapCovs = NULL, sigCovs = NULL,
    trapOperation = NULL, telemetry = NULL, rsfDF = NULL, type = "scr")
{
    if (any(is.null(caphist), is.null(traps)))
        stop("caphist and trap must be provided")
    if (!is.list(caphist))
        stop("caphist must be a list")
    n.sessions <- length(caphist)
    caphist.dimensions <- sapply(caphist, dim)
    if (nrow(caphist.dimensions) == 2)
        caphist.dimensions <- rbind(caphist.dimensions, 1)
    for (i in 1:n.sessions) {
        caphist[[i]] <- array(caphist[[i]], dim = caphist.dimensions[,
            i])
        all.zero <- apply(apply(caphist[[i]], c(1, 3), sum),
            1, sum)
        if (any(all.zero == 0)) {
            cat("At least one individual has an all-zero encounter history",
                fill = TRUE)
            cat("Make sure this is ok...", fill = TRUE)
        }
    }
    if (!is.null(indCovs)) {
        if (!is.list(indCovs))
            stop("indCovs must be a list")
        if (any(!sapply(indCovs, is.data.frame)))
            stop("indCovs must be a list of dataframes")
        if (length(indCovs) != length(caphist))
            stop("number of sessions in indCovs does not match caphist")
        check.dim <- sapply(indCovs, nrow)
        if (any(check.dim != caphist.dimensions[1, ]))
            stop("number of individuals in indCovs does not match caphist")
        if (!("rmv" %in% indCovs[[1]])) {
            for (i in 1:length(indCovs)) {
                indCovs[[i]]$removed <- dim(caphist[[i]])[3]
            }
        }
    }
    else {
        indCovs <- list()
        for (i in 1:length(caphist)) {
            indCovs[[i]] <- data.frame(removed = rep(dim(caphist[[i]])[3],
                dim(caphist[[i]])[1]))
        }
    }
    if (!is.list(traps))
        stop("traps must be a list")
    if (length(traps) != length(caphist))
        stop("number of sessions in traps does not match caphist")
    check.dim <- sapply(traps, nrow)
    if (!all(check.dim == caphist.dimensions[2, ]))
        stop("number of traps does not match caphist")
    if (!is.null(trapCovs)) {
        if (!is.list(trapCovs))
            stop("trapCovs must be a list")
        if (any(!sapply(trapCovs, is.list)))
            stop("trapCovs must be a list of lists")
        if (any(!unlist(sapply(trapCovs, function(x) sapply(x,
            is.data.frame)))))
            stop("trapCovs must be a list of dataframes")
        if (length(trapCovs) != length(caphist))
            stop("number of sessions in trapCovs does not match caphist")
        check.dim <- lapply(trapCovs, function(x) sapply(x, nrow))
        for (i in 1:length(check.dim)) {
            if (!all(check.dim[[i]] == caphist.dimensions[2,
                i]))
                stop("number of traps does not match caphist")
        }
    }
    if (!is.null(sigCovs)) {
        if (nrow(sigCovs) != length(caphist))
            stop("number of rows in sigCovs does not match number of sessions")
        if (!"session" %in% colnames(sigCovs)) {
            sigCovs$session <- factor(1:n.sessions)
        }
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    else {
        sigCovs <- data.frame(session = factor(1:n.sessions))
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    if (!is.null(trapOperation)) {
        if (!is.list(trapOperation))
            stop("trapOperation must be a list")
        if (length(trapOperation) != length(caphist))
            stop("number of sessions in trapOperation does not match caphist")
        check.dim <- sapply(trapOperation, nrow)
        if (!all(check.dim == caphist.dimensions[2, ]))
            stop("number of traps does not match caphist")
    }
    max.dist <- NULL
    for (i in 1:length(caphist)) {
        for (j in 1:nrow(caphist[[i]])) {
            if (dim(caphist[[i]])[3] > 1) {
                where <- apply(caphist[[i]][j, , ], 1, sum) >
                  0
            }
            else {
                where <- caphist[[i]][j, , ] > 0
            }
            if (sum(where) > 1)
                max.dist <- c(max.dist, max(0, dist(traps[[i]][where,
                  c("X", "Y")]), na.rm = T))
        }
    }
    mmdm <- mean(max.dist[max.dist > 0], na.rm = T)
    mdm <- max(max.dist, na.rm = T)
    if (!is.null(telemetry)) {
        if (!is.list(telemetry$fixfreq))
            stop("telemetry$fixfreq must be a list")
        fixfreq.dimensions <- sapply(telemetry$fixfreq, dim)
        if (nrow(fixfreq.dimensions) == 2)
            fixfreq.dimensions <- rbind(fixfreq.dimensions, 1)
        if (!is.null(telemetry$indCovs)) {
            if (!is.list(telemetry$indCovs))
                stop("telemetry$indCovs must be a list")
            if (any(!sapply(telemetry$indCovs, is.data.frame)))
                stop("telemetry$indCovs must be a list of dataframes")
            if (length(telemetry$indCovs) != length(telemetry$fixfreq))
                stop("number of sessions in telemetry$indCovs does not match telemetry$fixfreq")
            check.dim <- sapply(telemetry$indCovs, nrow)
            if (any(check.dim != fixfreq.dimensions[1, ]))
                stop("number of individuals in telemetry$indCovs does not match telemetry$fixfreq")
            if (any(!names(indCovs[[1]]) %in% c(names(telemetry$indCovs[[1]]),
                "removed")))
                stop("indCovs do not match between capture and telemetry data")
        }
        if (!is.null(telemetry$cap.tel)) {
            if (!is.list(telemetry$cap.tel))
                stop("telemetry$indCovs must be a list")
            warning("make sure captured individuals w/ collars sorted first!")
        }
    }
    if (!is.null(rsfDF)) {
        library(FNN)
        rsfCovs <- names(rsfDF[[1]][, -c(1:2), drop = F])
        if (is.null(trapCovs)) {
            trapCovs <- list()
            length(trapCovs) <- n.sessions
            for (s in 1:n.sessions) {
                trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                  c("X", "Y")], traps[[s]][, c("X",
                  "Y")], 1)$nn.index)
                trapCovs[[s]] <- list()
                length(trapCovs[[s]]) <- caphist.dimensions[3,
                  s]
                for (k in 1:caphist.dimensions[3, s]) {
                  trapCovs[[s]][[k]] <- data.frame(rsfDF[[s]][trap.grid,
                    rsfCovs])
                  names(trapCovs[[s]][[k]]) <- rsfCovs
                }
            }
        }
        else {
            for (s in 1:n.sessions) {
                if (any(!rsfCovs %in% trapCovs[[s]][[1]])) {
                  miss.rsfCovs <- rsfCovs[which(!rsfCovs %in%
                    trapCovs[[s]][[1]])]
                  trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                    c("X", "Y")], traps[[s]][, c("X",
                    "Y")], 1)$nn.index)
                  for (k in 1:caphist.dimensions[3, s]) {
                    newtrapCovs <- data.frame(rsfDF[[s]][trap.grid,
                      miss.rsfCovs])
                    names(newtrapCovs) <- miss.rsfCovs
                    trapCovs[[s]][[k]] <- data.frame(trapCovs[[s]][[k]],
                      newtrapCovs)
                  }
                }
            }
        }
    }
    scrFrame <- list(caphist = caphist, traps = traps, indCovs = indCovs,
        trapCovs = trapCovs, sigCovs = sigCovs, trapOperation = trapOperation,
        occasions = caphist.dimensions[3, ], type = type, mmdm = mmdm,
        mdm = mdm, telemetry = telemetry)
    class(scrFrame) <- "scrFrame"
    return(scrFrame)
}



e2dist <- function (x, y) {  # Function from scrbook package to calculate the distance between locations in 2 matrices.
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

e2dist.2 <- function (x, y) {  # Function from scrbook package to calculate the squared distance between locations in 2 matrices.
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- (x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

