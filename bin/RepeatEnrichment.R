RepeatEnrichment <- function(repcnt.tbl, lib.sizes, conditions, targets) {
# Based on the count table from diffRepeats, calculate the repeat element 
# enrichment.
# Args:
#   repcnt.tbl: repeat element count table.
#   lib.sizes: named vector representing library sizes.
#   conditions: a factor describing experimental conditons: Control, Treatment.
#   targets: a factor describing target: must be Target or Input.
#       character vectors will be converted to factors.
    
    # Check input arguments.
    stopifnot(is.data.frame(repcnt.tbl))
    stopifnot(identical(colnames(repcnt.tbl)[1:3], 
                        c("Name", "Type", "Origin")))
    stopifnot(is.integer(lib.sizes) && !is.null(names(lib.sizes)))
    stopifnot(colnames(repcnt.tbl)[4:ncol(repcnt.tbl)] %in% names(lib.sizes))
    stopifnot(is.factor(conditions) || is.character(conditions))
    stopifnot(is.factor(targets) || is.character(targets))
    conditions <- as.factor(conditions)
    targets <- as.factor(targets)
    stopifnot(nlevels(targets) == 2 && 
                  levels(targets) %in% c("Target", "Input"))
    
    # Step1: normalize by library sizes.
    cnt.tbl <- repcnt.tbl[, 4:ncol(repcnt.tbl)]
    ordered.libsz <- lib.sizes[colnames(cnt.tbl)]
    norm.tbl <- t((t(cnt.tbl) + 1) / ordered.libsz * 1e6)  # pseudo-count: 1
    tar.idx <- which(targets == "Target")
    norm.tbl[, tar.idx][cnt.tbl[, tar.idx] == 0] <- 0
    
    # Step2: calculate enrichment ratios.
    enr.tbl <- sapply(levels(conditions), function(con) {
        # browser()
        inp.idx <- which(conditions == con & targets == "Input")
        if(length(inp.idx) > 1) {
            inp.avg <- colMeans(norm.tbl[, inp.idx])
        } else {
            inp.avg <- norm.tbl[, inp.idx]
        }
        tar.idx <- which(conditions == con & targets == "Target")
        
        norm.tbl[, tar.idx] / inp.avg
    })

    enr.names <- sapply(levels(conditions), function(con) {
        tar.idx <- which(conditions == con & targets == "Target")
        colnames(norm.tbl)[tar.idx]
    })
    enr.names <- paste("Enr(", enr.names, ")", sep="")
    colnames(enr.tbl) <- enr.names

    data.frame(repcnt.tbl, enr.tbl, check.names=F)
}
















