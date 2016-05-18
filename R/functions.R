############################################################
## generally useful helper functions
############################################################

## helper function to create a cluser of
## cores with the doSNOW Package
#' Create Cluster
#'
#' initiate a cluster for .parallel call
#' @param noCores numeric - number of cores
#' @param logfile character - name of file to write logs to
#' @param export named list - items to export to each cluster
#' @param lib named list - libraries to export to each cluster
#' @keywords
#' @export
#' @examples
createCluster <- function(noCores, logfile = "/dev/null",
                        export = NULL, lib = NULL) {
    require(doSNOW)
    require(plyr)
    cl <- makeCluster(noCores, type = "SOCK", outfile = logfile)
    if(!is.null(export)) clusterExport(cl, export)
    if(!is.null(lib)) {
        l_ply(lib, function(dum) {
            clusterExport(cl, "dum", envir = environment())
            clusterEvalQ(cl, library(dum, character.only = TRUE))
        })
    }
    registerDoSNOW(cl)
    return(cl)
}

#' Multi-Grepl
#'
#' Grepl multiple patterns in a single string
#' @param patterns vector of strings to search for
#' @param x character string in which to search for patterns
#' @param ignore.case boolean - is case important
#' @param strict integer - how many patterns should match for a TRUE
#' @param fuzzy boolean - allow fuzzy matching
#' @param perl boolean - allow perl regex
#' @keywords
#' @export
#' @examples
mgrepl <- function(patterns, x, ignore.case = FALSE,
                      strict = 0, fuzzy = FALSE, perl = FALSE){
    require(plyr)
    f <- if (fuzzy){ agrepl } else { grepl }
    bool_df <- ldply(patterns, function(pattern){
        f(pattern, x, ignore.case = ignore.case, perl = perl)
    })
    colSums(bool_df) > strict
}

#' Multi-Grep
#'
#' Grep multiple patterns in a single string
#' @param SEE_mgrepl
#' @param value boolean - return positions or values
#' @keywords
#' @export
#' @examples
mgrep <- function(patterns, x, ignore.case = FALSE,
                     strict = 0, fuzzy = FALSE, value = FALSE, perl = FALSE){

    i <- mgrepl(patterns, x, ignore.case = ignore.case,
              strict = strict, fuzzy = fuzzy, perl = perl)

    if(value){
        x[i]
    } else {
        seq_along(x)[i]        
    }
}

#' Trim
#'
#' Remove leading/following blank space
#' @param x string to trim
#' @keywords
#' @export
#' @examples
trim <- function(x){ gsub("^\\s+|\\s+$", "", x)}

#' String Intersection
#'
#' Find all the words in common between two strings
#' @param s1 the first string
#' @param s2 the second string
#' @keywords
#' @export
#' @examples
stringIntersect <- function(s1, s2){
    require(plyr)
    s <- strsplit(gsub("[^[:alnum:]]", " ", c(s1,s2)), " +")
    s <- llply(s, function(si){
        si[which(si != "")]
    })
    intersect(tolower(s[[1]]), tolower(s[[2]]))
}

#' String containment
#'
#' See if string2 contains all of string1 (only alpha numeric)
#' @param s1 character string
#' @param s2 character string
#' @keywords
#' @export
#' @examples
stringContains <- function(s1, s2){
    s <- strsplit(gsub("[^[:alnum:]]", " ", c(s1,s2)), " ")
    all(tolower(s[[1]]) %in% tolower(s[[2]]))
}

#' Remove Parentheses
#'
#' Remove any parenthetical content from a character string
#' @param s character string to be culled of parentheses
#' @keywords
#' @export
#' @examples
removeParens <- function(s){
    gsub("(?=\\().*?(?<=\\))", "", s, perl=T)
}

#' Fuzzy string-set containment
#'
#' Match one set of strings within another set of strings
#' @param s1 vector of patterns
#' @param s2 vector of strings to be searched
#' @keywords
#' @export
#' @examples
aIn <- function(s1, s2, ...){

    if (all(is.na(s1) & all(is.na(s2)))){ return(TRUE) }
    
    s1 <- na.omit(s1)
    s2 <- na.omit(s2)
    sapply(s1, function(x){
        any(agrepl(x, s2, ...))
    })
}

#' 'Random' data frame
#'
#' Create a random dataframe of uniformly sampled integers
#' less than or equal to the product of the dimensions given
#' @param rows - integer number of rows
#' @param cols - integer number of cols
#' @keywords
#' @export
#' @examples
RandDF <- function(rows, cols){
    reps <- 1:ceiling(cols / 26)
    column.names <- c()
    for (i in reps){
        adenda <- sapply(letters, function(letter){
            paste(rep(letter, i), collapse = "")
        }, USE.NAMES = FALSE)
        
        column.names <- c(column.names, adenda)
    }
    data <- sample(1:(rows * cols), rows * cols,
                   replace = TRUE)
    tmp <- data.frame(matrix(data, nrow = rows, ncol = cols))
    colnames(tmp) <- column.names[1:cols]
    tmp
}

#' Segment a vector
#'
#' Chop a vector up into equally sized parts
#' @param vec the vector to butcher
#' @param num.segs number of segments to partition into
#' @param ordered boolean - should the vector be ordered?
#' @keywords
#' @export
#' @examples
SegmentVec <- function(vec, num.segs, ordered = TRUE){
    if(ordered){
        vec <- sort(vec, na.last = FALSE)
    }
    segs <- sort(rep(1:num.segs,length = length(vec)))
    segments <- lapply(1:num.segs, function(i) vec[which(segs == i)] )
    return(segments)
}

#' Insert values
#' 
#' Extend a vector in internal positions with given values
#' @param vector vector to be extended
#' @param values values to add to vector
#' @param indices locations in vector to insert values
#' @keywords
#' @export
#' @examples
insert.vals <- function(vector, values, indices){
    if(length(values) != length(indices)){
        values <- rep(values, length(indices))
    }
    indices <- sort(indices)
    index <- indices[1]
    new_value <- values[1]
    if(index == 1){
        vector <- c(new_value, vector)
    } else if(index >= length(vector)){
        vector <- c(vector, new_value)
    } else {
        vector <- c(
            vector[1:(index - 1)],
            values[1],
            vector[index:length(vector)]
        )
    }
    if(length(indices) != 1){
        indices <- indices[2:length(indices)]
        values <- values[2:length(values)]
        vector <- insert.vals(vector, values, indices)
    }
    vector
}

#' Multiple Gsub
#'
#' Replace multiple patterns within a character string
#' @param pattern vector of patterns to replace
#' @param replacement vector of replacements, synced to pattern
#' @param x character string to modify
#' @keywords
#' @export
#' @examples
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}

#' Symmetric Difference
#'
#' Determine the symmetric difference between two sets.
#' @param A set
#' @param B set
#' @keywords
#' @export
#' @examples
symdiff <- function(A, B){
    c(setdiff(A, B), setdiff(B, A))
}

#' Inquire User
#'
#' Ask the user a yes/no question
#' @param yesNoQuestion ...
#' @param answer set a default answer for testing
#' @keywords
#' @export
#' @examples
inquire <- function(yesNoQuestion, answer = NULL){
    if(is.null(answer)){
        ans <- readline(yesNoQuestion)
        while(ans != "n" & ans != "y"){
            ans <- readline("please answer y/n: ")
        }
    } else {
        return(answer)
    }
    return(ans)
}

#' Normalize Column
#'
#' Normalize a column within a dataframe
#' @param df data.frame containing the column
#' @param colname character string of the column's name
#' @keywords
#' @export
#' @examples
normalizeColumn <- function(df, colname){
    x <- df[, colname]
    new_x <- (x - min(x, na.rm = TRUE)) /
        (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    df[,paste("Normalized", colname, sep = "")] <- new_x
    df
}

#' Accumulative summing
#'
#' count number of values since the last TRUE value in a given vector
#' @param vector boolean vector
#' @keywords
#' @export
#' @examples  c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE,
#'              FALSE, FALSE, FALSE, TRUE) -> c(1, 0, 0, 1, 2, 0, 1, 2, 3, 0)
sumToIndex <- function(vector){
    vec_length <- length(vector)
    new_vector <- rep(0,vec_length)
    valued_indices <- which(vector == TRUE)
    j <- 1
    for (i in 1:vec_length){
        if(i %in% valued_indices){
            new_vector[i] <- 0
            j <- 1
        } else {
            new_vector[i] <- j
            j <- j + 1
        }
    }
    new_vector
}

############################################################
############################################################
