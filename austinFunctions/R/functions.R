################################################################################
################################################################################
## austin helper functions
data(sysdata, envir=environment())

#' Get Min Edges From Source
#'
#' Get Minimum Edge Length from Source Node
#' @param node numeric - node number
#' @keywords
#' @export
#' @examples
getMinEdgesFromSource <- function(node){
    if(!(node %in% edges$DestinationNodeInd)){
        return(0)
    }
    sources <- (edges %>%
                filter(DestinationNodeInd == node) %>%
                distinct(SourceNodeInd))$SourceNodeInd
    print(sources)
    min(sapply(sources, function(s){
        1 + getMinEdgesFromSource(s)
    }))
}

## returns the maximum number of edges between the
## node and Lake Austin.
#' Get Max Edges From Source
#'
#' Get Maximum Edge Length From Source Node
#' @param node numeric - node number
#' @keywords
#' @export
#' @examples
getMaxEdgesFromSource <- function(node){
    if(!(node %in% edges$DestinationNodeInd)){
        return(0)
    }
    sources <- (edges %>%
                filter(DestinationNodeInd == node) %>%
                distinct(SourceNodeInd))$SourceNodeInd
    print(sources)
    max(sapply(sources, function(s){
        1 + getMaxEdgesFromSource(s)
    }))
}

#' Get Destination
#'
#' Get all nodes leading from a node
#' @param node numeric - node number
#' @keywords
#' @export
#' @examples
getDestination <- function(node){
    ret <- edges %>%
        filter(SourceNodeInd == node) %>%
        left_join(nodes, by = c("DestinationNodeInd" = "nodeInd")) %>%
        filter(Type == "Pressure Zone")
    ifelse(empty(ret), NA, ret$DestinationNodeInd)
}

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

#' Query Austin Database
#'
#' send a sql query to the Austin database
#' @param sql chracter - sql query to be sent
#' @keywords
#' @export
#' @examples
queryAustinDB <- function(sql){
    require(RODBC)
    source("/Volumes/Yoni/Passwords/passwords.R")
    con <- try(odbcConnect("cbnetworkbuilder",
                           uid = "Yoni_temp",
                           pwd = credentials$passwords$Yoni_temp))
    if(identical(class(con), "try-error")){
        warning("Password was incorrect, please check password and try again")
    }
    sqlQuery(con, sql)
}

## take the derivative of a 'function' of the form y = f(x)
## inputted as x, fy
#' d/dx
#'
#' take the derivative of a function represented as two index-aligned
#' vectors. returns the derivative functions represented in the same
#' format.
#' @param x numeric vector - sequence of domain points
#' @param fx numeric
#' @keywords
#' @export
#' @examples
d_dx <- function(x, fx){
    list(df = diff(fx) / diff(x),
         dx = rowMeans(embed(x, 2))) ## embed is cool
}



#' Plot Distribution
#'
#' plot and save a density distribution
#' @param df - data frame containing data
#' @param descriptors - the grouping factors whose values will
#'                      determine the file save location
#' @param value_name - the name of the column containing data to be plotted
#' @param location - directory in which to save the plot
#' @keywords
#' @export
#' @examples
plotDistribution <- function(df, descriptors, value_name, location){
    require(ggplot2)
    descriptors <- unique(df[,descriptors])
    qplot(df[, value_name], geom = "density")
    ggsave(filename = paste0(location, "/",
                             paste(descriptors, collapse = "_"), ".pdf"))
    return(df)
}

#' Perform K Means
#'
#' perform K means clustering on multidplyr/plyr-grouped data
#' return clustering results in the form of a data frame
#' @param df - dataframe containing the data to be clustered
#' @param k - number of clusters
#' @param .nstart - number of random points from which to start clustering
#' @keywords
#' @export
#' @examples

performKmeans <- function(df, k, .nstart){
    require(helperFunctions)
    tsColumns <- grep("H\\.", colnames(df), value = TRUE)
    tsData <- df[, tsColumns]
    ## account for NA values
    tsData <- cleanTimeSeries(tsData)
    cl <- kmeans(tsData, centers = k, nstart = .nstart)
    fitted.cl <- as.data.frame(fitted(cl))
    clusters <- fitted.cl %>% distinct()
    ## create a 'group' variable - it will indicated clusters more clearly
    clusters$group <- 1:nrow(clusters)
    fitted.cl <- dplyr::left_join(fitted.cl, clusters, by = colnames(fitted.cl))
    ## account for any erroneous rows:
    omitted <- attr(tsData, "omitted")
    if(!is.null(omitted)){
        df$group <- insert.vals(fitted.cl$group, NA, omitted)
    } else {
        df$group <- fitted.cl$group
    }
    df
}


#' Normalize Column
#'
#' normalize a column in a dataframe
#' @param df - dataframe
#' @param colname - name of column to normalize
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


#' Clean Time Series
#'
#' clean NAs from each row of a timeSeries data.frame.
#' if a row is too brutally NA'd to be cleaned, remove
#' it and set the 'omitted' attributed of the returned
#' timeSeries data.frame to reflect the omission.
#' @param timeSeries_df - dataframe with time series to clean
#' @keywords
#' @export
#' @examples

cleanTimeSeries <- function(timeSeries_df){
    ## this will carry the indices of any omitted rows
    omitted <- NULL
    timeSeries_df <- ldply(1:nrow(timeSeries_df), function(i){
        timeSeries <- try(smoothNAs(timeSeries_df[i, ]))
        if(identical(class(timeSeries), "try-error")){
            ## go up 4 levels to find 'omitted':
            ## (call -> assign -> temporary function ->
            ##                            ldply -> cleanTimeSeries)
            omitted <- get("omitted", envir = parent.frame(n = 4))
            omitted <- c(omitted, i)
            assign("omitted", omitted, envir = parent.frame(n = 4))
            return(NULL)
        }
        timeSeries
    })
    attr(timeSeries_df, "omitted") <- omitted
    return(timeSeries_df)
}


#' Smooth NAs
#'
#' clear NA values from a timeseries vector by converting
#' them to the average of the non-na values nearest them
#' in the series
#' @param timeSeries - vector timeseries to be smoothed over
#' @keywords
#' @export
#' @examples
smoothNAs <- function(timeSeries){
    if(sum(is.na(timeSeries)) > length(timeSeries) / 2){
        stop("The time series is unsalvagable: too many NA's")
    }
    i <- which(is.na(timeSeries))
    j <- which(!is.na(timeSeries))
    ## loop through 
    for (val in i){
        m <- min(j[j > val])
        M <- max(j[j < val])
        m <- ifelse(is.finite(m), m, NA)
        M <- ifelse(is.finite(M), M, NA)
        ## unlist the timeSeries in case it is a data.frame
        new_val <- mean(unlist(timeSeries)[c(m, M)], na.rm = TRUE)
        timeSeries[val] <- new_val
    }
    timeSeries
}


#' Get Week
#'
#' determine the week of the month the
#' given date resides in
#' @param date - vector of dates
#' @keywords
#' @export
#' @examples
getWeek <- function(date){
    f <- Vectorize(function(date){
        weeks <- sort(rep(1:5, 7))
        day <- as.numeric(format(date, format = "%d"))
        weeks[day]
    })
    f(date)
}

#' Get Restriction
#'
#' Get the usage restriction in effect on a given date
#' @param
#' @keywords
#' @export
#' @examples
getRestriction <- function(date, restriction){
    ret <- h2oRestrictions %>%
        dplyr::select_("startDate", "endDate", restriction) %>%
        dplyr::distinct_("startDate", "endDate", restriction) %>%
        dplyr::filter(endDate > date & startDate <= date)
    ret[1, restriction]
}

## given a 24 hour timeseries, find the approximate area below
## the 'curve' between the specified hours (given as a list of
## hour pair vectors)
#' Function Name
#'
#' Function Description
#' @param series
#' @param hours
#' @keywords
#' @export
#' @examples
areaBelowTimeSeries <- function(series, hours){
    stopifnot(length(hours) >= 1)
    unlist(lapply(hours, function(pair){
        pair <- pair + 1
        sum(series[pair[1]:pair[2]])
    }))
}

#' Get Usage During (a period)
#'
#' to be used with dplyr:
#' given a data.frame with 'Hour' as a variable
#' that ranges uniquely across 0-23 and 'Value' as the
#' quantity at each time, find the percent
#' of the total area between the given hour pairs 
#' (given as list of vector pairs)
#' @param df group segment of a data frame
#' @param hours - vector of hours, each sequential pair defines a period
#' @keywords
#' @export
#' @examples

getUsageDuring <- function(df, hours){
    ts <- df$Value
    if(length(ts) != 24){
        ## instead of trying to extrapolate, etc.
        ## I'm just going to convert the values
        ## to NA and eventually remove them.
        areas <- rep(NA, length(hours))
    } else {
        areas <- areaBelowTimeSeries(ts, hours = hours)
    }
    tmp <- df %>%
        dplyr::select(-Hour, -Value) %>%
        dplyr::distinct()
    for (i in 1:length(hours)) {
        new_col_name <- paste0("period", paste(hours[[i]], collapse = "to"))
        tmp[, new_col_name] <- areas[i]
    }
    tmp
}


#' Sum To Index
#'
#' count number of values since the
#' last TRUE value in a given vector i.e.
#' c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE)
#' -> c(1, 0, 0, 1, 2, 0, 1, 2, 3, 0)
#' note: this is biased towards 1's...
#' @param vector boolean vector
#' @keywords
#' @export
#' @examples
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


#' Sum Dataframe Column To Index
#'
#' dplyr wrapper for sumToIndex
#' @param df - dataframe containing column to 'sum to index'
#' @param colname - name of the column within df to 'sum to index'
#' @param new_colname - name for the new 'summed' column
#' @keywords
#' @export
#' @examples
sumDfColToIndex <- function(df, colname, new_colname){
    df[, new_colname] <- sumToIndex(unlist(df[, colname]))
    df
}


#' Most Frequent
#'
#' dplyr wrapper to get most frequent value of a variable
#' in a grouped table
#' @param df dataframe containing column to freq
#' @param colname name of the column of interest
#' @keywords
#' @export
#' @examples
mostFrequent <- function(df, colname){
    vals <- table(unlist(df[, colname]))
    i <- which.max(vals)
    names(vals)[i]
}

#' Inquire Reload
#'
#' function to call the reload statement:
#' @param yesNoQuestion - character yes or no question
#' @param answer - supply an answer to skip the process
#' @keywords
#' @export
#' @examples
inquireReload <- function(yesNoQuestion, answer = NULL){
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


#' Summarize Dataframe
#'
#' quantiles, means, variances for a dataframe's variable
#' @param df dataframe
#' @param column_name name of column to summarize
#' @param probs quantiles to summarize by
#' @keywords
#' @export
#' @examples
summarizeDF <- function(df, column_name, probs = seq(0, 1, .25)){
    if(empty(df)){
        data.frame(Q1 = NA,
                   median = NA,
                   Q3 = NA,
                   mean = NA,
                   variance = NA,
                   min = NA,
                   max = NA
                   )
    } else {
        data <- df[, column_name]
        quantiles <- unname(quantile(data, probs = probs, na.rm = TRUE))
        data.frame(Q1 = quantiles[2],
                   median = quantiles[3],
                   Q3 = quantiles[4],
                   mean = mean(data, na.rm = TRUE),
                   variance = var(data, na.rm = TRUE),
                   min = min(data, na.rm = TRUE),
                   max = max(data, na.rm = TRUE)
                   )
    }
}

#' calcPeriodAggPumpData
#'
#' Calculate period-aggregated pump flow data
#' @param pumpFlowData data frame of the data to be aggregated
#' @keywords
#' @export
#' @examples
calcPeriodAggPumpData <- function(pumpFlowData = NULL){
    if(!is.null(pumpFlowData)){
        ## for each node, restriction, ymd-date, and
        ## each hour within that date, get the sum of
        ## all the timeseries' values. 
        periodAggPumpFlow <- pumpFlowData %>%
            ## divide by 24 to get MG/H from MG/D and multiply
            dplyr::mutate(Value = Value / 24) %>%
            dplyr::group_by(max_temp, RestPeriod, Restriction, Name,
                     NetworkNodeInd, TSDate, Hour) %>%
            dplyr::summarise(Value = sum(Value, na.rm = TRUE)) %>%
            ungroup() %>%
            ## now for each node, restriction, ymd-date,
            ## order by hour and pass to custom function getUsageDuring
            ## to get the area under the curve for the stated
            ## hour pairs
            group_by(max_temp, RestPeriod, Restriction, Name,
                     NetworkNodeInd, TSDate) %>%
            dplyr::arrange(Hour) %>%
            do(getUsageDuring(., hours = list(c(0, 23), c(0, 5),
                                              c(6, 21), c(6, 13),
                                              c(21, 21), c(22, 23),
                                              c(14, 19)))) %>%
            ungroup() %>%
            dplyr::mutate(offpeak = period22to23 + period0to5,
                          midpeakSummer = period6to13 + period21to21) %>%
            dplyr::select(-period22to23, -period0to5,
                          -period6to13, -period20to21) %>%
            ## give everyone appropriate names
            dplyr::rename(midpeak = period6to21,
                          peak = period14to20,
                          total = period0to23)
        ## Aggregate up to the full network
        periodAggPumpFlowTot <- periodAggPumpFlow %>%
            dplyr::select(max_temp, RestPeriod, Restriction, Name,
                          NetworkNodeInd, TSDate, total,
                          midpeak, midpeakSummer, peak, offpeak) %>%
            dplyr::group_by(max_temp, RestPeriod, Restriction, TSDate) %>%
            dplyr::summarise(total = sum(total, na.rm = TRUE),
                             midpeak = sum(midpeak, na.rm = TRUE),
                             midpeakSummer = sum(midpeakSummer, na.rm = TRUE),
                             peak = sum(peak, na.rm = TRUE),
                             offpeak = sum(offpeak, na.rm = TRUE)) %>%
            dplyr::ungroup()
        save(periodAggPumpFlow, periodAggPumpFlowTot,
             file = "/Volumes/Yoni/AustinR/periodAggPumpFlow.rda")
    } else {
        load("/Volumes/Yoni/AustinR/periodAggPumpFlow.rda")    
    }
    assign(x = "periodAggPumpFlow", periodAggPumpFlow,
           envir = globalenv())
    assign(x = "periodAggPumpFlowTot", periodAggPumpFlowTot,
           envir = globalenv())
}

#' Load Austin Energy
#'
#' loads the austin network energy (kwh) data
#' @param dataset - the type of data to be retrieved
#' @keywords
#' @export
#' @examples
load_austinData <- function(dataset = c("flow", "energy", "usage"))
{
    dataset <- "usage"
    ## passwords are stored on the external hard drive
    source("/Volumes/Yoni/Passwords/passwords.R")
    con <- try(odbcConnect("cbnetworkbuilder",
                           uid = "Yoni_temp",
                           pwd = credentials$passwords$Yoni_temp))
    if(identical(class(con), "try-error")){
        warning(paste0("Password was incorrect, please ",
                       "check password and try again"))
    }
    ## Statement to extract all mgd pump time series data.
    conditions <- switch(dataset,
                              flow = list(where = " WHERE Attributes = 'Flow'",
                                          assetInd = " 1"),
                              energy = list(where = " WHERE Unit = 'kwh'",
                                            assetInd = " 1"),
                              usage = list(where = paste0(" WHERE Unit = 'MGD'",
                                                          " AND Name LIKE",
                                                          " '%Usage%'"),
                                           assetInd = " 3"))
    sql <- paste0("SELECT TSDate, TSTime, Value, NetworkNodeInd,
                         TimeSeriesData.TimeSeriesInd AS TimeSeriesInd
                  FROM   TimeSeriesData INNER JOIN (
                         SELECT * FROM dbo.TimeSeriesToNetworkNodes
                         WHERE NetworkNodeInd IN (
                               SELECT Ind FROM dbo.NetworkNodes
                               WHERE NetworkInd = 40 AND AssetInd =",
                         conditions[["assetInd"]], ")) tn
                         ON TimeSeriesData.TimeSeriesInd = tn.TimeSeriesInd
                  WHERE  TimeSeriesData.TimeSeriesInd IN (
                         SELECT Ind FROM dbo.TimeSeries",
                         conditions[["where"]], ")")
    austinData <- sqlQuery(con, sql) %>%
        dplyr::mutate(
            TSDate = as.Date(TSDate, format = "%Y-%m-%d"),
            ymd = as.Date(TSDate, format = "%Y-%m-%d")) %>%
        dplyr::mutate(
            Year = lubridate::year(TSDate),
            Month = lubridate::month(TSDate),
            Day = lubridate::day(TSDate),
            Hour = lubridate::hour(TSDate)) %>%
        dplyr::mutate(jday = base::julian(ymd),
                  yday = lubridate::yday(ymd),
                  weekday = weekdays(ymd))
    austinData$weekday <- factor(pumpEnergyData$weekday,
                                   levels = c("Monday", "Tuesday", "Wednesday",
                                              "Thursday", "Friday", "Saturday",
                                              "Sunday"))
    austinData <- left_join(austinData, load_austinNetwork()[["nodes"]],
                         by = c("NetworkNodeInd" = "nodeInd"))
    ## close the connection
    close(con)
    return(austinData)
}

#' load Austin Network
#'
#' load austin network data, return a list with nodes and edges.
#' @param NONE
#' @keywords
#' @export
#' @examples
load_austinNetwork <- function()
{
    ## passwords are stored on the external hard drive
    source("/Volumes/Yoni/Passwords/passwords.R")

    con <- try(odbcConnect("cbnetworkbuilder",
                           uid = "Yoni_temp",
                           pwd = credentials$passwords$Yoni_temp))

    if(identical(class(con), "try-error")){
        warning("Password was incorrect, please check password and try again")
    }

    ## Statement to extract the network edges
    sql <- "SELECT * FROM NetworkEdges
        WHERE DestinationNodeInd IN (
              SELECT Ind FROM NetworkNodes
              WHERE NetworkInd = 40)
        OR SourceNodeInd IN (
              SELECT Ind FROM NetworkNodes
              WHERE NetworkInd = 40)"
    edges <- sqlQuery(con, sql)

    ## Statement to extract the network nodes
    sql <- "SELECT NetworkNodes.Ind as nodeInd, Name, Type FROM
               NetworkNodes LEFT JOIN Assets
               ON NetworkNodes.AssetInd = Assets.Ind
               WHERE NetworkInd = 40"

    nodes <- sqlQuery(con, sql)
    ## close the connection
    close(con)
    return(list("nodes" = nodes, "edges" = edges))
}

