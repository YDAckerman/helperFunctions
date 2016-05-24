################################################################################
################################################################################
## data.frame of all the water use policies
stage <- c(1,2,1,2,1,1,2,2)

periodLabel <- c("A", "B", "C", "D", "E", "E" ,"F", "F")

startDate <- c("5-1-2009", "8-24-2009", "11-20-2009", "9-6-2011",
              "7-16-2012", "7-16-2012", "9-4-2012", "9-4-2012")

endDate <- c("8-24-2009", "11-20-2009","9-6-2011","7-16-2012",
            "9-4-2012","9-4-2012", NA, NA)

perWeek <- c(2,1,2,1,2,2,1,1)

hoursPerDay <- c(15,15,15,15,10,15,10,15)

hourSpecification <- c("12-10, 7-12","12-10, 7-12", "12-10, 7-12","12-10, 7-12",
                      "12-5, 7-12", "12-10, 7-12", "12-5, 7-12", "12-10, 7-12")

type <- c(NA,NA,NA,NA, "automatic", "hose-end","automatic", "hose-end")

notes <- c("Mandatory Restrictions per Code (May through September)",
           "Suppy Trigger",
           "Permanent Mandatory Year Round Water Restrinctions started",
           "Supply Trigger", "Supply Trigger", "Supply Trigger",
           "Supply Trigger", "Supply Trigger")

h2oRestrictions <- data.frame(
    stage = stage,
    periodLabel = periodLabel,
    startDate = as.Date(startDate, "%m-%d-%Y"),
    endDate = as.Date(endDate, "%m-%d-%Y"),
    perWeek = perWeek,
    hoursPerDay = hoursPerDay,
    hourSpecification = hourSpecification,
    type = type,
    notes = notes
)

rm(periodLabel, startDate, endDate, perWeek, hoursPerDay,
   hourSpecification, type, notes, stage)
