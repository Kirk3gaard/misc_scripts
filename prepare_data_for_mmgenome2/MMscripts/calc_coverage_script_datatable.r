args = commandArgs(trailingOnly=TRUE)
if (!length(args) == 3) {
  stop("Three arguments must be supplied\n(input depth file, and output coverage file, and the number of threads)", 
       call. = FALSE)
}
if(!require("data.table")) {
  install.packages("data.table")
  suppressMessages(require("data.table"))
}
setDTthreads(args[3])
cov <- fread(args[1], header = F, col.names = c("scaffold", "position", "coverage"))
output <- cov[,
              .(Average.coverage = round(sum(coverage)/max(position),digits = 2),
                Reference.length = max(position)),
              by = scaffold]
fwrite(output, file = args[2])