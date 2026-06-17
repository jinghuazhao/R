set.seed(0)
knitr::opts_chunk$set(
  out.extra = 'style="display:block; margin: auto"',
  fig.align = "center",
  fig.height = 8,
  fig.path = "web/",
  fig.width = 8,
  collapse = TRUE,
  comment = "#>",
  dev = "CairoPNG")

pkgs <- c("httr","httpuv","jsonlite","plumber","Rsamtools")
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
missing <- pkgs[!vapply(pkgs, has_pkg, logical(1))]
if (length(missing)) {
  warning("Missing packages: ", paste(missing, collapse = ", "))
}
installed <- pkgs[vapply(pkgs, has_pkg, logical(1))]
invisible(suppressMessages(lapply(installed, require, character.only = TRUE)))

# httpuv::startServer("0.0.0.0", 8000, list(
#   call = function(req) {
#     list(
#       status = 200L,
#       headers = list("Content-Type" = "text/plain"),
#       body = "Hello, world!"
#     )
#   })
# )

check_url_availability <- function(url) {
# Send a GET request to the URL
  response <- tryCatch({
    httr::GET(url)
  }, error = function(e) {
# If there was an error, return FALSE
    return(NULL)
  })
# If response is NULL or status code is not 200, return FALSE
  if (is.null(response) || httr::status_code(response) != 200) {
    message(paste("Warning: The URL", url, "is not accessible. Please check the link."))
    return(FALSE)
  }
# If status code is 200, the URL is available
  message(paste("The URL", url, "is accessible."))
  return(TRUE)
}

url_to_check <- "http://www.zotero.org/styles/nature-genetics"
is_available <- check_url_availability(url_to_check)

if (is_available) {
  message("Using CSL as usual.")
} else {
  message("Using fallback to local CSL file instead.")
}

# get_data <- function(filename, region)
# {
# # query_result <- seqminer::tabix.read(filename, region)
#   tbx <- Rsamtools::TabixFile(filename)
#   query_result <- Rsamtools::scanTabix(tbx, param = region)[[1]]
#   if (length(query_result) == 0) {
#     return(data.frame())
#   }
#   hdr <- c("Chromosome", "Position",
#            "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq",
#            "Effect", "StdErr", "logP",
#            "Direction", "HetISq", "HetChiSq", "HetDf", "logHetP", "N")
#   df <- read.table(text = paste(query_result, collapse = "\n"), sep = "\t", col.names = hdr,
#                    stringsAsFactors = FALSE)
#   return(df)
# }
# 
# plbr <- plumber::Plumber$new()
# plbr$handle("GET", "/tests", function(req, res) {
#   protein <- req$args$protein
#   region <- req$args$region
#   if (is.null(protein) || is.null(region)) {
#     res$status <- 400
#     return(list(error = "Both 'protein' and 'region' must be provided"))
#   }
#   filename <- file.path("tests",paste0(protein,"-1.tbl.gz"))
#   print(filename)
#   if (!file.exists(filename)) {
#     res$status <- 404
#     return(list(error = paste("File for", protein, "not found")))
#   }
#   data <- get_data(filename, region)
#   json_data <- jsonlite::toJSON(data, dataframe = "rows", na = "null")
#   res$setHeader("Content-Type", "application/json")
#   return(json_data)
# })
# options(width = 200)
# filename <- file.path("tests","IL.18R1-1.tbl.gz")
# region <- "2:102700000-103800000"
# data <- get_data(filename, region)
# head(data,1)
# plbr$run(port = 8001)

# tmp <- tempfile()
# curl::curl_download("http://localhost:8001/tests?protein=IL.18R1&region=2:102700000-103800000", tmp)
# df <- jsonlite::fromJSON(readLines(tmp)) |>
#       jsonlite::fromJSON(flatten=TRUE) |>
#       as.data.frame()
# dim(df)

# dir.create("content/assets", recursive = TRUE)
# dir.create("content/lib", recursive = TRUE)
# s <- httpuv::startServer(
#   host = "0.0.0.0",
#   port = 5000,
#   app = list(
#     call = function(req) {
#       list(
#         status = 200L,
#         headers = list(
#           'Content-Type' = 'text/html',
#           'Access-Control-Allow-Origin' = '*',
#           'Access-Control-Allow-Methods' = 'GET, POST, OPTIONS',
#           'Access-Control-Allow-Headers' = 'Content-Type'
#         ),
#         body = "Hello world!"
#       )
#     },
#     staticPaths = list(
#       "/assets" = "content/assets/", # Static assets
#       "/lib" = httpuv::staticPath("content/lib", indexhtml = FALSE),
#       "/lib/dynamic" = httpuv::excludeStaticPath()
#     ),
#     staticPathOptions = httpuv::staticPathOptions(indexhtml = TRUE)
#   )
# )
# cat("Server running at http://0.0.0.0:5000\n")
# s$stop()
