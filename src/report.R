# output folder 
create_output_folder <- function(output_folder, verbose=FALSE){
  if (output_folder=="tmp") {
    # temporal_dir <- tempdir()
    # If the output folder is set to the default "tmp", create a unique folder name with a prefix "metabopathia_report"
    output_folder <- file.path("tmp",
                               paste0("metabopathia_report",length(list.files("tmp", pattern = "metabopathia_report")) + 1))
  }else{
    output_folder <- file.path(output_folder,
                               paste0("metabopathia_report",length(list.files(output_folder, pattern = "metabopathia_report")) + 1))
  }
  if (!file.exists(output_folder)){
    # If the output folder doesn't exist, create it
    dir.create(output_folder, showWarnings = T,recursive = T)
    if(verbose) message("Output folder created:", file.path(getwd(), output_folder))# maybe remove the getwd()
  }else if(verbose) warning(paste("The specified output folder '", output_folder, "' already exists. Results may be overwritten.", sep = ""))
  return(output_folder)
}
# Function to write a status value with a description to a file in the output folder
status <- function(value,descrip, output_folder) {
  # Construct the path to the status.txt file in the output folder
  status_file_path <- file.path(output_folder, "status.txt")
  # Write the value and description to the status.txt file, appending to existing content
  write(c(value,descrip), file = status_file_path, append = T, ncolumns = 2, sep = "% : ")
  if(verbose) message("Status updated to ", value, "% : ", descrip)
}

is_port_in_use <- function(port) {
  con <- tryCatch({
    socketConnection(host = "127.0.0.1", port = port, server = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(con)) {
    close(con)
    return(TRUE)
  } else {
    return(FALSE)
  }
}

serve_report <-function (output_folder, port = 4000, browser= T, daemon = T) 
{
  servr::httd(paste0(output_folder, "/pathway-viewer"), port = port, 
              browser = browser, daemon = daemon)
  cat("Open a web browser and go to URL http://127.0.0.1:", 
      port, "\n", sep = "")
}
