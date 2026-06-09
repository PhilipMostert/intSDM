# Set the geodata_path to the user data directory
#
# This function is used in the testthat tests to ensure that the geodata_path
# is set to the user data directory, or a temporary directory. This is
# important, to avoid using other folders, such as the git repository itself.
#
# @param force logical; if `TRUE and geodata_path("user_data_dir") is
# not working, create a temporary directory instead. Default `TRUE`
# @param envir The environment in which to unlink the temporary directory,
# if created.
#
# @returns logical; If successful, returns `TRUE`, otherwise `FALSE`
#
# @details Can be improved by explicitly checking if the resulting geodata_path
# exists and is writable.
#
# @examples
# skip_if_not(local_testthat_geodata_path())
local_testthat_geodata_path <- function(force = TRUE, envir = parent.frame()) {
  res <- tryCatch(
    geodata::geodata_path("user_data_dir", persistent = FALSE)
  )
  if (!inherits(res, "error")) {
    return(TRUE)
  }
  if (!force) {
    return(FALSE)
  }
  geodata_temp_path <- tempfile("geodata_user_data_dir_")
  withr::defer(unlink(geodata_temp_path, recursive = TRUE), envir = envir)
  dir.create(geodata_temp_path, recursive = TRUE, showWarnings = FALSE)
  geodata::geodata_path(geodata_temp_path, persistent = FALSE)
  return(TRUE)
}

