#' Helper function to replace slashes and backslashes depending on the OS type.
#'
#' @keywords internal
#'
#' @param pathstring A string of a filepath containing slashes or backslashes.
#'
#' @return A OS-dependent filepath string.

fileplatform <- function(pathstring) {
    interm <- pathstring
    ostype <- Sys.info()['sysname']
    if (ostype == 'Windows') {
        interm <- stringr::str_replace_all(interm, "/", "\\\\")
    } else {
        interm <- stringr::str_replace_all(interm, "\\\\", "/")
    }
    return(interm)
}
