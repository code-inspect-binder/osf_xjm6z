message(".Rprofile loaded: Starting initialization...")

if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Step 1: Confirmed interactive RStudio session. Waiting before next step...")

  later::later(function() {
    message("Step 2: Waiting before setting working directory...")

    later::later(function() {
      try({
        setwd("xjm6z_src/Annotated Analysis Script")
        message("Working directory set to: ", getwd())
      }, silent = FALSE)
    }, delay = 2)  # Second delay before setwd

  }, delay = 2)  # First delay after detecting RStudio

} else {
  message("interactive() is FALSE. Retrying in 3 seconds...")

  later::later(function() {
    if (interactive() && Sys.getenv("RSTUDIO") == "1") {
      message("Retried and now in RStudio interactive mode.")

      later::later(function() {
        try({
          setwd("xjm6z_src/Annotated Analysis Script")
          message("Working directory set to: ", getwd())
        }, silent = FALSE)
      }, delay = 2)
    } else {
      message("Still not in interactive RStudio session after retry.")
    }
  }, delay = 3)
}
