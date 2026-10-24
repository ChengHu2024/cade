# ============================================================================
# Package-level setup
# ============================================================================

## ggplot2 uses non-standard evaluation: inside aes(), column names are not
## evaluated in the calling environment, but R CMD check cannot see the data
## frame they refer to and reports them as undefined globals. Declaring them
## here silences the note without changing behaviour.
utils::globalVariables(c("comp", "gene", "label", "obs", "pred", "value"))
