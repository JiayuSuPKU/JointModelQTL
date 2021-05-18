load_models <- function() {
  stan_dirs <- dir(
    file.path("src", "stan_models"),
    full.names = TRUE,
    recursive = FALSE
  )

  stan_models <- lapply(stan_dirs, function(d) {
    cat(paste("Loading stan models in", d, "\n"))

    stan_files <- dir(d, full.names = TRUE, pattern = "\\.stan$")
    model_names <- sub("\\.stan$", "", basename(stan_files))

    models <- lapply(stan_files, function(f) {
      rstan::stan_model(file = f, obfuscate_model_name = FALSE)
    })

    names(models) <- model_names

    return(models)
  })

  names(stan_models) <- basename(stan_dirs)

  return(stan_models)
}
