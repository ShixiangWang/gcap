local({
  install = local({
    if (!requireNamespace("ASCAT", quietly = TRUE)) {
      message("ASCAT is not available, installing it...")

      remotes::install_github("ShixiangWang/ascat@v3-for-gcap-v1", subdir = "ASCAT")
    }
  })
})
