.onAttach <- function(libname, pkgname) {
  version <- packageDescription(pkgname, fields = "Version")

  msg <- paste0(pkgname, " version ", version, "\n- Project URL at https://github.com/ShixiangWang/gcap")
  packageStartupMessage(msg)
}
