.onAttach <- function(libname, pkgname) {
  packageStartupMessage("A R package to test ecological patches \nWARNING: This package is under development.")
}

# .onUnload <- function (libpath) {
#   library.dynam.unload('pat', libpath)
# }

