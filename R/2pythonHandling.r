pythonLoc <- virtualenv_starter()

use_python(pythonLoc)
pyBuiltins <- import_builtins()
install_pyFK <- function(envname="r-funcsKishore", reinstall = T) {
    if (reinstall && virtualenv_exists(envname)) {
        virtualenv_remove(envname)
    }
    if (!virtualenv_exists(envname))
        virtualenv_create(envname)
    
    virtualenv_install(envname = envname, ignore_installed = F, 
        packages = python_packages)
    use_virtualenv(envname)
}
install_pyFK()
# netx <- import("networkx",  delay_load = F)
# use_virtualenv("r-funcsKishore")


# os <- .Platform$OS.type
# # if (os == "windows")
# #     pythonLoc <- Sys.which("python") %>% unname
# # else {
# #    pythonLoc <- Sys.which("python3") %>% unname
# # }


# use_virtualenv()


