#!/rds/user/jhz22/hpc-work/bin/Rscript --vanilla

library(tools)
pkg <- commandArgs(trailingOnly = TRUE)
print(pkg)
package_native_routine_registration_skeleton(pkg,con="package_native_routine_registration_skeleton.c")
