pathmix <- function (iop=3, datfile="data", jobfile="job", profile="prolix", terfile="summary") {
   if (iop==1) {
      cat("Testing ALMINI with documentation examples (to ALMTEST.OUT and ALMTEST.PLX)\n\n")
      .Fortran("almtest")
   }
   else if (iop==2) {
      cat("Testing GEMINI with documentation examples (to GEMTEST.OUT and GEMTEST.PLX)\n\n")
      .Fortran("gemtest")
   }
   else {
      if (!file.exists(datfile)|!file.exists(jobfile)) stop("data or job file does not exist")
      if (file.exists(profile)) unlink(profile)
      if (file.exists(terfile)) unlink(terfile)
      if (iop==3) .Fortran("path3a", filnam1=as.character(jobfile), filnam2=as.character(datfile),
                           filnam3=as.character(terfile), filnam4=as.character(profile))
      else if (iop==4) .Fortran("path3b", filnam1=as.character(jobfile), filnam2=as.character(datfile),
                                filnam3=as.character(terfile), filnam4=as.character(profile))
  }
  invisible()
}
