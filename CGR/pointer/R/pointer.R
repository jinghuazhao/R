# 18/3/2005
control.pointer <- function(sexlink=FALSE,
                   split=FALSE,
                   mating.type='00,01,02,10,11,12,20,21,22',
                   pointer.selection='0,4,5,6',
                   pointer.degree='1,2,3,4,5',
                   ascertainment='C')
{
   if (sexlink) pointer.degree <- '1,2,3,4'
   list(sexlink=sexlink,split=split,mating.type=mating.type,
        pointer.selection=pointer.selection,ascertainment=ascertainment)
}

# 15/3/2005
pointer <- function(datfile="poidat",jobfile="poijob", profile="poipro", 
                    terfile="poiter", control=control.pointer())
{
## house-keeping
   if (!file.exists(datfile)|!file.exists(jobfile)) stop("check if data/job file exists")
   unlink(profile)
   unlink(terfile)
#  in case there are leftovers from previous session
# wild card does not work
   junkfile <- c("INTERMED.1","INTERMEDFILE","INTERMED.2","INTERMED.3","INTERMED.ALP")
   unlink(junkfile)
  .Fortran("nucfama",datfile=as.character(datfile), jobfile=as.character(jobfile), 
            profile=as.character(profile), terfile=as.character(terfile), PACKAGE="pointer")
  if (control$sexlink) isex="S"
  else isex="A"
  if (control$split) {
    .Fortran("emx", profile=as.character(profile), terfile=as.character(terfile),
              isex=as.character(isex),
              fm1=as.character(control$mating.type),
              fm2=as.character(control$pointer.selection),
              fm3=as.character(control$pointer.degree),
              fm4=as.character(control$ascertainment), PACKAGE="pointer")
     file.rename("INTERMED.2", "INTERMEDFILE")
    .Fortran("pointr", profile=as.character(profile), terfile=as.character(terfile),
             isex=as.character(isex), PACKAGE="pointer")
     file.copy("INTERMED.3", "INTERMEDFILE",overwrite=TRUE)
    .Fortran("pointr", profile=as.character(profile), terfile=as.character(terfile),
             isex=as.character(isex), PACKAGE="pointer")
  } else {
     file.rename("INTERMED.1", "INTERMEDFILE")
    .Fortran("pointr",jobfile=as.character(jobfile), profile=as.character(profile),
              terfile=as.character(terfile), PACKAGE="pointer")
  }
  unlink(junkfile)
}
