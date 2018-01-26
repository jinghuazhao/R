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
