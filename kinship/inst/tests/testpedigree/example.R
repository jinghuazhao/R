postscript('example.ps')

dtest <- read.table("example.dat",
                    col.names=c("upn","momid","dadid","sex","affect"))
 
ptemp <- pedigree(id=dtest$upn,dadid=dtest$dadid,momid=dtest$momid,
		 sex=dtest$sex,affected=dtest$affect)


## default plot
plot(ptemp)
title('Original Plot: 1 Affection Status')


## create some 'extra' affected status values
## - sometimes want to look at who has diseases X, Y, or Z

set.seed(24)
dtest$affect2 <- ifelse(runif(50)>.8,2,1) 
dtest$affect3 <- ifelse(runif(50)>.75,2,1)
dtest$affect4 <- sample(1:2,50,replace=T)

plot(ptemp,affected=cbind(dtest$affect,dtest$affect2))
title('Add 2 Affection Status values')

plot(ptemp,affected=cbind(dtest$affect,dtest$affect2,dtest$affect3))
title('Add 3 Affection Status values')

plot(ptemp,affected=cbind(dtest$affect,dtest$affect2,
                          dtest$affect3,dtest$affect4))
title('Add 4 Affection Status values')



## sometimes you may want to identify a person as the founder

ind.col <- rep(1,50)
ind.col[dtest$upn==11] <- 2

plot(ptemp,col=ind.col)
title('Use color to identify the founders')

## often you need to indicate those who are dead

dtest$dead <- ifelse(runif(50)>.7,1,0)
plot(ptemp,status=dtest$dead)
title('Indicate Alive/Dead status')


## you can use any types of labels you want (such as to indicate marker values)

plot(ptemp, id=rep(paste('age','bmi',sep='\n'),50))
title('Use labels to show off features of individuals')

## This shows how you can move people around (I'd eventually like this to be
## point & click)

  # gives position within this family
  plot(ptemp,id=1:50,packed=F)
  title('Use labels to find id numbers')

  # see the positioning by siblings
  # first column indicates what order do siblings get plotted
  # second column indicates information about other relationships
  print(ptemp$hints) 

  #switch positioning of siblings 9 & 11
  print(ptemp$hints[c(9,11),1])
  ptemp$hints[9,1] <- 8
  ptemp$hints[11,1] <- 7
  print(ptemp$hints[c(9,11),1])
  plot(ptemp,id=1:50,packed=F)
  title('Switched positions of siblings 9 & 11')

  #switch positioning of husband/wife (4,2)
  print(ptemp$hints[2,2])
  ptemp$hints[2,2] <- 4
  print(ptemp$hints[2,2])
  plot(ptemp,id=1:50,packed=F)
  title('Switched positioning of husband/wife (4 & 2)')



#create someone with unknown gender, terminated birth

dtest$sex2 <- dtest$sex
dtest$sex2[dtest$upn==41] <- 3 #unknown
dtest$sex2[dtest$upn==53] <- 4
 
ptemp2 <- pedigree(id=dtest$upn,momid=dtest$momid,dadid=dtest$dadid,
		 sex=dtest$sex2,affected=dtest$affect)

plot(ptemp2)
title('Add in unknown gender (diamond), terminated pregnancy (triangle)')

dev.off()
