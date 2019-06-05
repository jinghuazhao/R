# 4-12-2013 MRC-Epid JHZ

load pedigree 51.ped
load phenotypes 51.phen
model new
trait qt
polygenic -screen
trait qt
covariate sex
polygenic -screen -fix sex
