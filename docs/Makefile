HOME_R       := $(HOME)/R
SOFTWARE_R   := /rds/project/rds-4o5vpvAowP0/software/R

CHANGELOG    := $(HOME_R)/gap/ChangeLog
GAP_CHECK    := $(HOME_R)/gap.Rcheck/gap-manual.pdf
SHINYGAP     := $(SOFTWARE_R)/gap/doc/shinygap.html
GAP_EXAMPLES := $(SOFTWARE_R)/gap.examples/doc/gap.examples.pdf
KINSHIP      := $(SOFTWARE_R)/kinship/doc/kinship.pdf
PAN          := $(SOFTWARE_R)/pan/doc/pan.pdf
TDTHAP       := $(SOFTWARE_R)/tdthap/doc/tdthap.pdf

vignettes/ChangeLog.txt: $(CHANGELOG)
	cp $< $@

vignettes/gap-manual.pdf: $(GAP_CHECK)
	cp $< vignettes/

vignettes/shinygap.html: $(SHINYGAP)
	cp $< vignettes/

vignettes/gap.examples.pdf: $(GAP_EXAMPLES)
	cp $< vignettes/

vignettes/kinship.pdf: $(KINSHIP)
	cp $< vignettes/

vignettes/pan.pdf: $(PAN)
	cp $< vignettes/

vignettes/tdthap.pdf: $(TDTHAP)
	cp $< vignettes/
