/*bootstrap*/
proc surveyselect data=gwas out=gwas_boot seed = 1234567 method = urs
                  samplerate = 1 rep = &rep;
run;

/*delete-d jacknife*/
proc surveyselect data=gwas out=gwas_jack seed = 7654321 method = srs
                  sampsize = 6000 rep = &rep;
run;

data gwas; 
     set gwas; 
     replicate=0;
run;
data x.gwas_jack;
     set gwas gwas_jack;
run;
