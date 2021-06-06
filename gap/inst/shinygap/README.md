This is an initial attempt to enable easy calculation/visualization of three study designs.

## Family-based study

This is a call to `gap::fbsize()`. Note in particular that `alpha` parameter actually contains three elements and the 3rd of which is for ASP+TDT.

## Population-based study

This is a call to `gap::pbsize()`.

## Case-cohort study

This is a call to `gap::ccsize()`. Note that its `poper` argument indcates power (TRUE) or sample size (FALSE) calculation.

To set the default parameters, some compromises need to be made, e.g., Kp=1e-5~0.4, MAF=1e-3~0.8, alpha=5e-8~0.05, beta=0.01~0.4. The slider input provides upper bounds of a particular parameter.

It is somewhat surprising that although the gap vignette gap.Rmd gives good reasults, `fbsize()` does not handle `gamma` well, nor does ``ccsize()` with sample size calculation.
