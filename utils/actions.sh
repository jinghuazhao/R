# The canonical actions are inappropriate for a cluster of R packages here.

Rscript -e '
  library(usethis)
  use_github_actions()
  use_github_actions_badge(name = "R-CMD-check.yaml", repo_spec = NULL)
  use_github_action_check_release(save_as = "R-CMD-check.yaml", ref = NULL, ignore = TRUE, open = FALSE)
'

<!-- badges: start -->
[![R-CMD-check](https://github.com/jinghuazhao/R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jinghuazhao/R/actions/workflows/R-CMD-ch>
<!-- badges: end -->
