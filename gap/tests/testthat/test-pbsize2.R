test_that("pbsize2 agrees with tscc for a one-stage additive design", {
  grr <- c(1.1, 1.3, 1.5, 2.0, 3.0)

  power_tscc <- vapply(grr, function(rr) {
    gap::tscc(
      model = "additive",
      GRR = rr,
      p1 = 0.1,
      n1 = 100,
      n2 = 100,
      M = 1,
      alpha.genome = 0.05,
      pi.samples = 1,
      pi.markers = 1,
      K = 0.1
    )$power[1]
  }, numeric(1))

  power_pbsize2 <- vapply(grr, function(rr) {
    gap::pbsize2(
      N = 200,
      fc = 0.5,
      alpha = 0.05,
      gamma = rr,
      p = 0.1,
      kp = 0.1,
      model = "additive"
    )
  }, numeric(1))

  expect_equal(power_pbsize2, power_tscc, tolerance = 1e-3)
})
