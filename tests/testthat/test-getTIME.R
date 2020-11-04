test_that("Get the TIME phenotype and TIMEscore", {
  load('data.Rdata')
  rr <- getTIME(data2)
  expect_identical(getTIME(data2),rr)
})
