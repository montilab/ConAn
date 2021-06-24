test_that("lower_tri works", {
  mat <- matrix(1:9, nrow = 3, ncol = 3)
  expect_is(mat, "matrix")
  mat_lowertri_diag <- lower_tri(mat, diag = TRUE)
  expect_is(mat_lowertri_diag, "integer")
  mat_lowertri <- lower_tri(mat, diag = FALSE)
  expect_equal(sum(mat_lowertri_diag), 26)
  expect_is(mat_lowertri, "integer")
  expect_equal(sum(mat_lowertri), 11)
})

test_that("remove_na works", {
  vec <- c("test", "test1", NA, "test2")
  expect_is(vec, "character")
  vecNA <- vec[!is.na(vec)]
  rmNA <- remove_na(vec)
  expect_equal(length(vec), 4)
  expect_equal(length(vecNA), 3)
  expect_equal(rmNA, vecNA)
})

test_that("subtract_BG works", {
  mat <- matrix(1:9, nrow = 3, ncol = 3)
  bgmat <- mat - 1
  expect_is(mat, "matrix")
  expect_is(bgmat, "matrix")
  subtract <- mat-bgmat
  subtractfn <- subtract_bg(mat, bgmat)
  expect_equal(sum(subtract), 9)
  expect_equal(subtract, subtractfn)
})

test_that("utils:erase_mods works", {
  sampledata <- readRDS(file.path(system.file("extdata", package="ConAn"), "sampledata.rds"))
  mod_list <- readRDS(file.path(system.file("extdata", package="ConAn"), "samplemods.rds"))
  mat <- cor(t(exprs(sampledata)))
  testfunc <- erase_mods(mat, mod_list)
  expect_is(testfunc, "matrix")
  expect_true(is.na(testfunc[1,1]))
})
