test_that("connectivity:modular_differential_connectivity works", {
  sampledata <- readRDS(file.path(system.file("extdata", package="ConAn"), "sampledata.rds"))
  mod_list <- readRDS(file.path(system.file("extdata", package="ConAn"), "samplemods.rds"))
  r_edat <- t(exprs(sampledata))[1:10,]
  t_edat <- t(exprs(sampledata))[11:20,]
  testfunc <- modular_differential_connectivity(r_edat, t_edat, 1, 0.9, "difference")
  expect_is(testfunc, "numeric")
})

test_that("connectivity:lower_tri_erase_mods_cor works", {
  sampledata <- readRDS(file.path(system.file("extdata", package="ConAn"), "sampledata.rds"))
  mod_list <- readRDS(file.path(system.file("extdata", package="ConAn"), "samplemods.rds"))
  erased <- lower_tri_erase_mods_cor(t(exprs(sampledata)), mod_list)
  expect_is(erased, "numeric")
})

test_that("connectivity:atanh_lower_tri_cor works", {
  sampledata <- readRDS(file.path(system.file("extdata", package="ConAn"), "sampledata.rds"))
  edat <- t(exprs(sampledata))
  testfunc <- atanh_lower_tri_cor(edat)
  expect_is(testfunc, "numeric")
  expect_equal(length(testfunc), ((ncol(t(exprs(sampledata)))^2) - ncol(t(exprs(sampledata))))/2)
})

test_that("connectivity:bg_corrected_atanh_lower_tri_cor works", {
  sampledata <- readRDS(file.path(system.file("extdata", package="ConAn"), "sampledata.rds"))
  edat <- t(exprs(sampledata))
  testfunc <- bg_corrected_atanh_lower_tri_cor(edat, 1)
  expect_is(testfunc, "numeric")
  expect_equal(length(testfunc), ((ncol(t(exprs(sampledata)))^2) - ncol(t(exprs(sampledata))))/2)
  expect_equal(testfunc, atanh_lower_tri_cor(edat))
  testfunc <- bg_corrected_atanh_lower_tri_cor(edat, 0)
  expect_equal(sum(testfunc!=0), 0)
})

test_that("connectivity:atanh_lower_tri_erase_mods_cor works", {
  sampledata <- readRDS(file.path(system.file("extdata", package="ConAn"), "sampledata.rds"))
  mod_list <- readRDS(file.path(system.file("extdata", package="ConAn"), "samplemods.rds"))
  edat <- t(exprs(sampledata))
  testfunc <- atanh_lower_tri_erase_mods_cor(edat, mod_list)
  expect_is(testfunc, "numeric")
})
