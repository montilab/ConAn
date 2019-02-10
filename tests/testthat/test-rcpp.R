# Randomly generated data
ncols = 1000
nrows = 500 
set.seed(1)
edat <- matrix(runif(ncols*nrows), ncol=ncols)

genes <- c(1:ncols)
modlist <- list("a"=c(854,709,829,592,356,13,266,144,497,174),
                "b"=c(283,737,937,551,264,842,201,648),
                "c"=c(407,92,41,983,133,142,686,939),
                "d"=c(584,482,368,241,764,972),
                "e"=c(996,918,620,275))

matz <- modlist.to.matzindex(modlist, genes)

# Tests
test_that("Pearson correlation is working", {

    # Pearson correlation matrix
    cmat.1 <- stats::cor(edat, method="pearson")
    cmat.2 <- ConAn::pcor(edat)
    
    expect_equal(cmat.1, cmat.2)
})

test_that("Erasure of matrix values is working", {

    m <- t(matrix(c(2,4,4,5,2,4,5,6,8,5,4,3,5,6,7,5,4,4,3,2,2), ncol=3))
    cmat <- ConAn::pcor(m)
    
    genes <- c(1:7)
    modlist <- list("a"=c(2,3), "b"=c(5,6,7))
    matz <- modlist.to.matzindex(modlist, genes)
    
    # Erase values from matrix
    ecmat <- ConAn::erase_vals(cmat, matz)

    expect_true(table(is.na(ecmat[modlist$a, modlist$a])) == 4)
    expect_true(table(is.na(ecmat[modlist$b, modlist$b])) == 9)
})

test_that("Connectivity vector is working", {
    
    cv.1 <- edat %>%
            stats::cor() %>%
            get_lower_tri(diag=FALSE) %>%
            atanh()
    
    cv.2 <- ConAn::cv(edat)
   
    expect_equal(cv.1, cv.2)
})

test_that("Background connectivity vector is working", {
    
    bgcv.1 <- edat %>%
              stats::cor() %>%
              erase_mods(mod_list=modlist) %>%
              get_lower_tri(diag=FALSE) %>%
              remove_na() %>%
              atanh()
    
    bgcv.2 <- ConAn::bgcv(edat, matz)
   
    expect_equal(bgcv.1, bgcv.2)
})

test_that("Mean connectivity is working", {
    
    mc.1 <- edat %>%
            stats::cor() %>%
            get_lower_tri(diag=FALSE) %>%
            atanh() %>%
            mean()
    
    mc.2 <- ConAn::mc(edat)
   
    expect_equal(mc.1, mc.2)
})

test_that("Background mean connectivity is working", {
    
    bgmc.1 <- edat %>%
              stats::cor() %>%
              erase_mods(mod_list=modlist) %>%
              get_lower_tri(diag=FALSE) %>%
              remove_na() %>%
              atanh() %>%
              mean()
    
    bgmc.2 <- ConAn::bgmc(edat, matz)
   
    expect_equal(bgmc.1, bgmc.2)
})

test_that("Modular differential connectivity is", {
    
    edat.r <- matrix(runif(ncols*nrows), ncol=ncols) 
    edat.t <- matrix(runif(ncols*nrows), ncol=ncols)     
    
    bg.r <- 0.035
    bg.t <- 0.025
    
    cv.r <- edat.r %>%
            stats::cor() %>%
            get_lower_tri(diag=FALSE) %>%
            atanh()    
    
    cv.t <- edat.t %>%
            stats::cor() %>%
            get_lower_tri(diag=FALSE) %>%
            atanh()
    
    mc.r <- mean( tanh(cv.r - bg.r)^2 )
    mc.t <- mean( tanh(cv.t - bg.t)^2 )
    
    # Fraction
    mdc.1 <- mc.t/mc.r
    
    mdc.2 <- .Call("S_mdc", R_xr=edat.r, R_xt=edat.t, R_bgr=bg.r, R_bgt=bg.t, PACKAGE="ConAn")

    expect_equal(mdc.1, mdc.2)
})
