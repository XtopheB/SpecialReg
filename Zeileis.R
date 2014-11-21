set.seed(1090)
dat <- as.data.frame(matrix(round(runif(21), digits = 2), ncol = 7))
colnames(dat) <- c("y1", "y2", "y3", "x1", "x2", "x3", "x4")
for(i in c(2, 6:7)) dat[[i]] <- factor(dat[[i]] < 0.5,
                                          + labels = c("a", "b"))
dat$y2[1] <- NA
dat

is.data.frame(dat)

ivcoef <- function(formula, data, subset, na.action, ...)
  {
     mf <- match.call(expand.dots = FALSE)
     m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
     mf <- mf[c(1, m)]
    
    f <- Formula(formula)
     mf[[1]] <- as.name("model.frame")
     mf$formula <- f
     mf <- eval(mf, parent.frame())
    
     y <- model.response(mf)
     x <- model.matrix(f, data = mf, rhs = 1)
     z <- model.matrix(f, data = mf, rhs = 2)
     xz <- as.matrix(lm.fit(z, x)$fitted.values)
     lm.fit(xz, y)$coefficients
     }

ivcoef(log(y1) ~ x1 | x2, data = dat)
ivcoef(i_tap~prix_2008 + i_under18  |itap_2008, data = data.work)
specialreg3(i_tap~prix_2008 + i_under18  |itap_2008, data = data.work)
