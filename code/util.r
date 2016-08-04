estimatePopulationVariance = function(M, var, Mbar, maxit = 1e4) {
    sigma2hati = (M - Mbar)^2 - var
    sigma2mu0 = 0
    sigma2mu = mean((M - Mbar)^2)
    t = 0
    while(abs(sigma2mu - sigma2mu0) > 1e-5 & t < maxit) {
        t = t + 1
        sigma2mu0 = sigma2mu
        weight = (2*(var + sigma2mu))^{-1}
        sigma2mu = sum(weight*sigma2hati)/sum(weight)
    }
    sigma2mu
}
