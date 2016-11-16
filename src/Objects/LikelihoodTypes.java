package Objects;

/**
 * Enumerates the different parameter types.
 * 
 * @author Paul J Newcombe
 */
public enum LikelihoodTypes {
    // PropSds/WhichUpdated Ref:
    // 0: Logistic
    // 1: Weibull
    // 2: Gaussian
    // 3: Marginal statistics from a Gaussian
    LOGISTIC, CLOGLOG, WEIBULL, GAUSSIAN, GAUSSIAN_CONJ,
    JAM_MCMC, JAM, JAMv2, ROCAUC, ROCAUC_ANCHOR, COX,
    CASECOHORT_PRENTICE, CASECOHORT_BARLOW;
    
}
