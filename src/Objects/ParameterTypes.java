package Objects;

/**
 * Enumerates the different parameter types.
 * 
 * @author Paul J Newcombe
 */
public enum ParameterTypes {
    // PropSds/WhichUpdated Ref:
    // 0: Alpha
    // 1: Betas (when remain in model)
    // 2: Study Intercepts
    // 3: Between cluster SD
    // 4: Beta prior SD
    // 5: Weibull k
    // 6: Gaussian residual
    // 7: Beta add
    // 8: Beta swap
    ALPHA, BETAS, CLUSTER_INTERCEPTS, BETWEEN_CLUSTER_SD, BETA_PRIOR_SD,
    WEIBULL_SCALE, GAUSSIAN_RESIDUAL, BETA_ADD, BETA_SWAP;
    
}
