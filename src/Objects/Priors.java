package Objects;

import umontreal.iro.lecuyer.probdist.GammaDist;
import umontreal.iro.lecuyer.probdist.UniformDist;

/**
 * Contains all prior distributions in use.
 * 
 * @author Paul J Newcombe
 */
public class Priors {
    public UniformDist betweenClusterPrecisionUniformPrior;
    public UniformDist[] hierarchicalCovariatePriorSd_UniformPrior;
    public GammaDist[] hierarchicalCovariatePriorPrecision_GammaPrior;
    public UniformDist gaussianResidualUniformPrior;
    public GammaDist betaPrecisionConjugateGammaPrior;
    public GammaDist betweenClusterPrecisionGammaPrior;
    public GammaDist weibullScalePrior;
    public GammaDist gaussianPrecisionGammaPrior;
    public GammaDist dirichletConcentrationGammaPrior;
    public UniformDist dirichletConcentrationUniformPrior;

    // Contructor method has name as class
    // NewLikeData data1 = new NewLikeData(m, sg, sa)
    public Priors(Arguments arguments, Data data) {
        if (data.whichLikelihoodType==LikelihoodTypes.JAM.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            /***
             * For conjugate Gaussian models this is the prior on the variable selection
             * coefficient, tau.
             * 
             * tau ~ InvGamma(1/2, n/2) [From P. 586 of Bottolo and Richardson, section 2.2]
             * 
             * so: (1/tau) ~ Gamma(1/2, 2/n)
             * under Java Gamma parameterisation (parameterisation 2) this becomes:
             * (1/tau) ~ Gamma(1/2, n/2)
             * so will have 1/betaPriorSd^2 ~ Gamma(1/2, n/2)
             */
            betaPrecisionConjugateGammaPrior = new GammaDist(0.5,(double) (data.tau/2));            
        }
        
        dirichletConcentrationGammaPrior = new GammaDist(
                arguments.dirichletConcentrationGammaPriorHyperparameter1,
                arguments.dirichletConcentrationGammaPriorHyperparameter2);
        
        /**
         * Hierarchical covariate priors
         */
        hierarchicalCovariatePriorSd_UniformPrior = new UniformDist[data.numberOfHierarchicalCovariatePriorPartitions];
        hierarchicalCovariatePriorPrecision_GammaPrior = new GammaDist[data.numberOfHierarchicalCovariatePriorPartitions];
        for (int c=0; c<data.numberOfHierarchicalCovariatePriorPartitions; c++) {
            if (data.hierarchicalCovariatePriorPartitionFamilies[c] == HierarchicalCovariatePriorTypes.UNIFORM.ordinal()) {
                hierarchicalCovariatePriorSd_UniformPrior[c] = new UniformDist( // SMMR Paper UniformDist(0,2);
                        data.hierarchicalCovariatePriorSd_UniformHyperparameter1[c],
                        data.hierarchicalCovariatePriorSd_UniformHyperparameter2[c]);
            } else if (data.hierarchicalCovariatePriorPartitionFamilies[c] == HierarchicalCovariatePriorTypes.GAMMA.ordinal()) {
                hierarchicalCovariatePriorPrecision_GammaPrior[c] = new GammaDist(
                        data.hierarchicalCovariatePriorPrecision_GammaHyperparameter1[c],
                        data.hierarchicalCovariatePriorPrecision_GammaHyperparameter2[c]);                
            }
        }
        
        /**
         * Gaussian residual prior
         */
        gaussianResidualUniformPrior = new UniformDist(
                arguments.gaussianResidualUniformPriorHyperparameter1,
                arguments.gaussianResidualUniformPriorHyperparameter2);
        
        /**
         * Weibull scale prior
         */
        weibullScalePrior = new GammaDist(arguments.weibullScaleGammaPriorHyperparameter1,
                arguments.weibullScaleGammaPriorHyperparameter2);
        /**
         * Variance ~ InvGamma iff Precision ~ Gamma
         * Variance ~ InvGamma(a,b) 
         * -> Precision ~ Gamma(a,1/b)
         * -> Precision ~ Gamma(a,b) [Java parameterisation is upside down]
         */
        gaussianPrecisionGammaPrior = new GammaDist(
                arguments.gaussianResidualVarianceInvGammaPrior_a,
                arguments.gaussianResidualVarianceInvGammaPrior_b);
    }

}
