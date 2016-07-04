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
    public UniformDist betaSigmaUniformPrior;
    public UniformDist gaussianResidualUniformPrior;
    public GammaDist betaPrecisionGammaPrior;
    public GammaDist betweenClusterPrecisionGammaPrior;
    public GammaDist weibullScalePrior;
    public GammaDist gaussianPrecisionGammaPrior;
    public GammaDist dirichletConcentrationGammaPrior;
    public UniformDist dirichletConcentrationUniformPrior;

    // Contructor method has name as class
    // NewLikeData data1 = new NewLikeData(m, sg, sa)
    public Priors(Arguments arguments, Data data) {
        if (data.whichLikelihoodType==LikelihoodTypes.JAM.ordinal()) {
            /***
             * For JAM this is used as the prior on the variable selection
             * coefficient, tau.
             * 
             * tau ~ InvGamma(1/2, n/2) [From P. 586 of Bottolo and Richardson, section 2.2]
             * 
             * so: (1/tau) ~ Gamma(1/2, 2/n)
             * under Java Gamma parameterisation (parameterisation 2) this becomes:
             * (1/tau) ~ Gamma(1/2, n/2)
             * so will have 1/betaPriorSd^2 ~ Gamma(1/2, n/2)
             */
            betaPrecisionGammaPrior = new GammaDist(0.5,(double) (data.tau/2));            
        } else {
            betaPrecisionGammaPrior = new GammaDist(arguments.betaPrecisionGammaPriorHyperparameter1,
                    arguments.betaPrecisionGammaPriorHyperparameter2);
    // SMMR Paper        betaPrecisionUniformPrior = new UniformDist(0,2);
        }
        dirichletConcentrationGammaPrior = new GammaDist(
                arguments.dirichletConcentrationGammaPriorHyperparameter1,
                arguments.dirichletConcentrationGammaPriorHyperparameter2);
        betaSigmaUniformPrior = new UniformDist(
                arguments.betaSigmaUniformPriorHyperparameter1,
                arguments.betaSigmaUniformPriorHyperparameter2);
        gaussianResidualUniformPrior = new UniformDist(
                arguments.gaussianResidualUniformPriorHyperparameter1,
                arguments.gaussianResidualUniformPriorHyperparameter2);
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
