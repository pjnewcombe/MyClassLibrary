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
    public GammaDist gaussianResidualPrecisionPrior;
    public GammaDist scaledChiSquareGaussianResidualPrecisionPrior;

    // Contructor method has name as class
    // NewLikeData data1 = new NewLikeData(m, sg, sa)
    public Priors(Arguments arguments, Data data) {
        betaPrecisionGammaPrior = new GammaDist(arguments.betaPrecisionGammaPriorHyperparameter1,
                arguments.betaPrecisionGammaPriorHyperparameter2);
// SMMR Paper        betaPrecisionUniformPrior = new UniformDist(0,2);
        betaSigmaUniformPrior = new UniformDist(
                arguments.betaSigmaUniformPriorHyperparameter1,
                arguments.betaSigmaUniformPriorHyperparameter2);
        gaussianResidualUniformPrior = new UniformDist(
                arguments.gaussianResidualUniformPriorHyperparameter1,
                arguments.gaussianResidualUniformPriorHyperparameter2);
        betweenClusterPrecisionUniformPrior = new UniformDist(
                arguments.betweenClusterPrecisionUniformPriorHyperparameter1,
                arguments.betweenClusterPrecisionUniformPriorHyperparameter2);
        betweenClusterPrecisionGammaPrior = new GammaDist(
                arguments.betweenClusterPrecisionGammaPriorHyperparameter1,
                arguments.betweenClusterPrecisionGammaPriorHyperparameter2);
        weibullScalePrior = new GammaDist(arguments.weibullScaleGammaPriorHyperparameter1,
                arguments.weibullScaleGammaPriorHyperparameter2);
        gaussianResidualPrecisionPrior = new GammaDist(
                arguments.gaussianResidualPrecisionGammaPriorHyperparameter1,
                arguments.gaussianResidualPrecisionGammaPriorHyperparameter2);
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
            scaledChiSquareGaussianResidualPrecisionPrior = new GammaDist(
                    (data.gaussianResidualVarianceEstimateN/2),
                    data.gaussianResidualVarianceEstimateN/
                            (2*data.gaussianResidualVarianceEstimate)
            );            
        }
    }

}
