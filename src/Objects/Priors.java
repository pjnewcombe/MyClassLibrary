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
    public UniformDist betaPrecisionUniformPrior;
    public GammaDist betaPrecisionGammaPrior;
    public GammaDist betweenClusterPrecisionGammaPrior;
    public GammaDist weibullScalePrior;
    public GammaDist gaussianResidualPrecisionPrior;

    // Contructor method has name as class
    // NewLikeData data1 = new NewLikeData(m, sg, sa)
    public Priors(Arguments arguments) {
        betaPrecisionGammaPrior = new GammaDist(arguments.betaPrecisionGammaPriorHyperparameter1,
                arguments.betaPrecisionGammaPriorHyperparameter2);
        betaPrecisionUniformPrior = new UniformDist(0,2);
        betweenClusterPrecisionUniformPrior = new UniformDist(arguments.betweenClusterPrecisionUniformPriorHyperparameter1,
                arguments.betweenClusterPrecisionUniformPriorHyperparameter2);
        betweenClusterPrecisionGammaPrior = new GammaDist(
                arguments.betweenClusterPrecisionGammaPriorHyperparameter1,
                arguments.betweenClusterPrecisionGammaPriorHyperparameter2);
        weibullScalePrior = new GammaDist(arguments.weibullScaleGammaPriorHyperparameter1,
                arguments.weibullScaleGammaPriorHyperparameter2);
        gaussianResidualPrecisionPrior = new GammaDist(0.001, 0.001);
    }

}
