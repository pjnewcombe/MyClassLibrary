package Methods;

/**
 * A collection of generic methods associated with regression modeling.
 * 
 * @author Paul J Newcombe
 */

import Objects.Arguments;
import Objects.Data;
import Objects.IterationValues;
import Objects.LikelihoodTypes;
import java.io.BufferedWriter;
import java.io.IOException;

public class RegressionMethods {


    /***
     * Final values of the iteration (stored in RegressionIterationValues object 'curr')
     * are written to the results text file. The first values written indicate
     * whether or not each parameter was present for the iteration. Note that
     * only marker OR parameters may of been absent for the iteration.
     *
     * @param curr
     *          RegressionIterationValues Object: Current likelihoodFamily and parameter values
     * @param buffer
     *          BufferedWriter Object: Corresponds to results text file
     */
    public static void writeToResultsFile (
            Arguments arguments,
            Data data,
            IterationValues curr,
            BufferedWriter buffer
            )
            throws IOException {
            // Write values
            if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
                buffer.write(curr.logWeibullScale +" ");   // Weibull k                                
            } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()) {
                buffer.write(curr.dirichletConcentration +" ");   // Weibull k                                
                buffer.write(curr.auc +" ");   // Weibull k                                
            } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
                buffer.write(curr.auc +" ");   // Weibull k                                
            } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                    data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
                buffer.write(curr.logGaussianResidual +" ");   // Gaussian residual                
            }
            buffer.write(curr.alpha+" ");
            for (int v=0; v<data.totalNumberOfCovariates; v++) {
                //System.out.println("TEST");
                buffer.write(curr.betas.get(v, 0) +" ");
            }
            for (int c=0; c<data.numberOfHierarchicalCovariatePriorPartitions; c++) {
                buffer.write(curr.logBetaPriorSds[c] +" ");   // between study var                
            }
            buffer.write(curr.logLikelihood+" ");   // log-Likelihood

            buffer.newLine();
    }


}
