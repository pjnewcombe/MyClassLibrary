package Methods;

/**
 * A collection of generic methods associated with regression modeling.
 * 
 * @author Paul J Newcombe
 */

import Objects.Arguments;
import Objects.Data;
import Objects.IterationValues;
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
            if (data.survivalAnalysis==1) {
                buffer.write(curr.logWeibullScale +" ");   // between study var                                
            }
            buffer.write(curr.alpha+" ");
            for (int v=0; v<data.totalNumberOfCovariates; v++) {
                buffer.write(curr.betas.get(v, 0) +" ");
            }
            for (int c=0; c<data.numberOfModelSpacePartitions; c++) {
                buffer.write(curr.logBetaPriorSd[c] +" ");   // between study var                
            }
            if (data.numberOfClusters > 0) {
                buffer.write(curr.logBetweenClusterSd +" ");   // between study var
                if (arguments.recordClusterIntercepts == 1) {
                    for (int r=0; r<data.numberOfClusters; r++) {
                        buffer.write(curr.clusterIntercepts[r]+" ");
                    }
                }
            }
            buffer.write(curr.logLikelihood+" ");   // log-Likelihood

            buffer.newLine();
    }


}
