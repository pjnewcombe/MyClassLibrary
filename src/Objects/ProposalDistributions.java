package Objects;

/**
 * Contains information relating to the proposal distribution SDs, and their
 * adaption.
 * 
 * @author Paul J Newcombe
 */
public class ProposalDistributions {
    // static means each instaniation will have fixed param values
    private int adaptionInterval;
    private int adaptionLength;
    private int[] adapting;
    public int nSds;  // Number of Sds total
    public int[] acceptanceRateNumerator;
    public int[] acceptanceRateDenominator;
    public double[] acceptanceRate;
    public double[] sds;
   
    // Contructor method has name as class
    // NewLikeData data1 = new NewLikeData(m, sg, sa)
    public ProposalDistributions(Arguments arguments, Data data) {
        nSds = 8;
        acceptanceRateNumerator = new int[nSds];
        acceptanceRateDenominator = new int[nSds];
        acceptanceRate = new double[nSds];
        sds = new double[nSds];
        adapting = new int[nSds];
        sds[ParameterTypes.ALPHA.ordinal()] = arguments.proposalDistributionSdForAlpha;
        adapting[ParameterTypes.ALPHA.ordinal()] = 1;
        sds[ParameterTypes.BETAS.ordinal()] = arguments.proposalDistributionSdForBetas;
        adapting[ParameterTypes.BETAS.ordinal()] = 1;
        if (data.numberOfClusters>0) {
            sds[ParameterTypes.CLUSTER_INTERCEPTS.ordinal()] = arguments.proposalDistributionSdForClusterIntercepts;
            sds[ParameterTypes.BETWEEN_CLUSTER_SD.ordinal()] = arguments.proposalDistributionSdForLogBetweenClusterSd;
            adapting[ParameterTypes.CLUSTER_INTERCEPTS.ordinal()] = 1;
            adapting[ParameterTypes.BETWEEN_CLUSTER_SD.ordinal()] = 1;            
        }
        if (data.numberOfModelSpacePartitions>0) {
            sds[ParameterTypes.BETA_PRIOR_SD.ordinal()] = arguments.proposalDistributionSdForBetaPriorSd;
            adapting[ParameterTypes.BETA_PRIOR_SD.ordinal()] = 1;            
        }
        if (data.survivalAnalysis==1) {
            sds[ParameterTypes.WEIBULL_SCALE.ordinal()] = arguments.proposalDistributionSdForLogWeibullScale;
            adapting[ParameterTypes.WEIBULL_SCALE.ordinal()] = 1;
        }
        sds[ParameterTypes.BETA_ADD.ordinal()] = arguments.proposalDistributionSdForAddingBeta;
        sds[ParameterTypes.BETA_SWAP.ordinal()] = arguments.proposalDistributionSdForSwappedInBeta;
        adaptionInterval = arguments.adaptionBinSize;
        adaptionLength = arguments.adaptionLength;
    }
    
    public void adapt(Data data, IterationValues Its, int i) {
        if (i<adaptionLength) {
            for (int j=0; j<nSds; j++) {
                if (adapting[j]==1) {
                    // NOTE: Order of following steps must remain for same results
                    // Reset numerator and denominator at end of bin
                    if (acceptanceRateDenominator[j]==adaptionInterval) {
                        acceptanceRateNumerator[j]=0;
                        acceptanceRateDenominator[j]=0;
                    }
                    // Add 1 to denominator of parameter attempted to update
                    if (Its.whichParameterTypeUpdated ==j) {
                        acceptanceRateDenominator[j]=acceptanceRateDenominator[j]+1;
                    }
                    // Perform adaption
                    if (acceptanceRateDenominator[j]==adaptionInterval) { 
                        acceptanceRate[j] = (double)acceptanceRateNumerator[j]
                                /(double)acceptanceRateDenominator[j];
                        if (acceptanceRate[j]>0.42) {
                            sds[j]=sds[j]*1.01;
                        }
                        if (acceptanceRate[j]<0.42) {
                            sds[j]=sds[j]*0.99;
                        }
                    }
                    // Add 1 to numerator of parameter updated if proposalAccepted
                    if (Its.proposalAccepted==1&&Its.whichParameterTypeUpdated==j) {
                        acceptanceRateNumerator[j]=acceptanceRateNumerator[j]+1;
                    }
                }
            }            
        } else if (i==adaptionLength) {
            System.out.println("End of proposal SD adaption...");
            for (int v=0; v<nSds; v++) {
                if (adapting[v]==1) {
                    System.out.println(ParameterTypes.values()[v]+" acceptance rate: "+acceptanceRate[v]);
                    System.out.println(ParameterTypes.values()[v]+" final proposal SD: "+sds[v]);                            
                }
            }
            // SET ADDITION AND SWAP PROPOSAL SDS TO FRACTION OF LOGOR SD
            sds[ParameterTypes.BETA_ADD.ordinal()] = (double)(sds[ParameterTypes.BETAS.ordinal()]);
            // swap SD should be smaller than addition SD
            sds[ParameterTypes.BETA_SWAP.ordinal()] = (double)(sds[ParameterTypes.BETAS.ordinal()]/(data.totalNumberOfCovariates-data.numberOfCovariatesToFixInModel));            
        }
    }
}
