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
    /**
     * Total number of proposal distributions
     */
    public int numberOfProposalDistributions;  // Number of Sds total
    /**
     * Numerators for the acceptance rate of each proposal distribution
     */
    public int[] acceptanceRateNumerators;
    /**
     * Denominators for the acceptance rate of each proposal distribution
     */
    public int[] acceptanceRateDenominators;
    /**
     * Proposal distribution acceptance rates
     */
    public double[] acceptanceRates;
    /**
     * Proposal distribution SDs
     */
    public double[] proposalDistributionSds;
   
    // Contructor method has name as class
    // NewLikeData data1 = new NewLikeData(m, sg, sa)
    public ProposalDistributions(Arguments arguments, Data data) {
        numberOfProposalDistributions = 8;
        acceptanceRateNumerators = new int[numberOfProposalDistributions];
        acceptanceRateDenominators = new int[numberOfProposalDistributions];
        acceptanceRates = new double[numberOfProposalDistributions];
        proposalDistributionSds = new double[numberOfProposalDistributions];
        adapting = new int[numberOfProposalDistributions];
        proposalDistributionSds[ParameterTypes.ALPHA.ordinal()] = arguments.proposalDistributionSdForAlpha;
        adapting[ParameterTypes.ALPHA.ordinal()] = 1;
        proposalDistributionSds[ParameterTypes.BETAS.ordinal()] = arguments.proposalDistributionSdForBetas;
        adapting[ParameterTypes.BETAS.ordinal()] = 1;
        if (data.numberOfClusters>0) {
            proposalDistributionSds[ParameterTypes.CLUSTER_INTERCEPTS.ordinal()] = arguments.proposalDistributionSdForClusterIntercepts;
            proposalDistributionSds[ParameterTypes.BETWEEN_CLUSTER_SD.ordinal()] = arguments.proposalDistributionSdForLogBetweenClusterSd;
            adapting[ParameterTypes.CLUSTER_INTERCEPTS.ordinal()] = 1;
            adapting[ParameterTypes.BETWEEN_CLUSTER_SD.ordinal()] = 1;            
        }
        if (data.numberOfUnknownBetaPriors>0) {
            proposalDistributionSds[ParameterTypes.BETA_PRIOR_SD.ordinal()] = arguments.proposalDistributionSdForBetaPriorSd;
            adapting[ParameterTypes.BETA_PRIOR_SD.ordinal()] = 1;            
        }
        if (data.survivalAnalysis==1) {
            proposalDistributionSds[ParameterTypes.WEIBULL_SCALE.ordinal()] = arguments.proposalDistributionSdForLogWeibullScale;
            adapting[ParameterTypes.WEIBULL_SCALE.ordinal()] = 1;
        }
        proposalDistributionSds[ParameterTypes.BETA_ADD.ordinal()] = arguments.proposalDistributionSdForAddingBeta;
        proposalDistributionSds[ParameterTypes.BETA_SWAP.ordinal()] = arguments.proposalDistributionSdForSwappedInBeta;
        adaptionInterval = arguments.adaptionBinSize;
        adaptionLength = arguments.adaptionLength;
    }
    
    public void adapt(Data data, IterationValues Its, int i) {
        if (i<adaptionLength) {
            for (int j=0; j<numberOfProposalDistributions; j++) {
                if (adapting[j]==1) {
                    // NOTE: Order of following steps must remain for same results
                    // Reset numerator and denominator at end of bin
                    if (acceptanceRateDenominators[j]==adaptionInterval) {
                        acceptanceRateNumerators[j]=0;
                        acceptanceRateDenominators[j]=0;
                    }
                    // Add 1 to denominator of parameter attempted to update
                    if (Its.whichParameterTypeUpdated ==j) {
                        acceptanceRateDenominators[j]=acceptanceRateDenominators[j]+1;
                    }
                    // Perform adaption
                    if (acceptanceRateDenominators[j]==adaptionInterval) { 
                        acceptanceRates[j] = (double)acceptanceRateNumerators[j]
                                /(double)acceptanceRateDenominators[j];
                        if (acceptanceRates[j]>0.42) {
                            proposalDistributionSds[j]=proposalDistributionSds[j]*1.01;
                        }
                        if (acceptanceRates[j]<0.42) {
                            proposalDistributionSds[j]=proposalDistributionSds[j]*0.99;
                        }
                    }
                    // Add 1 to numerator of parameter updated if proposalAccepted
                    if (Its.proposalAccepted==1&&Its.whichParameterTypeUpdated==j) {
                        acceptanceRateNumerators[j]=acceptanceRateNumerators[j]+1;
                    }
                }
            }            
        } else if (i==adaptionLength) {
            System.out.println("Proposal distribution adaption phase complete;");
            for (int v=0; v<numberOfProposalDistributions; v++) {
                if (adapting[v]==1) {
                    System.out.println(ParameterTypes.values()[v]+" final acceptance rate: "+acceptanceRates[v]+" value: "+proposalDistributionSds[v]);
                    //System.out.println(ParameterTypes.values()[v]+" final proposal SD: "+proposalDistributionSds[v]);                            
                }
            }
            // SET ADDITION AND SWAP PROPOSAL SDS TO FRACTION OF LOGOR SD
            proposalDistributionSds[ParameterTypes.BETA_ADD.ordinal()] = (double)(proposalDistributionSds[ParameterTypes.BETAS.ordinal()]);
            // swap SD should be smaller than addition SD
            proposalDistributionSds[ParameterTypes.BETA_SWAP.ordinal()] = (double)(proposalDistributionSds[ParameterTypes.BETAS.ordinal()]/(data.totalNumberOfCovariates-data.numberOfCovariatesToFixInModel));            
        }
    }
}
