package Objects;

import Methods.GeneralMaths;
import Jama.Matrix;
import Methods.GeneralMethods;
import java.util.Random;

/**
 * This class contains a model and parameter values, in addition to it's
 * likelihood and prior support.
 * 
 * @author Paul J Newcombe
 */
public class IterationValues {
    /**
     * 
     * Model information. Below are key relevant bits of information on the
     * model setup.
     * 
     */
    
    /**
     * The total number of covariates.
     */
    private final int totalNumberOfCovariates;
    /**
     * The total number of covariates fixed in the model.
     */
    private final int numberOfCovariatesToFixInModel;
    /**
     * The total number of clusters (in the case of random intercepts).
     */
    private final int numberOfClusters;
    
    /**
     * 
     * Modeling parameters. This includes the model parameters, information
     * describing the selected model, and information describing the beta
     * priors.
     * 
     */
    
    /**
     * Model intercept.
     */
    public double alpha;
    /**
     * Model covariate parameters (eg log-ORs for logistic regression).
     */
    public Matrix betas;
    /**
     * Vector of clusters specific random intercepts.
     */
    public double[] clusterIntercepts;
    /**
     * Between clusters standard deviation (logged).
     */
    public double logBetweenClusterSd;
    /**
     * Weibull scale parameter.
     */
    public double weibullScale;
    /**
     * Log-Weibull scale parameter.
     */
    public double logWeibullScale;
    /**
     * Gaussian residual parameter.
     */
    public double gaussianResidual;
    /**
     * Log-Gaussian residual parameter.
     */
    public double logGaussianResidual;
    /**
     * Vector indicating which covariates are included in the model
     * (length = total number of covariates being searched over).
     */
    public int[] model;
    /**
     * Total number of covariates included in the model.
     */
    public int modelDimension;
    /**
     * Total number of covariates included in each component of the model space
     * (length = number of model space components).
     */
    public int[] modelSpacePartitionDimensions;
    /**
     * Unknown common prior SD across the betas (length = number of model space
     * components).
     */
    public double[] betaPriorSds;
    /**
     * Log-unknown common prior SD across the betas (length = number of model space
     * components).
     */
    public double[] logBetaPriorSd;    
    
    
    /**
     * 
     * 
     * Below are quantities which are calculated as part of the likelihood
     * calculations.
     * 
     */
    
    /**
     * Product of the covariate matrix and coefficient vector.
     */
    public Matrix Xbeta;        
    /**
     * Product of the covariate matrix and coefficient vector, by block (for
     * GuassianMarg).
     */
    public Matrix[] XbetaBlocks;        
    /**
     * Contains the random intercepts.
     */
    public Matrix alphaJama;
    /**
     * GaussianMarg: Likelihood score term - saved so that minor changes can be
     * made
     **/
    public double likelihoodScoreTerm;
    /**
     * GaussianMarg: Block -specific likelihood score term - saved so that minor
     * changes can be made
     **/
    public double[] likelihoodScoreTermsBlocks;
    /**
     * Y - X*Beta.
     */ 
    private Matrix yMinXbeta;
    /**
     * Prior support (logged).
     */
    public double logPrior;
    /**
     * Log-Likelihood.
     */
    public double logLikelihood;
    /**
     * Acceptance probability for this model and parameter values in comparison
     * to the current state.
     */
    public double acceptanceProbability;
    
    /**
     * 
     * Book-keeping. The following variables are used internally for various
     * book-keeping tasks. For example some of these are required to calculate
     * the reversible jump move probabilities.
     * 
     */
    
    /**
     * Which blocks had betas updated for the proposal.
     */
    private boolean[] whichBlocksUpdated;
    /**
     * Flags whether to accept the model and parameters stored in this object
     * as the new model and parameters.
     */
    public int proposalAccepted;
    /**
     * Indicates which parameter was updated according to the dictionary
     * (@link Objects.ParameterTypes).
     */
    public int whichParameterTypeUpdated;
    /**
     * Indicates which model search move is performed at this iteration.
     * 0: Removal, 1: Addition, 2: Swap, 3: Null.
     */
    public int whichMove;
    /**
     * In the case of a null move, and beta update, indicates which beta was
     * updated.
     */
    public int whichBetaUpdated;
    /**
     * In the case of a null move, and beta update, indicates which block the
     * updated beta was in.
     */
    public int whichBetaUpdatedBlock;
    /**
     * In the case of a removal move, indicates which beta was removed (for
     * GuassianMarg).
     */
    public int whichBetaRemoved;
    /**
     * In the case of a removal move, indicates which block the removed beta
     * was in (for GuassianMarg).
     */
    public int whichBetaRemovedBlock;
    /**
     * In the case of an addition move, indicates which beta was added.
     */
    public int whichBetaAdded;
    /**
     * In the case of an addition move, indicates which block the beta was
     * added to (for GaussianMarg).
     */
    public int whichBetaAddedBlock;

    /**
     * Constructor; used to initiate an analysis with a model (usually the null)
     * and some initial parameter values
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param priors {@link Objects.Priors} class object, containing information
     * on the prior distributions
     */
    public IterationValues(
            Arguments arguments,
            Data data,
            Priors priors) {
        totalNumberOfCovariates = data.totalNumberOfCovariates;
        numberOfCovariatesToFixInModel = data.numberOfCovariatesToFixInModel;
        numberOfClusters = data.numberOfClusters;
        model = new int[data.totalNumberOfCovariates]; // indicates markers present in current model
        // Parameters that are updated
        betas = new Matrix(data.totalNumberOfCovariates,1);
        clusterIntercepts = new double[data.numberOfClusters]; // log-ORs of current model
        likelihoodScoreTermsBlocks = new double[data.nBlocks];
        Xbeta = new Matrix(0,0); // Needs to be initialised, if not used, for setTo
        XbetaBlocks = new Matrix[data.nBlocks];
        whichBlocksUpdated = new boolean[data.nBlocks];
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
            alpha = 0;
        } else {
            alpha = arguments.initialAlpha;
        }
        for (int v=0; v<data.totalNumberOfCovariates; v++) {
            if ((arguments.useSaturatedInitialModel == 1)|(v < data.numberOfCovariatesToFixInModel)) {
                model[v] = 1;
            } else {
                model[v] = 0;
            }
            betas.set(v, 0, arguments.initialBetas);
        }
        modelDimension = GeneralMethods.countPresVars(data.totalNumberOfCovariates, model);
        modelSpacePartitionDimensions = GeneralMethods.countPresVarsComps(
                arguments.numberOfModelSpacePriorPartitions, data.modelSpacePartitionIndices, model);
        logBetweenClusterSd = Math.log(arguments.initialBetweenClusterSd);
        logBetaPriorSd = new double[data.numberOfUnknownBetaPriors];
        betaPriorSds = new double[data.numberOfUnknownBetaPriors];
        for (int c=0; c<data.numberOfUnknownBetaPriors; c++) {
            logBetaPriorSd[c] = Math.log(arguments.initialBetaPriorSd);
            betaPriorSds[c] = Math.exp(logBetaPriorSd[c]);
        }
        weibullScale = arguments.initialWeibullScale;
        logWeibullScale = Math.log(weibullScale);
        gaussianResidual = arguments.initialGaussianResidual;
        logGaussianResidual = Math.log(gaussianResidual);
        
        
        calcLogLike(arguments, data);
        calcLogPrior(arguments, data, priors);
    }
    
    /**
     * Sets this object equal to another object also of the class
     * {@link Objects.IterationValues}
     * 
     * @param its Another instance of {@link Objects.IterationValues} class to set all
     * parameters equal to.
     */
    public void setTo(IterationValues its) {
        // Vectors
        clusterIntercepts = its.clusterIntercepts.clone();
        model = its.model.clone();
        modelSpacePartitionDimensions = its.modelSpacePartitionDimensions.clone();
        logBetaPriorSd = its.logBetaPriorSd.clone();
        betaPriorSds = its.betaPriorSds.clone();
        likelihoodScoreTermsBlocks = its.likelihoodScoreTermsBlocks.clone();
        
        // Single values
        alpha = its.alpha;
        modelDimension = its.modelDimension;
        logBetweenClusterSd = its.logBetweenClusterSd;
        weibullScale = its.weibullScale;
        logWeibullScale = its.logWeibullScale;
        gaussianResidual = its.gaussianResidual;
        logGaussianResidual = its.logGaussianResidual;
        likelihoodScoreTerm = its.likelihoodScoreTerm;
        logLikelihood = its.logLikelihood;
        logPrior = its.logPrior;
        
        // Jama objects
        alphaJama = its.alphaJama;  // Copy does not work for some reason
        betas = its.betas.copy();
        Xbeta = its.Xbeta.copy();
        XbetaBlocks = its.XbetaBlocks.clone();
    }
    
    /**
     * Generates a new model and set of parameter values, from the current
     * state and proposal distributions, and stores in the place of the current
     * model and parameter values.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param curr {@link Objects.IterationValues} class object containing the current
     * model and parameter values
     * @param priors {@link Objects.Priors} class object, containing information
     * on the prior distributions
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param randomDraws A random number generator object
     */    
    public void update(
            Arguments arguments,
            Data data,
            IterationValues curr,
            Priors priors,
            ProposalDistributions propsds,
            Random randomDraws) {
        
        // Start by choosing model selection move
        updateModel(arguments, data, randomDraws);

        double parameterDraw = 0;
        
        // ProposalDistributions/WhichUpdated Ref:
        // 0: Alpha
        // 1: Betas (when remain in model)
        // 2: Study Intercepts
        // 3: Between clusters SD
        // 4: Beta prior SD
        // 5: Weibull k
        // 6: Beta (when added)
        // 7: Beta (when swapped in)

        if (whichMove == 1) {
            updateBetaAdd(propsds, randomDraws);
        } else if (whichMove == 2) {
            updateBetaSwap(propsds, randomDraws);
        } else if (whichMove == 3) {
            // NULL MOVE (split by whether there are random intercepts)
            if (data.numberOfClusters > 0) {
                // NULL MOVE - Region intercepts ------------------------------
                double paramTypeDraw = randomDraws.nextInt( (4+data.numberOfUnknownBetaPriors+data.nExtraParametersBeyondLinPred) );
                if (paramTypeDraw == 0) {
                    if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                        // Alpha is fixed for the approximate meta-analysis
                        // methods, hence update the gaussianResidual instead
                        updateGaussianResidual(propsds, randomDraws);                        
                    } else {
                        updateAlpha(propsds, randomDraws);
                    }
                } else if (paramTypeDraw == 1) {
                    if (modelDimension>0) {
                        updateBetas(propsds, randomDraws, data);
                    } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                        updateGaussianResidual(propsds, randomDraws);                        
                    } else {
                        updateAlpha(propsds, randomDraws);                        
                    }
                } else if (paramTypeDraw == 2) {
                    updateBetweenClusterSd(propsds, randomDraws);
                } else if (paramTypeDraw == 3) {
                    updateClusterIntercepts(propsds, randomDraws);
                } else if (paramTypeDraw == 4) {
                    if (data.numberOfUnknownBetaPriors>0) {
                        updateBetaPriorSd(propsds, randomDraws, data.numberOfUnknownBetaPriors);                        
                    } else if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
                        updateWeibullK(propsds, randomDraws);                        
                    } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                            data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                        updateGaussianResidual(propsds, randomDraws);                        
                    }
                } else if (paramTypeDraw == 5) {
                    if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
                        updateWeibullK(propsds, randomDraws);
                    } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                            data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                        updateGaussianResidual(propsds, randomDraws);                                                
                    }
                }
            } else {
            // NULL MOVE for when no clusters intercepts -----------------------
                double paramTypeDraw = randomDraws.nextInt((2+data.numberOfUnknownBetaPriors+data.nExtraParametersBeyondLinPred));
                if (paramTypeDraw == 0) {
                    if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                        updateGaussianResidual(propsds, randomDraws);                        
                    } else {
                        updateAlpha(propsds, randomDraws);
                    }                    
                } else if (paramTypeDraw == 1) {
                    if (modelDimension>0) {
                        updateBetas(propsds, randomDraws, data);
                    } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                        updateGaussianResidual(propsds, randomDraws);                        
                    } else {
                        updateAlpha(propsds, randomDraws);                        
                    }
                } else if (paramTypeDraw == 2) {
                    if (data.numberOfUnknownBetaPriors>0) {
                        updateBetaPriorSd(propsds, randomDraws, data.numberOfUnknownBetaPriors);                        
                    } else if (data.nExtraParametersBeyondLinPred==1) {
                        if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
                            updateWeibullK(propsds, randomDraws);
                        } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                            data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                            updateGaussianResidual(propsds, randomDraws);                                                
                        }
                    }
                } else if (paramTypeDraw == 3) {
                    // Weibull k or Gaussian residual
                    if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
                        updateWeibullK(propsds, randomDraws);
                    } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                            data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                        updateGaussianResidual(propsds, randomDraws);                                                
                    }
                }
            }
        }
        
        calcLogPrior(arguments, data, priors);
        if (data.numberOfClusters>0) {
            calcLogLike(arguments, data);
        } else {
            calcLogLike_IncrementFromProposal(arguments, data, curr);
        }
        
        acceptanceProbability(arguments, data, curr, propsds);
    }
    
    /**
     * Updates the intercept parameter {@link Objects.IterationValues#alpha}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     */
    public void updateAlpha(
            ProposalDistributions propsds,
            Random r) {
        whichParameterTypeUpdated = ParameterTypes.ALPHA.ordinal();
        alpha = GeneralMaths.normalDraw(
                alpha,  propsds.proposalDistributionSds[ParameterTypes.ALPHA.ordinal()], r);
    }
    
    /**
     * Updates a currently included beta parameter in 
     * {@link Objects.IterationValues#betas}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     * @param data The data object (to determine the likelihood, and for block
     * information)
     */
    public void updateBetas(
            ProposalDistributions propsds,
            Random r,
            Data data) {
        whichParameterTypeUpdated = ParameterTypes.BETAS.ordinal();
        int markUpdate = r.nextInt(modelDimension);
        int count = 0;
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
            /**
             * For GaussianMarg, start by looping over the fixed effects
             */
            for (int m=0; m<numberOfCovariatesToFixInModel; m++) {
                if (count==markUpdate) {
                    double newbeta = GeneralMaths.normalDraw(
                                betas.get(m, 0), 
                                propsds.proposalDistributionSds[ParameterTypes.BETAS.ordinal()], r);
                    betas.set(m, 0, newbeta);
                    whichBetaUpdated = m;
                }
                count++;
            }
            /**
             * For GaussianMarg, loop through the blocks
             */
            for (int b=0; b<data.nBlocks; b++) {
                for (int m=data.blockIndices[b]; m<data.blockIndices[b+1]; m++) {
                    if (model[m]==1) {
                        if (count==markUpdate) {
                            double newbeta = GeneralMaths.normalDraw(
                                        betas.get(m, 0), 
                                        propsds.proposalDistributionSds[ParameterTypes.BETAS.ordinal()], r);
                            betas.set(m, 0, newbeta);
                            whichBetaUpdated = m;
                            whichBetaUpdatedBlock = b;
                        }
                        count++;
                    }                        
                }
            }                
        } else {
            for (int m=0; m<totalNumberOfCovariates; m++) {
                if (model[m] == 1) {
                    if (count==markUpdate) {
                        double newbeta = GeneralMaths.normalDraw(
                                    betas.get(m, 0), 
                                    propsds.proposalDistributionSds[ParameterTypes.BETAS.ordinal()], r);
                        betas.set(m, 0, newbeta);
                        whichBetaUpdated = m;
                    }
                    count++;
                }
            }
        }
        
    }
    
    /**
     * Updates the beta parameter vector {@link Objects.IterationValues#betas} according
     * to an addition move (a beta previously set to 0 is assigned a value).
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     */
    public void updateBetaAdd(
            ProposalDistributions propsds,
            Random r) {
        whichParameterTypeUpdated = ParameterTypes.BETA_ADD.ordinal();
        double newbeta = GeneralMaths.normalDraw(
                0,
                propsds.proposalDistributionSds[ParameterTypes.BETA_ADD.ordinal()],
                r);
        betas.set(whichBetaAdded,0,newbeta);
    }

    /**
     * Updates the beta parameter vector {@link Objects.IterationValues#betas} according
     * to a swap move (a beta previously set to 0 is assigned a value, and a 
     * previously included beta is set to 0).
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     */
    public void updateBetaSwap(
            ProposalDistributions propsds,
            Random r) {
        whichParameterTypeUpdated = ParameterTypes.BETA_SWAP.ordinal();
        double newbeta = GeneralMaths.normalDraw(
                0,
                propsds.proposalDistributionSds[ParameterTypes.BETA_SWAP.ordinal()],
                r);
        betas.set(whichBetaAdded,0,newbeta);
    }
    
    /**
     * Updates the between clusters SD parameter 
     * {@link Objects.IterationValues#logBetweenClusterSd}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     */
    public void updateBetweenClusterSd(
            ProposalDistributions propsds,
            Random r) {
        whichParameterTypeUpdated = ParameterTypes.BETWEEN_CLUSTER_SD.ordinal();
        logBetweenClusterSd = GeneralMaths.normalDraw(
                logBetweenClusterSd,
                propsds.proposalDistributionSds[ParameterTypes.BETWEEN_CLUSTER_SD.ordinal()], r);
    }
    
    /**
     * Updates the unknown prior SD common to the betas 
     * {@link Objects.IterationValues#betaPriorSds}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     * @param nSds The number of model space components
     */
    public void updateBetaPriorSd(
            ProposalDistributions propsds,
            Random r,
            int nSds) {
        whichParameterTypeUpdated = ParameterTypes.BETA_PRIOR_SD.ordinal();
        int sdUpdate = 0;
        if (nSds>1) {   // Pick one of the SDs at random to update
            sdUpdate = r.nextInt(nSds);
        }
        logBetaPriorSd[sdUpdate] = GeneralMaths.normalDraw(
                logBetaPriorSd[sdUpdate],
                propsds.proposalDistributionSds[ParameterTypes.BETA_PRIOR_SD.ordinal()], r);
        betaPriorSds[sdUpdate] = Math.exp(logBetaPriorSd[sdUpdate]);
    }
    
     /**
     * Updates the matrix product XBeta efficiently, with minimal computation,
     * by considering only betas which have been changed
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param curr {@link Objects.IterationValues} class object, containing the
     * current parameter values/states
     */
    public final void updateXBetaEfficiently(
            Arguments arguments,
            Data data,
            IterationValues curr) {
        /**
         * Efficiently modify XBeta by an increment to reflect removal of a
         * parameter
         */
        if (whichMove == 0 | whichMove == 2) {
            Xbeta = Xbeta.minus(
                    data.dataJama.getMatrix(0, data.numberOfIndividuals-1,
                            whichBetaRemoved, whichBetaRemoved)
                    .times(
                    curr.betas.get(whichBetaRemoved,0)
                    ) );                
        }
        
        /**
         * Efficiently modify XBeta by an increment to reflect addition of a
         * parameter
         */
        if (whichMove == 1 | whichMove == 2) {
            Xbeta = Xbeta.plus(
                    data.dataJama.getMatrix(0, data.numberOfIndividuals-1,
                            whichBetaAdded, whichBetaAdded)
                    .times(
                    (betas.get(whichBetaAdded,0)-
                    curr.betas.get(whichBetaAdded,0))
                    )
                    );
        }
        
        /**
         * Efficiently modify XBeta by an increment to reflect a beta update
         */
        if (whichMove == 3 & whichParameterTypeUpdated == ParameterTypes.BETAS.ordinal()) {
            Xbeta = Xbeta.plus(
                    data.dataJama.getMatrix(0, data.numberOfIndividuals-1,
                            whichBetaUpdated, whichBetaUpdated)
                    .times(
                    (betas.get(whichBetaUpdated, 0)-
                    curr.betas.get(whichBetaUpdated, 0))
                    )
                    );
        }
        
        /**
         * Efficiently modify XBeta to reflect an intercept update
         */
        if (whichMove == 3 & whichParameterTypeUpdated == ParameterTypes.ALPHA.ordinal()) {
            Matrix alphaDelta = new Matrix(data.numberOfIndividuals,1);
            if (data.numberOfClusters >0 ) {
                for (int i=0; i<data.numberOfIndividuals; i++) {
                    alphaDelta.set(i, 0,
                            (
                            (alpha-curr.alpha)
                            + (clusterIntercepts[(data.clusters[i]-1)]-
                            curr.clusterIntercepts[(data.clusters[i]-1)])
                            )
                            );
                }
            } else {
                for (int i=0; i<data.numberOfIndividuals; i++) {
                    alphaDelta.set(i, 0, (alpha-curr.alpha) );
                }
            }
            Xbeta = Xbeta.plus(alphaDelta);
        }
    }
    
     /**
     * Updates the matrix product XBeta efficiently, with minimal computation,
     * by considering only betas which have been changed - this version is for
     * block independence of the covariate matrix
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param curr {@link Objects.IterationValues} class object, containing the
     * current parameter values/states
     */
    public final void updateXBetaBlocksEfficiently(
            Arguments arguments,
            Data data,
            IterationValues curr) {
        /**
         * Efficiently modify XBeta by an increment to reflect removal of a
         * parameter
         */
        if (whichMove == 0 | whichMove == 2) {
            XbetaBlocks[whichBetaRemovedBlock]
                    = XbetaBlocks[whichBetaRemovedBlock]
                            .minus((data.dataJamaBlocks[whichBetaRemovedBlock]
                                    .getMatrix(0,
                                            data.blockSizes[whichBetaRemovedBlock]-1,
                                            (whichBetaRemoved-data.cumulativeBlockSizes[whichBetaRemovedBlock]),
                                            (whichBetaRemoved-data.cumulativeBlockSizes[whichBetaRemovedBlock])))
                                    .times(curr.betas.get(whichBetaRemoved,0)));
        }
        
        /**
         * Efficiently modify XBeta by an increment to reflect addition of a
         * parameter
         */
        if (whichMove == 1 | whichMove == 2) {
            XbetaBlocks[whichBetaAddedBlock]
                    = XbetaBlocks[whichBetaAddedBlock]
                            .plus((data.dataJamaBlocks[whichBetaAddedBlock]
                                    .getMatrix(0,
                                            data.blockSizes[whichBetaAddedBlock]-1,
                                            (whichBetaAdded-data.cumulativeBlockSizes[whichBetaAddedBlock]),
                                            (whichBetaAdded-data.cumulativeBlockSizes[whichBetaAddedBlock])))
                                    .times(betas.get(whichBetaAdded,0)
                                            -curr.betas.get(whichBetaAdded,0)));
        }
        
        /**
         * Efficiently modify XBeta by an increment to reflect a beta update
         */
        if (whichMove == 3 & whichParameterTypeUpdated == ParameterTypes.BETAS.ordinal()) {
            XbetaBlocks[whichBetaUpdatedBlock]
                    = XbetaBlocks[whichBetaUpdatedBlock]
                            .plus((data.dataJamaBlocks[whichBetaUpdatedBlock]
                                    .getMatrix(0,
                                            data.blockSizes[whichBetaUpdatedBlock]-1,
                                            (whichBetaUpdated-data.cumulativeBlockSizes[whichBetaUpdatedBlock]),
                                            (whichBetaUpdated-data.cumulativeBlockSizes[whichBetaUpdatedBlock])))
                                    .times(betas.get(whichBetaUpdated,0)
                                            -curr.betas.get(whichBetaUpdated,0)));
        }        
    }
    
    
    /**
     * Updates the Weibull scale parameter {@link Objects.IterationValues#weibullScale}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     */
    public void updateWeibullK(
            ProposalDistributions propsds,
            Random r) {
        whichParameterTypeUpdated = ParameterTypes.WEIBULL_SCALE.ordinal();
        logWeibullScale = GeneralMaths.normalDraw(
                logWeibullScale,
                propsds.proposalDistributionSds[ParameterTypes.WEIBULL_SCALE.ordinal()], r);
        weibullScale = Math.exp(logWeibullScale);
    }

    /**
     * Updates the Gaussian residual parameter {@link Objects.IterationValues#gaussianResidual}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     */
    public void updateGaussianResidual(
            ProposalDistributions propsds,
            Random r) {
        whichParameterTypeUpdated = ParameterTypes.GAUSSIAN_RESIDUAL.ordinal();
        logGaussianResidual = GeneralMaths.normalDraw(
                logGaussianResidual,
                propsds.proposalDistributionSds[ParameterTypes.GAUSSIAN_RESIDUAL.ordinal()], r);
        gaussianResidual = Math.exp(logGaussianResidual);
    }
    
    /**
     * Updates the clusters intercepts {@link Objects.IterationValues#clusterIntercepts}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     */
    public void updateClusterIntercepts(
            ProposalDistributions propsds,
            Random r) {
        whichParameterTypeUpdated = ParameterTypes.CLUSTER_INTERCEPTS.ordinal();
        for (int i=0; i<numberOfClusters; i++) {
            clusterIntercepts[i] = GeneralMaths.normalDraw(
                    clusterIntercepts[i],
                    propsds.proposalDistributionSds[ParameterTypes.CLUSTER_INTERCEPTS.ordinal()], r);
            
        }
    }

    /**
     * Updates the {@link Objects.IterationValues#model} by first choosing a move,
     * then choosing which (if any) of the betas to add/remove
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param randomDraws Random number generator object
     */
    public void updateModel(
            Arguments arguments,
            Data data,
            Random randomDraws) {
            // Count number variables included in current model to determine
            // whichMove
        whichMove = GeneralMethods.chooseMove( (totalNumberOfCovariates-numberOfCovariatesToFixInModel),
                    (modelDimension-numberOfCovariatesToFixInModel),
                    arguments.moveProbabilities,randomDraws);
        if (arguments.useReversibleJump == 0) {whichMove = 3;}
        
        int presentMarkN = modelDimension - numberOfCovariatesToFixInModel;
        int markRem;
        int markAdd;

        if ((whichMove==0)||(whichMove==2)) {
            // If whichMove is removal (0) or swap (2) choose at random a marker
            // to remove and place a 1 in corresponding position in presentMark
            markRem = randomDraws.nextInt(presentMarkN);
            int count = 0;
            if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                /**
                 * Loop over blocks for GaussianMarg
                 */
                for (int b=0; b<data.nBlocks; b++) {
                    for (int m=data.blockIndices[b]; m<data.blockIndices[b+1]; m++) {
                        if (model[m]==1) {
                            if (count==markRem) {
                                whichBetaRemoved = m;
                                whichBetaRemovedBlock = b;
                                model[m] = 0;
                                betas.set(m, 0, 0);
                            }
                            count++;
                        }                        
                    }
                }
            } else {
                for (int m=numberOfCovariatesToFixInModel; m<totalNumberOfCovariates; m++) {
                    if (model[m]==1) {
                        if (count==markRem) {
                            whichBetaRemoved = m;
                            model[m] = 0;
                            betas.set(m, 0, 0);
                        }
                        count++;
                    }
                }                
            }
        }

        int missMarkN = totalNumberOfCovariates - numberOfCovariatesToFixInModel - presentMarkN;
        if ((whichMove==1)||(whichMove==2)) {
            // If whichMove is addition (1) or swap (2) choose at random a marker
            // to add and place a 1 in corresponding position in presentMark
            markAdd = randomDraws.nextInt(missMarkN);
            int count = 0;
            if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                /**
                 * Loop over blocks for GaussianMarg
                 */
                for (int b=0; b<data.nBlocks; b++) {
                    for (int m=data.blockIndices[b]; m<data.blockIndices[b+1]; m++) {
                        if (model[m]==0) {
                            if (count==markAdd) {
                                whichBetaAdded = m;
                                whichBetaAddedBlock = b;
                                model[m] = 1;
                            }
                            count++;
                        }                        
                    }
                }                
            } else {
                for (int m=numberOfCovariatesToFixInModel; m<totalNumberOfCovariates; m++) {
                    if (model[m]==0) {
                        if (count==markAdd) {
                            whichBetaAdded = m;
                            model[m] = 1;
                        }
                        count++;
                    }
                }                
            }
        }        

        modelDimension = GeneralMethods.countPresVars(totalNumberOfCovariates, model);
        modelSpacePartitionDimensions = GeneralMethods.countPresVarsComps(
                arguments.numberOfModelSpacePriorPartitions, data.modelSpacePartitionIndices, model);
        
    }

    /**
     * Calculates the log-likelihood according to the stored model and parameter
     * values, and a dataset, and stores in {@link Objects.IterationValues#logLikelihood}.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     */
    public final void calcLogLike(
            Arguments arguments,
            Data data) {

        /***
         * Start by calculating the linear predictor; alpha + Xbeta
         */
        // Alpha
        Matrix alphaMat = new Matrix(data.numberOfIndividuals,1);
        if (data.numberOfClusters >0 ) {
            for (int i=0; i<data.numberOfIndividuals; i++) {
                alphaMat.set(i, 0, (alpha + clusterIntercepts[(data.clusters[i]-1)]) );
            }
        } else {
            for (int i=0; i<data.numberOfIndividuals; i++) {
                alphaMat.set(i, 0, alpha);
            }
        }
        // Xbeta

        /***
         * Use the linear predictor to calculate the log-likelihood
         */
        logLikelihood = 0;        
        if (data.whichLikelihoodType==LikelihoodTypes.LOGISTIC.ordinal()) {
            // Logistic likelihood contributions
            Xbeta = data.dataJama.times(betas);
            Xbeta.plusEquals(alphaMat);
            double p;
            for (int i=0; i<data.numberOfIndividuals; i++) {
                p = (double) (Math.exp( Xbeta.get(i, 0) )
                        / (1+Math.exp(Xbeta.get(i, 0)) ));
                if (data.outcomes[i] == 1) {
                    logLikelihood = logLikelihood + Math.log(p);
                } else {
                    logLikelihood = logLikelihood + Math.log(1-p);
                }
            }            
        } else if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
            // Weibull likelihood contributions
            Xbeta = data.dataJama.times(betas);
            Xbeta.plusEquals(alphaMat);
            for (int i=0; i<data.numberOfIndividuals; i++) {
                // Log-Hazard function if event occurred at the end of follow up
                if (data.outcomes[i] == 1) {
                    logLikelihood = logLikelihood
                            + Math.log(weibullScale)
                            + weibullScale*Xbeta.get(i, 0)
                            + (weibullScale-1)*Math.log(data.times[i]);
                }
                // Log-Survival function
                logLikelihood = logLikelihood
                        -Math.pow(
                        (double)(data.times[i]*Math.exp( Xbeta.get(i, 0) )),
                        weibullScale);
            }            
        } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()) {
            // Gaussian likelihood contibution (use Weibullk as sigma
            Xbeta = data.dataJama.times(betas);
            Xbeta.plusEquals(alphaMat);
            yMinXbeta = data.continuousOutcomesJama.minus(Xbeta);            
            logLikelihood = logLikelihood
                    - (double) ((yMinXbeta.transpose().times(yMinXbeta)).get(0, 0)/(2*gaussianResidual*gaussianResidual))
                    - data.numberOfIndividuals*logGaussianResidual;
        } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
            /**
             * Calculate the likelihood score term within each block, summing
             * during the loop
             */
            likelihoodScoreTerm = 0;
            for (int b=0; b<data.nBlocks; b++) {
                XbetaBlocks[b] = data.dataJamaBlocks[b].times(betas.getMatrix(
                                data.blockIndices[b],
                                (data.blockIndices[(b+1)]-1), 0, 0));                
                yMinXbeta = data.continuousOutcomesJama.getMatrix(
                        data.blockIndices[b],
                        (data.blockIndices[(b+1)]-1),0,0)
                        .minus(XbetaBlocks[b]);            
                likelihoodScoreTermsBlocks[b] = (yMinXbeta.transpose()
                        .times(data.xTxInvBlocks[b])
                        .times(yMinXbeta)).get(0,0);
                likelihoodScoreTerm=likelihoodScoreTerm+likelihoodScoreTermsBlocks[b];
            }
            logLikelihood = logLikelihood
                    - (double) (likelihoodScoreTerm/(2*gaussianResidual*gaussianResidual))
                    - data.numberOfIndividuals*logGaussianResidual;            
        }
    }

    /**
     * Efficiently calculates the log-likelihood by modifying the previously
     * stored value according to the specific changes to the model and parameter
     * values since then, and stores in {@link Objects.IterationValues#logLikelihood}.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param curr {@link Objects.IterationValues} class object containing the previous
     * model and parameter values, to use as a starting point
     */
    public final void calcLogLike_IncrementFromProposal(
            Arguments arguments,
            Data data,
            IterationValues curr) {
        
        /***
         * Calculate log-likelihood, using Xbeta efficiently updated above
         **/        
        logLikelihood = 0;
        
        /**
         * Logistic regression
         */
        if (data.whichLikelihoodType==LikelihoodTypes.LOGISTIC.ordinal()) {
            // Logistic likelihood contributions
            updateXBetaEfficiently(arguments, data, curr);
            double p;
            for (int i=0; i<data.numberOfIndividuals; i++) {
                p = (double) (Math.exp( Xbeta.get(i, 0) )
                        / (1+Math.exp(Xbeta.get(i, 0)) ));
                if (data.outcomes[i] == 1) {
                    logLikelihood = logLikelihood + Math.log(p);
                } else {
                    logLikelihood = logLikelihood + Math.log(1-p);
                }
            }            
        }
        
        /**
         * Weibull model
         */
        if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
            // Weibull likelihood contributions
            updateXBetaEfficiently(arguments, data, curr);
            for (int i=0; i<data.numberOfIndividuals; i++) {
                logLikelihood = logLikelihood 
                        -(Math.exp( Xbeta.get(i, 0) )*Math.pow(data.times[i], weibullScale));
                if (data.outcomes[i] == 1) {
                    logLikelihood = logLikelihood
                            + Math.log(weibullScale)
                            + ((weibullScale-1)*Math.log(data.times[i]))
                            + Xbeta.get(i, 0);
                }
            }
        }
        
        /**
         * Gaussian linear regression
         */
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()) {
            updateXBetaEfficiently(arguments, data, curr);
            yMinXbeta = data.continuousOutcomesJama.minus(Xbeta);            
            logLikelihood = logLikelihood
                    -(double) ((yMinXbeta.transpose().times(yMinXbeta)).get(0, 0)/(2*gaussianResidual*gaussianResidual))
                    -data.numberOfIndividuals*logGaussianResidual;
        }
        
        /**
         * Regression of Gaussian marginal statistics
         */
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
            if ((whichMove<3)|
                    (whichParameterTypeUpdated==ParameterTypes.BETAS.ordinal())) {
                updateXBetaBlocksEfficiently(arguments, data, curr);
                /**
                 * Determine which blocks were updated, and update score for
                 * those only
                 */
                whichBlocksUpdated = GeneralMethods.determineUpdatedBlocks(
                        data.nBlocks,
                        data.blockIndices,
                        whichMove,
                        whichBetaRemoved,
                        whichBetaAdded,
                        whichBetaUpdated);                
                for (int b=0; b<data.nBlocks; b++) {
                    if (whichBlocksUpdated[b]) {
                        // Calculate score term for block only
                        yMinXbeta = data.continuousOutcomesJama.getMatrix(
                                data.blockIndices[b],
                                (data.blockIndices[(b+1)]-1),
                                0,
                                0).minus(XbetaBlocks[b]) ;            
                        likelihoodScoreTermsBlocks[b] = (yMinXbeta.transpose()
                                .times(data.xTxInvBlocks[b])
                                .times(yMinXbeta)
                                ).get(0, 0);
                        // Update total score with the difference
                        likelihoodScoreTerm=likelihoodScoreTerm+
                                (likelihoodScoreTermsBlocks[b]
                                -curr.likelihoodScoreTermsBlocks[b]);
                    }
                }                
            }
            /**
             * Calculate the log-likelihood, using the score term from above
             * Note that this covers the case for a gaussianResidual update
             **/
            logLikelihood = logLikelihood
                    - (double) (likelihoodScoreTerm/(2*gaussianResidual*gaussianResidual))
                    - data.numberOfIndividuals*logGaussianResidual;            
        }
    }

    /**
     * Calculates the log-prior support of the stored model and parameter
     * values according to the prior distributions, and stores in 
     * {@link Objects.IterationValues#logPrior}.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param priors {@link Objects.Priors} class object, containing information
     * on the prior distributions
     */
    public final void calcLogPrior(
            Arguments arguments,
            Data data,
            Priors priors) {

        /*
         * Initiate log-prior.
         */
        logPrior = 0;

        /***
         * Global intercept, alpha
         ***/
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
            logPrior = logPrior + GeneralMaths.logNormDens(
                alpha,
                0,
                1000
                );            
        } else {
            logPrior = logPrior + GeneralMaths.logNormDens(
                alpha,
                arguments.alphaPriorMu,
                arguments.alphaPriorSd
                );
        }
        
        /***
         * Study specific intercepts
         ***/        
        if (data.numberOfClusters >0) {
            for (int r=0; r<data.numberOfClusters; r++) {
                // dont need all H, since Hth is implied by rest
                logPrior = logPrior + GeneralMaths.logNormDens(
                        clusterIntercepts[r],
                        0,
                        Math.exp(logBetweenClusterSd)
                        );
            }

            if (arguments.betweenClusterSdPriorFamily == 0) {
                // uniform prior for SD parameters
                logPrior = logPrior +
                        Math.log(priors.betweenClusterPrecisionUniformPrior.density(Math.exp(logBetweenClusterSd)));
            } else if (arguments.betweenClusterSdPriorFamily==1) {
                // uniform prior for SD parameters
                logPrior = logPrior +
                        Math.log(priors.betweenClusterPrecisionGammaPrior.density(
                        (double)1/(Math.exp(logBetweenClusterSd)*Math.exp(logBetweenClusterSd))
                        )
                        );
            }
        }
        
        /***
         * Covariates with informative priors
         ***/        
        for (int v = 0; v < data.numberOfCovariatesWithInformativePriors; v++) {
            if (model[v]==1) {
                logPrior = logPrior
                        + GeneralMaths.logNormDens(
                        betas.get(v, 0),
                        data.betaPriorMus[v] ,
                        data.betaPriorSds[v] );
            }
        }
        
        /***
         * Covariates with unknown priors (and common SD)
         ***/        
        for (int c=0; c<data.numberOfUnknownBetaPriors; c++) {
            for (int v = data.commonBetaPriorPartitionIndices[c]; v < data.commonBetaPriorPartitionIndices[c+1]; v++) {
                if (model[v]==1) {
                    logPrior = logPrior
                            + GeneralMaths.logNormDens(
                            betas.get(v, 0),
                            data.betaPriorMus[v] ,
                            betaPriorSds[c] );
                }
            }
        }
        
        /***
         * Covariate prior precision(s)
         ***/        
        if (data.numberOfUnknownBetaPriors > 0) {        
            // Beta prior SD hyperparameter
            if (arguments.betweenClusterSdPriorFamily == 0) {
                // uniform prior for SD parameters
                for (int c=0; c<data.numberOfUnknownBetaPriors; c++) {
                    logPrior = logPrior +
                            Math.log(priors.betweenClusterPrecisionUniformPrior.density(betaPriorSds[c]));
                }
            } else if (arguments.betweenClusterSdPriorFamily==1) {
                for (int c=0; c<data.numberOfUnknownBetaPriors; c++) {
                    logPrior = logPrior +
                            Math.log(priors.betaPrecisionUniformPrior.density(betaPriorSds[c]));
                }
            }
        }
        
        /***
         * Weibull scale parameter k
         ***/
        if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
            // Normal prior on the log N(0,10e6) is used in the paper
            logPrior = logPrior
                    + GeneralMaths.logNormDens(logWeibullScale, 0, 1000);            
        }
        
        /***
         * Gaussian residual precision
         ***/
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
//            logPrior = logPrior
//                    + Math.log(priors.gaussianResidualPrecisionPrior.density(
//                            (double) (1/(gaussianResidual*gaussianResidual)) ));
            // Jeffrey's prior: P(sigma^2)=1/(sigma^2)
            logPrior = logPrior - (double) (2*Math.log(gaussianResidual) );
        }
    }
    
    /**
     * Calculates the Metropolis Hastings Reversible Jump acceptance probability
     * for the model and parameter values in this object vs the previously
     * saved model and parameter values, and stores in the field
     * {@link Objects.IterationValues#acceptanceProbability}.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param curr {@link Objects.IterationValues} class object containing the current
     * model and parameter values
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     */
    public void acceptanceProbability(
            Arguments arguments,
            Data data,
            IterationValues curr,
            ProposalDistributions propsds
            ) {
            // log-P('prop')
            double logNumerator = logLikelihood+logPrior;
            // log-P('curr')
            double logDenominator = curr.logLikelihood+curr.logPrior;
            if (whichMove==0) {
                // If removal happened logNumerator needs added to it the
                // log-probability of adding the removed marker, and drawing
                // the corresponding removed log-OR
                logNumerator = logNumerator +
                        GeneralMaths.logNormDens(curr.betas.get(whichBetaRemoved, 0),
                        0, propsds.proposalDistributionSds[ParameterTypes.BETA_ADD.ordinal()]);
                logNumerator = logNumerator
                        +Math.log(arguments.moveProbabilities[1]
                        -arguments.moveProbabilities[0])
                        -Math.log(data.totalNumberOfCovariates-modelDimension);
                logDenominator = logDenominator
                        +Math.log(arguments.moveProbabilities[0])
                        - Math.log(modelDimension-data.numberOfCovariatesToFixInModel +1);
            }
            if (whichMove==1) {
                // If addition happened logDenominator needs added to it the
                // log-probability of adding the added marker, and drawing
                // the corresponding new log-OR
                logDenominator = logDenominator +
                        GeneralMaths.logNormDens(betas.get(whichBetaAdded, 0), 0,
                        propsds.proposalDistributionSds[ParameterTypes.BETA_ADD.ordinal()]);
                logNumerator = logNumerator
                        +Math.log(arguments.moveProbabilities[0])
                        - Math.log(modelDimension-data.numberOfCovariatesToFixInModel);
                logDenominator = logDenominator
                        + Math.log(arguments.moveProbabilities[1]
                        -arguments.moveProbabilities[0])
                        -Math.log(data.totalNumberOfCovariates-modelDimension+1);
            }
            if (whichMove==2) {
                // If swap happened logNumerator needs added to it the
                // log-probability of drawing the removed log-OR and
                // logDenominator needs added to it the log-probability of
                // drawing the added log-OR. NOTE: The probabilties of removing
                // a random marker and adding a random marker cancel since both
                // 'curr' and 'prop' have the same number of markers present.
                logNumerator = logNumerator +
                        GeneralMaths.logNormDens(curr.betas.get(whichBetaRemoved, 0),
                        0, propsds.proposalDistributionSds[ParameterTypes.BETA_SWAP.ordinal()]);
                logDenominator = logDenominator +
                        GeneralMaths.logNormDens(betas.get(whichBetaAdded, 0),
                        0, propsds.proposalDistributionSds[ParameterTypes.BETA_SWAP.ordinal()]);
            }

            // Model Space prior is accounted for, which is not accounted
            // for in method 'calcLogPrior'. A truncated Poisson Prior for
            // is used for the model space.
            double logModelPriorRatio = 0;
            if (arguments.modelSpacePriorFamily==0) {
                // Poisson model space prior
                for (int c=0; c<arguments.numberOfModelSpacePriorPartitions; c++) {
                    logModelPriorRatio = logModelPriorRatio +
                            GeneralMethods.logModelPriorRatio(
                            curr.modelSpacePartitionDimensions[c],
                            modelSpacePartitionDimensions[c],
                            data.modelSpacePoissonPriorMeans[c]);
                } 
            } else if (arguments.modelSpacePriorFamily==1) {
                // Beta-binomial model space prior
                for (int c=0; c<arguments.numberOfModelSpacePriorPartitions; c++) {
                    logModelPriorRatio = logModelPriorRatio +
                            GeneralMethods.logModelPriorRatioBetaBin(
                            curr.modelSpacePartitionDimensions[c],
                            modelSpacePartitionDimensions[c],
                            data.modelSpacePartitionSizes[c],
                            arguments.modelSpaceBetaBinomialPriorHyperparameterA[c],
                            arguments.modelSpaceBetaBinomialPriorHyperparameterB[c]);
                }                
            }

            acceptanceProbability = Math.min(1,
                    Math.exp(logNumerator
                    -logDenominator+logModelPriorRatio));
    }
    
}

