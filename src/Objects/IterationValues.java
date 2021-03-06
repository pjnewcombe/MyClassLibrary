package Objects;

import Methods.GeneralMaths;
import Jama.Matrix;
import Methods.GeneralMethods;
import java.util.Arrays;
import java.util.Comparator;
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
     * 
     * Modeling parameters. This includes the model parameters, information
     * describing the selected model, and information describing the beta
     * priors.
     * 
     */
    
    /**
     * General: Model intercept.
     */
    public double alpha;
    /**
     * General: Model covariate parameters (eg log-ORs for logistic regression).
     */
    public Matrix betas;
    /**
     * Normalised weights for use in the ROCAUC model.
     */
    public Matrix rocAucWeights;
    /**
     * Weibull model: Weibull scale parameter.
     */
    public double weibullScale;
    /**
     * Weibull model: Log-Weibull scale parameter.
     */
    public double logWeibullScale;
    /**
     * ROC_AUC: Dirichlet concentration parameter.
     */
    public double dirichletConcentration;
    /**
     * ROC_AUC: log-Dirichlet concentration parameter.
     */
    public double logDirichletConcentration;    
    /**
     * Guassian model: Residual parameter.
     */
    public double gaussianResidual;
    /**
     * Guassian model: Log-residual parameter.
     */
    public double logGaussianResidual;
    /**
     * General: Model. Vector indicating which covariates are included in the
     * model (length = total number of covariates being searched over).
     */
    public int[] model;
    /**
     * General: Model dimension. The total number of covariates included in the
     * model.
     */
    public int modelDimension;
    /**
     * General: Model dimension by partition. Total number of covariates included
     * in each component of the model space (length = number of model space
     * components).
     */
    public int[] modelSpacePartitionDimensions;
    /**
     * General: Prior SD for the betas. Assumed common across all betas in
     * a model space partition (length = number of model space components).
     */
    public double[] betaPriorSds;
    /**
     * General: log of the prior SD for the betas. Log-unknown common prior SD
     * across the betas (length = number of model space components).
     */
    public double[] logBetaPriorSds;    
    
    /**
     * 
     * Prior, likelihood and acceptance ratio.
     * 
     */
    
    /**
     * Prior support (logged).
     */
    public double logPrior;
    /**
     * Log-Likelihood.
     */
    public double logLikelihood;
    /**
     * Ratio of priors for current vis proposed models.
     */
    public double logModelPriorRatio;
    /**
     * Acceptance probability for this model and parameter values in comparison
     * to the current state.
     */
    public double acceptanceProbability;
    
    /**
     * 
     * Partial likelihood calculations. Below are quantities which are
     * calculated and stored during the likelihood calculations.
     * 
     */
    
    /**
     * Product of the covariate matrix and coefficient vector.
     */
    public Matrix Xbeta;        
    /**
     * Product of the covariate matrix and coefficient vector. Stored by block
     * for GuassianMarg.
     */
    public Matrix[] XbetaBlocks;        
    /**
     * Block specific likelihood terms to be stored. For GaussianMarg these are
     * the block-specific likelihood score terms.
     **/
    public double[] likelihoodTermBlocks;
    /**
     * GaussianMarg: Y - X*Beta.
     */ 
    private Matrix yMinXbeta;
    /**
     * GaussianMarg: A part of the likelihood calculation to be stored, allowing minor updates
     * iteration to iteration. For GaussianMarg this is the likelihood score
     * term.
     **/
    public double likelihoodTerm;
    /**
     * GaussianMargCong: Block specific matrix products. This is the
     * (LInvt)'X_g(X_g'X_g)^-1X_gLInvt
     * which are useful when only tau is updated
     **/
    public double[] S_Gamma_MatrixProductBlocks;
    /**
     * Cox regression: Vector of cumulative risk set Xbs. This is the
     **/
    public double[] cumulativeExpXbSumsOverRiskSets;
    
    /**
     * 
     * Extra quantities for the ROC likelihood calculation.
     * 
     */

    /**
     * Vector of indices. For use with the ROC likelihood.
     */
    private final Integer[] indices;
    /**
     * ROC AUC.
     */
    public double auc;
    /**
     * Current x-axis tick. Used in calculating the ROC AUC.
     */
    private double xNow;
    /**
     * Next x-axis tick. Used in calculating the ROC AUC.
     */
    private double xNext = 0;
    /**
     * Current y-axis tick. Used in calculating the ROC AUC.
     */
    private double yNow;
    /**
     * Next y-axis tick. Used in calculating the ROC AUC.
     */
    private double yNext;

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
     * Indicates what the probability of removing a covariate was for this
     * instance of IterationValues.
     * 0: Removal, 1: Addition, 2: Swap, 3: Null.
     */
    public double probRemove;
    /**
     * Indicates what the probability of adding a covariate was for this
     * instance of IterationValues.
     */
    public double probAdd;
    /**
     * Indicates what the probability of a swap move was for this
     * instance of IterationValues.
     */
    public double probSwap;
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
        /**
         * General. Setup parameters which are always used.
         */
        
        totalNumberOfCovariates = data.totalNumberOfCovariates;
        numberOfCovariatesToFixInModel = data.numberOfCovariatesToFixInModel;
        
        /**
         * Alpha, model and betas.
         */
        
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
            alpha = 0;
        } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
            alpha = 0;
        } else {
            alpha = arguments.initialAlpha;
        }
        model = new int[data.totalNumberOfCovariates];
        betas = new Matrix(data.totalNumberOfCovariates,1);
        for (int v=0; v<data.totalNumberOfCovariates; v++) {
            if ((data.initialModelOption == 1)|(v < data.numberOfCovariatesToFixInModel)) {
                model[v] = 1;
            } else if (data.initialModelOption == 0) {
                model[v] = 0;
            } else if (data.initialModelOption == 2) {
                model[v] = data.initialModel[v];                
            }
            betas.set(v, 0, arguments.initialBetas);
        }
        modelDimension = GeneralMethods.countPresVars(data.totalNumberOfCovariates, model);
        if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()) {
            // Initiate the weights
            for (int v=0; v<data.totalNumberOfCovariates; v++) {
                if (model[v]==1) {
                    betas.set(v, 0, (double) 1/modelDimension);
                }
            }
        } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
            betas.set(0, 0, 1); // First covariate (which is fixed in model) is set to 1
        } 
        modelSpacePartitionDimensions = GeneralMethods.countPresVarsComps(
                arguments.numberOfModelSpacePriorPartitions, data.modelSpacePartitionIndices, model);
        updateMoveProbabilities(arguments);
        
        /**
         * Model space prior parameters.
         */
        
        betaPriorSds = new double[data.numberOfHierarchicalCovariatePriorPartitions];
        logBetaPriorSds = new double[data.numberOfHierarchicalCovariatePriorPartitions];
        for (int c=0; c<data.numberOfHierarchicalCovariatePriorPartitions; c++) {
            betaPriorSds[c] = data.initialBetaPriorSds[c];
            if (
                    data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()|
                    data.whichLikelihoodType==LikelihoodTypes.JAM.ordinal()|
                    data.whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()
                    ) {
                /**
                 * Model selection coefficient
                 */
                betaPriorSds[c] = Math.sqrt(data.tau);
            }
            logBetaPriorSds[c] = Math.log(betaPriorSds[c]);
        }
        
        /**
         * Weibull scale parameter.
         */
        
        weibullScale = arguments.initialWeibullScale;
        logWeibullScale = Math.log(weibullScale);
        
        /**
         * Dirichlet concentration parameter
         */
        
        dirichletConcentration = arguments.initialDirichletConcentration;
        logDirichletConcentration = Math.log(dirichletConcentration);
        
        /**
         * Gaussian residual.
         */
        if (data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
            /**
             * Set according to the informative prior. NECESSARY!
             */
            gaussianResidual = Math.sqrt(arguments.gaussianResidualVarianceInvGammaPrior_a/arguments.gaussianResidualVarianceInvGammaPrior_b);
        } else {
            gaussianResidual = arguments.initialGaussianResidual;        
        }
        logGaussianResidual = Math.log(gaussianResidual);
        
        /**
         * Terms used in the likelihood calculation.
         */
        
        likelihoodTermBlocks = new double[data.nBlocks];
        S_Gamma_MatrixProductBlocks = new double[data.nBlocks];
        cumulativeExpXbSumsOverRiskSets = new double[data.numberOfIndividuals];
        Xbeta = new Matrix(0,0); // Needs to be initialised, if not used, for setTo
        XbetaBlocks = new Matrix[data.nBlocks];
        whichBlocksUpdated = new boolean[data.nBlocks];
        
        /*
        * Extra terms used for ROC likelihood.
        */
        indices = new Integer[data.numberOfIndividuals];
        for (int i=0; i<data.numberOfIndividuals; i++) {
            indices[i] = i;
        }
        
        /**
         * Initiate the prior and likelihood.
         */
        
        if (data.whichLikelihoodType==LikelihoodTypes.JAM.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            conjugate_calcLogLike(arguments, data);
            if (data.modelTau==1) {
                conjugate_calcLogPrior(arguments, data, priors);                
            }
        } else {
            calcLogLike(arguments, data);
            calcLogPrior(arguments, data, priors);            
        }
    }
    
    /**
     * Sets this object equal to another object also of the class
     * {@link Objects.IterationValues}
     * 
     * @param its Another instance of {@link Objects.IterationValues} class to set all
     * parameters equal to.
     */
    public void setTo(IterationValues its) {
        /**
         * Vector array parameters.
         */
        model = its.model.clone();
        modelSpacePartitionDimensions = its.modelSpacePartitionDimensions.clone();
        logBetaPriorSds = its.logBetaPriorSds.clone();
        betaPriorSds = its.betaPriorSds.clone();
        likelihoodTermBlocks = its.likelihoodTermBlocks.clone();
        S_Gamma_MatrixProductBlocks = its.S_Gamma_MatrixProductBlocks.clone();
        
        /**
         * Single parameters.
         */
        alpha = its.alpha;
        weibullScale = its.weibullScale;
        logWeibullScale = its.logWeibullScale;
        dirichletConcentration = its.dirichletConcentration;
        logDirichletConcentration = its.logDirichletConcentration;
        gaussianResidual = its.gaussianResidual;
        logGaussianResidual = its.logGaussianResidual;
        auc = its.auc;
        
        /**
         * Likelihood and prior.
         */
        likelihoodTerm = its.likelihoodTerm;
        logLikelihood = its.logLikelihood;
        logPrior = its.logPrior;
        
        /**
         * Model dimension and move probabilities.
         */
        modelDimension = its.modelDimension;
        probRemove = its.probRemove;
        probAdd = its.probAdd;
        probSwap = its.probSwap;
        
        /**
         * Matrices used in the likelihood calculation.
         */
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
            // NULL MOVE for when no clusters intercepts -----------------------
            double paramTypeDraw = randomDraws.nextInt((2+data.numberOfHierarchicalCovariatePriorPartitions+data.nExtraParametersBeyondLinPred));
            if (paramTypeDraw == 0) { // Update the intercept (except for marginal)
                if (data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
                    // No intercept so update the residual
                    updateGaussianResidual(propsds, randomDraws);                        
                } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()|
                        data.whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()|
                        data.whichLikelihoodType==LikelihoodTypes.COX.ordinal()|
                        data.whichLikelihoodType==LikelihoodTypes.CASECOHORT_BARLOW.ordinal()|
                        data.whichLikelihoodType==LikelihoodTypes.CASECOHORT_PRENTICE.ordinal()) {
                    // No intercept in these likelihoods so update the betas
                    if (modelDimension>0) {
                        updateBetas(propsds, randomDraws, data);
                    }
                } else {
                    updateAlpha(propsds, randomDraws);
                }                    
            } else if (paramTypeDraw == 1) { // Update betas
                if (modelDimension>0) {
                    updateBetas(propsds, randomDraws, data);
                } else if (data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
                    updateGaussianResidual(propsds, randomDraws);                        
                } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()|
                        data.whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
                    // NO POSSIBLE FURTHER ACTION IF MODEL DIM 0 (no intercept)
                } else {
                    updateAlpha(propsds, randomDraws);                        
                }
            } else if (paramTypeDraw == 2) { // Extra parameter -- weibull k, residual -- or hierarchical beta SD
                if (data.nExtraParametersBeyondLinPred==1) {
                    if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
                        updateWeibullK(propsds, randomDraws);
                    } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                        data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
                        updateGaussianResidual(propsds, randomDraws);                                                
                    } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()) {
                        updateDirichletConcentration(propsds, arguments.logDirichletConcentrationMinimum, randomDraws);
                   }
                } else {
                    updateBetaPriorSd(propsds, randomDraws, data.numberOfHierarchicalCovariatePriorPartitions);                    
                }
            } else if (paramTypeDraw >= 3) { // Only occurs if there is an extra param and hierarchical SD(s) or there is more than 1 hierarchical SD
                updateBetaPriorSd(propsds, randomDraws, data.numberOfHierarchicalCovariatePriorPartitions);                    
            }
        }
        
        if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()) {
            // Normalise betas
            betas = GeneralMethods.normaliseAbsoluteBetasToSumToOne(betas, model);
        }
        
        calcLogPrior(arguments, data, priors);            
        calcLogLike_IncrementFromProposal(arguments, data, curr);
        acceptanceProbability(arguments, data, curr, propsds);
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
    public final void conjugate_calcLogLike(
            Arguments arguments,
            Data data) {

        /***
         * Use the linear predictor to calculate the log-likelihood
         */
        logLikelihood = 0;
        
        /**
         * 
         * IPD conjugate Gaussian model likelihood calculation.
         * 
         */
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            if (modelDimension==0) {
                logLikelihood = 
                        -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                        *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+data.YtY);                
            } else {
                /**
                 * Get list of included covariate(s). Used to extract X_gamma.
                 */
                int[] covInds = GeneralMethods.getVarIndices(
                        modelDimension,model);
                /**
                 * Extract X_g and X_g'X_g.
                 */
                Matrix X_g = data.X.getMatrix(
                        0,
                        data.numberOfIndividuals-1,
                        covInds);
                Matrix X_gTX_g = data.XtX.getMatrix(
                        covInds,
                        covInds);                
                if (data.useGPrior==0) {
                    /**
                     * For the independence prior we add tau to each element
                     * of the X_gamma'X_gamma diaganol before taking the
                     * inverse and determinant below.
                     */
                    for (int v=0; v<modelDimension; v++) {
                        X_gTX_g.set(v, v, (X_gTX_g.get(v, v)+data.tau) );
                    }
                }

                /**
                 * Calculate the complex matrix product in S_g.
                 */
                double S_Gamma_MatrixProduct = (
                        data.Y.transpose().times(
                                X_g.times(
                                        X_gTX_g.inverse()
                                ).times(
                                        X_g.transpose().times(data.Y)
                                ))
                        ).get(0, 0);

                /**
                 * Calculate S_gamma, and the loglikelihood contribution from
                 * the block.
                 */
                if (data.useGPrior==1) {
                    double S_gamma = data.YtY
                            - (data.tau/(1+data.tau))*S_Gamma_MatrixProduct;
                    logLikelihood = 
                            -((double)modelDimension/2)*Math.log(data.tau+1)
                            -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                            *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                                
                } else if (data.useGPrior==0) {
                    double S_gamma = data.YtY
                            - S_Gamma_MatrixProduct;
                    logLikelihood = 
                            -((double)modelDimension/2)*Math.log(data.tau)
                            -0.5*Math.log(X_gTX_g.det())
                            -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                            *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                    
                }                
            }
        }
        
        /**
         * 
         * JAM model marginal (over beta and sigma) likelihood calculation.
         * 
         */
        
        if (data.whichLikelihoodType==LikelihoodTypes.JAM.ordinal()) {
            int[] modelDimByBlock = GeneralMethods.countPresVarsComps(
                    data.nBlocks,
                    data.blockIndices,
                    model);
            for (int b=0; b<data.nBlocks; b++) {
                if (modelDimByBlock[b]==0) {
                    /**
                     * MUCH simpler loglikelihood contribution from the block.
                     * Should be the same whether using independence or g-prior?
                     */
                    likelihoodTermBlocks[b] = 
                            -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.blockSizes[b]-1)/2)
                            *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+data.YtY_Blocks[b]);                
                } else if (modelDimByBlock[b] > 0) {
                    /**
                     * Get list of included covariate(s). Used to extract X_gamma.
                     */
                    int[] covIndsBlock = GeneralMethods.getVarIndices(
                            modelDimByBlock[b],
                            Arrays.copyOfRange(
                                    model,
                                    data.blockIndices[b],
                                    data.blockIndices[b+1])
                    );
                    /**
                     * Extract X_g and X_g'X_g.
                     */
                    Matrix X_g = data.X_Blocks[b].getMatrix(
                            0,
                            data.blockSizes[b]-1,
                            covIndsBlock);
                    Matrix X_gTX_g = data.XtX_Blocks[b].getMatrix(
                            covIndsBlock,
                            covIndsBlock);
                    if (data.useGPrior==0) {
                        /**
                         * For the independence prior we add tau to each element
                         * of the X_gamma'X_gamma diaganol before taking the
                         * inverse and determinant below.
                         */
                        for (int v=0; v<modelDimByBlock[b]; v++) {
                            X_gTX_g.set(v, v, (X_gTX_g.get(v, v)+data.tau) );
                        }
                    }
                    /**
                     * Calculate the complex matrix product in S_g.
                     */
                    S_Gamma_MatrixProductBlocks[b] = (
                            data.Y_Blocks[b].transpose().times(
                                    X_g.times(
                                            X_gTX_g.inverse()
                                    ).times(
                                            X_g.transpose().times(data.Y_Blocks[b])
                                    ))
                            ).get(0, 0);
                    
                    /**
                     * Calculate S_gamma, and the loglikelihood contribution from
                     * the block.
                     */
                    if (data.useGPrior==1) {
                        double S_gamma = data.YtY_Blocks[b]
                                - (data.tau/(1+data.tau))*S_Gamma_MatrixProductBlocks[b];
                        likelihoodTermBlocks[b] = 
                                -((double)modelDimByBlock[b]/2)*Math.log(data.tau+1)
                                -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.blockSizes[b]-1)/2)
                                *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                                
                    } else if (data.useGPrior==0) {
                        double S_gamma = data.YtY_Blocks[b]
                                - S_Gamma_MatrixProductBlocks[b];
                        likelihoodTermBlocks[b] = 
                                -((double)modelDimByBlock[b]/2)*Math.log(data.tau)
                                -0.5*Math.log(X_gTX_g.det())
                                -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.blockSizes[b]-1)/2)
                                *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                    
                    }
                }

                /**
                 * Add block contribution to likelihood total.
                 */

                logLikelihood = logLikelihood + likelihoodTermBlocks[b];
            }            
        } // End of marginal model loop

        /**
         * 
         * JAM V2 model marginal (over beta and sigma) likelihood calculation.
         * 
         */
        
        if (data.whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()) {
            int[] modelDimByBlock = GeneralMethods.countPresVarsComps(
                    data.nBlocks,
                    data.blockIndices,
                    model);
            for (int b=0; b<data.nBlocks; b++) {
                if (modelDimByBlock[b]==0) {
                    /**
                     * MUCH simpler loglikelihood contribution from the block.
                     * Should be the same whether using independence or g-prior?
                     */
                    likelihoodTermBlocks[b] = 
                            -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                            *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+data.YtY_Blocks[b]);                
                } else if (modelDimByBlock[b] > 0) {
                    /**
                     * Get list of included covariate(s). Used to extract X_gamma.
                     */
                    int[] covIndsBlock = GeneralMethods.getVarIndices(
                            modelDimByBlock[b],
                            Arrays.copyOfRange(
                                    model,
                                    data.blockIndices[b],
                                    data.blockIndices[b+1])
                    );
                    /**
                     * Extract X'y_g and X_g'X_g.
                     */
                    Matrix Xy_g = data.Y_Blocks[b].getMatrix(
                            covIndsBlock,
                            0, 0);
                    Matrix X_gTX_g = data.XtX_Blocks[b].getMatrix(
                            covIndsBlock,
                            covIndsBlock);
                    if (data.useGPrior==0) {
                        /**
                         * For the independence prior we add tau to each element
                         * of the X_gamma'X_gamma diaganol before taking the
                         * inverse and determinant below.
                         */
                        for (int v=0; v<modelDimByBlock[b]; v++) {
                            X_gTX_g.set(v, v, (X_gTX_g.get(v, v)+data.tau) );
                        }
                    }
                    /**
                     * Calculate the matrix product in S_g.
                     */
                    S_Gamma_MatrixProductBlocks[b] = (
                            Xy_g.transpose().times(
                                            X_gTX_g.inverse()
                                    ).times(Xy_g)
                            ).get(0, 0);
                    
                    /**
                     * Calculate S_gamma, and the loglikelihood contribution from
                     * the block.
                     */
                    if (data.useGPrior==1) {
                        double S_gamma = data.YtY_Blocks[b]
                                - (data.tau/(1+data.tau))*S_Gamma_MatrixProductBlocks[b];
                        likelihoodTermBlocks[b] = 
                                -((double)modelDimByBlock[b]/2)*Math.log(data.tau+1)
                                -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                                *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                                
                    } else if (data.useGPrior==0) {
                        double S_gamma = data.YtY_Blocks[b]
                                - S_Gamma_MatrixProductBlocks[b];
                        likelihoodTermBlocks[b] = 
                                -((double)modelDimByBlock[b]/2)*Math.log(data.tau)
                                -0.5*Math.log(X_gTX_g.det())
                                -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                                *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                    
                    }
                }

                /**
                 * Add block contribution to likelihood total.
                 */

                logLikelihood = logLikelihood + likelihoodTermBlocks[b];
            }            
            
        } // End of marginal model loop        
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
    public final void conjugate_calcLogLike_IncrementFromProposal(
            Arguments arguments,
            Data data,
            IterationValues curr) {
        
        
        /**
         * IPD conjugate Gaussian model likelihood calculation.
         */
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            if (whichMove < 3) {
                conjugate_calcLogLike(arguments, data);
            } else if (whichMove==3) {
                conjugate_calcLogLike(arguments, data);                
            }            
        }
        
        /**
         * JAM marginal (over beta and sigma) model likelihood calculation.
         */
        if (data.whichLikelihoodType==LikelihoodTypes.JAM.ordinal()) {
            /**
             * Get model dimension by block and re-calculate tau. This will be used
             * whether moveType is 3 (null) or less than 3 (add, delete or swap).
             */        
            int[] modelDimByBlock = GeneralMethods.countPresVarsComps(
                    data.nBlocks,
                    data.blockIndices,
                    model);

            if (whichMove<3) {
                /**
                 * Determine which blocks were updated. and update score for
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
                        if (modelDimByBlock[b]==0) {
                            /**
                             * MUCH simpler loglikelihood contribution from the block.
                             * I think this is the same for the independence and
                             * g-prior.
                             */
                            likelihoodTermBlocks[b] = 
                                    -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.blockSizes[b]-1)/2)
                                    *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+data.YtY_Blocks[b]);                
                        } else if (modelDimByBlock[b] > 0) {

                            /**
                             * Get list of included covariate(s). Used to extract X_gamma.
                             */
                            int[] covIndsBlock = GeneralMethods.getVarIndices(
                                    modelDimByBlock[b],
                                    Arrays.copyOfRange(
                                            model,
                                            data.blockIndices[b],
                                            data.blockIndices[b+1])
                            );
                            /**
                             * Extract X_g and X_g'X_g.
                             */
                            Matrix X_g = data.X_Blocks[b].getMatrix(
                                    0,
                                    data.blockSizes[b]-1,
                                    covIndsBlock);
                            Matrix X_gTX_g = data.XtX_Blocks[b].getMatrix(
                                    covIndsBlock,
                                    covIndsBlock);
                            if (data.useGPrior==0) {
                                /**
                                 * For the independence prior we add tau to each element
                                 * of X_gamma'X_gamma before taking the inverse
                                 */
                                for (int v=0; v<modelDimByBlock[b]; v++) {
                                    X_gTX_g.set(v, v, (X_gTX_g.get(v, v)+data.tau) );
                                }
                            }
                            /**
                             * Calculate the complex matrix product in S_g.
                             */
                            S_Gamma_MatrixProductBlocks[b] = (
                                    data.Y_Blocks[b].transpose().times(
                                            X_g.times(
                                                    X_gTX_g.inverse()
                                            ).times(
                                                    X_g.transpose().times(data.Y_Blocks[b])
                                            ))
                                    ).get(0, 0);
                            /**
                             * Calculate S_gamma, and the loglikelihood contribution from
                             * the block. This is all a function of tau.
                             */
                            if (data.useGPrior==1) {
                                double S_gamma = data.YtY_Blocks[b]
                                        - (data.tau/(1+data.tau))*S_Gamma_MatrixProductBlocks[b];
                                likelihoodTermBlocks[b] = 
                                        -((double)modelDimByBlock[b]/2)*Math.log(data.tau+1)
                                        -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.blockSizes[b]-1)/2)
                                        *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);
                            } else if (data.useGPrior==0) {
                                double S_gamma = data.YtY_Blocks[b]
                                        - S_Gamma_MatrixProductBlocks[b];
                                likelihoodTermBlocks[b] = 
                                        -((double)modelDimByBlock[b]/2)*Math.log(data.tau)
                                        -0.5*Math.log(X_gTX_g.det())
                                        -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.blockSizes[b]-1)/2)
                                        *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                            
                            }
                        }

                        /**
                         * Update total likelihood with the difference.
                         */
                        logLikelihood=logLikelihood+
                                (likelihoodTermBlocks[b]
                                -curr.likelihoodTermBlocks[b]);
                    }
                }
            } else if (whichMove==3) {
                /**
                 * If tau is updated must recalculate likelihood in each block. BUT:
                 * all the matrix multiplication is already stored.
                 * 
                 * NB: For blocks of 0 dimension, tau does not enter the likelihood
                 * calculation
                 */
                logLikelihood = 0;
                for (int b=0; b<data.nBlocks; b++) {
                    if (modelDimByBlock[b] > 0) {
                        /**
                         * Calculate S_gamma, and the loglikelihood contribution from
                         * the block. This is all a function of tau.
                         */
                        double S_gamma = data.YtY_Blocks[b]
                                - (data.tau/(1+data.tau))*S_Gamma_MatrixProductBlocks[b];
                        likelihoodTermBlocks[b] = 
                                -((double)modelDimByBlock[b]/2)
                                *Math.log(data.tau+1)
                                -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.blockSizes[b]-1)/2)
                                *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                                
                    }

                    /**
                     * Add block contribution to likelihood total.
                     */
                    logLikelihood = logLikelihood + likelihoodTermBlocks[b];
                }
            }            
        } // End of marginal JAM model loop
        
        /**
         * JAM V2 marginal (over beta and sigma) model likelihood calculation.
         */
        if (data.whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()) {
            /**
             * Get model dimension by block and re-calculate tau. This will be used
             * whether moveType is 3 (null) or less than 3 (add, delete or swap).
             */        
            int[] modelDimByBlock = GeneralMethods.countPresVarsComps(
                    data.nBlocks,
                    data.blockIndices,
                    model);

            if (whichMove<3) {
                /**
                 * Determine which blocks were updated. and update score for
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
                        if (modelDimByBlock[b]==0) {
                            /**
                             * MUCH simpler loglikelihood contribution from the block.
                             * I think this is the same for the independence and
                             * g-prior.
                             */
                            likelihoodTermBlocks[b] = 
                                    -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                                    *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+data.YtY_Blocks[b]);                
                        } else if (modelDimByBlock[b] > 0) {

                            /**
                             * Get list of included covariate(s). Used to extract X_gamma.
                             */
                            int[] covIndsBlock = GeneralMethods.getVarIndices(
                                    modelDimByBlock[b],
                                    Arrays.copyOfRange(
                                            model,
                                            data.blockIndices[b],
                                            data.blockIndices[b+1])
                            );
                            /**
                             * Extract X'y_g and X_g'X_g.
                             */
                            Matrix Xy_g = data.Y_Blocks[b].getMatrix(
                                    covIndsBlock,
                                    0, 0);
                            Matrix X_gTX_g = data.XtX_Blocks[b].getMatrix(
                                    covIndsBlock,
                                    covIndsBlock);
                            if (data.useGPrior==0) {
                                /**
                                 * For the independence prior we add tau to each element
                                 * of X_gamma'X_gamma before taking the inverse
                                 */
                                for (int v=0; v<modelDimByBlock[b]; v++) {
                                    X_gTX_g.set(v, v, (X_gTX_g.get(v, v)+data.tau) );
                                }
                            }
                            /**
                             * Calculate the complex matrix product in S_g.
                             */
                            S_Gamma_MatrixProductBlocks[b] = (
                                    Xy_g.transpose().times(
                                                    X_gTX_g.inverse()
                                            ).times(Xy_g)
                                    ).get(0, 0);
                            /**
                             * Calculate S_gamma, and the loglikelihood contribution from
                             * the block. This is all a function of tau.
                             */
                            if (data.useGPrior==1) {
                                double S_gamma = data.YtY_Blocks[b]
                                        - (data.tau/(1+data.tau))*S_Gamma_MatrixProductBlocks[b];
                                likelihoodTermBlocks[b] = 
                                        -((double)modelDimByBlock[b]/2)*Math.log(data.tau+1)
                                        -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                                        *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);
                            } else if (data.useGPrior==0) {
                                double S_gamma = data.YtY_Blocks[b]
                                        - S_Gamma_MatrixProductBlocks[b];
                                likelihoodTermBlocks[b] = 
                                        -((double)modelDimByBlock[b]/2)*Math.log(data.tau)
                                        -0.5*Math.log(X_gTX_g.det())
                                        -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                                        *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                            
                            }
                        }

                        /**
                         * Update total likelihood with the difference.
                         */
                        logLikelihood=logLikelihood+
                                (likelihoodTermBlocks[b]
                                -curr.likelihoodTermBlocks[b]);
                    }
                }
            } else if (whichMove==3) {
                /**
                 * If tau is updated must recalculate likelihood in each block. BUT:
                 * all the matrix multiplication is already stored.
                 * 
                 * NB: For blocks of 0 dimension, tau does not enter the likelihood
                 * calculation
                 */
                logLikelihood = 0;
                for (int b=0; b<data.nBlocks; b++) {
                    if (modelDimByBlock[b] > 0) {
                        /**
                         * Calculate S_gamma, and the loglikelihood contribution from
                         * the block. This is all a function of tau.
                         */
                        double S_gamma = data.YtY_Blocks[b]
                                - (data.tau/(1+data.tau))*S_Gamma_MatrixProductBlocks[b];
                        likelihoodTermBlocks[b] = 
                                -((double)modelDimByBlock[b]/2)
                                *Math.log(data.tau+1)
                                -((2*arguments.gaussianResidualVarianceInvGammaPrior_a+data.numberOfIndividuals-1)/2)
                                *Math.log(2*arguments.gaussianResidualVarianceInvGammaPrior_b+S_gamma);                                
                    }

                    /**
                     * Add block contribution to likelihood total.
                     */
                    logLikelihood = logLikelihood + likelihoodTermBlocks[b];
                }
            }
            
        } // End of marginal model loop        
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
    public final void conjugate_calcLogPrior(
            Arguments arguments,
            Data data,
            Priors priors) {

        /*
         * Initiate log-prior.
         */
        logPrior = 0;
        
        /***
         * Prior contribution of g-prior tau.
         ***/
        
        /**
         * NB: Because this prior is in terms of priorVarianceEstimateN, which
         * is also used for the sigma inverse gamma hyper priors and can take 
         * very large values leave it COMMENTED OUT. It
         * can result in NaN's with very varianceEstimateN's meaning the program
         * fails (this is called from the IterationValues Constructor),
         */
        logPrior = logPrior +
                Math.log(priors.betaPrecisionConjugateGammaPrior.density(
                        (double) (1/data.tau)));
    }
    
    
    
    /**
     * FOR USE WITH THE G-PRIOR: Generates a new model from the current state
     * and proposal distributions, and stores in the place of the current model.
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
    public void conjugate_update(
            Arguments arguments,
            Data data,
            IterationValues curr,
            Priors priors,
            ProposalDistributions propsds,
            Random randomDraws) {
        
        /**
         * 
         * Update the model (move type is chosen too.
         * 
         */
        updateModel(arguments, data, randomDraws);
        
        if (whichMove == 3) {
            /**
             * NULL move. Only occurs if modelTau = 1.
             */
            updateBetaPriorSd(propsds, randomDraws, data.numberOfHierarchicalCovariatePriorPartitions);
            // Recalculate prior for tau = betaPriorSd^2
            conjugate_calcLogPrior(arguments, data, priors);
        } else if (whichMove == 0) {
            betas.set(whichBetaRemoved,0,0);
        }  else if (whichMove == 1) {
            betas.set(whichBetaAdded,0,1);
        }  else if (whichMove == 2) {
            betas.set(whichBetaAdded,0,1);
            betas.set(whichBetaRemoved,0,0);
        }

        conjugate_calcLogLike_IncrementFromProposal(arguments, data, curr);
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
        if (data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
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
        } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
            /**
             * Must start at m=1. There is always a covariate present at the
             * beginning which is left fixed at 1.
             */
            for (int m=1; m<totalNumberOfCovariates; m++) {
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
        } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()) {
            /**
             * Update betas on the log-scale. Flip the sign with prob
             * 1/3.
             */
            for (int m=0; m<totalNumberOfCovariates; m++) {
                if (model[m] == 1) {
                    if (count==markUpdate) {
                        double newbeta = GeneralMaths.normalDraw(
                                    Math.log(Math.abs(betas.get(m, 0))),
                                    propsds.proposalDistributionSds[ParameterTypes.BETAS.ordinal()], r);
                        int flipSign = r.nextInt(3);
                        if (flipSign==0) {
                            // Flip sign with 1/3 chance
                            betas.set(m, 0, -Math.signum(betas.get(m, 0))*Math.exp(newbeta));                                                        
                        } else {
                            betas.set(m, 0, Math.signum(betas.get(m, 0))*Math.exp(newbeta));                                                        
                        }
                        whichBetaUpdated = m;
                    }
                    count++;
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
        logBetaPriorSds[sdUpdate] = GeneralMaths.normalDraw(logBetaPriorSds[sdUpdate],
                propsds.proposalDistributionSds[ParameterTypes.BETA_PRIOR_SD.ordinal()], r);
        betaPriorSds[sdUpdate] = Math.exp(logBetaPriorSds[sdUpdate]);
    }

    /**
     * Updates the unknown prior SD for the Dirichlet Alpha concentration parameter 
     * {@link Objects.IterationValues#betaPriorSds}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param logMinValue Specify a minimum for the log-Dirichlet concentration
     * parameter.
     * @param r Random number generator object
     */
    public void updateDirichletConcentration(
            ProposalDistributions propsds,
            double logMinValue,
            Random r) {
        whichParameterTypeUpdated = ParameterTypes.DIRICHLET_CONCENTRATION.ordinal();
        
        boolean acceptableValueDrawn = false;
        while (!acceptableValueDrawn) {
            double logDirichletConcentrationSample = GeneralMaths.normalDraw(logDirichletConcentration,
                    propsds.proposalDistributionSds[ParameterTypes.DIRICHLET_CONCENTRATION.ordinal()], r);
            if (logDirichletConcentrationSample>=logMinValue) {
                acceptableValueDrawn = true;
                logDirichletConcentration = logDirichletConcentrationSample;
                dirichletConcentration = Math.exp(logDirichletConcentration);
            }
        }
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
                    data.X.getMatrix(0, data.numberOfIndividuals-1,
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
                    data.X.getMatrix(0, data.numberOfIndividuals-1,
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
                    data.X.getMatrix(0, data.numberOfIndividuals-1,
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
            for (int i=0; i<data.numberOfIndividuals; i++) {
                alphaDelta.set(i, 0, (alpha-curr.alpha) );
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
                            .minus((data.X_Blocks[whichBetaRemovedBlock]
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
                            .plus((data.X_Blocks[whichBetaAddedBlock]
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
                            .plus((data.X_Blocks[whichBetaUpdatedBlock]
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

        /**
         * Randomly draw move according to probabilities associated with the
         * current dimension.
         */
        whichMove = GeneralMethods.chooseMove(
                probRemove, probAdd, probSwap, randomDraws);

        if (arguments.useReversibleJump == 0) {whichMove = 3;}
        
        int presentMarkN = modelDimension - numberOfCovariatesToFixInModel;
        int missMarkN = totalNumberOfCovariates - numberOfCovariatesToFixInModel - presentMarkN;
        int markRem;
        int markAdd;
        
        if (whichMove==0) {
            /**
             * If move is removal, choose a marker to remove.
             */
            markRem = randomDraws.nextInt(presentMarkN);
            int count = 0;
            if (data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
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
        } else if (whichMove==1) {
            /**
             * If move is addition, choose a marker to add.
             */
            markAdd = randomDraws.nextInt(missMarkN);
            int count = 0;
            if (data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
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
        } else if (whichMove==2) {
            /**
             * If move is swap, choose marker to remove and another to add.
             * MUST HANDLE SWAP SEPEARTELY. Or else re-add a removed marker
             * after updating the model index to 0.
             */
            markRem = randomDraws.nextInt(presentMarkN);
            markAdd = randomDraws.nextInt(missMarkN);
            int countRemove = 0;
            int countAdd = 0;
            if (data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
                /**
                 * Loop over blocks for GaussianMarg
                 */
                for (int b=0; b<data.nBlocks; b++) {
                    for (int m=data.blockIndices[b]; m<data.blockIndices[b+1]; m++) {
                        if (model[m]==1) {
                            if (countRemove==markRem) {
                                whichBetaRemoved = m;
                                whichBetaRemovedBlock = b;
                                model[m] = 0;
                                betas.set(m, 0, 0);
                            }
                            countRemove++;
                        } else if (model[m]==0) {
                            if (countAdd==markAdd) {
                                whichBetaAdded = m;
                                whichBetaAddedBlock = b;
                                model[m] = 1;
                            }
                            countAdd++;                            
                        }          
                    }
                }
            } else {
                for (int m=numberOfCovariatesToFixInModel; m<totalNumberOfCovariates; m++) {
                    if (model[m]==1) {
                        if (countRemove==markRem) {
                            whichBetaRemoved = m;
                            model[m] = 0;
                            betas.set(m, 0, 0);
                        }
                        countRemove++;
                    } else if (model[m]==0) {
                        if (countAdd==markAdd) {
                            whichBetaAdded = m;
                            model[m] = 1;
                        }
                        countAdd++;                        
                    }
                }                
            }            
        }

        modelDimension = GeneralMethods.countPresVars(totalNumberOfCovariates, model);
        modelSpacePartitionDimensions = GeneralMethods.countPresVarsComps(
                arguments.numberOfModelSpacePriorPartitions, data.modelSpacePartitionIndices, model);
        /**
         * Update move probabilities according to the new dimension. Only
         * necessary if a non-null move.
         */
        if (whichMove != 3) {
            updateMoveProbabilities(arguments);            
        }
    }
    
    /**
     * Updates the move probabilities. For use after a change (and corresponding
     * recalculation of) the model dimension.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     */
    public final void updateMoveProbabilities(
            Arguments arguments) {
        /* 
        * NB: modelDimension INCLUDES the fixed covariates
        */
        if ( (modelDimension-numberOfCovariatesToFixInModel)==0) {
            /**
             * Null model. Can only add (and null)
             */
            probRemove = 0;
            probAdd = 1-arguments.probNull;
            probSwap = 0;
        } else if (modelDimension == totalNumberOfCovariates) {
            /**
             * Saturated model. Can only remove (and null)
             */
            probRemove = 1-arguments.probNull;
            probAdd = 0;
            probSwap = 0;
        } else if ( (modelDimension-numberOfCovariatesToFixInModel) ==
                arguments.maximumModelDimension) {
            /**
             * Non-saturated but maximum dimension. Can still swap.
             */
            probRemove = arguments.probRemove+arguments.probAdd;
            probAdd = 0;
            probSwap = arguments.probSwap;
        } else {
            /**
             * All moves possible.
             */
            probRemove = arguments.probRemove;
            probAdd = arguments.probAdd;
            probSwap = arguments.probSwap;
        }
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
        Matrix alphaMat = new Matrix(data.numberOfIndividuals,1);
        for (int i=0; i<data.numberOfIndividuals; i++) {
            alphaMat.set(i, 0, alpha);
        }

        /***
         * Use the linear predictor to calculate the log-likelihood
         */
        logLikelihood = 0;        
        if (data.whichLikelihoodType==LikelihoodTypes.LOGISTIC.ordinal()) {
            // Logistic likelihood contributions
            Xbeta = data.X.times(betas);
            Xbeta.plusEquals(alphaMat);
            double p;
            for (int i=0; i<data.numberOfIndividuals; i++) {
                p = (double) (Math.exp( Xbeta.get(i, 0) )
                        / (1+Math.exp(Xbeta.get(i, 0)) ));
                if (data.binaryOutcomes[i] == 1) {
                    logLikelihood = logLikelihood + Math.log(p);
                } else {
                    logLikelihood = logLikelihood + Math.log(1-p);
                }
            }
        } else if (data.whichLikelihoodType==LikelihoodTypes.CLOGLOG.ordinal()) {
            // Complementary log-log likelihood contributions
            Xbeta = data.X.times(betas);
            Xbeta.plusEquals(alphaMat);
            double p;
            for (int i=0; i<data.numberOfIndividuals; i++) {
                p = (double) 1 - Math.exp(-Math.exp(Xbeta.get(i, 0)));
                if (data.binaryOutcomes[i] == 1) {
                    logLikelihood = logLikelihood + Math.log(p);
                } else {
                    logLikelihood = logLikelihood + Math.log(1-p);
                }
            }            
        } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
            
            // Predicted risk order, i.e. ranks
            Xbeta = data.X.times(betas);
            
            /**
             * !! Copied and pasted from the incremental likelihood function. !!
             * Works as follows:
             * 1) Gets the ordering of people by REVERSE predicted risk (ie X x Beta)
             * 2) This is the order of the points along the ROC curve; each
             *    corresponds to a risk threshold, from bottom left 
             *    (low sensitivity because v high threshold) to top right 
             *    (high sensitivity because a low threshold)
             * 3) Y is proportion of cases detected (sensitivity)
             *    X is proportion of controls called cases (FPR)
             * 4) Ideally, from origin, we would loop through all cases, and
             * and then through all controls; so y units are in terms of cases
             * and x units are in terms of controls
             * 5) As we loop through:
             *   If there is a case, we move up one unit on Y (sens increases).
             *   If there is a control we move one unit on X, FPR increase
             */
            if (data.minSensitivityOrYval==0) {
                // Regular AUC, or X-truncated calculation
                // Start at bottom left; reverse risk order
                // Want all cases then all controls
                Arrays.sort(indices, new Comparator<Integer>() {
                    @Override public int compare(final Integer o1, final Integer o2) {
                        if (Xbeta.get(o1, 0) > Xbeta.get(o2, 0)) {
                            return -1;
                        } else if (Xbeta.get(o1, 0) == Xbeta.get(o2, 0)) {
                            return 0;
                        } else {
                            return 1;                        
                        }
                    }
                });
            } else {
                // Y-truncated calculation
                // Start at top right; increasing risk order
                // Want all controls then all cases
                // x and y units are flipped below, and corresponding encoding
                Arrays.sort(indices, new Comparator<Integer>() {
                    @Override public int compare(final Integer o1, final Integer o2) {
                        if (Xbeta.get(o1, 0) < Xbeta.get(o2, 0)) {
                            return -1;
                        } else if (Xbeta.get(o1, 0) == Xbeta.get(o2, 0)) {
                            return 0;
                        } else {
                            return 1;                        
                        }
                    }
                });                
            }
            auc=0;
            xNow = 0;
            xNext = 0;
            yNow = 0;
            yNext = 0;
            int xTick = 0; // which x-tick are we at - useful for truncated AUC
            int person = 0; // which person are we at
            do {
                if (data.binaryOutcomes[indices[person]]==1) {
                    // As we reduce through risk threhsolds, and continue
                    // to find cases, we keep moving up the Y-axis
                    yNext = yNow + data.rocYunit;                    
                } else {
                    // otheriwse, if we come across a control, we nudge along 
                    // the x-axis
                    xNext = xNow + data.rocXunit;
                    xTick++;
                    if (xTick > data.lastXtickNumberForTruncatedRoc) { 
                        // Truncate when we reach the truncation point
                        xNext = data.maxFprOrXval;
                    }
                }
                auc = auc + (double) ((xNext-xNow)*(yNext+yNow)/2);
                xNow=xNext;
                yNow=yNext;
                person++;
            } while (
                    (xNow < data.maxFprOrXval) & person<data.numberOfIndividuals);
            
            
            // Loglikelihood
            logLikelihood = data.aucMultiplier*(data.maxAuc-auc);
            
            // Scale partial AUC to the same range as the regular AUC
            //auc = 0.5 + (auc - data.minAuc)*0.5/(data.maxAuc - data.minAuc);
            // Take log of AUC and scale by the logistic log-likelihood range
            //logLikelihood = data.aucMultiplier*Math.log(auc);
            
        } else if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
            // Weibull likelihood contributions
            Xbeta = data.X.times(betas);
            Xbeta.plusEquals(alphaMat);
            for (int i=0; i<data.numberOfIndividuals; i++) {
                // Log-Hazard function if event occurred at the end of follow up
                if (data.binaryOutcomes[i] == 1) {
                    logLikelihood = logLikelihood
                            + Math.log(weibullScale)
                            + weibullScale*Xbeta.get(i, 0)
                            + (weibullScale-1)*Math.log(data.survivalTimes[i]);
                }
                // Log-Survival function
                logLikelihood = logLikelihood
                        -Math.pow(
                        (double)(data.survivalTimes[i]*Math.exp( Xbeta.get(i, 0) )),
                        weibullScale);
            }            
        } else if (data.whichLikelihoodType==LikelihoodTypes.COX.ordinal()) {
            // Cox partial likelihood contributions
            Xbeta = data.X.times(betas);
            for (int i=0; i<data.numberOfIndividuals; i++) {
                if (i==0) {
                    cumulativeExpXbSumsOverRiskSets[i] = Math.exp(Xbeta.get(i, 0));
                } else {
                    cumulativeExpXbSumsOverRiskSets[i] = cumulativeExpXbSumsOverRiskSets[(i-1)] + Math.exp(Xbeta.get(i, 0));
                }
                if (data.binaryOutcomes[i] == 1) {
                    // Cox likelihood contributions, for cases
                    // Numerator: exp(Xb) for the case
                    // Denominator: sum[exp(Xb)] over anyone at risk the time the case dies
                    logLikelihood = logLikelihood
                            + Xbeta.get(i, 0) - Math.log(cumulativeExpXbSumsOverRiskSets[i]);
                }
            } 
        } else if (data.whichLikelihoodType==LikelihoodTypes.CASECOHORT_PRENTICE.ordinal()) {
            // Case-cohort pseudo likelihood contributions (0/1 Prentice weighting)
            Xbeta = data.X.times(betas);
            double eXbeta_i;            
            for (int i=0; i<data.numberOfIndividuals; i++) {
                // Data are ordered in reverse survival time - the longest
                // survior is index 0 (who must be in the sub-cohort as a 
                // survivor), at which point they are the only one in the risk
                // set.
                eXbeta_i = Math.exp(Xbeta.get(i, 0));
                if (i==0) {
                    cumulativeExpXbSumsOverRiskSets[i] = eXbeta_i;
                } else {
                    cumulativeExpXbSumsOverRiskSets[i] = cumulativeExpXbSumsOverRiskSets[(i-1)] + eXbeta_i;
                }
                if (data.binaryOutcomes[i] == 1) {
                    // Cox type likelihood contributions for cases, whether or
                    // not they are in the sub-cohort.
                    // Numerator: exp(Xb) for the case
                    // Denominator: sum[exp(Xb)] over anyone at risk the time 
                    // the case dies (including the case for a non-subcohort 
                    // case)
                    logLikelihood = logLikelihood
                            + Xbeta.get(i, 0) - Math.log(cumulativeExpXbSumsOverRiskSets[i]);
                    if (data.subcohort[i]==0) {
                        // Remove non-subcohort cases from further riskets
                        // This is the Prentice 0/1 weighting.
                        cumulativeExpXbSumsOverRiskSets[i] =cumulativeExpXbSumsOverRiskSets[i] - eXbeta_i;
                    }
                }
            }
            logLikelihood = data.casecohortPseudoLikelihoodMultiplier*logLikelihood;
        } else if (data.whichLikelihoodType==LikelihoodTypes.CASECOHORT_BARLOW.ordinal()) {
            // Barlow weighting tries to estimate denominator if whole cohort was used
            Xbeta = data.X.times(betas);
            double eXbeta_i;
            for (int i=0; i<data.numberOfIndividuals; i++) {
                // Data are ordered in reverse survival time - the longest
                // survior is index 0 (who must be in the sub-cohort as a 
                // survivor), at which point they are the only one in the risk
                // set.
                eXbeta_i = Math.exp(Xbeta.get(i, 0));
                if (i==0) {
                    cumulativeExpXbSumsOverRiskSets[i] = eXbeta_i;
                } else {
                    cumulativeExpXbSumsOverRiskSets[i] = cumulativeExpXbSumsOverRiskSets[(i-1)] + eXbeta_i;
                }
                if (data.binaryOutcomes[i] == 1) {
                    // Cox type likelihood contributions for cases, whether or
                    // not they are in the sub-cohort.
                    // Numerator: exp(Xb) for the case
                    // Denominator: sum[exp(Xb)] over anyone at risk the time 
                    cumulativeExpXbSumsOverRiskSets[i] -= eXbeta_i; // Temporarily remove from risk set
                    logLikelihood = logLikelihood
                            + Xbeta.get(i, 0) - Math.log(cumulativeExpXbSumsOverRiskSets[i]*data.barlowMultiplier + eXbeta_i);
                    if (data.subcohort[i]==1) {
                        // Add subcohort cases back for further riskets
                        cumulativeExpXbSumsOverRiskSets[i] += eXbeta_i;
                    }
                }
            }
            logLikelihood = data.casecohortPseudoLikelihoodMultiplier*logLikelihood;            
        } else if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()) {
            // Gaussian likelihood contibution (use Weibullk as sigma
            Xbeta = data.X.times(betas);
            Xbeta.plusEquals(alphaMat);
            yMinXbeta = data.Y.minus(Xbeta);            
            logLikelihood = logLikelihood
                    - (double) ((yMinXbeta.transpose().times(yMinXbeta)).get(0, 0)/(2*gaussianResidual*gaussianResidual))
                    - data.numberOfIndividuals*logGaussianResidual;
        } else if (data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
            /**
             * Calculate the likelihood score term within each block, summing
             * during the loop
             */
            likelihoodTerm = 0;
            for (int b=0; b<data.nBlocks; b++) {
                XbetaBlocks[b] = data.X_Blocks[b].times(betas.getMatrix(
                                data.blockIndices[b],
                                (data.blockIndices[(b+1)]-1), 0, 0));                
                yMinXbeta = data.Y.getMatrix(
                        data.blockIndices[b],
                        (data.blockIndices[(b+1)]-1),0,0)
                        .minus(XbetaBlocks[b]);            
                likelihoodTermBlocks[b] = (yMinXbeta.transpose()
                        .times(data.InverseCovarianceMatrix_Blocks[b])
                        .times(yMinXbeta)).get(0,0);
                likelihoodTerm=likelihoodTerm+likelihoodTermBlocks[b];
            }
            logLikelihood = logLikelihood
                    - (double) (likelihoodTerm/(2*gaussianResidual*gaussianResidual))
                    - (double) (data.totalNumberOfCovariates*logGaussianResidual);  
        }
        System.out.println("Initial log-likelihood "+logLikelihood);
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
                if (data.binaryOutcomes[i] == 1) {
                    logLikelihood = logLikelihood + Math.log(p);
                } else {
                    logLikelihood = logLikelihood + Math.log(1-p);
                }
            }
        }

        /**
         * Complementary log-log regression
         */
        if (data.whichLikelihoodType==LikelihoodTypes.CLOGLOG.ordinal()) {
            // Complementary log-log likelihood contributions
            updateXBetaEfficiently(arguments, data, curr);
            double p;
            for (int i=0; i<data.numberOfIndividuals; i++) {
                p = (double) 1 - Math.exp(-Math.exp(Xbeta.get(i, 0)));
                if (data.binaryOutcomes[i] == 1) {
                    logLikelihood = logLikelihood + Math.log(p);
                } else {
                    logLikelihood = logLikelihood + Math.log(1-p);
                }
            }
        }
        
        if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()|
                data.whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {

            
            if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()) {
                Xbeta = data.X.times(betas); // Cannot update efficiently - one beta changes, they all change
                // Can probably write a new efficient updater
            } else if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
                updateXBetaEfficiently(arguments, data, curr);
            }
            
            /**
             * Works as follows.
             * 1) Gets the ordering of people by REVERSE predicted risk (ie X x Beta)
             * 2) This is the order of the points along the ROC curve; each
             *    corresponds to a risk threshold, from bottom left 
             *    (low sensitivity because v high threshold) to top right 
             *    (high sensitivity because a low threshold)
             * 3) Y is proportion of cases detected (sensitivity)
             *    X is proportion of controls called cases (FPR)
             * 4) Ideally, from origin, we would loop through all cases, and
             * and then through all controls; so y units are in terms of cases
             * and x units are in terms of controls
             * 5) As we loop through:
             *   If there is a case, we move up one unit on Y (sens increases).
             *   If there is a control we move one unit on X, FPR increase
             */
            if (data.minSensitivityOrYval==0) {
                // Regular AUC, or X-truncated calculation
                // Start at bottom left; reverse risk order
                // Want all cases then all controls
                Arrays.sort(indices, new Comparator<Integer>() {
                    @Override public int compare(final Integer o1, final Integer o2) {
                        if (Xbeta.get(o1, 0) > Xbeta.get(o2, 0)) {
                            return -1;
                        } else if (Xbeta.get(o1, 0) == Xbeta.get(o2, 0)) {
                            return 0;
                        } else {
                            return 1;                        
                        }
                    }
                });
            } else {
                // Y-truncated calculation
                // Start at top right; increasing risk order
                // Want all controls then all cases
                // x and y units are flipped below, and corresponding encoding
                Arrays.sort(indices, new Comparator<Integer>() {
                    @Override public int compare(final Integer o1, final Integer o2) {
                        if (Xbeta.get(o1, 0) < Xbeta.get(o2, 0)) {
                            return -1;
                        } else if (Xbeta.get(o1, 0) == Xbeta.get(o2, 0)) {
                            return 0;
                        } else {
                            return 1;                        
                        }
                    }
                });                
            }
            
            // Calculate ROC AUC
            auc=0;
            xNow = 0;
            xNext = 0;
            yNow = 0;
            yNext = 0;
            int xTick = 0; // which x-tick are we at - useful for truncated AUC
            int person = 0; // which person are we at
            do {
                if (data.binaryOutcomes[indices[person]]==1) {
                    // As we reduce through risk threhsolds, and continue
                    // to find cases, we keep moving up the Y-axis
                    yNext = yNow + data.rocYunit;                    
                } else {
                    // otheriwse, if we come across a control, we nudge along 
                    // the x-axis
                    xNext = xNow + data.rocXunit;
                    xTick++;
                    if (xTick > data.lastXtickNumberForTruncatedRoc) { 
                        // Truncate when we reach the truncation point
                        xNext = data.maxFprOrXval;
                    }
                }
                auc = auc + (double) ((xNext-xNow)*(yNext+yNow)/2);
                xNow=xNext;
                yNow=yNext;
                person++;
            } while (xNow < data.maxFprOrXval& person<data.numberOfIndividuals);
                        
            // Loglikelihood
            logLikelihood = data.aucMultiplier*(data.maxAuc-auc);
            
            // Scale partial AUC to the same range as the regular AUC
            //auc = 0.5 + (auc - data.minAuc)*0.5/(data.maxAuc - data.minAuc);
            // Take log of AUC and scale by the logistic log-likelihood range
            //logLikelihood = data.aucMultiplier*Math.log(auc);
        }                
        
        /**
         * Weibull model
         */
        if (data.whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
            // Weibull likelihood contributions
            updateXBetaEfficiently(arguments, data, curr);
            for (int i=0; i<data.numberOfIndividuals; i++) {
                logLikelihood = logLikelihood 
                        -(Math.exp( Xbeta.get(i, 0) )*Math.pow(data.survivalTimes[i], weibullScale));
                if (data.binaryOutcomes[i] == 1) {
                    logLikelihood = logLikelihood
                            + Math.log(weibullScale)
                            + ((weibullScale-1)*Math.log(data.survivalTimes[i]))
                            + Xbeta.get(i, 0);
                }
            }
        }

        /**
         * Cox model
         */
        if (data.whichLikelihoodType==LikelihoodTypes.COX.ordinal()) {
            updateXBetaEfficiently(arguments, data, curr);
            for (int i=0; i<data.numberOfIndividuals; i++) {
                if (i==0) {
                    cumulativeExpXbSumsOverRiskSets[i] = Math.exp(Xbeta.get(i, 0));
                } else {
                    cumulativeExpXbSumsOverRiskSets[i] = cumulativeExpXbSumsOverRiskSets[(i-1)] + Math.exp(Xbeta.get(i, 0));
                }
                if (data.binaryOutcomes[i] == 1) {
                    // Cox likelihood contributions, for cases
                    // Numerator: exp(Xb) for the case
                    // Denominator: sum[exp(Xb)] over anyone at risk the time the case dies
                    logLikelihood = logLikelihood
                            + Xbeta.get(i, 0) - Math.log(cumulativeExpXbSumsOverRiskSets[i]);
                }
            }            
        }
        
        /**
         * Case-Cohort model with Prentice 0/1 weighting
         */
        if (data.whichLikelihoodType==LikelihoodTypes.CASECOHORT_PRENTICE.ordinal()) {
            updateXBetaEfficiently(arguments, data, curr);
            double eXbeta_i;
            for (int i=0; i<data.numberOfIndividuals; i++) {
                // Data are ordered in reverse survival time - the longest
                // survior is index 0 (who must be in the sub-cohort as a 
                // survivor), at which point they are the only one in the risk
                // set.
                eXbeta_i = Math.exp(Xbeta.get(i, 0));
                if (i==0) {
                    cumulativeExpXbSumsOverRiskSets[i] = eXbeta_i;
                } else {
                    cumulativeExpXbSumsOverRiskSets[i] = cumulativeExpXbSumsOverRiskSets[(i-1)] + eXbeta_i;
                }
                if (data.binaryOutcomes[i] == 1) {
                    // Cox type likelihood contributions for cases, whether or
                    // not they are in the sub-cohort.
                    // Numerator: exp(Xb) for the case
                    // Denominator: sum[exp(Xb)] over anyone at risk the time 
                    // the case dies (including the case for a non-subcohort 
                    // case)
                    logLikelihood = logLikelihood
                            + Xbeta.get(i, 0) - Math.log(cumulativeExpXbSumsOverRiskSets[i]);
                    if (data.subcohort[i]==0) {
                        // Remove non-subcohort cases from further riskets
                        cumulativeExpXbSumsOverRiskSets[i] = cumulativeExpXbSumsOverRiskSets[i] - eXbeta_i;
                    }
                }
            }
            logLikelihood = data.casecohortPseudoLikelihoodMultiplier*logLikelihood;            
       }

        /**
         * Case-Cohort model with Barlow weighting
         */
        if (data.whichLikelihoodType==LikelihoodTypes.CASECOHORT_BARLOW.ordinal()) {
            updateXBetaEfficiently(arguments, data, curr);            
            double eXbeta_i;
            for (int i=0; i<data.numberOfIndividuals; i++) {
                // Data are ordered in reverse survival time - the longest
                // survior is index 0 (who must be in the sub-cohort as a 
                // survivor), at which point they are the only one in the risk
                // set.
                eXbeta_i = Math.exp(Xbeta.get(i, 0));
                if (i==0) {
                    cumulativeExpXbSumsOverRiskSets[i] = eXbeta_i;
                } else {
                    cumulativeExpXbSumsOverRiskSets[i] = cumulativeExpXbSumsOverRiskSets[(i-1)] + eXbeta_i;
                }
                if (data.binaryOutcomes[i] == 1) {
                    // Cox type likelihood contributions for cases, whether or
                    // not they are in the sub-cohort.
                    // Numerator: exp(Xb) for the case
                    // Denominator: sum[exp(Xb)] over anyone at risk the time 
                    cumulativeExpXbSumsOverRiskSets[i] -= eXbeta_i; // Temporarily remove from risk set
                    logLikelihood = logLikelihood
                            + Xbeta.get(i, 0) - Math.log(cumulativeExpXbSumsOverRiskSets[i]*data.barlowMultiplier + eXbeta_i);
                    if (data.subcohort[i]==1) {
                        // Add subcohort cases back for further riskets
                        cumulativeExpXbSumsOverRiskSets[i] += eXbeta_i;
                    }
                }
            }
            logLikelihood = data.casecohortPseudoLikelihoodMultiplier*logLikelihood;            
       }
        
        /**
         * Gaussian linear regression
         */
        if (data.whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()) {
            updateXBetaEfficiently(arguments, data, curr);
            yMinXbeta = data.Y.minus(Xbeta);            
            logLikelihood = logLikelihood
                    -(double) ((yMinXbeta.transpose().times(yMinXbeta)).get(0, 0)/(2*gaussianResidual*gaussianResidual))
                    -data.numberOfIndividuals*logGaussianResidual;
        }
        
        /**
         * Regression of Gaussian marginal statistics
         */
        if (data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
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
                        yMinXbeta = data.Y.getMatrix(
                                data.blockIndices[b],
                                (data.blockIndices[(b+1)]-1),
                                0,
                                0).minus(XbetaBlocks[b]) ;            
                        likelihoodTermBlocks[b] = (yMinXbeta.transpose()
                                .times(data.InverseCovarianceMatrix_Blocks[b])
                                .times(yMinXbeta)
                                ).get(0, 0);
                        // Update total score with the difference
                        likelihoodTerm=likelihoodTerm+
                                (likelihoodTermBlocks[b]
                                -curr.likelihoodTermBlocks[b]);
                    }
                }                
            }
            /**
             * Calculate the log-likelihood, using the score term from above
             * Note that this covers the case for a gaussianResidual update
             **/
            logLikelihood = logLikelihood
                    - (double) (likelihoodTerm/(2*gaussianResidual*gaussianResidual))
                    - (double) (data.totalNumberOfCovariates*logGaussianResidual);            
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
                data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
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
        for (int v=data.numberOfCovariatesWithInformativePriors; v<data.totalNumberOfCovariates; v++) {
            if (model[v]==1) {
                logPrior = logPrior
                        + GeneralMaths.logNormDens(
                        betas.get(v, 0),
                        data.betaPriorMus[v] ,
                        betaPriorSds[data.hierarchicalCovariatePriorPartitionPicker[v]] );
            }
        }
        
        /**
         * Dirichlet prior contribution for the absolute betas
         */
        if (data.whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()) {
            /**
             * Dirichlet concentration parameter is the sum of all the alphas.
             * Setting every alpha to the concentration parameter divided by
             * modelDimension ensures this is the case, and allows the same
             * interpretation when the dimension changes
             */
            for (int v = 0; v < data.totalNumberOfCovariates; v++) {
                if (model[v]==1) {
                    logPrior += (dirichletConcentration-1)*Math.log(
                            Math.abs(betas.get(v, 0)));
                }
            }
            logPrior += Math.log(priors.dirichletConcentrationGammaPrior.density(dirichletConcentration)); //Inverse-Gamma
        }
        
        /***
         * Covariate prior precision(s)
         ***/
        if (data.numberOfHierarchicalCovariatePriorPartitions > 0) {
            for (int c=0; c<data.numberOfHierarchicalCovariatePriorPartitions; c++) {
                if (data.hierarchicalCovariatePriorPartitionFamilies[c]==HierarchicalCovariatePriorTypes.UNIFORM.ordinal()) {
                    logPrior = logPrior +
                            Math.log(priors.hierarchicalCovariatePriorSd_UniformPrior[c].density(betaPriorSds[c]));
                } else if (data.hierarchicalCovariatePriorPartitionFamilies[c]==HierarchicalCovariatePriorTypes.GAMMA.ordinal()) {
                    logPrior = logPrior +
                            Math.log(priors.hierarchicalCovariatePriorPrecision_GammaPrior[c].density(
                                    (double) (1/(betaPriorSds[c]*betaPriorSds[c])) ));                    
                } else if (data.hierarchicalCovariatePriorPartitionFamilies[c]==HierarchicalCovariatePriorTypes.JEFFREYS.ordinal()) {
                    // Jeffreys prior: P(sigma^2)=1/(sigma^2)
                    logPrior = logPrior - (double) (2*Math.log(betaPriorSds[c]));                    
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
                data.whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
            if (arguments.gaussianResidualPriorFamily==0) {
                // Jeffrey's prior: P(sigma^2)=1/(sigma^2)
                logPrior = logPrior - (double) (2*Math.log(gaussianResidual) );
            } else if (arguments.gaussianResidualPriorFamily==1) {
                // Uniform prior
                logPrior = logPrior +
                        Math.log(priors.gaussianResidualUniformPrior.density(gaussianResidual));
            } else if (arguments.gaussianResidualPriorFamily==2) {
                // Gamma on the Gaussian precision, i.e. inverse gamma on the
                // Gaussian variance.
                logPrior = logPrior
                        + Math.log(priors.gaussianPrecisionGammaPrior.density(
                                (double) (1/(gaussianResidual*gaussianResidual)) ));   
            }
        }
        
    }
    
    /**
     * Calculates the prior ratio of the proposed vs current model. For use in
     * the Reversible Jump acceptance probability calculation.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param curr {@link Objects.IterationValues} class object containing the current
     * model and parameter values
     */
    public void calcLogModelPriorRatio(
            Arguments arguments,
            Data data,
            IterationValues curr
            ) {
            // Model Space prior is accounted for, which is not accounted
            // for in method 'calcLogPrior'. A truncated Poisson Prior for
            // is used for the model space.
            logModelPriorRatio = 0;
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
                            GeneralMethods.logModelPriorRatioBetaBin_Bottolo_SpecificModels(
                            curr.modelSpacePartitionDimensions[c],
                            modelSpacePartitionDimensions[c],
                            data.modelSpacePartitionSizes[c],
                            arguments.modelSpaceBetaBinomialPriorHyperparameterA[c],
                            arguments.modelSpaceBetaBinomialPriorHyperparameterB[c]);
                }
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
                if (data.whichLikelihoodType!=
                        LikelihoodTypes.JAM.ordinal()&
                        data.whichLikelihoodType!=
                        LikelihoodTypes.JAMv2.ordinal()&
                        data.whichLikelihoodType!=
                        LikelihoodTypes.GAUSSIAN_CONJ.ordinal()
                        ) {
                    // Only do this when not intergating betas out
                    logNumerator = logNumerator +
                            GeneralMaths.logNormDens(curr.betas.get(whichBetaRemoved, 0),
                            0, propsds.proposalDistributionSds[ParameterTypes.BETA_ADD.ordinal()]);                    
                }
                logNumerator = logNumerator
                        +Math.log(probAdd) // Prob prop -> curr with new prob
                        -Math.log(data.totalNumberOfCovariates-modelDimension);
                logDenominator = logDenominator
                        +Math.log(curr.probRemove) // Prob curr -> prob with old prob
                        - Math.log(modelDimension-data.numberOfCovariatesToFixInModel +1);
            }
            if (whichMove==1) {
                // If addition happened logDenominator needs added to it the
                // log-probability of adding the added marker, and drawing
                // the corresponding new log-OR
                if (data.whichLikelihoodType!=
                        LikelihoodTypes.JAM.ordinal()&
                        data.whichLikelihoodType!=
                        LikelihoodTypes.JAMv2.ordinal()&
                        data.whichLikelihoodType!=
                        LikelihoodTypes.GAUSSIAN_CONJ.ordinal()
                        ) {
                    // Only do this when not intergating betas out
                    logDenominator = logDenominator +
                            GeneralMaths.logNormDens(betas.get(whichBetaAdded, 0), 0,
                            propsds.proposalDistributionSds[ParameterTypes.BETA_ADD.ordinal()]);                    
                }
                logNumerator = logNumerator
                        +Math.log(probRemove) // Prob prop -> curr under new removal prob
                        - Math.log(modelDimension-data.numberOfCovariatesToFixInModel);
                logDenominator = logDenominator
                        + Math.log(curr.probAdd) // Prob curr -> prop under old addition
                        -Math.log(data.totalNumberOfCovariates-modelDimension+1);
            }
            if (data.whichLikelihoodType!=
                    LikelihoodTypes.JAM.ordinal()&
                    data.whichLikelihoodType!=
                    LikelihoodTypes.JAMv2.ordinal()&
                    data.whichLikelihoodType!=
                    LikelihoodTypes.GAUSSIAN_CONJ.ordinal()
                    ) {
                /**
                 * ONLY consists of normal density calculations. So do not need
                 * to do when integrating out the betas
                 */
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
            }
            
            calcLogModelPriorRatio(arguments,data,curr);
            acceptanceProbability = Math.min(1,
                    Math.exp(logNumerator
                    -logDenominator+logModelPriorRatio));
    }    
}
