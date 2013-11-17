package Objects;

import Methods.GeneralMaths;
import Jama.Matrix;
import Methods.GeneralMethods;
import java.util.Random;

/**
 * This class contains a likelihoodFamily and parameter values, in addition to it's
 * likelihood and prior support.
 * 
 * @author Paul J Newcombe
 */
public class IterationValues {
    // PRIVATE
    /**
     * The total number of covariates.
     */
    private int totalNumberOfCovariates;
    /**
     * The total number of covariates fixed in the likelihoodFamily.
     */
    private int numberOfCovariatesToFixInModel;
    /**
     * The total number of clusters (in the case of random intercepts).
     */
    private int numberOfClusters;
    /**
     * Vector indicating which covariates are included in the likelihoodFamily
     * (length = total number of covariates being searched over).
     */
    // PUBLIC
    public int[] model;
    /**
     * Total number of covariates included in the likelihoodFamily.
     */
    public int modelDimension;
    /**
     * Total number of covariates included in each component of the likelihoodFamily space
     * (length = number of likelihoodFamily space components).
     */
    public int[] modelSpacePartitionDimensions;
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
     * Unknown common prior SD across the betas (length = number of likelihoodFamily space
     * components).
     */
    public double[] betaPriorSds;
    /**
     * Log-unknown common prior SD across the betas (length = number of likelihoodFamily space
     * components).
     */
    public double[] logBetaPriorSd;
    /**
     * Weibull scale parameter.
     */
    public double weibullScale;
    /**
     * Log-Weibull scale parameter.
     */
    public double logWeibullScale;
    /**
     * Product of the covariate matrix and coefficient vector.
     */
    public Matrix Xbeta;
    /**
     * Contains the random intercepts.
     */
    public Matrix alphaJama;
    
    // Prior, likelihood etc
    /**
     * Prior support (logged).
     */
    public double logPrior;
    /**
     * Log-Likelihood.
     */
    public double logLikelihood;
    /**
     * Acceptance probability for this likelihoodFamily and parameter values in comparison
     * to the current state.
     */
    public double acceptanceProbability;
    /**
     * Flags whether to accept the likelihoodFamily and parameters stored in this object
     * as the new likelihoodFamily and parameters.
     */
    public int proposalAccepted;
    /**
     * Indicates which parameter was updated according to the dictionary
     * (@link Objects.ParameterTypes).
     */
    public int whichParameterTypeUpdated;
    /**
     * Indicates which likelihoodFamily search move is performed at this iteration.
     * 0: Removal, 1: Addition, 2: Swap, 3: Null.
     */
    public int whichMove;
    /**
     * In the case of a null move, and beta update, indicates which beta was
     * updated.
     */
    public int whichBetaUpdated;
    /**
     * In the case of a removal move, indicates which beta was removed.
     */
    public int whichBetaRemoved;
    /**
     * In the case of an addition move, indicates which beta was added.
     */
    public int whichBetaAdded;

    /**
     * Constructor; used to initiate an analysis with a likelihoodFamily (usually the null)
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
        model = new int[data.totalNumberOfCovariates]; // indicates markers present in current likelihoodFamily
        // Parameters that are updated
        betas = new Matrix(data.totalNumberOfCovariates,1);
        clusterIntercepts = new double[data.numberOfClusters]; // log-ORs of current likelihoodFamily
        alpha = arguments.initialAlpha;
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
        logBetaPriorSd = new double[data.numberOfModelSpacePartitions];
        betaPriorSds = new double[data.numberOfModelSpacePartitions];
        for (int c=0; c<data.numberOfModelSpacePartitions; c++) {
            logBetaPriorSd[c] = Math.log(arguments.initialBetaPriorSd);
            betaPriorSds[c] = Math.exp(logBetaPriorSd[c]);
        }
        weibullScale = arguments.initialWeibullScale;
        logWeibullScale = Math.log(weibullScale);
        
        
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
        
        // Single values
        alpha = its.alpha;
        modelDimension = its.modelDimension;
        logBetweenClusterSd = its.logBetweenClusterSd;
        weibullScale = its.weibullScale;
        logWeibullScale = its.logWeibullScale;
        logLikelihood = its.logLikelihood;
        logPrior = its.logPrior;
        
        // Jama objects
        alphaJama = its.alphaJama;  // Copy does not work for some reason
        betas = its.betas.copy();
        Xbeta = its.Xbeta.copy();
    }
    
    /**
     * Generates a new likelihoodFamily and set of parameter values, from the current
     * state and proposal distributions, and stores in the place of the current
     * likelihoodFamily and parameter values.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param curr {@link Objects.IterationValues} class object containing the current
     * likelihoodFamily and parameter values
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
        
        // Start by choosing likelihoodFamily selection move
        updateModel(arguments, data, randomDraws);

        double parameterDraw = 0;
        
        // ProposalDistributions/WhichUpdated Ref:
        // 0: Alpha
        // 1: Betas (when remain in likelihoodFamily)
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
                double paramTypeDraw = randomDraws.nextInt( (4+data.numberOfModelSpacePartitions+data.survivalAnalysis) );
                if (paramTypeDraw == 0) {
                    updateAlpha(propsds, randomDraws);
                } else if (paramTypeDraw == 1) {
                    updateBetas(propsds, randomDraws);
                } else if (paramTypeDraw == 2) {
                    updateBetweenClusterSd(propsds, randomDraws);
                } else if (paramTypeDraw == 3) {
                    updateClusterIntercepts(propsds, randomDraws);
                } else if (paramTypeDraw == 4) {
                    if (data.numberOfModelSpacePartitions>0) {
                        updateBetaPriorSd(propsds, randomDraws, data.numberOfModelSpacePartitions);                        
                    } else if (data.survivalAnalysis==1) {
                        updateWeibullK(propsds, randomDraws);                        
                    }
                } else if (paramTypeDraw == 5) {
                    updateWeibullK(propsds, randomDraws);
                }
            } else {
            // NULL MOVE for when no clusters intercepts -----------------------
                double paramTypeDraw = randomDraws.nextInt((2+data.numberOfModelSpacePartitions+data.survivalAnalysis));
                if (paramTypeDraw == 0) {
                    updateAlpha(propsds, randomDraws);
                } else if (paramTypeDraw == 1) {
                    updateBetas(propsds, randomDraws);
                } else if (paramTypeDraw == 2) {
                    if (data.numberOfModelSpacePartitions>0) {
                        updateBetaPriorSd(propsds, randomDraws, data.numberOfModelSpacePartitions);                        
                    } else if (data.survivalAnalysis==1) {
                        updateWeibullK(propsds, randomDraws);                        
                    }
                } else if (paramTypeDraw == 3) {
                    updateWeibullK(propsds, randomDraws);
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
                alpha,  propsds.sds[ParameterTypes.ALPHA.ordinal()], r);
    }
    
    /**
     * Updates a currently included beta parameter in 
     * {@link Objects.IterationValues#betas}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     */
    public void updateBetas(
            ProposalDistributions propsds,
            Random r) {
        if (modelDimension>0) {
            whichParameterTypeUpdated = ParameterTypes.BETAS.ordinal();
            int markUpdate = r.nextInt(modelDimension);
            int count = 0;
            for (int m=0; m<totalNumberOfCovariates; m++) {
                if (model[m] == 1) {
                    if (count==markUpdate) {
                        double newbeta = GeneralMaths.normalDraw(
                                    betas.get(m, 0), 
                                    propsds.sds[ParameterTypes.BETAS.ordinal()], r);
                        betas.set(m, 0, newbeta);
                        whichBetaUpdated = m;
                    }
                    count++;
                }
            }                    
        } else {
            updateAlpha(propsds, r);
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
                propsds.sds[ParameterTypes.BETA_ADD.ordinal()],
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
                propsds.sds[ParameterTypes.BETA_SWAP.ordinal()],
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
                propsds.sds[ParameterTypes.BETWEEN_CLUSTER_SD.ordinal()], r);
    }
    
    /**
     * Updates the unknown prior SD common to the betas 
     * {@link Objects.IterationValues#betaPriorSds}
     * 
     * @param propsds {@link Objects.ProposalDistributions} class object, containing all
     * proposal distribution SDs
     * @param r Random number generator object
     * @param nSds The number of likelihoodFamily space components
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
                propsds.sds[ParameterTypes.BETA_PRIOR_SD.ordinal()], r);
        betaPriorSds[sdUpdate] = Math.exp(logBetaPriorSd[sdUpdate]);
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
                propsds.sds[ParameterTypes.WEIBULL_SCALE.ordinal()], r);
        weibullScale = Math.exp(logWeibullScale);
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
                    propsds.sds[ParameterTypes.CLUSTER_INTERCEPTS.ordinal()], r);
            
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
            // Count number variables included in current likelihoodFamily to determine
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

        int missMarkN = totalNumberOfCovariates - numberOfCovariatesToFixInModel - presentMarkN;
        if ((whichMove==1)||(whichMove==2)) {
            // If whichMove is addition (1) or swap (2) choose at random a marker
            // to add and place a 1 in corresponding position in presentMark
            markAdd = randomDraws.nextInt(missMarkN);
            int count = 0;
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

        modelDimension = GeneralMethods.countPresVars(totalNumberOfCovariates, model);
        modelSpacePartitionDimensions = GeneralMethods.countPresVarsComps(
                arguments.numberOfModelSpacePriorPartitions, data.modelSpacePartitionIndices, model);
        
    }

    /**
     * Calculates the log-likelihood according to the stored likelihoodFamily and parameter
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

        //beta*X
        Xbeta = data.dataJama.times(betas);

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

        Xbeta.plusEquals(alphaMat);

        logLikelihood = 0;
        
        if (data.survivalAnalysis==0) {
            // Logistic likelihood contributions
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
        } else if (data.survivalAnalysis==1) {
            // Weibull likelihood contributions
//            for (int i=0; i<data.numberOfIndividuals; i++) {
//                logLikelihood = logLikelihood -Math.exp( Xbeta.get(i, 0) )
//                        *Math.pow(data.times[i],weibullScale);
//                if (data.outcomes[i] == 1) {
//                    logLikelihood = logLikelihood
//                            + Math.log(weibullScale)
//                            + (weibullScale-1)*Math.log(data.times[i])
//                            + Xbeta.get(i, 0);
//                }
//            }            
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
        }
    }

    /**
     * Efficiently calculates the log-likelihood by modifying the previously
     * stored value according to the specific changes to the likelihoodFamily and parameter
     * values since then, and stores in {@link Objects.IterationValues#logLikelihood}.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param curr {@link Objects.IterationValues} class object containing the previous
     * likelihoodFamily and parameter values, to use as a starting point
     */
    public final void calcLogLike_IncrementFromProposal(
            Arguments arguments,
            Data data,
            IterationValues curr) {

        if (whichMove == 0) {
            Xbeta = Xbeta.minus(
                    data.dataJamaColumns[whichBetaRemoved]
                    .times(
                    curr.betas.get(whichBetaRemoved,0)
                    )
                    );
        } else if (whichMove == 1) {
            Xbeta = Xbeta.plus(
                    data.dataJamaColumns[whichBetaAdded]
                    .times(
                    (betas.get(whichBetaAdded,0)-
                    curr.betas.get(whichBetaAdded,0))
                    )
                    );
        } else if (whichMove == 2) { //Swap - 2 changes
            Xbeta = Xbeta.minus(
                    data.dataJamaColumns[whichBetaRemoved]
                    .times(
                    curr.betas.get(whichBetaRemoved,0)
                    )
                    );
            Xbeta = Xbeta.plus(
                    data.dataJamaColumns[whichBetaAdded]
                    .times(
                    (betas.get(whichBetaAdded,0)-
                    curr.betas.get(whichBetaAdded,0))
                    )
                    );
        }
        else if (whichMove == 3) { // Null
            if (whichParameterTypeUpdated == 1) { // Modify by a logOR
                Xbeta = Xbeta.plus(
                        data.dataJamaColumns[whichBetaUpdated]
                        .times(
                        (betas.get(whichBetaUpdated, 0)-
                        curr.betas.get(whichBetaUpdated, 0))
                        )
                        );
            } else {
//            else {    // Modify by alpha
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

        logLikelihood = 0;
        
        if (data.survivalAnalysis==0) {
            // Logistic likelihood contributions
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
        } else if (data.survivalAnalysis==1) {
            // Weibull likelihood contributions
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
    }

    /**
     * Calculates the log-prior support of the stored likelihoodFamily and parameter
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
         * Calculate log-prior.
         */
        logPrior = 0;

        // Global intercept
        logPrior = logPrior + GeneralMaths.logNormDens(
                alpha,
                arguments.alphaPriorMu,
                arguments.alphaPriorSd
                );        
        // Study specific intercepts
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

        // Only included variables make a normal contribution for their log-ORs
        // depending on presence of the causal variant presentMarkN-1 of these
        // may be converted to truncated normal priors, by multiplying the pdf
        // by 2
        
        // Beta contributions where fixed priors are provided
        for (int v = 0; v < data.numberOfCovariatesWithInformativePriors; v++) {
            if (model[v]==1) {
                logPrior = logPrior
                        + GeneralMaths.logNormDens(
                        betas.get(v, 0),
                        data.betaPriorMus[v] ,
                        data.betaPriorSds[v] );
            }
        }
        // Beta contributions conditional on Beta Prior SD which has its own
        // hyper prior
        for (int c=0; c<data.numberOfModelSpacePartitions; c++) {
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
        
        // Beta prior sds (if they have hyperpiors)
        if (data.numberOfModelSpacePartitions > 0) {        
            // Beta prior SD hyperparameter
            if (arguments.betweenClusterSdPriorFamily == 0) {
                // uniform prior for SD parameters
                for (int c=0; c<data.numberOfModelSpacePartitions; c++) {
                    logPrior = logPrior +
                            Math.log(priors.betweenClusterPrecisionUniformPrior.density(betaPriorSds[c]));
                }
            } else if (arguments.betweenClusterSdPriorFamily==1) {
                // My home made gamma contribution (the constant multiplier cancels)
                for (int c=0; c<data.numberOfModelSpacePartitions; c++) {
                    logPrior = logPrior
                            -2*(arguments.betaPrecisionGammaPriorHyperparameter1+1)*logBetaPriorSd[c]
                            -(double)arguments.betaPrecisionGammaPriorHyperparameter2/(betaPriorSds[c]*betaPriorSds[c]);                    
                }
            }
        }
        
        // Weibull shape parameter
        if (data.survivalAnalysis==1) {
            logPrior = logPrior
                    +priors.weibullScalePrior.density(weibullScale);
        }
    }
    
    /**
     * Calculates the Metropolis Hastings Reversible Jump acceptance probability
     * for the likelihoodFamily and parameter values in this object vs the previously
     * saved likelihoodFamily and parameter values, and stores in the field
     * {@link Objects.IterationValues#acceptanceProbability}.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments
     * @param data {@link Objects.Data} class object, containing all the data to
     * be analysed
     * @param curr {@link Objects.IterationValues} class object containing the current
     * likelihoodFamily and parameter values
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
                        0, propsds.sds[ParameterTypes.BETA_ADD.ordinal()]);
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
                        propsds.sds[ParameterTypes.BETA_ADD.ordinal()]);
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
                        0, propsds.sds[ParameterTypes.BETA_SWAP.ordinal()]);
                logDenominator = logDenominator +
                        GeneralMaths.logNormDens(betas.get(whichBetaAdded, 0),
                        0, propsds.sds[ParameterTypes.BETA_SWAP.ordinal()]);
            }

            // Model Space prior is accounted for, which is not accounted
            // for in method 'calcLogPrior'. A truncated Poisson Prior for
            // is used for the likelihoodFamily space.
            double logModelPriorRatio = 0;
            if (arguments.modelSpacePriorFamily==0) {
                // Poisson likelihoodFamily space prior
                for (int c=0; c<arguments.numberOfModelSpacePriorPartitions; c++) {
                    logModelPriorRatio = logModelPriorRatio +
                            GeneralMethods.logModelPriorRatio(
                            curr.modelSpacePartitionDimensions[c],
                            modelSpacePartitionDimensions[c],
                            data.modelSpacePoissonPriorMeans[c]);
                } 
            } else if (arguments.modelSpacePriorFamily==1) {
                // Beta-binomial likelihoodFamily space prior
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

