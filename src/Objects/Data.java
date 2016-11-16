package Objects;

import Jama.Matrix;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * A class to contain all data to be analysed.
 * 
 * @author Paul J Newcombe
 */
public class Data {
    
    /**
     * Repeatedly used to read into.
     */
    private String dataNameInFile;
    
    /**
     * 
     * Basic data. These variables contain the data, and key aspects such as the
     * number of covariates and individuals.
     * 
     */
    
    /**
     * Number of individuals in the dataset.
     */
    public int numberOfIndividuals;
    /**
     * Total number of covariates.
     */
    public int totalNumberOfCovariates;
    /**
     * String of variable names (to include in the results file).
     */
    public String[] covariateNames;
    /**
     * Data matrix in the JAMA format.
     */
    public Matrix X;
    /**
     * Data matrix in the JAMA format, in LD blocks.
     */
    public Matrix[] X_Blocks;
    /**
     * Vector of binary outcomes for each individual.
     */
    public int[] binaryOutcomes;
    /**
     * Quantities required for the ROC AUC calculation
     */
    public int nCases;
    public int nControls;
    public double prevalence;
    public double aucMultiplier;
    public double rocXunit;
    public double rocYunit;
    public double maxFprOrXval;
    public int lastXtickNumberForTruncatedRoc;
    public double minSensitivityOrYval;
    public double minAuc;
    public double maxAuc;
    
    /**
     * Vector of continuous binaryOutcomes for each individual.
     */
    public Matrix Y;
    /**
     * Vector of survival survivalTimes for each individual (if it is survival data).
     */
    public double[] survivalTimes;    
    /**
     * Vector of sub-cohort membership for each individual (if case-cohort data).
     */
    public int[] subcohort;
    public double barlowMultiplier;
    public double casecohortPseudoLikelihoodMultiplier;
    
    /**
     * 
     * User options. These are options set by the user, for example how many
     * covariates should be fixed in the model.
     * 
     */
    
    /**
     * Likelihood type. Indicates the likelihood according to (@link Objects.LikelihoodTypes).
     */
    public int whichLikelihoodType;
    /**
     * Number of covariates to fix in the model (eg confounders).
     */
    public int numberOfCovariatesToFixInModel;    
    /**
     * User provided initial model option. 0: Null model, 1: Saturated model
     * 2: Provided initial model
     */
    public int initialModelOption;
    /**
     * User provided initial model. Used if option 2 set for initialModelOption.
     */
    public int[] initialModel;
    /**
     * How many extra parameters beyond linear predictor. (e.g. 1 for Weibull).
     */
    public int nExtraParametersBeyondLinPred;
    
    /**
     * 
     * Model space partitions. These variables describe the
     * model space, including how it is partitioned.
     * 
     */
    
    /**
     * Number of covariates in each model space component.
     */
    public int[] modelSpacePartitionSizes;
    /**
     * Covariate indices to partition among the different model space components.
     */
    public int[] modelSpacePartitionIndices;
    /**
     * Vector of prior means for the different model space components.
     */
    public double[] modelSpacePoissonPriorMeans;
    
    /**
     * 
     * Hierarchical covariate prior partitions. Information defining the
     * partitions of covariates (i.e. betas) each of which is ascribed a common
     * prior with unknown precision.
     * 
     * NB: These partitions can be different to the modelSpacePartitions, which
     * relate instead to priors on dimension.
     * 
     */
    
    /**
     * Number of covariate partitions that will use a common prior with unknown
     * SD.
     */
    public int numberOfHierarchicalCovariatePriorPartitions;
    /**
     * Integer vector containing an index for which hierarchical covariate
     * prior partition each covariate belongs to.
     */
    public int[] hierarchicalCovariatePriorPartitionPicker;

    /**
     * Integer vector containing an index for the prior family of each
     * hierarchical covariate prior (1=Uniform, 2=Gamma).
     */
    public int[] hierarchicalCovariatePriorPartitionFamilies;
    
    /**
     * Vector of doubles containing the first Uniform hyper-parameter of the
     * beta SD prior for each hierarchical partition.
     */
    public double[] hierarchicalCovariatePriorSd_UniformHyperparameter1;
    
    /**
     * Vector of doubles containing the second Uniform hyper-parameter of the
     * beta SD prior for each hierarchical partition.
     */
    public double[] hierarchicalCovariatePriorSd_UniformHyperparameter2;

    /**
     * Vector of doubles containing the first Gamma hyper-parameter of the
     * beta SD prior for each hierarchical partition.
     */
    public double[] hierarchicalCovariatePriorPrecision_GammaHyperparameter1;

    /**
     * Vector of doubles containing the second Gamma hyper-parameter of the
     * beta SD prior for each hierarchical partition.
     */
    public double[] hierarchicalCovariatePriorPrecision_GammaHyperparameter2;
    
    /**
     * Initial values for the unknown prior SDs common to all betas in a
     * model space component {@link Objects.IterationValues#logBetaPriorSd}.
     */
    public double[] initialBetaPriorSds; // init SD value of beta hyper prior
    
    /**
     * 
     * Covariate prior setup. Information defining the covariate prior setup,
     * such as the prior mean and standard deviations for the beta covariates.
     * 
     */
    
    /**
     * Number of covariates with informative priors provided. Counted from the
     * beginning.
     */
    public int numberOfCovariatesWithInformativePriors;
    /**
     * Vector of normal prior means for the different betas.
     */
    public double[] betaPriorMus;         // Vector of prior means
    /**
     * Vector of normal prior SDs for the different betas.
     */
    public double[] betaPriorSds;         // Vector of prior SDs
    /**
     * Vector of Dirichlet Alphas for the ROCAUC likelihood.
     */
    public double[] priorDirichletWeights;         // Vector of prior means
    
    /**
     * 
     * Information describing block independence of the covariates. This is for
     * the marginal statistics models.
     * 
     */
    
    /**
     * Number of blocks (for Gaussian marginal).
     */
    public int nBlocks;         
    /**
     * Vector of block indices (for Gaussian marginal).
     */
    public int[] blockIndices;         
    /**
     * Vector of block sizes (for Gaussian marginal).
     */
    public int[] blockSizes;         
    /**
     * Vector of cumulative block sizes (for Gaussian marginal).
     */
    public int[] cumulativeBlockSizes;         
    
    /**
     * 
     * Additional quantities for marginal statistics models.
     * 
     */
    
    /**
     * Under the conjugate model: Whether to use a g-prior (1) or independent
     * priors (0).
     */
    public int useGPrior;
    /**
     * Under the conjugate model: Value to use for tau. If provide as 0, it will
     * be modeled.
     */
    public double tau;
    /**
     * Under the conjugate model: Whether or not to model tau (0/1).
     */
    public int modelTau;
    /**
     * Under the conjugate model: Whether or not to write the posterior scores
     * for all single SNP models at the top of the results file.
     */
    public int enumerateUpToDim;
    
    /**
     * 
     * Data transforms useful for likelihood calculations.
     * 
     */
    
    /**
     * List of inverted X'X matrices, by block. Useful for the (non-conjugate)
     * Guassian marginal effects analysis model.
     */
    public Matrix[] InverseCovarianceMatrix_Blocks;
    /**
     * List of X'X matrices, by block. Useful for the conjugate marginal
     * statistics model.
     */
    public Matrix[] XtX_Blocks;
    public Matrix XtX;
    /**
     * List of outcomes, by block.
     */
    public Matrix[] Y_Blocks;
    /**
     * List of (L^-1*t)'(L^-1*t), by block. Useful for the conjugate marginal
     * statistics model since these are repeatedly used.
     */
    public double[] YtY_Blocks;
    public double YtY;


    /**
     * Constructor function - populates the object by reading data from a text
     * file, whose path is stored in the {@link Objects.Arguments} object.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments, including the location of the data file.
     * @throws FileNotFoundException 
     */
    public Data(Arguments arguments) throws FileNotFoundException {
        
        /**
         * 
         * Read in basic information at the top of the data file.
         * 
         */
        Scanner dataScan = new Scanner(new File(arguments.pathToDataFile));
        String likelihoodFamily = dataScan.next();             // which likelihoodFamily type
        if (likelihoodFamily.equals("Logistic")) {
            nExtraParametersBeyondLinPred=0;
            whichLikelihoodType=LikelihoodTypes.LOGISTIC.ordinal();
        } else if (likelihoodFamily.equals("CLogLog")) {
            nExtraParametersBeyondLinPred=0;
            whichLikelihoodType=LikelihoodTypes.CLOGLOG.ordinal();
        } else if (likelihoodFamily.equals("RocAUC")) {
            nExtraParametersBeyondLinPred=1; // For beta prior SD
            whichLikelihoodType=LikelihoodTypes.ROCAUC.ordinal();
        } else if (likelihoodFamily.equals("RocAUC_Anchoring")) {
            nExtraParametersBeyondLinPred=0;
            whichLikelihoodType=LikelihoodTypes.ROCAUC_ANCHOR.ordinal();
        } else if (likelihoodFamily.equals("Weibull")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.WEIBULL.ordinal();
        } else if (likelihoodFamily.equals("Gaussian")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.GAUSSIAN.ordinal();
        } else if (likelihoodFamily.equals("GaussianConj")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.GAUSSIAN_CONJ.ordinal();
        } else if (likelihoodFamily.equals("JAM_MCMC")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.JAM_MCMC.ordinal();
        } else if (likelihoodFamily.equals("JAM")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.JAM.ordinal();            
        } else if (likelihoodFamily.equals("JAMv2")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.JAMv2.ordinal();            
        } else if (likelihoodFamily.equals("Cox")) {
            nExtraParametersBeyondLinPred=0;
            whichLikelihoodType=LikelihoodTypes.COX.ordinal();            
        } else if (likelihoodFamily.equals("CaseCohort_Prentice")) {
            nExtraParametersBeyondLinPred=0;
            whichLikelihoodType=LikelihoodTypes.CASECOHORT_PRENTICE.ordinal();            
        } else if (likelihoodFamily.equals("CaseCohort_Barlow")) {
            nExtraParametersBeyondLinPred=0;
            whichLikelihoodType=LikelihoodTypes.CASECOHORT_BARLOW.ordinal();            
        }
        
        dataNameInFile = dataScan.next();
        totalNumberOfCovariates = dataScan.nextInt();
        
        covariateNames = new String[totalNumberOfCovariates];
        for (int v=0; v<totalNumberOfCovariates; v++) {
            covariateNames[v] = dataScan.next();
        }
        
        dataNameInFile = dataScan.next();
        numberOfCovariatesToFixInModel = dataScan.nextInt();         // Number of vars which have no RJ (at begininng)
        
        dataNameInFile = dataScan.next();
        numberOfIndividuals = dataScan.nextInt();
        
        /**
         * 
         * Read in information Dirichlet prior information (for ROC model).
         * 
         */
        
        if (whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()) {
            priorDirichletWeights = new double[totalNumberOfCovariates];
            for (int v=0; v<totalNumberOfCovariates; v++) {
                priorDirichletWeights[v]  = dataScan.nextDouble();
            }
        }
        
        /**
         * 
         * Read in information regarding the covariate effect priors.
         * 
         */
        
        if (whichLikelihoodType!=LikelihoodTypes.ROCAUC.ordinal()) {
            
            // Must be initiated regardless of whether fixed
            betaPriorMus = new double[totalNumberOfCovariates];
            betaPriorSds = new double[totalNumberOfCovariates];
            
            dataNameInFile = dataScan.next(); // "NumberOfCovariatesWithInformativePriors"
            numberOfCovariatesWithInformativePriors = dataScan.nextInt();  // Number of variables from the beginning
            
            if (numberOfCovariatesWithInformativePriors>0) {
                dataNameInFile = dataScan.next(); // "FixedPriors"                
                for (int v=0; v<numberOfCovariatesWithInformativePriors; v++) {
                    betaPriorMus[v]  = dataScan.nextDouble();
                    betaPriorSds[v]  = dataScan.nextDouble();
                }
            }
        }
        
        /**
         * 
         * Write information regarding the hierarchical beta prior partitions.
         * 
         */
                
        dataNameInFile = dataScan.next();
        numberOfHierarchicalCovariatePriorPartitions = dataScan.nextInt();
                
        if (numberOfHierarchicalCovariatePriorPartitions>0) {
            // Partition picker
            hierarchicalCovariatePriorPartitionPicker = new int[totalNumberOfCovariates];
            dataNameInFile = dataScan.next();
            for (int v=numberOfCovariatesWithInformativePriors; v<totalNumberOfCovariates; v++) {
                hierarchicalCovariatePriorPartitionPicker[v]  = (dataScan.nextInt()-1);
            }
            // Partition-specific hierarchical prior parameters
            hierarchicalCovariatePriorPartitionFamilies = new int[numberOfHierarchicalCovariatePriorPartitions];
            hierarchicalCovariatePriorSd_UniformHyperparameter1 = new double[numberOfHierarchicalCovariatePriorPartitions];
            hierarchicalCovariatePriorSd_UniformHyperparameter2 = new double[numberOfHierarchicalCovariatePriorPartitions];
            hierarchicalCovariatePriorPrecision_GammaHyperparameter1 = new double[numberOfHierarchicalCovariatePriorPartitions];
            hierarchicalCovariatePriorPrecision_GammaHyperparameter2 = new double[numberOfHierarchicalCovariatePriorPartitions];
            initialBetaPriorSds = new double[numberOfHierarchicalCovariatePriorPartitions];
            // Hyper-priors
            dataNameInFile = dataScan.next();
            for (int c=0; c<numberOfHierarchicalCovariatePriorPartitions; c++) {
                String hierarchicalPriorFamily = dataScan.next();
                if (hierarchicalPriorFamily.equals("Uniform")) {
                    hierarchicalCovariatePriorPartitionFamilies[c] = HierarchicalCovariatePriorTypes.UNIFORM.ordinal();
                    hierarchicalCovariatePriorSd_UniformHyperparameter1[c] = dataScan.nextDouble();
                    hierarchicalCovariatePriorSd_UniformHyperparameter2[c] = dataScan.nextDouble();
                } else if (hierarchicalPriorFamily.equals("Gamma")) {
                    hierarchicalCovariatePriorPartitionFamilies[c] = HierarchicalCovariatePriorTypes.GAMMA.ordinal();
                    hierarchicalCovariatePriorPrecision_GammaHyperparameter1[c] = dataScan.nextDouble();
                    hierarchicalCovariatePriorPrecision_GammaHyperparameter2[c] = dataScan.nextDouble();                    
                } else if (hierarchicalPriorFamily.equals("Jeffreys")) { // No hyper-parameters
                    hierarchicalCovariatePriorPartitionFamilies[c] = HierarchicalCovariatePriorTypes.JEFFREYS.ordinal();
                }
                initialBetaPriorSds[c] = dataScan.nextDouble(); // Common for all families
            }
        }
                
        /**
         * 
         * Read in the initial model (if provided).
         * 
         */
        
        dataNameInFile = dataScan.next();
        initialModelOption = dataScan.nextInt();
        
        if (arguments.useReversibleJump==0) {
            /**
             * Saturated model forced if Reversible Jump is disabled.
             */
            initialModelOption = 1;
        }
        if (initialModelOption == 2) {
            /**
             * Read in user specified initial model.
             */
            initialModel = new int[totalNumberOfCovariates];
            dataNameInFile = dataScan.next(); //"InitalModel"
            for (int v=0; v<totalNumberOfCovariates; v++) {
                initialModel[v] = dataScan.nextInt();
            }
        }
        
        /**
         * 
         * Read in conjugate-only modeling options.
         * 
         */
        if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()|
                whichLikelihoodType==LikelihoodTypes.JAM.ordinal()|
                whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()
                ) {
            dataNameInFile = dataScan.next();
            useGPrior = dataScan.nextInt();
            
            dataNameInFile = dataScan.next();
            tau = dataScan.nextDouble();
            
            dataNameInFile = dataScan.next();
            modelTau = dataScan.nextInt();
            
            dataNameInFile = dataScan.next();
            enumerateUpToDim = dataScan.nextInt();            
        }
        if (whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()){
            dataNameInFile = dataScan.next();
            YtY = dataScan.nextDouble(); // An estimate of the trait variance is provided
        }
        
        /**
         * 
         * Read in covariate data. By LD block for marginal methods.
         * 
         */
                
        if (whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()|
                whichLikelihoodType==LikelihoodTypes.JAM.ordinal()|
                whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()
                ) {
            /**
             * Read in block information.
             */
            
            dataNameInFile = dataScan.next();
            nBlocks = dataScan.nextInt();
            
            blockIndices = new int[(nBlocks+1)];
            blockSizes = new int[nBlocks];
            cumulativeBlockSizes = new int[nBlocks];
            dataNameInFile = dataScan.next(); // "blockIndices"
            blockIndices[0] = (dataScan.nextInt()-1);
            blockIndices[0] = numberOfCovariatesToFixInModel;
            for (int b=0; b<nBlocks; b++) {
                blockIndices[(b+1)] = (dataScan.nextInt()-1);
                blockSizes[b] = blockIndices[(b+1)] - blockIndices[b];
                // Cumulative block sizes are one block down total
                // for indice adjustments
                if (b>0) {
                    cumulativeBlockSizes[b] = blockSizes[b-1]
                            +cumulativeBlockSizes[b-1];                    
                }
            }
            /**
             * Read in X'X / L^-1*X'X, block by block.
             */
            if (whichLikelihoodType==LikelihoodTypes.JAM.ordinal()|
                    whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()                    
                    ) {
                X_Blocks = new Matrix[nBlocks];
                for (int b=0; b<nBlocks; b++) {
                    X_Blocks[b] = new Matrix(blockSizes[b], blockSizes[b]);
                    for (int v1=0; v1<blockSizes[b]; v1++) {
                        for (int v2=0; v2<blockSizes[b]; v2++) {
                            X_Blocks[b].set(v1, v2,
                                    dataScan.nextDouble());
                        }
                    }
                }
            } else if (whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()) {
                XtX_Blocks = new Matrix[nBlocks];
                for (int b=0; b<nBlocks; b++) {
                    XtX_Blocks[b] = new Matrix(blockSizes[b], blockSizes[b]);
                    for (int v1=0; v1<blockSizes[b]; v1++) {
                        for (int v2=0; v2<blockSizes[b]; v2++) {
                            XtX_Blocks[b].set(v1, v2,
                                    dataScan.nextDouble());
                        }
                    }
                }                
            }
        } else {
            /**
             * IPD covariate matrix.
             */
            X = new Matrix(numberOfIndividuals,totalNumberOfCovariates);
            for (int i=0; i<numberOfIndividuals; i++) {
                for (int v=0; v<totalNumberOfCovariates; v++) {
                    X.set(i, v, dataScan.nextDouble());
                }
            }
        }
        
        /**
         * 
         * Read in outcome vector.
         * 
         */
        
        binaryOutcomes = new int[numberOfIndividuals];
        if (whichLikelihoodType==LikelihoodTypes.LOGISTIC.ordinal()|
                whichLikelihoodType==LikelihoodTypes.CLOGLOG.ordinal()|
                whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()|
                whichLikelihoodType==LikelihoodTypes.COX.ordinal()|
                whichLikelihoodType==LikelihoodTypes.CASECOHORT_PRENTICE.ordinal()|
                whichLikelihoodType==LikelihoodTypes.CASECOHORT_BARLOW.ordinal()|
                whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()|
                whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
            for (int i=0; i<numberOfIndividuals; i++) {
                binaryOutcomes[i]  = dataScan.nextInt();
            }
            /**
             * Read in survival times, if a Weibull model.
             */
            if (whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
                survivalTimes = new double[numberOfIndividuals];
                for (int i=0; i<numberOfIndividuals; i++) {
                    survivalTimes[i] = dataScan.nextDouble();
                }
            }            
            /**
             * Read in sub-cohort membership, if Case-Cohort model.
             */
            if (whichLikelihoodType==LikelihoodTypes.CASECOHORT_PRENTICE.ordinal()|
                    whichLikelihoodType==LikelihoodTypes.CASECOHORT_BARLOW.ordinal()) {
                subcohort = new int[numberOfIndividuals];
                for (int i=0; i<numberOfIndividuals; i++) {
                    subcohort[i] = dataScan.nextInt();
                }
                casecohortPseudoLikelihoodMultiplier = dataScan.nextDouble();
            }
            if (whichLikelihoodType==LikelihoodTypes.CASECOHORT_BARLOW.ordinal()) {
                barlowMultiplier = (double)(1/dataScan.nextDouble()); // Recprical of subcohort sampling fraction
            }
        } else if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            Y = new Matrix(numberOfIndividuals,1);
            for (int i=0; i<numberOfIndividuals; i++) {
                Y.set(i, 0, dataScan.nextDouble());
            }
        } else if (whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()|
                whichLikelihoodType==LikelihoodTypes.JAM.ordinal()|
                whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()
                ) {
            Y = new Matrix(totalNumberOfCovariates,1);
            for (int i=0; i<totalNumberOfCovariates; i++) {
                Y.set(i, 0, dataScan.nextDouble());
            }
        }
        
       /**
        * 
        * Store various calculations which will be repeatedly used in the
        * MCMC. These are calculated once now for efficiency. Different 
        * calculations are required for the different models.
        * 
        **/
        
        if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            /**
             * Calculate X'X once. Sub-matrices used each iteration.
             */
            XtX = (X.transpose()).times(X);            
            /**
             * Calculate (L^-1t)'(L^-1t) for each block. Repeatedly used in
             * the likelihood calculation.
             */
            YtY = ((Y.transpose()).times(Y)).get(0,0);
        } else if (whichLikelihoodType==LikelihoodTypes.JAM.ordinal()) {
            /**
             * Calculate X'X for each block. Can use sub-matrices each
             * iteration.
             */
            XtX_Blocks = new Matrix[nBlocks];
            for (int b=0; b<nBlocks; b++) {
                XtX_Blocks[b] = (X_Blocks[b].transpose()).times(X_Blocks[b]);
            }
            /**
             * Calculate (L^-1t)'(L^-1t) for each block. Repeatedly used in
             * the likelihood calculation.
             */
            Y_Blocks = new Matrix[nBlocks];
            YtY_Blocks = new double[nBlocks];
            for (int b=0; b<nBlocks; b++) {
                Y_Blocks[b] = Y.getMatrix(
                        blockIndices[b],
                        (blockIndices[b+1]-1),
                        0, 0);
                YtY_Blocks[b] = (Y_Blocks[b].transpose().times(
                                Y_Blocks[b])).get(0, 0);
            }
        } else if (whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()) {
            Y_Blocks = new Matrix[nBlocks];
            YtY_Blocks = new double[nBlocks];
            for (int b=0; b<nBlocks; b++) {
                Y_Blocks[b] = Y.getMatrix(
                        blockIndices[b],
                        (blockIndices[b+1]-1),
                        0, 0);
                YtY_Blocks[b] = YtY; // Set equal to the estimated trait variance * N
            }
        }  else if (whichLikelihoodType==LikelihoodTypes.JAM_MCMC.ordinal()) {
            /**
             * Calculate inverse of X'X for each block. These are repeatedly
             * re-used.
             */
            System.out.println("Taking inverse of xTx...");
            InverseCovarianceMatrix_Blocks = new Matrix[nBlocks];
            for (int b=0; b<nBlocks; b++) {
                InverseCovarianceMatrix_Blocks[b] = X_Blocks[b].inverse();
                System.out.println("  block "+(b+1)+"/"+nBlocks+" inverted");
            }
            System.out.println("...inverse calculated");   
        } else if (whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()|
                whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
            /**
             * Read in maximum acceptable FPR or minimum acceptable sensitvity.
             */
            maxFprOrXval = dataScan.nextDouble();
            minSensitivityOrYval = dataScan.nextDouble();
            
            /**
             * Set quantities for the AUC calculation.
             */
            // Calculate number of controls and cases; required for axis units
            nControls = 0;
            nCases = 0;
            for (int i=0; i<numberOfIndividuals; i++) {
                if (binaryOutcomes[i]==1) {
                    nCases = nCases + 1;
                } else if (binaryOutcomes[i]==0) {
                    nControls = nControls + 1;                    
                }
            }
            // Calculate X and Y tick units
            rocXunit = (double) 1/nControls; // X-axis is specificity; determined by controls
            rocYunit = (double) 1/nCases; // Y-axis is sensitivity; determined by cases  
            // For a truncated ROC curve calculate final X-tick as f(no. controls)
            lastXtickNumberForTruncatedRoc = (int) Math.round(maxFprOrXval*nControls); // Final control on truncated axis
            // If truncating by sensitivity, i.e. y-axis, it is better to
            // invert the curve, and use the same x-axis truncating code
            if (minSensitivityOrYval > 0) {
                // Re-calculate the X and Y tick units for the inverted
                // ROC curve
                rocXunit = (double) 1/nCases; 
                rocYunit = (double) 1/nControls; 
                // Flip the disease coding
                for (int i = 0; i<numberOfIndividuals; i++) {
                    binaryOutcomes[i] = 1 - binaryOutcomes[i];
                }
                // Set the X-trunctation point for the inverted ROC curve
                maxFprOrXval = 1 - minSensitivityOrYval;
                lastXtickNumberForTruncatedRoc = (int) Math.round(maxFprOrXval*nCases); // Final control on truncated axis
            }
            
            /**
             * Calculate weighting for the ROC relative to the prior.
             */
            prevalence = (double) nCases / numberOfIndividuals;
            if (maxFprOrXval < 1) {
                minAuc = (maxFprOrXval*maxFprOrXval)/2;
                maxAuc = maxFprOrXval;
            } else {
                minAuc = 0.5;
                maxAuc = 1;
            }
            aucMultiplier = (double) (arguments.aucMultiplierWeight/(maxAuc-minAuc))*(nCases*Math.log(prevalence)+nControls*Math.log(1-prevalence));
            //NEW: for logging the AUC: aucMultiplier = (double) arguments.aucMultiplierWeight * (nCases*Math.log(prevalence)+nControls*Math.log(1-prevalence))/Math.log(0.5);
        }
        
        /**
         * 
         * Setup information describing the model space. 
         * Could move this to the arguments class? 
         * Here because historically was used for beta prior partitions.
         * 
         */
        
        modelSpacePoissonPriorMeans = new double[arguments.numberOfModelSpacePriorPartitions];
        modelSpacePartitionSizes = new int[arguments.numberOfModelSpacePriorPartitions];
        modelSpacePartitionIndices = new int[(arguments.numberOfModelSpacePriorPartitions+1)];
        modelSpacePartitionIndices[0] = numberOfCovariatesToFixInModel;
        modelSpacePartitionIndices[arguments.numberOfModelSpacePriorPartitions] = totalNumberOfCovariates;
        int compSizeTotal = 0;
        if (arguments.numberOfModelSpacePriorPartitions>1) {     // Only in data if >1 components
            for (int i=0; i<(arguments.numberOfModelSpacePriorPartitions-1); i++) {  //1st & last split already done
                if (arguments.modelSpacePriorFamily==0) { // Poisson
                    modelSpacePoissonPriorMeans[i] = arguments.modelSpacePoissonPriorRate[i]
                            *arguments.modelSpacePriorPartitionSizes[i];
                } else if (arguments.modelSpacePriorFamily==1) { // Beta-binomial
                    modelSpacePartitionSizes[i] = arguments.modelSpacePriorPartitionSizes[i];                    
                }
                modelSpacePartitionIndices[i+1]=modelSpacePartitionIndices[i]
                        +arguments.modelSpacePriorPartitionSizes[i];
                compSizeTotal = compSizeTotal + arguments.modelSpacePriorPartitionSizes[i];
            }
            int finalCompSize = totalNumberOfCovariates-numberOfCovariatesToFixInModel-compSizeTotal;
            if (arguments.modelSpacePriorFamily==0) { // Poisson
                modelSpacePoissonPriorMeans[arguments.numberOfModelSpacePriorPartitions-1] = 
                        arguments.modelSpacePoissonPriorRate[arguments.numberOfModelSpacePriorPartitions-1]
                        *finalCompSize;
            } else if (arguments.modelSpacePriorFamily==1) { // Beta-binomial
                modelSpacePartitionSizes[arguments.numberOfModelSpacePriorPartitions-1] = finalCompSize;                
            }
        } else {
            if (arguments.modelSpacePriorFamily==0) { // Poisson
                modelSpacePoissonPriorMeans[0] = arguments.modelSpacePoissonPriorRate[0]*
                        (totalNumberOfCovariates-numberOfCovariatesToFixInModel);                
            } else if (arguments.modelSpacePriorFamily==1) { // Beta-binomial
                modelSpacePartitionSizes[0] = (totalNumberOfCovariates-numberOfCovariatesToFixInModel);                
            }
        }        
               
        /**
         * 
         * Feedback basic information about the dataset just read in, to the
         * console.
         * 
         */
        
        System.out.println("------------");
        System.out.println("--- DATA ---");
        System.out.println("------------");
        System.out.println("Data read from "+arguments.pathToDataFile);
        System.out.println("Likelihood: "+likelihoodFamily);
        if (whichLikelihoodType==LikelihoodTypes.ROCAUC.ordinal()|
                whichLikelihoodType==LikelihoodTypes.ROCAUC_ANCHOR.ordinal()) {
            if (minSensitivityOrYval>0) {
                System.out.println("Truncated to minimum sensitivity "+minSensitivityOrYval);                                
            } else {
                System.out.println("Truncated to maximum FPR "+maxFprOrXval);                
            }
        }
        if (
                whichLikelihoodType==LikelihoodTypes.JAM.ordinal()|
                whichLikelihoodType==LikelihoodTypes.JAMv2.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            if (useGPrior==1) {
                System.out.println("G-prior in use");
            } else if (useGPrior==0) {
                System.out.println("Independent priors in use");                
            }
            if (modelTau==0) {
                System.out.println("Tau set to "+tau);                
            } else {
                System.out.println("Tau is being modelled");                
            }
        }
        System.out.println(numberOfIndividuals+" individuals");
        System.out.println(totalNumberOfCovariates+" covariates");
        
        System.out.println("");
        System.out.println("--------------");
        System.out.println("--- PRIORS ---");
        System.out.println("--------------");
        /*
        * Covariate priors
        */
        if (numberOfCovariatesWithInformativePriors==totalNumberOfCovariates) {
            System.out.println("Fixed priors on effect sizes have been"
                    + " provided for all covariates");
        } else {
            System.out.println("Fixed effect priors are in use for "
                    +numberOfCovariatesWithInformativePriors+ " covariates");
            System.out.println("The other covariates are partitioned into "
                    + numberOfHierarchicalCovariatePriorPartitions+" group(s)"
                    + " within each of which a common prior with unknown"
                    + " variance is used for the effects");
        }
        /*
        * Model space priors
        */
        if (arguments.useReversibleJump==0) {
            // No Model selection
            System.out.println("Model selection is DISABLED - all covariates"
                    + " are fixed in the model");
        } else {
            // Model selection will be used
            System.out.println(numberOfCovariatesToFixInModel
                    +" covariates will be fixed in the model");
            if (arguments.numberOfModelSpacePriorPartitions==1) {
                if (arguments.modelSpacePriorFamily==0) {
                    System.out.println("Model selection will be performed for "
                            +(totalNumberOfCovariates-numberOfCovariatesToFixInModel)
                            + " covariates, using a Poisson prior on the"
                            + " proportion to include");
                } else if (arguments.modelSpacePriorFamily==1) {
                    System.out.println("Model selection will be performed for "
                            +(totalNumberOfCovariates-numberOfCovariatesToFixInModel)
                            + " covariates, using a Beta-binomial prior on the"
                            + " proportion to include");
                }
            } else if (arguments.numberOfModelSpacePriorPartitions>1) {
                System.out.println("Model selection will be performed for "
                        +(totalNumberOfCovariates-numberOfCovariatesToFixInModel)
                        +" covariates, split across "
                        +arguments.numberOfModelSpacePriorPartitions
                        +" partitions of differing apriori support");
            }            
        }
    }
}
