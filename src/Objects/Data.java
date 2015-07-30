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
     * Vector of continuous binaryOutcomes for each individual.
     */
    public Matrix Y;
    /**
     * Vector of survival survivalTimes for each individual (if it is survival data).
     */
    public double[] survivalTimes;    
    
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
     * Number of covariates in each model space component using the same
     * common prior with unknown SD.
     */
    public int[] hierarchicalCovariatePriorPartitionSizes;
    /**
     * Covariate indices which partition the covariates into the different
     * model space components using the same common prior with unknown SD over
     * the betas.
     */
    public int[] hierarchicalCovariatePriorPartitionIndices;
        
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
     * Residual variance inverse-gamma hyper-parameter 1.
     */
    public double sigma2_invGamma_a;
    /**
     * Residual variance inverse-gamma hyper-parameter 1.
     */
    public double sigma2_invGamma_b;
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
     * Under the conjugate model: Initial proposal SD for adapting tau 
     * proposals. This can be specified since is very different to the normal
     * use of betaPriorSd.
     */
    public double tauInitialProposalSd;
    /**
     * Under the conjugate model: Whether or not to write the posterior scores
     * for all single SNP models at the top of the results file.
     */
    public int allModelScoresUpToDim;
    
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
        } else if (likelihoodFamily.equals("Weibull")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.WEIBULL.ordinal();
        } else if (likelihoodFamily.equals("Gaussian")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.GAUSSIAN.ordinal();
        } else if (likelihoodFamily.equals("GaussianConj")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.GAUSSIAN_CONJ.ordinal();
        } else if (likelihoodFamily.equals("GaussianMarg")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal();
        } else if (likelihoodFamily.equals("GaussianMargConj")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.GAUSSIAN_MARGINAL_CONJ.ordinal();            
        }
        totalNumberOfCovariates = dataScan.nextInt();             // no. variables (including interaction terms)
        covariateNames = new String[totalNumberOfCovariates];
        for (int v=0; v<totalNumberOfCovariates; v++) {
            covariateNames[v] = dataScan.next();
        }
        numberOfCovariatesToFixInModel = dataScan.nextInt();         // Number of vars which have no RJ (at begininng)
        numberOfIndividuals = dataScan.nextInt();             // No. Studies
        int numberOfClustersDummy = dataScan.nextInt();             // No. clusters
        
        /**
         * 
         * Setup information describing the model space. Required below.
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
         * Read in covariate data. By LD block for marginal methods.
         * 
         */
        
        /**
         * Residual prior info for all Gaussian models.
         */
        if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL_CONJ.ordinal()
                ) {
            /**
             * Gamma prior parameters for the residual variance. Calculate the
             * hyper-parameters of an informative inverse gamma of an
             * informative inverse gamma using the inverse-gamma/inverse-chi
             * squared equivalence. Note that java uses parameterisation 2 for
             * gammas
             */
            sigma2_invGamma_a = dataScan.nextDouble();
            sigma2_invGamma_b = dataScan.nextDouble();            
        }
        
        /**
         * Modeling options for conjugate Gaussian models.
         */
        if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL_CONJ.ordinal()
                ) {
            useGPrior = dataScan.nextInt();
            tau = dataScan.nextDouble();
            modelTau = dataScan.nextInt();
            tauInitialProposalSd = dataScan.nextDouble();
            allModelScoresUpToDim = dataScan.nextInt();            
        }
        
        
        /**
         * Read in covariate information. Very different for marginal models.
         */
        if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL_CONJ.ordinal()) {
            /**
             * Read in block information.
             */
            
            nBlocks = dataScan.nextInt();
            blockIndices = new int[(nBlocks+1)];
            blockSizes = new int[nBlocks];
            cumulativeBlockSizes = new int[nBlocks];
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
        } else {
            /**
             * Full IPD regression models. No blocks.
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
        * Store various calculations which will be re-used, for efficiency.
        * 
        **/
        
        if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            /**
             * Calculate X'X once. Can use sub-matrices each
             * iteration.
             */
            XtX = (X.transpose()).times(X);            
        } else if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
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
        } else if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL_CONJ.ordinal()) {
            /**
             * Calculate X'X for each block. Can use sub-matrices each
             * iteration.
             */
            XtX_Blocks = new Matrix[nBlocks];
            for (int b=0; b<nBlocks; b++) {
                XtX_Blocks[b] = (X_Blocks[b].transpose()).times(X_Blocks[b]);
            }
         }

        /**
         * 
         * Read in outcome vector.
         * 
         */
        
        binaryOutcomes = new int[numberOfIndividuals];
        if (whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()|
                whichLikelihoodType==LikelihoodTypes.LOGISTIC.ordinal()) {
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
        } else if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            Y = new Matrix(numberOfIndividuals,1);
            for (int i=0; i<numberOfIndividuals; i++) {
                Y.set(i, 0, dataScan.nextDouble());
            }
        } else if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL_CONJ.ordinal()) {
            Y = new Matrix(totalNumberOfCovariates,1);
            for (int i=0; i<totalNumberOfCovariates; i++) {
                Y.set(i, 0, dataScan.nextDouble());
            }
        }
        
       /**
        * 
        * Store various calculations which will be re-used, for efficiency.
        * 
        **/
        
        if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            /**
             * Calculate (L^-1t)'(L^-1t) for each block. Repeatedly used in
             * the likelihood calculation.
             */
            YtY = ((Y.transpose()).times(Y)).get(0,0);
        } else if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL_CONJ.ordinal()) {
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
        }        
        
        
        /**
         * 
         * Read in information regarding the prior setup.
         * 
         */
        
        betaPriorMus = new double[totalNumberOfCovariates];
        betaPriorSds = new double[totalNumberOfCovariates];
        numberOfCovariatesWithInformativePriors = dataScan.nextInt();  // Number of variables from the beginning
        for (int v=0; v<numberOfCovariatesWithInformativePriors; v++) {
            betaPriorMus[v]  = dataScan.nextDouble();
            betaPriorSds[v]  = dataScan.nextDouble();
        }            
        numberOfHierarchicalCovariatePriorPartitions = dataScan.nextInt();  // Could be 0
        if (numberOfHierarchicalCovariatePriorPartitions>0) {
            numberOfHierarchicalCovariatePriorPartitions=1; // FORCE TO 1 FOR NOW
        }        
        hierarchicalCovariatePriorPartitionIndices = new int[(numberOfHierarchicalCovariatePriorPartitions+1)];
        hierarchicalCovariatePriorPartitionIndices[0] = numberOfCovariatesWithInformativePriors; // could be 0
        hierarchicalCovariatePriorPartitionIndices[numberOfHierarchicalCovariatePriorPartitions] = totalNumberOfCovariates;
        if (numberOfHierarchicalCovariatePriorPartitions>1) {     // Only in data if >1 components
            hierarchicalCovariatePriorPartitionSizes = new int[(numberOfHierarchicalCovariatePriorPartitions-1)];
            for (int i=0; i<(numberOfHierarchicalCovariatePriorPartitions-1); i++) {
                hierarchicalCovariatePriorPartitionSizes[i] = dataScan.nextInt(); //Model Space Poisson Means
                hierarchicalCovariatePriorPartitionIndices[i+1]=hierarchicalCovariatePriorPartitionIndices[i]+hierarchicalCovariatePriorPartitionSizes[i];
            }
        }
        
        /**
         * 
         * Initial model.
         * 
         */
        
        initialModelOption = dataScan.nextInt();
        if (arguments.useReversibleJump==0) {
            /**
             * Saturated model used if Reversible Jump is disabled.
             */
            initialModelOption = 1;
        }
        if (initialModelOption==2) {
            /**
             * Read in user specified initial model.
             */
            initialModel = new int[totalNumberOfCovariates];
            for (int v=0; v<totalNumberOfCovariates; v++) {
                initialModel[v] = dataScan.nextInt();
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
        if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL_CONJ.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_CONJ.ordinal()) {
            if (useGPrior==1) {
                System.out.println("G-prior in use.");
            } else if (useGPrior==0) {
                System.out.println("Independent priors in use.");                
            }
            if (modelTau==0) {
                System.out.println("Tau set to "+tau);                
            } else {
                System.out.println("Tau is being modelled.");                
            }
        }
        System.out.println(numberOfIndividuals+" individuals.");
        System.out.println(totalNumberOfCovariates+" covariates");
        
        System.out.println("");
        System.out.println("--------------");
        System.out.println("--- PRIORS ---");
        System.out.println("--------------");
        if (arguments.useReversibleJump==0) {
            // No Model selection
            System.out.println("Model selection is DISABLED - all covariates"
                    + " are fixed in the model");
            System.out.println("Priors have been provided on the effect of"
                    + " all covariates");
            
        } else {
            // Model selection will be used
            System.out.println(numberOfCovariatesToFixInModel
                    +" covariates will be fixed in the model");
            if (arguments.numberOfModelSpacePriorPartitions>1) {
                System.out.println("Model selection will be performed for "
                        +(totalNumberOfCovariates-numberOfCovariatesToFixInModel)
                        +" covariates, split across "
                        +arguments.numberOfModelSpacePriorPartitions
                        +" partitions of differing apriori support");
                // Model space
                if (arguments.modelSpacePriorFamily==0) {
                    System.out.println("Within each partition, Poisson priors"
                            + " will be used on the number of covariates to"
                            + " include");
                } else if (arguments.modelSpacePriorFamily==1) {
                    System.out.println("Within each partition, Beta-binomial"
                            + " priors will be used for the number of"
                            + " covariates to include");            
                }
                // Covariates
                if (numberOfCovariatesWithInformativePriors==numberOfCovariatesToFixInModel) {
                    System.out.println("Common priors with unknown SDs will be"
                            + " used across covariate effects within each model"
                            + " space partition");
                } else {
                    System.out.println("Fixed priors on effect sizes have been"
                            + " provided for all covariates");
                }
            } else if (arguments.numberOfModelSpacePriorPartitions==1) {
                System.out.println("Model selection will be performed for "
                        +(totalNumberOfCovariates-numberOfCovariatesToFixInModel)
                        +" covariates");
                // Model space
                if (arguments.modelSpacePriorFamily==0) {
                    System.out.println("A Poisson prior on the number of"
                            + " covariates to include will be used");
                } else if (arguments.modelSpacePriorFamily==1) {
                    System.out.println("A Beta-binomial prior on the number of"
                            + " covariates to include will be used");            
                }
                // Covariate effects
                if (numberOfCovariatesWithInformativePriors==numberOfCovariatesToFixInModel) {
                    System.out.println("A common prior with unknown SD will be"
                            + " used across the effects of covariates under"
                            + " model selection");
                } else {
                    System.out.println("Fixed priors on effect sizes have been"
                            + " provided for all covariates");
                }
            }
            
        }
    }
}
