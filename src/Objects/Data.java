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
     * The data. These variables contain the data, and key aspects such as the
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
     * Number of clusters if fitting random intercepts.
     */
    public int numberOfClusters;
    /**
     * Data matrix in the JAMA format.
     */
    public Matrix dataJama;
    /**
     * Data matrix in the JAMA format, in blocks.
     */
    public Matrix[] dataJamaBlocks; // for GuassianMarg
    /**
     * Vector of clusters assignments for each individual.
     */
    public int[] clusters;         // Vector of clusters assignments
    /**
     * Vector of outcomes for each individual.
     */
    public int[] outcomes;         // Binary outcomes vector (logistic model)
    /**
     * Vector of continuous outcomes for each individual.
     */
    public Matrix continuousOutcomesJama;
    /**
     * Vector of survival times for each individual (if it is survival data).
     */
    public double[] times;         // Event/death times vector (survival model)    
    
    /**
     * 
     * Information on the model. These variables contain key information on the
     * modeling setup such as which likelihood, how many covariates to fix in.
     * 
     */
    
    /**
     * Number of covariates to fix in the model (eg confounders).
     */
    public int numberOfCovariatesToFixInModel;    
    /**
     * Indicates the likelihood type according to the dictionary
     * (@link Objects.LikelihoodTypes).
     */
    public int whichLikelihoodType;
    /**
     * How many extra parameters beyond logistic regression.
     */
    public int nExtraParametersBeyondLinPred;
    /**
     * Whether to fit random intercepts (0/1).
     */
    public int useRandomIntercepts;
    
    /**
     * 
     * Information describing the model space. These variables describe the
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
     * Information on the prior setup. These variables contain information on
     * the prior setup, such as the prior mean and standard deviations for the
     * beta covariates.
     * 
     */
    
    /**
     * Number of covariates (from the beginning) which informative priors have
     * been provided for.
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
     * Number of covariate partitions that will use a common prior with unknown
     * SD.
     */
    public int numberOfUnknownBetaPriors;
    /**
     * Number of covariates in each model space component using the same
     * common prior with unknown SD.
     */
    public int[] commonBetaPriorPartitionSizes;
    /**
     * Covariate indices which partition the covariates into the different
     * model space components using the same common prior with unknown SD over
     * the betas.
     */
    public int[] commonBetaPriorPartitionIndices;
    /**
     * Inverted xTx matrix for the Guassian marginal effects analysis.
     */
    public Matrix xTxInvForGPrior;    
    
    /**
     * 
     * Information describing block independence of the covariates. This is for
     * the Gaussian marginal statistics model.
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
     * List of inverted xTx matrices, within blocks, for the Guassian marginal
     * effects analysis.
     */
    public Matrix[] xTxInvBlocks;


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
         * Read in basic information at the top of the data file
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
        } else if (likelihoodFamily.equals("GaussianMarg")) {
            nExtraParametersBeyondLinPred=1;
            whichLikelihoodType=LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal();
        }
        totalNumberOfCovariates = dataScan.nextInt();             // no. variables (including interaction terms)
        covariateNames = new String[totalNumberOfCovariates];
        for (int v=0; v<totalNumberOfCovariates; v++) {
            covariateNames[v] = dataScan.next();
        }
        numberOfCovariatesToFixInModel = dataScan.nextInt();         // Number of vars which have no RJ (at begininng)
        numberOfIndividuals = dataScan.nextInt();             // No. Studies
        numberOfClusters = dataScan.nextInt();             // No. clusters
        if (numberOfClusters>0) {
            useRandomIntercepts = 1;
        } else {
            useRandomIntercepts = 0;
        }
        
        /**
         * Setup basic information describing the model space - do now, since
         * needed below
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
         * Read in covariate data
         */
        if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
            /**
             * Read in block indices
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
                if (b>0) {
                    cumulativeBlockSizes[b] = blockSizes[b]
                            +cumulativeBlockSizes[b-1];                    
                }
            }
            /**
             * Read in covariate data, block by block
             */
            dataJamaBlocks = new Matrix[nBlocks];
            for (int b=0; b<nBlocks; b++) {
                dataJamaBlocks[b] = new Matrix(blockSizes[b], blockSizes[b]);
                for (int v1=0; v1<blockSizes[b]; v1++) {
                    for (int v2=0; v2<blockSizes[b]; v2++) {
                        dataJamaBlocks[b].set(v1, v2,
                                dataScan.nextDouble());
                    }
                }
            }            
            /**
             * Calculate inverse of xTx, block by block
             */
            System.out.println("Taking inverse of xTx...");
            xTxInvBlocks = new Matrix[nBlocks];
            for (int b=0; b<nBlocks; b++) {
                xTxInvBlocks[b] = dataJamaBlocks[b].inverse();
                System.out.println("  block "+(b+1)+"/"+nBlocks+" inverted");
            }
            System.out.println("...inverse calculated");
        } else {
            /**
             * No blocks - simpler to read in data
             */
            dataJama = new Matrix(numberOfIndividuals,totalNumberOfCovariates);
            for (int i=0; i<numberOfIndividuals; i++) {
                for (int v=0; v<totalNumberOfCovariates; v++) {
                    dataJama.set(i, v, dataScan.nextDouble());
                }
            }
        }

        /**
         * Read in cluster labels for each individual
         */
        clusters = new int[numberOfIndividuals];        
        if (numberOfClusters >0) {
            for (int i=0; i<numberOfIndividuals; i++) {
                clusters[i]  = dataScan.nextInt();
            }
        }
        
        /**
         * Read in outcome vector
         */
        outcomes = new int[numberOfIndividuals];
        if (whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()|
                whichLikelihoodType==LikelihoodTypes.LOGISTIC.ordinal()) {
            for (int i=0; i<numberOfIndividuals; i++) {
                outcomes[i]  = dataScan.nextInt();
            }
        } else if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN.ordinal()|
                whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
            continuousOutcomesJama = new Matrix(numberOfIndividuals,1);
            for (int i=0; i<numberOfIndividuals; i++) {
                continuousOutcomesJama.set(i, 0, dataScan.nextDouble());
            }        
        }
        
        /**
         * Read in survival times, if a Weibull model
         */
        times = new double[numberOfIndividuals];
        if (whichLikelihoodType==LikelihoodTypes.WEIBULL.ordinal()) {
            for (int i=0; i<numberOfIndividuals; i++) {
                times[i] = dataScan.nextDouble();
            }
        }
        
        /**
         * Read in information regarding the prior setup
         */
        betaPriorMus = new double[totalNumberOfCovariates];
        betaPriorSds = new double[totalNumberOfCovariates];
        numberOfCovariatesWithInformativePriors = dataScan.nextInt();  // Number of variables from the beginning
        for (int v=0; v<numberOfCovariatesWithInformativePriors; v++) {
            betaPriorMus[v]  = dataScan.nextDouble();
            betaPriorSds[v]  = dataScan.nextDouble();
        }            
        numberOfUnknownBetaPriors = dataScan.nextInt();  // Could be 0
        if (numberOfUnknownBetaPriors>0) {
            numberOfUnknownBetaPriors=1; // FORCE TO 1 FOR NOW
        }        
        commonBetaPriorPartitionIndices = new int[(numberOfUnknownBetaPriors+1)];
        commonBetaPriorPartitionIndices[0] = numberOfCovariatesWithInformativePriors; // could be 0
        commonBetaPriorPartitionIndices[numberOfUnknownBetaPriors] = totalNumberOfCovariates;
        if (numberOfUnknownBetaPriors>1) {     // Only in data if >1 components
            commonBetaPriorPartitionSizes = new int[(numberOfUnknownBetaPriors-1)];
            for (int i=0; i<(numberOfUnknownBetaPriors-1); i++) {
                commonBetaPriorPartitionSizes[i] = dataScan.nextInt(); //Model Space Poisson Means
                commonBetaPriorPartitionIndices[i+1]=commonBetaPriorPartitionIndices[i]+commonBetaPriorPartitionSizes[i];
            }
        }
        
        /**
         * If a G-prior is being used, set up the inverted matrix
         */
        if (arguments.useGPrior==1) {
            if (whichLikelihoodType==LikelihoodTypes.GAUSSIAN_MARGINAL.ordinal()) {
                // Already calculated above
            } else {
                xTxInvForGPrior = ( dataJama.getMatrix(
                        0,
                        (numberOfIndividuals-1),
                        numberOfUnknownBetaPriors,
                        (totalNumberOfCovariates-1)).transpose()
                        .times(dataJama.getMatrix(
                                0,
                                (numberOfIndividuals-1),
                                numberOfUnknownBetaPriors,
                                (totalNumberOfCovariates-1))) ).inverse();                
            }
        }
                
        /**
         * Feedback basic information about the dataset just read in, to the
         * console
         */
        System.out.println("------------");
        System.out.println("--- DATA ---");
        System.out.println("------------");
        System.out.println("Data read from "+arguments.pathToDataFile);
        System.out.println("Likelihood: "+likelihoodFamily);
        if(numberOfClusters>0) {
            System.out.println(numberOfIndividuals+" individuals from "+numberOfClusters+" clusters " +
                    "- random intercepts will be used");
        } else if (numberOfClusters == 0) {
            System.out.println(numberOfIndividuals+" individuals from a single cluster");
        }
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
                    System.out.println("Common priors with unknown Uniform[0,2] SDs will be"
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
                    System.out.println("A common prior with unknown Uniform[0,2] SD will be"
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
