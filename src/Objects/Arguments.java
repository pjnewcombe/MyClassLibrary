package Objects;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * A class to contain all user specified arguments via both the command line,
 * and via an Arguments text file.
 * 
 * @author Paul J Newcombe
 */
public class Arguments {
    // Command line arguments -----------------------------------
    private String pathToArgumentsFile;
    private String argumentNameInFile;
    /**
     * Path to a text file containing the data to analyse.
     */
    public String pathToDataFile;
    /**
     * Desired path for the text file which will contain the results.
     */
    public String pathToResultsFile;
    /**
     * Number of numberOfIterations.
     */
    public int numberOfIterations;
    /**
     * Number of numberOfIterations to discard as burn in.
     */
    public int burnInLength;
    /**
     * Every nth iteration to save to the results file.
     */
    public int thinningInterval;
    /**
     * Random number seed used in the MCMC sampling.
     */
    public int whichSeed;
    /**
     * Number of model space prior components.
     */
    public int numberOfModelSpacePriorPartitions;
    /**
     * Prior family distribution to use over the model space (0=Poisson,
     * 1=Beta-Binomial).
     */
    public int modelSpacePriorFamily; // 0: Poisson, 1: Beta-binomial
    /**
     * Vector of number of covariates in each model space component.
     */
    public int[] modelSpacePriorPartitionSizes;
    /**
     * Vector of prior Poisson rates for the different model space components.
     */
    public double[] modelSpacePoissonPriorRate;
    /**
     * Vector of Beta-Binomial hyper-parameters a for the different model space 
     * components.
     */
    public double[] modelSpaceBetaBinomialPriorHyperparameterA;
    /**
     * Vector of Beta-Binomial hyper-parameters b for the different model space 
     * components.
     */
    public double[] modelSpaceBetaBinomialPriorHyperparameterB;
    
    // ARGUMENTS FILE --------------------------------------------
    // static means each instaniation will have fixed param values
    /**
     * Every nth iteration to write to the results file.
     */
    public int consoleOutputInterval;  // Number of Sds total
    /**
     * Whether to store random intercepts in results file.
     */
    public int recordClusterIntercepts;
    
    // Initial Values
    /**
     * Use saturated initial model (0/1).
     */
    public int useSaturatedInitialModel; // whether to use saturated initial model
    /**
     * Initial intercept value.
     */
    public double initialAlpha; // init intercept value
    /**
     * Initial values for the betas {@link Objects.IterationValues#betas}.
     */
    public double initialBetas; // init logOR values
    /**
     * Initial values for the unknown prior SDs common to all betas in a
     * model space component {@link Objects.IterationValues#logBetaPriorSd}.
     */
    public double initialBetaPriorSd; // init SD value of beta hyper prior
    /**
     * Initial value for the between cluster SD when fitting random intercepts
     * {@link Objects.IterationValues#logBetweenClusterSd}.
     */
    public double initialBetweenClusterSd; // init between-cluster SD value
    /**
     * Initial value for the Weibull scale parameter {@link Objects.IterationValues#weibullScale}.
     */
    public double initialWeibullScale; // init Weibull scale value
    /**
     * Whether to use alternative initial values (0/1).
     */
    public int useAlternativeInitialValues; // whether to use alternative initial values
    /**
     * Whether to use model selection (0/1).
     */
    public int useReversibleJump; // indicates whether to use Model Prior
    // Prior Params
    /**
     * Normal prior mean for the intercept parameter
     * {@link Objects.IterationValues#alpha}.
     */
    public double alphaPriorMu;
    /**
     * Normal prior SD for the intercept parameter
     * {@link Objects.IterationValues#alpha}.
     */
    public double alphaPriorSd;
    /**
     * Normal prior mean for the beta parameters
     * {@link Objects.IterationValues#betas}.
     */
    public double betaPriorMu; // log-OR prior mean
    /**
     * Normal prior SD for the beta parameters
     * {@link Objects.IterationValues#betas}.
     */
    public double betaPriorSd; // log-OR prior sd
    /**
     * Hyperparameter 1 of the Gamma prior on the unknown precision used for
     * the beta's priors.
     */
    public double betaPrecisionGammaPriorHyperparameter1;
    /**
     * Hyperparameter 2 of the Gamma prior on the unknown precision used for
     * the beta's priors.
     */
    public double betaPrecisionGammaPriorHyperparameter2;
    /**
     * Which prior family to use for the between-cluster SD (0=uniform, 1=gamma).
     */
    public int betweenClusterSdPriorFamily; // whether to use uniform or gamma sd prior
    /**
     * Uniform prior hyper-parameter 1 for the between cluster SD.
     */
    public double betweenClusterPrecisionUniformPriorHyperparameter1;
    /**
     * Uniform prior hyper-parameter 2 for the between cluster SD.
     */
    public double betweenClusterPrecisionUniformPriorHyperparameter2;
    /**
     * Gamma prior hyper-parameter 1 for the between cluster SD.
     */
    public double betweenClusterPrecisionGammaPriorHyperparameter1;
    /**
     * Gamma prior hyper-parameter 2 for the between cluster SD.
     */
    public double betweenClusterPrecisionGammaPriorHyperparameter2;
    /**
     * Hyper-parameter 1 for a Gamma prior on the Weibull scale parameter.
     */
    public double weibullScaleGammaPriorHyperparameter1;
    /**
     * Hyper-parameter 2 for a Gamma prior on the Weibull scale parameter.
     */
    public double weibullScaleGammaPriorHyperparameter2;
    /**
     * Vector of move type probabilities with four elements corresponding to
     * (in order) removal, addition, swap and null moves.
     */
    public double[] moveProbabilities = new double[3]; //rem, add,
                                                              // swap, null
    // Proposal Sds And Parameters
    /**
     * SD of the proposal distribution for the intercept 
     * {@link Objects.IterationValues#alpha}.
     */
    public double proposalDistributionSdForAlpha;
    /**
     * SD of the proposal distribution for the betas 
     * {@link Objects.IterationValues#betas}.
     */
    public double proposalDistributionSdForBetas;
    /**
     * SD of the proposal distribution for the cluster intercepts 
     * {@link Objects.IterationValues#clusterIntercepts}.
     */
    public double proposalDistributionSdForClusterIntercepts;
    /**
     * SD of the proposal distribution for the between cluster SD 
     * {@link Objects.IterationValues#logBetweenClusterSd}.
     */
    public double proposalDistributionSdForLogBetweenClusterSd;
    /**
     * SD of the proposal distribution for a beta upon addition to the model.
     */
    public double proposalDistributionSdForAddingBeta;
    /**
     * SD of the proposal distribution for a beta upon addition to the model 
     * during a swap move.
     */
    public double proposalDistributionSdForSwappedInBeta;
    /**
     * SD of the proposal distribution for the unknown prior SD(s) for the betas
     * {@link Objects.IterationValues#logBetaPriorSd}.
     */
    public double proposalDistributionSdForBetaPriorSd;
    /**
     * SD of the proposal distribution for the Weibull scale parameter
     * {@link Objects.IterationValues#logWeibullScale}.
     */
    public double proposalDistributionSdForLogWeibullScale;
    /**
     * Every nth iteration to modify the proposal SDs, during adaption.
     */
    public int adaptionBinSize;
    /**
     * Number of numberOfIterations to perform adaption of proposal distribution SDs 
     * for.
     */
    public int adaptionLength;
   
    /**
     * Constructor - populates the various fields by reading in from an
     * Arguments text file, the path of which is included in args.
     * 
     * @param args Command line arguments passed to bayesglm.
     * 
     * @throws FileNotFoundException 
     * 
     * @author Paul J Newcombe
     */
    public Arguments(String[] args) throws FileNotFoundException {
        // Extract command line arguments ------------------------------
        pathToArgumentsFile = args[0];// data file
        pathToDataFile = args[1];// data file
        pathToResultsFile = args[2]; // results file
        numberOfIterations = Integer.parseInt(args[3]); // no itns
        burnInLength = Integer.parseInt(args[4]); // no burn in itns
        thinningInterval = Integer.parseInt(args[5]); //every 'thinningInterval'th iteration is saved
        consoleOutputInterval = Integer.parseInt(args[6]); // how often to output progress to terminal
        whichSeed = Integer.parseInt(args[7]); // no burn in itns
        useReversibleJump = Integer.parseInt(args[8]); // whether to use model selection
        useAlternativeInitialValues = Integer.parseInt(args[9]); // whether to use model selection
        numberOfModelSpacePriorPartitions = Integer.parseInt(args[10]); //Number of model space components
        modelSpacePriorFamily = Integer.parseInt(args[11]); //Number of model space components
        if (modelSpacePriorFamily==0) {
            // Poisson
            modelSpacePoissonPriorRate = new double[numberOfModelSpacePriorPartitions];
            for (int i=0; i<numberOfModelSpacePriorPartitions; i++) {
                modelSpacePoissonPriorRate[i] = Double.parseDouble(args[12+i]); //Model Space Poisson Means
            }
            if (numberOfModelSpacePriorPartitions > 1) {
                modelSpacePriorPartitionSizes = new int[numberOfModelSpacePriorPartitions-1];
                for (int i=0; i<(numberOfModelSpacePriorPartitions-1); i++) {
                    modelSpacePriorPartitionSizes[i] = Integer.parseInt(args[12+numberOfModelSpacePriorPartitions+i]); //Model Space Poisson Means
                }
            }
        } else if (modelSpacePriorFamily==1) {
            // Beta-binomial
            modelSpaceBetaBinomialPriorHyperparameterA = new double[numberOfModelSpacePriorPartitions];
            modelSpaceBetaBinomialPriorHyperparameterB = new double[numberOfModelSpacePriorPartitions];
            for (int i=0; i<numberOfModelSpacePriorPartitions; i++) {
                modelSpaceBetaBinomialPriorHyperparameterA[i] = Double.parseDouble(args[12+(i*2)]); //Model Space Poisson Means
                modelSpaceBetaBinomialPriorHyperparameterB[i] = Double.parseDouble(args[12+(i*2)+1]); //Model Space Poisson Means
            }
            if (numberOfModelSpacePriorPartitions > 1) {
                modelSpacePriorPartitionSizes = new int[numberOfModelSpacePriorPartitions-1];
                for (int i=0; i<(numberOfModelSpacePriorPartitions-1); i++) {
                    modelSpacePriorPartitionSizes[i] = Integer.parseInt(args[12+(2*numberOfModelSpacePriorPartitions)+i]); //Model Space Poisson Means
                }
            }
        }

        // Arguments file --------------------------------------------
        Scanner argScan = new Scanner(new File(pathToArgumentsFile));
        argumentNameInFile = argScan.next();
        recordClusterIntercepts = argScan.nextInt(); // how often to output progress

        // Prior Params
        argumentNameInFile = argScan.next();
        alphaPriorMu = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        alphaPriorSd = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betaPrecisionGammaPriorHyperparameter1 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betaPrecisionGammaPriorHyperparameter2 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betweenClusterSdPriorFamily = argScan.nextInt(); // Precision prior; 0: Unif 1: Gamma
        
        argumentNameInFile = argScan.next();        
        betweenClusterPrecisionUniformPriorHyperparameter1 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betweenClusterPrecisionUniformPriorHyperparameter2 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betweenClusterPrecisionGammaPriorHyperparameter1 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betweenClusterPrecisionGammaPriorHyperparameter2 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        weibullScaleGammaPriorHyperparameter1 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        weibullScaleGammaPriorHyperparameter2 = argScan.nextDouble();

        // Initial Values        
        argumentNameInFile = argScan.next();        
        useSaturatedInitialModel = argScan.nextInt(); // whether to use staurated initial model
        // If model selection is disabled force saturated initial model
        if (useReversibleJump==0) {
            useSaturatedInitialModel = 1;
        }
        
        argumentNameInFile = argScan.next();        
        initialAlpha = argScan.nextDouble(); // init logOR values
        
        argumentNameInFile = argScan.next();        
        initialBetas = argScan.nextDouble(); // init logOR values
        
        argumentNameInFile = argScan.next();        
        initialBetaPriorSd = argScan.nextDouble();

        argumentNameInFile = argScan.next();        
        initialBetweenClusterSd = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        initialWeibullScale = argScan.nextDouble();
        
        // Modelling Arguments
        
        argumentNameInFile = argScan.next();        
        adaptionBinSize = argScan.nextInt();
        
        argumentNameInFile = argScan.next();        
        adaptionLength = argScan.nextInt();
        
        argumentNameInFile = argScan.next();        
        for (int j=0;j<3;j++) {
            moveProbabilities[j] = argScan.nextDouble(); //rem, add, swap, null
        }

        // Proposal Sds
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForAlpha = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForBetas = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForBetaPriorSd = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForClusterIntercepts = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForLogBetweenClusterSd = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForLogWeibullScale = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForAddingBeta = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForSwappedInBeta = argScan.nextDouble();
        
        // Alternative initial values
        if (useAlternativeInitialValues==1) {
            // Initial Values        
            argumentNameInFile = argScan.next();        
            useSaturatedInitialModel = argScan.nextInt(); // whether to use staurated initial model
            // If model selection is disabled force saturated initial model
            if (useReversibleJump==0) {
                useSaturatedInitialModel = 1;
            }
            argumentNameInFile = argScan.next();        
            initialAlpha = argScan.nextDouble(); // init logOR values
            argumentNameInFile = argScan.next();        
            initialBetas = argScan.nextDouble(); // init logOR values
            argumentNameInFile = argScan.next();        
            initialBetaPriorSd = argScan.nextDouble();
            argumentNameInFile = argScan.next();        
            initialBetweenClusterSd = argScan.nextDouble();
            argumentNameInFile = argScan.next();        
            initialWeibullScale = argScan.nextDouble();            
        }
        
    }
    
}
