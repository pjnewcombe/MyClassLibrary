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
    /**
     * Repeatedly used to read into.
     */
    private String argumentNameInFile;

    /**
     * 
     * Command Line Arguments. These are the arguments passed by the command
     * line call.
     * 
     */
    
    /**
     * Path to the Arguments file.
     */
    public String pathToArgumentsFile;
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
     * Every nth iteration to write to the results file.
     */
    public int consoleOutputInterval;
    /**
     * Random number seed used in the MCMC sampling.
     */
    public int whichSeed;
    /**
     * Whether to use model selection (0/1).
     */
    public int useReversibleJump;
    /**
     * Whether to use alternative initial values (0/1).
     */
    public int useAlternativeInitialValues; // whether to use alternative initial values
    /**
     * Maximum number of covariates allowed.
     */
    public int maximumModelDimension;
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
    
    /**
     * 
     * Initial values.
     * 
     */
    
    /**
     * Initial intercept value.
     */
    public double initialAlpha;
    /**
     * Initial values for the betas {@link Objects.IterationValues#betas}.
     */
    public double initialBetas;
    /**
     * Initial values for the unknown prior SDs common to all betas in a
     * model space component {@link Objects.IterationValues#logBetaPriorSd}.
     */
    public double initialBetaPriorSd; // init SD value of beta hyper prior
    /**
     * Initial value for the Weibull scale parameter {@link Objects.IterationValues#weibullScale}.
     */
    public double initialWeibullScale; // init Weibull scale value
    /**
     * Initial value for the Dirichlet concentration parameter {@link Objects.IterationValues#weibullScale}.
     */
    public double initialDirichletConcentration; // init Dirichlet concentration value
    /**
     * Initial value for the Weibull scale parameter {@link Objects.IterationValues#weibullScale}.
     */
    public double initialGaussianResidual; // init Weibull scale value
    
    /**
     * Vector of move type probabilities with four elements corresponding to
     * (in order) removal, addition, swap and null moves.
     */
    public double probRemove;
    public double probAdd;
    public double probSwap;
    public double probNull;
    
    /**
     * 
     * Prior distributions.
     * 
     */
    
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
     * Which prior family to use for the beta precision (0=uniform, 1=gamma).
     */
    public int betaPrecisionPriorFamily; // whether to use uniform or gamma sd prior
    /**
     * Hyper-parameter 1 of the Gamma prior on the unknown precision used for
     * the beta priors.
     */
    public double betaSigmaUniformPriorHyperparameter1;
    /**
     * Hyper-parameter 2 of the Gamma prior on the unknown precision used for
     * the beta priors.
     */
    public double betaSigmaUniformPriorHyperparameter2;
    /**
     * Hyper-parameter 1 of the Gamma prior on the unknown precision used for
     * the beta priors.
     */
    public double betaPrecisionGammaPriorHyperparameter1;
    /**
     * Hyper-parameter 2 of the Gamma prior on the unknown precision used for
     * the beta priors.
     */
    public double betaPrecisionGammaPriorHyperparameter2;
    /**
     * Which prior family to use for the Gaussian residuals (0=Jeffreys, 1=Inverse Gamma).
     */
    public int gaussianResidualPriorFamily; // whether to use uniform or gamma sd prior
    /**
     * Uniform prior hyper-parameter 1 for Gaussian residual SD.
     */
    public double gaussianResidualUniformPriorHyperparameter1;
    /**
     * Uniform prior hyper-parameter 2 for Gaussian residual SD.
     */
    public double gaussianResidualUniformPriorHyperparameter2;
    /**
     * Inverse-Gamma prior hyper-parameter 1 for Gaussian residual SD.
     */
    public double gaussianResidualVarianceInvGammaPrior_a;
    /**
     * Inverse-Gamma prior hyper-parameter 2 for Gaussian residual SD.
     */
    public double gaussianResidualVarianceInvGammaPrior_b;
    
    /**
     * Hyper-parameter 1 for a Gamma prior on the Weibull scale parameter.
     */
    public double weibullScaleGammaPriorHyperparameter1;
    /**
     * Hyper-parameter 2 for a Gamma prior on the Weibull scale parameter.
     */
    public double weibullScaleGammaPriorHyperparameter2;
    /**
     * Hyper-parameter 1 for a Gamma prior on the Dirichlet concentration parameter.
     */
    public double dirichletConcentrationGammaPriorHyperparameter1;
    /**
     * Hyper-parameter 2 for a Gamma prior on the Dirichlet concentration parameter.
     */
    public double dirichletConcentrationGammaPriorHyperparameter2;
    /**
     * Lower-bound for truncated Gamma prior on the Dirichlet concentration parameter.
     */
    public double dirichletConcentrationMinimum;
    /**
     * Lower-bound for truncated Gamma prior on the Dirichlet concentration parameter.
     */
    public double logDirichletConcentrationMinimum;
    /**
     * Multiplier for AUC multiplier within the posterior.
     */
    public double aucMultiplierWeight;
    
    /**
     * 
     * Proposal distributions.
     * 
     */
    
    /**
     * Initial SD of the proposal distribution for the intercept 
     * {@link Objects.IterationValues#alpha}.
     */
    public double proposalDistributionSdForAlpha;
    /**
     * Initial SD of the proposal distribution for the betas 
     * {@link Objects.IterationValues#betas}.
     */
    public double proposalDistributionSdForBetas;
    /**
     * Initial SD of the proposal distribution for a beta upon addition to the model.
     */
    public double proposalDistributionSdForAddingBeta;
    /**
     * Initial SD of the proposal distribution for a beta upon addition to the model 
     * during a swap move.
     */
    public double proposalDistributionSdForSwappedInBeta;
    /**
     * Initial SD of the proposal distribution for the unknown prior SD(s) for the betas
     * {@link Objects.IterationValues#logBetaPriorSd}.
     */
    public double proposalDistributionSdForBetaPriorSd;
    /**
     * Initial SD of the proposal distribution for the Weibull scale parameter
     * {@link Objects.IterationValues#logWeibullScale}.
     */
    public double proposalDistributionSdForLogWeibullScale;
    /**
     * Initial SD of the proposal distribution for the Dirichlet concentration parameter
     */
    public double proposalDistributionSdForDirichletConcentration;    
    /**
     * Initial SD of the proposal distribution for the Gaussian residual parameter
     * {@link Objects.IterationValues#logGaussianResidual}.
     */
    public double proposalDistributionSdForLogGaussianResidual;
    /**
     * Initial SD of the proposal distribution for the tau variable selection
     * coefficient
     * {@link Objects.IterationValues#tau}.
     */
    public double proposalDistributionSdForTau;
    
    /**
     * 
     * Proposal distribution adaption.
     * 
     */
    
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
        
        /**
         * 
         * Process/extract command line arguments.
         * 
         */
        
        pathToArgumentsFile = args[0];
        pathToDataFile = args[1];
        pathToResultsFile = args[2];
        numberOfIterations = Integer.parseInt(args[3]);
        burnInLength = Integer.parseInt(args[4]);
        thinningInterval = Integer.parseInt(args[5]);
        consoleOutputInterval = Integer.parseInt(args[6]);
        whichSeed = Integer.parseInt(args[7]);
        useReversibleJump = Integer.parseInt(args[8]);
        int dummyUseGPrior= Integer.parseInt(args[9]);
        useAlternativeInitialValues = Integer.parseInt(args[10]);
        maximumModelDimension = Integer.parseInt(args[11]);
        numberOfModelSpacePriorPartitions = Integer.parseInt(args[12]);
        modelSpacePriorFamily = Integer.parseInt(args[13]);
        if (modelSpacePriorFamily==0) {
            // Poisson
            modelSpacePoissonPriorRate = new double[numberOfModelSpacePriorPartitions];
            for (int i=0; i<numberOfModelSpacePriorPartitions; i++) {
                modelSpacePoissonPriorRate[i] = Double.parseDouble(args[14+i]); //Model Space Poisson Means
            }
            if (numberOfModelSpacePriorPartitions > 1) {
                modelSpacePriorPartitionSizes = new int[numberOfModelSpacePriorPartitions-1];
                for (int i=0; i<(numberOfModelSpacePriorPartitions-1); i++) {
                    modelSpacePriorPartitionSizes[i] = Integer.parseInt(args[14+numberOfModelSpacePriorPartitions+i]); //Model Space Poisson Means
                }
            }
        } else if (modelSpacePriorFamily==1) {
            // Beta-binomial
            modelSpaceBetaBinomialPriorHyperparameterA = new double[numberOfModelSpacePriorPartitions];
            modelSpaceBetaBinomialPriorHyperparameterB = new double[numberOfModelSpacePriorPartitions];
            for (int i=0; i<numberOfModelSpacePriorPartitions; i++) {
                modelSpaceBetaBinomialPriorHyperparameterA[i] = Double.parseDouble(args[14+(i*2)]); //Model Space Poisson Means
                modelSpaceBetaBinomialPriorHyperparameterB[i] = Double.parseDouble(args[14+(i*2)+1]); //Model Space Poisson Means
            }
            if (numberOfModelSpacePriorPartitions > 1) {
                modelSpacePriorPartitionSizes = new int[numberOfModelSpacePriorPartitions-1];
                for (int i=0; i<(numberOfModelSpacePriorPartitions-1); i++) {
                    modelSpacePriorPartitionSizes[i] = Integer.parseInt(args[14+(2*numberOfModelSpacePriorPartitions)+i]); //Model Space Poisson Means
                }
            }
        }
        
        /**
         * 
         * Read/process Arguments text file. argumentNameInFile is repeatedly
         * re-used to read the character string labels contained in the file.
         * 
         */

        Scanner argScan = new Scanner(new File(pathToArgumentsFile));

        /**
         * Prior distributions.
         */
        
        argumentNameInFile = argScan.next();
        alphaPriorMu = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        alphaPriorSd = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betaPrecisionPriorFamily = argScan.nextInt();
        
        argumentNameInFile = argScan.next();        
        betaSigmaUniformPriorHyperparameter1 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betaSigmaUniformPriorHyperparameter2 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betaPrecisionGammaPriorHyperparameter1 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        betaPrecisionGammaPriorHyperparameter2 = argScan.nextDouble();

        argumentNameInFile = argScan.next();        
        weibullScaleGammaPriorHyperparameter1 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        weibullScaleGammaPriorHyperparameter2 = argScan.nextDouble();

        argumentNameInFile = argScan.next();        
        dirichletConcentrationGammaPriorHyperparameter1 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        dirichletConcentrationGammaPriorHyperparameter2 = argScan.nextDouble();

        argumentNameInFile = argScan.next();        
        dirichletConcentrationMinimum = argScan.nextDouble();
        logDirichletConcentrationMinimum = Math.log(dirichletConcentrationMinimum);

        argumentNameInFile = argScan.next();        
        gaussianResidualPriorFamily = argScan.nextInt();
        
        argumentNameInFile = argScan.next();        
        gaussianResidualUniformPriorHyperparameter1 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        gaussianResidualUniformPriorHyperparameter2 = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        gaussianResidualVarianceInvGammaPrior_a = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        gaussianResidualVarianceInvGammaPrior_b = argScan.nextDouble();

        argumentNameInFile = argScan.next();        
        aucMultiplierWeight = argScan.nextDouble();
        
        /**
         * Initial Values.
         */
        
        argumentNameInFile = argScan.next();        
        initialAlpha = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        initialBetas = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        initialBetaPriorSd = argScan.nextDouble();

        argumentNameInFile = argScan.next();        
        initialWeibullScale = argScan.nextDouble();

        argumentNameInFile = argScan.next();        
        initialDirichletConcentration = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        initialGaussianResidual = argScan.nextDouble();
        
        /**
         * Adaption of proposal distributions.
         */
        
        argumentNameInFile = argScan.next();        
        adaptionBinSize = argScan.nextInt();
        
        argumentNameInFile = argScan.next();        
        adaptionLength = argScan.nextInt();
        
        /**
         * Reversible jump move probabilities.
         */
        
        argumentNameInFile = argScan.next();        
        probRemove = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        probAdd = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        probSwap = argScan.nextDouble();
        
        probNull = 1-(probRemove+probAdd+probSwap);

        /**
         * Proposal distributions.
         */
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForAlpha = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForBetas = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForBetaPriorSd = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForLogWeibullScale = argScan.nextDouble();

        argumentNameInFile = argScan.next();        
        proposalDistributionSdForDirichletConcentration = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForLogGaussianResidual = argScan.nextDouble();

        argumentNameInFile = argScan.next();        
        proposalDistributionSdForTau = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForAddingBeta = argScan.nextDouble();
        
        argumentNameInFile = argScan.next();        
        proposalDistributionSdForSwappedInBeta = argScan.nextDouble();
        
        /**
         * Alternative initial values.
         */
        
        if (useAlternativeInitialValues==1) {
            // Initial Values
            argumentNameInFile = argScan.next();        
            initialAlpha = argScan.nextDouble();
            argumentNameInFile = argScan.next();        
            initialBetas = argScan.nextDouble();
            argumentNameInFile = argScan.next();        
            initialBetaPriorSd = argScan.nextDouble();
            argumentNameInFile = argScan.next();        
            argumentNameInFile = argScan.next();        
            initialWeibullScale = argScan.nextDouble();            
        }
        
    }
    
}
