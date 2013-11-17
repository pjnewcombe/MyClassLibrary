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
    // static means each instaniation will have fixed param values
    String likelihoodFamily;
    /**
     * Whether this is survival analysis data (0/1).
     */
    public int survivalAnalysis;
    /**
     * Whether to fit random intercepts (0/1).
     */
    public int useRandomIntercepts;
    /**
     * Total number of covariates.
     */
    public int totalNumberOfCovariates;
    /**
     * String of variable names (to include in the results file).
     */
    public String[] covariateNames;
    /**
     * Number of covariates to fix in the likelihoodFamily (eg confounders).
     */
    public int numberOfCovariatesToFixInModel;
    /**
     * Number of covariates in each likelihoodFamily space component.
     */
    public int[] modelSpacePartitionSizes;
    /**
     * Covariate indices to partition among the different likelihoodFamily space components.
     */
    public int[] modelSpacePartitionIndices;
    /**
     * Vector of prior means for the different likelihoodFamily space components.
     */
    public double[] modelSpacePoissonPriorMeans;
    /**
     * Number of individuals in the dataset.
     */
    public int numberOfIndividuals;
    /**
     * Number of clusters if fitting random intercepts.
     */
    public int numberOfClusters;
    /**
     * Matrix of covariates values (rows individuals, columns covariates).
     */
    public double[][] data;         // Covariate matrix
    /**
     * Data matrix in the JAMA format.
     */
    public Matrix dataJama;
    public Matrix[] dataJamaColumns; // Study hapFqs as a JAMA object
    /**
     * Vector of clusters assignments for each individual.
     */
    public int[] clusters;         // Vector of clusters assignments
    /**
     * Vector of outcomes for each individual.
     */
    public int[] outcomes;         // Binary outcomes vector (logistic likelihoodFamily)
    /**
     * Vector of survival times for each individual (if it is survival data).
     */
    public double[] times;         // Event/death times vector (survival likelihoodFamily)
    // Beta priors
    /**
     * Vector of normal prior means for the different betas.
     */
    public double[] betaPriorMus;         // Vector of prior means
    /**
     * Vector of normal prior SDs for the different betas.
     */
    public double[] betaPriorSds;         // Vector of prior SDs
    /**
     * Number of covariates (from the beginning) which informative priors have
     * been provided for.
     */
    public int numberOfCovariatesWithInformativePriors;
    /**
     * Number of likelihoodFamily space components, in which a common prior with unknown
     * SD will be used over betas.
     */
    public int numberOfModelSpacePartitions;
    /**
     * Number of covariates in each likelihoodFamily space component using the same
     * common prior with unknown SD.
     */
    public int[] commonBetaPriorPartitionSizes;
    /**
     * Covariate indices which partition the covariates into the different
     * likelihoodFamily space components using the same common prior with unknown SD over
     * the betas.
     */
    public int[] commonBetaPriorPartitionIndices;


    // Contructor method has name as class
    // NewLikeData data1 = new NewLikeData(m, sg, sa)
    /**
     * Constructor function - populates the object by reading data from a text
     * file, whose path is stored in the {@link Objects.Arguments} object.
     * 
     * @param arguments {@link Objects.Arguments} class object, containing all 
     * modeling arguments, including the location of the data file.
     * @throws FileNotFoundException 
     */
    public Data(Arguments arguments) throws FileNotFoundException {
        Scanner dataScan = new Scanner(new File(arguments.pathToDataFile));
        likelihoodFamily = dataScan.next();             // which likelihoodFamily type
        if (likelihoodFamily.equals("Weibull")) {
            survivalAnalysis=1;
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

        ///////////////////////////////////////////
        // Setup info for likelihoodFamily space components //
        ///////////////////////////////////////////
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
        
        data = new double[numberOfIndividuals][totalNumberOfCovariates];
        clusters = new int[numberOfIndividuals];
        outcomes = new int[numberOfIndividuals];
        times = new double[numberOfIndividuals];
        betaPriorMus = new double[totalNumberOfCovariates];
        betaPriorSds = new double[totalNumberOfCovariates];

        // READ IN Data

        dataJama = new Matrix(numberOfIndividuals,totalNumberOfCovariates);
        dataJamaColumns = new Matrix[totalNumberOfCovariates];
        for (int v=0; v<totalNumberOfCovariates; v++) {
            dataJamaColumns[v] = new Matrix(numberOfIndividuals,1);
        }

        for (int i=0; i<numberOfIndividuals; i++) {
            for (int v=0; v<totalNumberOfCovariates; v++) {
                data[i][v] = dataScan.nextDouble();
                dataJama.set(i, v, data[i][v]);
                dataJamaColumns[v].set(i, 0, data[i][v]);
            }
        }

        if (numberOfClusters >0) {
            for (int i=0; i<numberOfIndividuals; i++) {
                clusters[i]  = dataScan.nextInt();
            }
        }

        for (int i=0; i<numberOfIndividuals; i++) {
            outcomes[i]  = dataScan.nextInt();
        }
        
        if (survivalAnalysis==1) {
            for (int i=0; i<numberOfIndividuals; i++) {
                times[i] = dataScan.nextDouble();
            }
        }

        numberOfCovariatesWithInformativePriors = dataScan.nextInt();  // Number of variables from the beginning
        for (int v=0; v<numberOfCovariatesWithInformativePriors; v++) {
            betaPriorMus[v]  = dataScan.nextDouble();
            betaPriorSds[v]  = dataScan.nextDouble();
        }            
        numberOfModelSpacePartitions = dataScan.nextInt();  // Could be 0
        commonBetaPriorPartitionIndices = new int[(numberOfModelSpacePartitions+1)];
        commonBetaPriorPartitionIndices[0] = numberOfCovariatesWithInformativePriors; // could be 0
        commonBetaPriorPartitionIndices[numberOfModelSpacePartitions] = totalNumberOfCovariates;
        if (numberOfModelSpacePartitions>1) {     // Only in data if >1 components
            commonBetaPriorPartitionSizes = new int[(numberOfModelSpacePartitions-1)];
            for (int i=0; i<(numberOfModelSpacePartitions-1); i++) {
                commonBetaPriorPartitionSizes[i] = dataScan.nextInt(); //Model Space Poisson Means
                commonBetaPriorPartitionIndices[i+1]=commonBetaPriorPartitionIndices[i]+commonBetaPriorPartitionSizes[i];
            }
        }
        
        // INITIATION MESSAGE
        System.out.println("Data read from "+arguments.pathToDataFile);
        System.out.println(likelihoodFamily+" model");
        if(numberOfClusters>0) {
            System.out.println(numberOfIndividuals+" individuals from "+numberOfClusters+" clusters " +
                    "- random intercept will be used");
        } else if (numberOfClusters == 0) {
            System.out.println(numberOfIndividuals+" individuals from a single cluster");
        }
        if (arguments.modelSpacePriorFamily==0) {
            System.out.println("Poisson prior(s) on model space will be used");
        } else if (arguments.modelSpacePriorFamily==1) {
            System.out.println("Beta-binomial prior(s) on model space will be used");            
        }
        System.out.println(totalNumberOfCovariates+" variables total, "+arguments.numberOfModelSpacePriorPartitions+" model space " +
                "components, "+numberOfCovariatesToFixInModel+" excluded from RJ");
        if (numberOfModelSpacePartitions>0) {
            System.out.println("Fixed priors for "+numberOfCovariatesWithInformativePriors+" of the beta's have been provided");            
            System.out.println(numberOfModelSpacePartitions+" common hyperprior(s) for the remaining beta prior SDs will be used");
        } else {
            System.out.println("Fixed priors for all the beta's have been provided");            
        }
        
    }
}
