package Methods;

import Jama.Matrix;
import java.util.Random;
import static java.lang.Math.PI;

/**
 * Collection of generic mathematical functions.
 * 
 * @author Paul J Newcombe
 */
public class GeneralMaths {

    /***
     * Given a mean and Sd of a normal distribution, a value is drawn and 
     * returned. A value is drawn from a N(0,1) distribution, then converted to
     * a draw from a N(mean,Sd) distribution. 
     * 
     * @param mean Normal mean
     * @param sd Normal standard deviation
     * @return A single draw from a normal distribution
     */      
     public static double normalDraw(
             double mean, 
             double sd,
             Random randomDraws) {
         double draw = mean + sd*randomDraws.nextGaussian();
         return draw;     
     }
     
    /***
     * Given the mean and SD of a normal distribution, the log-Normal density of
     * 'value' is calculated and returned. 
     * 
     * @param value Value to calculate the log-normal density of
     * @param mu  Normal mean
     * @param sd  Normal standard deviation
     * @return Log-Normal density of 'value'
     */ 
    public static double logNormDens(
             double value, 
             double mu, 
             double sd) {
         double logDensity = (double) (-Math.pow((value-mu),2)
                 /(2*Math.pow(sd,2))) 
         -Math.log(sd*Math.sqrt(2*PI));
         return logDensity;     
     }
    
    /***
     * The Multivariate Normal log-density of JAMA Matrix object 'values' (a 
     * vector) is calculated and returned where the multivariate normal is 
     * defined by required JAMA Matrix objects 'means' and 'covMatrix'.
     * 
     * @param values
     *          JAMA Matrix object: Vector of values to calculate multivariate 
     *          normal log-density of
     * @param means
     *          JAMA Matrix object: Vector of means of multivariate normal 
     *          distribution
     * @param covMatrix
     *          JAMA Matrix object: Covariance matrix of multivariate normal 
     *          distribution
     * @return Multivariate Normal log-density of 'values'
     */ 
     public static double logMVNormDens(
             Matrix values,      // observed marginal counts
             Matrix means,    // log linear model means
             Matrix covMatrix      ) {
         // calculate determinant 
         double detCovMatrix = covMatrix.det();

         // calculate inverse of covariance matrix for likelihood calculation
         Matrix invCovMatrix = covMatrix.inverse();
         
         // calculates the log of the square likelihood contribution, then halfs 
         // it (equivalent to rooting likelihood), to avoid having to root the
         // inverse covariance matrix - which often is not possible
         double density =
                 (-values.minus(means).transpose()
                 .times(invCovMatrix)
                 .times(values.minus(means))
                 .get(0, 0)/2)-Math.log(Math.sqrt(detCovMatrix));
         return density;     
     }

    /***
     * Given an integer, this method calculates 2 to the power of that integer.
     * 
     * @param power Power that 2 should be raised to
     * @return 2^'power'
     */             
    public static int powerOf2(int power) {
        int x = 1;
        for (int i=0; i<power; i++) {x = x*2;}
        return x;
    }
    
    /***
     * Given a value the logit is calculated and returned.
     * 
     * @param value A value to calculate the logit of
     * @return The logit of 'value'
     */             
     public static double logit (
             double value
             ) {
         double logitValue = Math.log((double)(value/(1-value)));
         return logitValue;
     }

    /***
     * Given a value the inverse-logit is calculated and returned.
     * 
     * @param value A value to calculate the inverse-logit of
     * @return The inverse-logit of 'value'
     */             
     public static double invLogit (
             double value
             ) {
         double logitValue = (double)(Math.exp(value)/(1+Math.exp(value)));
         return logitValue;
     }


    /***
     * Given an integer the log-factorial is calculated and returned.
     * 
     * @param n Integer to calculate the log-factorial of
     * @return The log-factorial of 'n'
     */   
      public static double logFactorial(int n) {
            double logFact = 0;
            for (int i=1; i<=n; i++) {
                logFact = logFact + Math.log(i);
            }        
            return logFact;
        }
      

    /**
     * Given a vector of N control proportions and N ORs (one of which is 1 
     * corresponding to the proportion taken as baseline), Seaman's formula is 
     * applied to calculate and return the N corresponding case proportions.
     * 
     * @param logORs Vector of N log-ORs (one of which is 1 corresponding to 
     * proportion taken as baseline)
     * @param contProportions Vector of N control proportions 
     * @param N Total number of proportions 
     * @return Vector of N case proportions 
     */   
    public static double[] seamansFormula (
                double[] logORs,
                double[] contProportions,
                int N) {
        double[] caseProportions = new double[N];
        double denominator = 0; 
        for (int n=0; n<N; n++) { 
            // numerator of each probabiltiy is OR*case prob
            caseProportions[n] = Math.exp(logORs[n])*contProportions[n];
            // Denominator is total of all numerators as in Seaman's formula
            denominator = denominator + caseProportions[n];
        }            
        for (int n=0; n<N; n++) {
            caseProportions[n] = (double)caseProportions[n]/(double)denominator;
        }
        return caseProportions;
    }

    /***
     * Normalises a vector of probabilities, such that they sum to one.
     * 
     * @param props Vector of proportions to normalise
     * @param N Total number of proportions 
     * @return Normalised proportions 
     */   
    public static double[] normaliseProps (
                double[] props,
                int N) {

        double[] normalisedProportions = new double[N];

        double total = 0;
        for (int h=0; h<N; h++) {
            total = total + props[h];
        }
        for (int h=0; h<N; h++) {
            normalisedProportions[h] 
                = (double)(props[h]
                /total);                                            
        }        

        return normalisedProportions;
    }   
    

    

     
     
     

}
