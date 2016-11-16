package Methods;

import Jama.Matrix;
import java.util.Random;

/**
 * Collection of generic methods.
 * 
 * @author Paul J Newcombe
 */
public class GeneralMethods {    
    /**
     * Chooses a Reversible Jump move type when using the conjugate model.
     * 
     * @param probRemove Probability of removing a covariate
     * @param probAdd Probability of adding a covariate
     * @param probSwap Probability of swapping two covariates
     * @param randomDraws Random number generator object.
     * 
     * @return Indicates chosen move type (0=Removal, 1=Addition, 2=Swap, 3=Null)
     */
    public static int chooseMove(
            double probRemove,
            double probAdd,
            double probSwap,
            Random randomDraws) {
        int whichMove;
        double eventDraw = randomDraws.nextFloat();
        
        /**
         * 
         * Must be one of the below.
         * 
         */
        
        /**
         * Null move.
         */
        if (eventDraw>=(probRemove+probAdd+probSwap)) {whichMove =3;}
        /**
         * Swap move.
         */
        else if (eventDraw>=(probRemove+probAdd)) {whichMove = 2;} // Never assesed if this number is bigger than 2
        /**
         * Add move.
         */
        else if (eventDraw>=probRemove) {whichMove = 1;}
        /**
         * Remove move.
         */
        else {whichMove=0;}
        
        return whichMove;          
        
    }
    
    /**
     * Chooses a Reversible Jump move type.
     * 
     * @param M Total number of variables.
     * @param maxAllowedMarkers Specifies a maximum number of selected variables.
     * @param presentMarkN Number of variables currently included.
     * @param moveProbabilities Vector of probabilities with which to select the
     * different move types.
     * @param randomDraws Random number generator object.
     * 
     * @return Indicates chosen move type (0=Removal, 1=Addition, 2=Swap, 3=Null)
     */
    public static int chooseMove_OLD(
            int M,
            int maxAllowedMarkers,
            int presentMarkN,
            double[] moveProbabilities,
            Random randomDraws) {
        int whichMove;
        double eventDraw = randomDraws.nextFloat();
        if (eventDraw>=moveProbabilities[2]) {whichMove =3;}
        else if (eventDraw>=moveProbabilities[1]) {whichMove = 2;} // Never assesed if this number is bigger than 2
        else if (eventDraw>=moveProbabilities[0]) {whichMove = 1;}
        else {whichMove =0;}

        if (presentMarkN==M) { 
            // can not have add or swap
            if (whichMove==2) {whichMove = 3;}
            if (whichMove==1) {whichMove = 3;}
        } else if (presentMarkN==maxAllowedMarkers) { 
            // can not have add
            // If not limiting model space, maxAllowedMarkers will be set
            // greater than M
            if (whichMove==1) {whichMove = 3;}
        } else if (presentMarkN==0) { 
            // can not have rem or swap
            if (whichMove==2) {whichMove = 3;}
            if (whichMove==0) {whichMove = 3;}            
        }

        return whichMove;          
        
    }
        
    /***
     * Converts a random draw from a uniform distribution into a random draw
     * all positive integers up to a specified value.
     *
     * @param max Maximum index
     * @param unifDraw double between 0 and 1
     * @return whichIndex Selected index
     */
    public static int randomInteger(
            int max,
            double unifDraw) {
        
        int whichIndex = 0;
        for (int m=0; m<max; m++) {
            if ( (unifDraw>=((double)m/(double)max))&&(unifDraw
                    < ((double)(m+1)/(double)max)) ) {
                whichIndex = m;
            }
        }
        return whichIndex;

    }
      
    /**
     * Calculates the log-ratio of prior support for two models, under a
     * Poisson prior on model dimension. The results is used in the
     * calculation of the acceptance probability.
     * 
     * @param nVarsCurr Number of variables selected in the current model
     * @param nVarsProp Number of variables selected in the proposed model
     * @param modPriMean Mean of the Poisson prior on model dimension
     * @return Log-ratio of the prior support for both models. 
     */
    public static double logModelPriorRatio(int nVarsCurr, int nVarsProp,
          double modPriMean) {
      double logModelPriorRatio = 
              GeneralMaths.logFactorial(nVarsCurr) 
              - GeneralMaths.logFactorial(nVarsProp) 
              + Math.log(Math.pow(modPriMean, nVarsProp)) 
              - Math.log(Math.pow(modPriMean, nVarsCurr));
      return logModelPriorRatio ;
    }

  /**
   * Calculates the log-ratio of prior support for two models, under a
   * Beta-Binomial prior on model dimension (to be used in the
   * calculation of the acceptance probability), using the density defined in
   * Bottolo (6) ESS 2010 paper
   * 
   * This density is designed to be used when specifying a mean and variance
   * for the number of effects. In practice I found we can not make this
   * varianece small.
   * 
   * @param nVarsCurr Number of variables selected in the current model
   * @param nVarsProp Number of variables selected in the proposed model
   * @param modSpaceSize Total number of covariates
   * @param modPriBetaBinA Beta-Binomial prior hyper-parameter a
   * @param modPriBetaBinB Beta-Binomial prior hyper-parameter b
   * @return Log-ratio of the prior support for both models.
   */
  public static double OLD_logModelPriorRatioBetaBin_THIS_IS_WRONG(
          int nVarsCurr,
          int nVarsProp,
          int modSpaceSize,
          double modPriBetaBinA,
          double modPriBetaBinB) {
      
      double logModelPriorRatio = 0;
      // Work out difference in model size between iterations to take advantage
      // of cancelling factorials
      int sizeDiff = nVarsProp-nVarsCurr;
      if (sizeDiff==1) {
          // Log-ratio of beta-bimoial coefficients
          logModelPriorRatio =
                  Math.log(modSpaceSize + 1 - nVarsProp)
                  -Math.log(nVarsProp);
          // Log-ratio of Beta distibutions
          logModelPriorRatio = logModelPriorRatio                  
                  +Math.log(nVarsProp + modPriBetaBinA - 1)
                  -Math.log(modSpaceSize - nVarsCurr + modPriBetaBinB);
      } else if (sizeDiff==-1) {
          // Log-ratio of beta-bimoial coefficients
          logModelPriorRatio =
                  Math.log(nVarsCurr)
                  -Math.log(modSpaceSize + 1 - nVarsCurr);
          // Log-ratio of Beta distibutions
          logModelPriorRatio = logModelPriorRatio                  
                  -Math.log(nVarsCurr + modPriBetaBinA - 1)
                  +Math.log(modSpaceSize - nVarsProp + modPriBetaBinB);
      }
      return logModelPriorRatio ;
    }

  /**
   * Calculates the log-ratio of prior support for two model sizes, using the 
   * density given in Bottolo et al, 2010 eq (6). The density defined by
   * bottolo et al is the prior probability for a specific model. This is defined
   * in terms of the dimension of the model p_gamma:
   * 
   * P(gamma) = Beta(p_gamma + a, P - p_gamma + b)/ Beta(a,b)
   * P(gamma_prop)/p(gamma_curr) = 
   *    Beta(p_gamma_prop + a, P - p_gamma_prop + b)/ Beta(p_gamma_curr + a, P - p_gamma_curr + b)
   * 
   * The Beta(a,b) terms cancel when defining the ratios of these probabilities.
   * Note that:
   * 
   * Beta(x,y) = (x-1)!(y-1)!/(x+y-1)!
   * 
   * I worked through this - see PDF.
   * 
   * @param nVarsCurr Number of variables selected in the current model
   * @param nVarsProp Number of variables selected in the proposed model
   * @param modSpaceSize Total number of covariates
   * @param modPriBetaBinA Beta-Binomial prior hyper-parameter a
   * @param modPriBetaBinB Beta-Binomial prior hyper-parameter b
   * @return Log-ratio of the prior support for both models.
   */
  public static double logModelPriorRatioBetaBin_Bottolo_SpecificModels(
          int nVarsCurr,
          int nVarsProp,
          int modSpaceSize,
          double modPriBetaBinA,
          double modPriBetaBinB) {
      
      double logModelPriorRatio = 0;
      // Work out difference in model size between iterations to take advantage
      // of cancelling factorials
      int sizeDiff = nVarsProp-nVarsCurr;
      if (sizeDiff==1) {
          // Log-ratio of beta-binomial coefficients
          logModelPriorRatio =
                  Math.log(nVarsCurr + modPriBetaBinA)
                  -Math.log(modSpaceSize - nVarsCurr + modPriBetaBinB - 1);
      } else if (sizeDiff==-1) {
          logModelPriorRatio =
                  Math.log(modSpaceSize - nVarsCurr + modPriBetaBinB)
                  -Math.log(nVarsCurr + modPriBetaBinA - 1);
      }
      return logModelPriorRatio ;
    }
  
    /***
     * Counts the number of covariates included in a `model' described by
     * a vector of 0's and 1's.
     * 
     * @param M Total number of covariates.
     * @param model A vector of 0's and 1's indicating which covariates are present
     * 
     * @return The number of included variables. 
     */     
      public static int countPresVars(int M, int[] model) {
          int noPresentMarkers = 0;
          for (int m=0; m<M; m++ ) {
              if (model[m]==1) {noPresentMarkers++;}
          }
          return noPresentMarkers;
        }
      
    /***
     * Counts the number of covariates included in a `model' described by
     * a vector of 0's and 1's.
     * 
     * @param betas A vector of betas stored in a Jama Matrix object.
     * @param model A vector of 0's and 1's indicating which covariates are present
     * 
     * @return The number of included variables. 
     */     
      public static Matrix normaliseAbsoluteBetasToSumToOne(Matrix betas, int[] model) {
          // Calculate the normalising absolute total
          double sumOfAbsoluteValues = 0;
          for (int m=0; m<model.length; m++ ) {
              if (model[m]==1) {
                  sumOfAbsoluteValues += betas.get(m, 0)*Math.signum(betas.get(m, 0));}
          }
          // Divide through by the normalising total
          for (int m=0; m<model.length; m++ ) {
              if (model[m]==1) {
                  betas.set(m, 0, betas.get(m, 0)/sumOfAbsoluteValues);
              }
          }
          // Return
          return betas;
        }      

    /***
     * Converts the betas (log-scale for updating) to exponentiated weights,
     * whose absolute values sum to one, but for which the signs are preserved
     * 
     * @param betas A vector of betas stored in a Jama Matrix object.
     * @param model A vector of 0's and 1's indicating which covariates are present
     * 
     * @return A vector of weights in a Matrix object with absolute values 
     * summing to one. 
     */     
      public static Matrix EXTRA_betasToNormalisedWeights(Matrix betas, int[] model) {
          //Exponentiate while maintinaing the signs
          for (int m=0; m<model.length; m++ ) {
              if (model[m]==1) {
                  betas.set(m, 0, 
                          Math.signum(betas.get(m, 0))*
                                  Math.exp(Math.abs(betas.get(m, 0))) );
              }
          }
                System.out.println("betasExp "+betas.get(1,0));
          
          // Calculate the normalising constant - the sum of the absolute values
          double sumOfAbsoluteValues = 0;
          for (int m=0; m<model.length; m++ ) {
              if (model[m]==1) {
                  sumOfAbsoluteValues += Math.abs(betas.get(m, 0));}
          }
          
          // Divide through by the normalising constant
          for (int m=0; m<model.length; m++ ) {
              if (model[m]==1) {
                  betas.set(m, 0, betas.get(m, 0)/sumOfAbsoluteValues);
              }
          }
          System.out.println("betasNorm "+betas.get(1,0));
          
          // Return
          return betas;
        }      

    /***
     * Converts the normalised weights back to betas on the log-scale, while
     * preserving signs.
     * 
     * @param betas A vector of normalised weights stored in a Jama Matrix object.
     * @param model A vector of 0's and 1's indicating which covariates are present
     * 
     * @return A vector of betas on the log-scale
     */     
      public static Matrix EXTRA_normalisedWeightsToLogBetas(Matrix betas, int[] model) {
          //Take logs while maintinaing the signs
          for (int m=0; m<model.length; m++ ) {
              if (model[m]==1) {
                  betas.set(m, 0, 
                          Math.signum(betas.get(m, 0))*
                                  Math.log(Math.abs(betas.get(m, 0))) );
              }
          }
          
          // Return
          return betas;
        }      
      
    /***
     * Gets the number of covariates included in a `model' described by
     * a vector of 0's and 1's.
     * 
     * @param modelDimension Total number of included covariates.
     * @param model A vector of 0's and 1's indicating which covariates are present
     * 
     * @return The number of included variables. 
     */     
      public static int[] getVarIndices(int modelDimension, int[] model) {
          int[] varIndices = new int[modelDimension];
          int enteredVars = 0;
          int var=0;
          while (enteredVars < modelDimension) {
              if (model[var]==1) {
                  varIndices[enteredVars]=var;
                  enteredVars++;
              }
              var++;
          }
          return varIndices;
      }
      
    /**
     * Counts the number of covariates included in a `model' described by
     * a vector of 0's and 1's, by each different model space component.
     * 
     * @param Ncomps Number of model space components.
     * @param compSplits Covariate inidices which partition into the different
     * model space components.
     * @param model A vector of 0's and 1's indicating which covariates are present
     * 
     * @return The number of included variables in each component. 
     */     
      public static int[] countPresVarsComps(
              int Ncomps,
              int[] compSplits,
              int[] model) {
          int[] noPresentMarkers = new int[Ncomps];
          for (int c=0; c<Ncomps; c++) {
              noPresentMarkers[c]=0;
              for (int m=compSplits[c]; m<compSplits[c+1]; m++ ) {
                  if (model[m]==1) {noPresentMarkers[c]++;}
              }
          }
          return noPresentMarkers;
      }
      
    /**
     * Counts the number of covariates included in a `model' described by
     * a vector of 0's and 1's, by each different model space component.
     * 
     * @param Nblocks Number of blocks in the X matrix.
     * @param blockIndices Vector of indices for the different blocks
     * @param whichMove Move type - to determine if swap and therefore 2 blocks
     * may have changed.
     * @param whichBetaRemoved Beta removed
     * @param whichBetaAdded Beta added
     * @param whichBetaUpdated Beta added
     * 
     * @return A vector indicating which blocks have been updated
     */     
      public static boolean[] determineUpdatedBlocks(
              int Nblocks,
              int[] blockIndices,
              int whichMove,
              int whichBetaRemoved,
              int whichBetaAdded,
              int whichBetaUpdated) {
          boolean[] whichBlocksUpdated = new boolean[Nblocks];
          for (int b=0; b<Nblocks; b++) {
              if (whichMove==0|whichMove==2) {
                  if ((whichBetaRemoved>=blockIndices[b])&(whichBetaRemoved<blockIndices[(b+1)])) {
                      whichBlocksUpdated[b] = true;
                  }
              }
              if (whichMove==1|whichMove==2) {
                  if ((whichBetaAdded>=blockIndices[b])&(whichBetaAdded<blockIndices[(b+1)])) {
                      whichBlocksUpdated[b] = true;
                  }
              }
              if (whichMove==3) {
                  if ((whichBetaUpdated>=blockIndices[b])&(whichBetaUpdated<blockIndices[(b+1)])) {
                      whichBlocksUpdated[b] = true;
                  }                  
              }              
          }
          return whichBlocksUpdated;
      } 
}
