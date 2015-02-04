package Methods;

import java.util.Random;

/**
 * Collection of generic methods.
 * 
 * @author Paul J Newcombe
 */
public class GeneralMethods {    
    /**
     * Chooses a Reversible Jump move type.
     * 
     * @param M Total number of variables.
     * @param presentMarkN Number of variables currently included.
     * @param moveProbabilities Vector of probabilities with which to select the
     * different move types.
     * @param randomDraws Random number generator object.
     * 
     * @return Indicates chosen move type (0=Removal, 1=Addition, 2=Swap, 3=Null)
     */
    public static int chooseMove(
            int M,
            int maxAllowedMarkers,
            int presentMarkN,
            double[] moveProbabilities,
            Random randomDraws) {
        int whichMove;
        double eventDraw = randomDraws.nextFloat();
        if (eventDraw>=moveProbabilities[2]) {whichMove =3;}
        else if (eventDraw>=moveProbabilities[1]) {whichMove = 2;}
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
   * calculation of the acceptance probability).
   * 
   * @param nVarsCurr Number of variables selected in the current model
   * @param nVarsProp Number of variables selected in the proposed model
   * @param modSpaceSize Total number of covariates
   * @param modPriBetaBinA Beta-Binomial prior hyper-parameter a
   * @param modPriBetaBinB Beta-Binomial prior hyper-parameter b
   * @return Log-ratio of the prior support for both models.
   */
  public static double logModelPriorRatioBetaBin(
          int nVarsCurr,
          int nVarsProp,
          int modSpaceSize,
          double modPriBetaBinA,
          double modPriBetaBinB) {
      // BetaBinomial for proposal is:
      // Beta(nVarsCurr+modPriBetaA, modSpaceSize-nVarsCurr+modPriBetaB)
      // = (nVarsCurr+modPriBetaA-1)!(modSpaceSize-nVarsCurr+modPriBetaB-1)!
      //    /(modSpaceSize+modPriBetaA+modPriBetaB-1)!
      // Can take advantage of cancelling factorials
      //
      // Numerator - corresponds to proposed likelihood/prior
      // Denominator - corresponds to curren likelihood/prior
      double logModelPriorRatio = 0;
      // Work out difference in model size between iterations to take advantage
      // of cancelling factorials
      int sizeDiff = nVarsProp-nVarsCurr;
      if (sizeDiff==1) {
          // addition
          //(nVarsProp+modPriBetaA-1)
          // / (modSpaceSize-nVarsCurr+modPriBetaB-1)
          logModelPriorRatio =
                  Math.log((nVarsProp+modPriBetaBinA-1))
                  -Math.log((modSpaceSize-nVarsCurr+modPriBetaBinB-1));
      } else if (sizeDiff==-1) {
          // removal
          //(modSpaceSize-nVarsProp+modPriBetaB-1)
          // / (nVarsCurr+modPriBetaA-1)
          logModelPriorRatio =
                  Math.log((modSpaceSize-nVarsProp+modPriBetaBinB-1))
                  -Math.log((nVarsCurr+modPriBetaBinA-1));
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
