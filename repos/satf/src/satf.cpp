

#include <math.h>
#include <assert.h>

inline double SATF(double t, double lambda, double beta, double delta) {
   return( (t >= delta)*(lambda*(1-exp(-(1/beta)*(t-delta)))) );
}

inline int index(int row, int col, int nRows) {
   return( (col*nRows + row)-1 );
}

inline double get_param(int* rows, int rowsLen, double *dm,
                        int dmRows, int col, double* coefs)  {
   double param = 0.0;
   for(int i=0; i < rowsLen; i++) {
       int row = rows[i];
       double contrast = dm[index(row, col, dmRows)];
       if(contrast) {
          param += coefs[row-1]*contrast;
       }
   }
   return(param);
}

static double pnorm(double d, double mean)
{
    const double A1 = 0.31938153;
    const double A2 = -0.356563782;
    const double A3 = 1.781477937;
    const double A4 = -1.821255978;
    const double A5 = 1.330274429;
    const double RSQRT2PI = 0.39894228040143267793994605993438;
    d = d - mean;
    double K = 1.0 / (1.0 + 0.2316419 * fabs(d));
    double cnd = RSQRT2PI * exp(- 0.5 * d * d) *
          (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));
    if (d > 0)
        cnd = 1.0 - cnd;
    return(cnd);
}

double log_dnorm(double x, double mean, double sd ) {
 return( -log(sqrt(2*M_PI*pow(sd,2))) - pow((x-mean),2)/(2*pow(sd,2)) );
}

double dnorm(double x, double mean, double sd ) {
  return( ( 1/sqrt(2*M_PI*pow(sd,2)) )*exp( -pow((x-mean),2)/(2*pow(sd,2)) ) );
}




/* Approximation of the conditional bivariate normal density based on the following
   paper: http://www.dtic.mil/dtic/tr/fulltext/u2/a125033.pdf (approximation to X, 
   given Y < b, best when correlation is low)
   Background: http://www.nasa.gov/centers/dryden/pdf/87929main_H-1120.pdf (density of X, 
   given Y < b).  Auch: Continuous Multivariate Distributions, Models and Applications (Kotz et al.)
*/

double conditional_pnorm(double nboundary, double rho,
                         double lboundary, int lboundary_upper)
{
  int mirror = (lboundary_upper*2-1);
  rho = rho*mirror;
  lboundary = lboundary*mirror;
  double mu = -rho*dnorm(lboundary, 0, 1)/pnorm(lboundary, 0);
  double sigma = sqrt(1 + rho*lboundary*mu - pow(mu, 2) );
  return( pnorm((nboundary-mu)/sigma, 0) );
}

double log_dgamma(double x, double scale, double shape)
{
  if(scale<0) {
    scale = abs(scale);
    x = -x;
  }
   return( -lgamma(shape) - shape*log(scale) + (shape-1)*log(x) -x/scale );
}

void* duplicate(void *src, size_t size) {
  void *dst = malloc(size);
  memcpy (dst, src, size );
  return(dst);
}
#define DUPLICATE(SRC, SIZE, TYPE) (TYPE*) duplicate(SRC, SIZE*sizeof(TYPE))
#define ALLOC(SIZE, TYPE) (TYPE*) malloc(SIZE*sizeof(TYPE))

class AllArgs 
{
  public:
   double *dm; int dmRows; int dmCols;
   int *transformVec; double *paramsMinimum; double *paramsMaximum;
   int *contrastIDs; int maxContrastID;
   int *rowsLambda; int rowsLambdaLen;
   int *rowsBeta; int rowsBetaLen;
   int *rowsDelta; int rowsDeltaLen; 
   double *time;
   double *predictedCriterion;
   int responseVariable; int *responseYes; double* responseDprime;
   int *trialId; int *responseHits; int *nGrammatical; double *bGrammatical;
   int *responseFAs; int *nUngrammatical; double *bUngrammatical;
   double *hyperparams1; double *hyperparams2;
   bool fitByLL;
};


// TODO: Fix the interaction of maximim and minimum values with transformations here: Currently, it is a hack.

void coefs_transform(double* coefs, int* transform_vec, double* minimum, double* maximum, int n_coefs) {
  for(int i=0; i < n_coefs; i++) {
    if(transform_vec[i] != 0) {
      coefs[i] = pow(coefs[i],2)*transform_vec[i];
    }
    if(coefs[i] < minimum[i]) {
      coefs[i] = minimum[i];
    }
     else if(coefs[i] > maximum[i]) {
      coefs[i] = maximum[i];
    }
  }
}

double coefs_LL(double *coefs, AllArgs &args)
{
  double LL=0;
  if(!args.hyperparams1)
    return(0);
  
  for(int i=0; i < args.dmRows; i++)
  {
    double hyperparam1, hyperparam2, coef;
    // LL += log_dnorm(coefs[i], args.hyperparams1[i], args.hyperparams2[i]);
    hyperparam2 = args.hyperparams2[i];
    
    if(args.hyperparams1[i] < 0) {
      hyperparam1 = -1*args.hyperparams1[i];
      coef = -1*coefs[i];
    } else {
      hyperparam1 = args.hyperparams1[i];
      coef = coefs[i];
    }
    LL += log_dgamma(coef, hyperparam1, hyperparam2);
    printf("LL <%f, %f, %f> %f\\n", coef, hyperparam1, hyperparam2, LL);
  }
  return(LL);
}

inline double logodds2p(double lodds) { return( exp(lodds)/(1+exp(lodds)) ); }

#define MAX_CORRELATION .96

double compute_data_fit(double *coefs, int coefsLen, double *LLparams,
                        int fitCorrelation, AllArgs &args, bool fitByLL)
{
  double LL=0; 
  bool *paramsForID = ALLOC(args.maxContrastID, bool);
  double *lambdas = ALLOC(args.maxContrastID, double);
  double *betas = ALLOC(args.maxContrastID, double);
  double *deltas = ALLOC(args.maxContrastID, double);
  double last_dprime = 0;
  double withinTrialCorrelation = 0;
  double SSresidual = 0;
  double SStotal = 0;

  if(coefsLen > args.dmRows) {
    withinTrialCorrelation = logodds2p(coefs[args.dmRows]);
    if(withinTrialCorrelation > MAX_CORRELATION) {
      withinTrialCorrelation = MAX_CORRELATION;
    }
  }
  for(int i=0; i < args.maxContrastID; i++) {
    paramsForID[i] = FALSE;
  }
 
  coefs_transform(coefs, args.transformVec, args.paramsMinimum, args.paramsMaximum, args.dmRows);
  *LLparams = coefs_LL(coefs, args);
  LL += *LLparams;

  double meanDprime = 0;
 
  if(!fitByLL) {
    double sumDprime = 0;
    for(int col=0; col < args.dmCols; col++) {
      sumDprime += args.responseDprime[col];
    }
    meanDprime = sumDprime/args.dmCols;
  }

  for(int col=0; col < args.dmCols; col++)
  {
    int contrastID = args.contrastIDs[col]-1;
    if(!paramsForID[contrastID]) {
      lambdas[contrastID] = get_param(args.rowsLambda, args.rowsLambdaLen,
                                      args.dm, args.dmRows, col, coefs);
      betas[contrastID] = get_param(args.rowsBeta, args.rowsBetaLen,
                                    args.dm, args.dmRows, col, coefs);
      deltas[contrastID] = get_param(args.rowsDelta, args.rowsDeltaLen,
                                     args.dm, args.dmRows, col, coefs);
      paramsForID[contrastID] = TRUE;
    }
    double dprime = SATF(args.time[col], lambdas[contrastID],
                         betas[contrastID], deltas[contrastID]);
    // printf("<%f>,<%f>,<%f>,<%f>, dprime <%f>\\n", args.time[col], lambdas[contrastID],
    //                      betas[contrastID], deltas[contrastID], dprime);
    if(args.responseVariable==1) 
    {
      if(fitByLL)
      {
        bool dependent = FALSE;
        if(args.trialId && col > 0 && withinTrialCorrelation != 0) {
          assert(trialNoise != 0);
          dependent = (args.trialId[col] == args.trialId[col-1]);
        }
        if(dependent)
        {
          double curLL;
          double CriterionMinusPsi = args.predictedCriterion[col] - dprime;
          double LastCriterionMinusPsi = args.predictedCriterion[col-1] - last_dprime;

          int lboundary_upper = 1 - args.responseYes[col-1];
          double p_No = conditional_pnorm(CriterionMinusPsi, withinTrialCorrelation,
                                          LastCriterionMinusPsi, lboundary_upper);
          if(args.responseYes[col] == 1) curLL = log(1-p_No);
          else                           curLL = log(p_No);
          LL += curLL;
            
        } else {
          int response = args.responseYes[col];  
          double p_No = pnorm(args.predictedCriterion[col], dprime);
          double curLL = log( response - (2*response-1)*p_No );
          LL += curLL;
        }
      } else
      {
        SSresidual += pow(dprime-args.responseDprime[col],2);
        SStotal += pow(dprime-meanDprime,2);
      }
    }
    else
    {
      double p_No = pnorm(args.predictedCriterion[col], dprime);
      int k = args.responseHits[col];
      int n = args.nGrammatical[col];
      LL += args.bGrammatical[col] + k*log(1-p_No) + (n-k)*log(p_No);
      if(args.nUngrammatical[col] != 0 ) {
        double p_No = pnorm(args.predictedCriterion[col], 0);
        int k = args.responseFAs[col];
        int n = args.nUngrammatical[col];
        LL += args.bUngrammatical[col] + k*log(1-p_No) + (n-k)*log(p_No);
      }
    }
    last_dprime = dprime;
  }

  if(!fitByLL) {
    // n.data = args.dmCols
    // n.params = args.dmRows
    double adjR2numerator = SSresidual / (args.dmCols-args.dmRows);
    double adjR2denominator = SStotal / (args.dmCols-1);
    LL = (1-adjR2numerator/adjR2denominator);
  }
  
  free(paramsForID); free(lambdas); free(betas); free(deltas);
  return(LL);
}




.sig_predict_dprime <- signature(initialize="integer", 
                                 fixedparams="double", coefs="double", coefsLen="int",
                                 dm="double", dmRows="integer", dmCols="integer",
                                 transformVec="integer", paramsMinimum="double", paramsMaximum="double", 
                                 contrastIDs="integer", maxContrastID="integer",
                                 rowsLambda="integer", rowsLambdaLen="integer", 
                                 rowsBeta="integer", rowsBetaLen="integer", 
                                 rowsDelta="integer", rowsDeltaLen="integer", 
                                 time="double", predictedCriterion="double",
                                 responseVariable="integer", responseYes="int",
                                 responseDprime="double", trialId="integer",
                                 responseHits="integer", nGrammatical="integer",
                                 bGrammatical="double",
                                 responseFAs="integer", nUngrammatical="integer",
                                 bUngrammatical="double",
                                 hyperparams="integer", hyperparams1="double", hyperparams2="double",
                                 LL="double", LLparams="double",
                                 fitCorrelation="integer")
.code_predict_dprime <- '
  AllArgs *args = NULL;
  double *tst = NULL;
  if(*initialize == 1) {
    args = new AllArgs();
    assert(sizeof(AllArgs *) == sizeof(double));
    memcpy(fixedparams, &args, sizeof(double));
    args->dmRows = *dmRows; args->dmCols = *dmCols;
    args->dm = DUPLICATE(dm, (args->dmRows)*(args->dmCols), double);
    args->transformVec = DUPLICATE(transformVec, args->dmCols, int);
    args->paramsMinimum = DUPLICATE(paramsMinimum, args->dmCols, double);
    args->paramsMaximum = DUPLICATE(paramsMaximum, args->dmCols, double);
    args->contrastIDs = DUPLICATE(contrastIDs, args->dmCols, int);
    args->maxContrastID = *maxContrastID;
    args->rowsLambdaLen=*rowsLambdaLen;
    args->rowsLambda = DUPLICATE(rowsLambda, args->rowsLambdaLen, int);
    args->rowsBetaLen=*rowsBetaLen;
    args->rowsBeta = DUPLICATE(rowsBeta, args->rowsBetaLen, int);
    args->rowsDeltaLen=*rowsDeltaLen;
    args->rowsDelta = DUPLICATE(rowsDelta, args->rowsDeltaLen, int);
    args->time=DUPLICATE(time, args->dmCols, double);
    args->predictedCriterion=DUPLICATE(predictedCriterion, args->dmCols, double);
    
    args->responseVariable = *responseVariable;
    if(args->responseVariable==1)
    {
      if(responseYes) {
        args->responseYes=DUPLICATE(responseYes, args->dmCols, int);
        args->fitByLL = true;
      } else {
        args->responseDprime=DUPLICATE(responseDprime, args->dmCols, double);
        args->fitByLL = false;
      }
      if(trialId != NULL)
        args->trialId=DUPLICATE(trialId, args->dmCols, int);
      else 
        args->trialId=NULL;
    }
    else {
      args->fitByLL = true;
      args->responseHits=DUPLICATE(responseHits, args->dmCols, int);
      args->nGrammatical=DUPLICATE(nGrammatical, args->dmCols, int);
      args->bGrammatical=DUPLICATE(bGrammatical, args->dmCols, double);
      args->responseFAs=DUPLICATE(responseFAs, args->dmCols, int);
      args->nUngrammatical=DUPLICATE(nUngrammatical, args->dmCols, int);
      args->bUngrammatical=DUPLICATE(bUngrammatical, args->dmCols, double);
    }
    if(*hyperparams) {
      args->hyperparams1=DUPLICATE(hyperparams1, args->dmRows, double);
      args->hyperparams2=DUPLICATE(hyperparams2, args->dmRows, double);
    } else {
      args->hyperparams1=NULL;
      args->hyperparams2=NULL;
    }
    return;
  } else {
    args = *((AllArgs**)fixedparams);    
  }
  if( args == NULL) {
     printf("error 1\\n");
    *LL = 0;
    return;
  }
  if(*initialize == 2) {
    free(args->dm);
    free(args->transformVec); free(args->paramsMinimum); free(args->paramsMaximum);
    free(args->contrastIDs);
    free(args->rowsLambda); free(args->rowsBeta);
    free(args->rowsDelta); free(args->time);
    free(args->predictedCriterion);
    if(args->responseVariable==1) {
      free(args->responseYes);
      free(args->trialId);
    } else {
      free(args->responseHits); free(args->nGrammatical);
      free(args->responseFAs);  free(args->nUngrammatical);
      free(args->bGrammatical); free(args->bUngrammatical);
    }
    if(args->hyperparams1) {
      free(args->hyperparams1);
      free(args->hyperparams2);
    }
    delete args;
    return;
  }
  *LL = compute_data_fit(coefs, *coefsLen, LLparams, *fitCorrelation, *args, args->fitByLL);
'
.predict_dprime <- cfunction(.sig_predict_dprime, .code_predict_dprime,
                             includes=.includes_predict_dprime, convention=".C")

