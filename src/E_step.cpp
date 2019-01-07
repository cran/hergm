// Files to look to get MM version: ReciprocityModel for basic functions, 
// BinaryReciprocityModel.cpp for version written already, 
// MMBinaryReciprocityModel.cpp for more eleborate version
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

const double minTau = 1e-6;
double minPi =  1e-4;

void logMatrix(NumericMatrix& pi, NumericMatrix& logPi) {
  for (int k = 0; k < pi.nrow(); k++)
    for (int l = 0; l < pi.ncol(); l++)
      logPi(k, l) = log(pi(k, l));
}


void logTransposedMatrix(NumericMatrix& pi,NumericMatrix& logPi) {
  for ( int k = 0; k < pi.nrow(); k++)
    for (int l = 0; l < pi.ncol(); l++)
      logPi(l, k) = log(pi(k, l));
}

void normalizeVector(NumericVector& vector, double minValue) {
  for (int k = 0; k < vector.size(); k++)
    if (vector(k) < minValue)
      vector(k) = minValue;
    double sumAlpha = 0;
    for (int k = 0; k < vector.size(); k++)
      sumAlpha += vector(k);
    for (int k = 0; k < vector.size(); k++)
      vector(k) /= sumAlpha;
}

void sumDoubleMatrixByRow(const NumericMatrix & matrix,
                          NumericVector & vector) {
  for (int j = 0; j < matrix.ncol(); j++) {
    vector(j) = 0;
    for (int i = 0; i < matrix.nrow(); i++)
      vector(j) += matrix(i, j);
  }
}

void normalizeLogTau2Tau(NumericMatrix& tau, double minValue) {
  int numOfVertices = tau.nrow();
  int numOfClasses = tau.ncol();
  
  //const double logDoubleMax = log(DBL_MAX) - 1;
  const double logDoubleMax = 100;

  for (int i = 0; i < numOfVertices; i++) {
    // First it holds the max value
    double slidingValue = tau(i, 0);
    for (int k = 1; k < numOfClasses; k++)
      if (slidingValue < tau(i, k))
        slidingValue = tau(i, k);
      // Now it actually holds the sliding value
      slidingValue = logDoubleMax - slidingValue;
      for (int k = 0; k < numOfClasses; k++)
        tau(i, k) += slidingValue;
  }
  
  // normalize
  for (int i = 0; i < numOfVertices; i++) {
    double denominator = 0;
    for (int k = 0; k < numOfClasses; k++) {
      tau(i, k) = exp(tau(i, k));
      denominator += tau(i, k);
    }
    bool again = false;
    for (int k = 0; k < numOfClasses; k++) {
      tau(i, k) /= denominator;
      if (tau(i, k) < minValue) {
        tau(i, k) = minValue;
        again = true;
      }
    }
    if (again) {
      denominator = 0;
      for (int k = 0; k < numOfClasses; k++) {
        denominator += tau(i, k);
      }
      for (int k = 0; k < numOfClasses; k++)
        tau(i, k) /= denominator;
    }
  }
  
}

void normalizeTau(NumericMatrix& tau, double minValue) {
  int numOfVertices = tau.nrow();
  int numOfClasses = tau.ncol();
  // normalize
  for (int i = 0; i < numOfVertices; i++) {
    double denominator = 0;
    for (int k = 0; k < numOfClasses; k++) {
      denominator += tau(i, k);
    }
    bool again = false;
    for (int k = 0; k < numOfClasses; k++) {
      tau(i, k) /= denominator;
      if (tau(i, k) < minValue) {
        tau(i, k) = minValue;
        again = true;
      }
    }
    if (again) {
      denominator = 0;
      for (int k = 0; k < numOfClasses; k++) {
        denominator += tau(i, k);
      }
      for (int k = 0; k < numOfClasses; k++)
        tau(i, k) /= denominator;
    }
  }
}

void solveQP(const NumericMatrix& m, const NumericMatrix& s,
             NumericMatrix& tau, double precision) {
  
  if (m.nrow() <= 0)
    return;
  
  int n = m.ncol();
  
  
  bool * J_k = new bool[n]();
  bool * J_lambda_a = new bool[n]();
  
  bool * J_lambda_k = new bool[n]();
  bool * J_lambda_k_a = new bool[n]();
  bool * J_lambda_k_b = new bool[n]();
  
  double lambda_k;
  
  for (int rowIndex = 0; rowIndex < m.nrow(); rowIndex++) {
    // Step 1:
    for (int j = 0; j < n; j++) {
      J_k[j] = true;
      J_lambda_a[j] = false;
    }
    
    bool done = false;
    // we never loop over n times.
    int count = 0;
    while (!done) {
      count++;
      
      // Step 2:
      double value1 = 0, value2 = 0;
      for (int j = 0; j < n; j++)
        if (J_k[j]) {
          value1 += 1 / m(rowIndex, j);
          value2 += s(rowIndex, j) / m(rowIndex, j);
        }
        lambda_k = (value2 - 2) / value1;
        
        // Step3:
        for (int j = 0; j < n; j++) {
          if (J_k[j]) {
            if (lambda_k >= s(rowIndex, j)) {
              J_lambda_k[j] = false;
              J_lambda_k_a[j] = true;
              J_lambda_k_b[j] = false;
            } else if (lambda_k
                         > (-2 * m(rowIndex, j) + s(rowIndex, j))) {
              J_lambda_k[j] = true;
              J_lambda_k_a[j] = false;
              J_lambda_k_b[j] = false;
            } else {
              J_lambda_k[j] = false;
              J_lambda_k_a[j] = false;
              J_lambda_k_b[j] = true;
            }
          } else {
            J_lambda_k[j] = false;
            J_lambda_k_a[j] = false;
            J_lambda_k_b[j] = false;
          }
        }
        
        // Step 4:
        double delta = 0;
        double value3 = 0, value4 = 0;
        bool J_lambda_k_empty = true;
        for (int j = 0; j < n; j++) {
          delta += J_lambda_k_b[j];
          if (J_lambda_k[j]) {
            value3 += s(rowIndex, j) / m(rowIndex, j);
            value4 += 1 / m(rowIndex, j);
            J_lambda_k_empty = false;
          }
        }
        delta += (value3 / 2 - lambda_k * value4 / 2 - 1);
        
        if (fabs(delta) < precision || J_lambda_k_empty || count >= n) {
          for (int j = 0; j < n; j++) {
            if (J_lambda_a[j] || J_lambda_k_a[j])
              tau(rowIndex, j) = 0;
            else if (J_lambda_k_b[j])
              tau(rowIndex, j) = 1;
            else
              tau(rowIndex, j) = (s(rowIndex, j) - lambda_k) / (2
                                                                  * m(rowIndex, j));
          }
          done = true;
        } else if (delta > 0) {
          for (int j = 0; j < n; j++) {
            if (J_lambda_k_a[j]) {
              J_lambda_a[j] = true;
              J_k[j] = false;
            }
          }
        } else {
          for (int j = 0; j < n; j++) {
            if (J_lambda_k_b[j])
              tau(rowIndex, j) = 1;
            else
              tau(rowIndex, j) = 0;
          }
          done = true;
        }
    }
  }
  // Don't forget to delete the bool arrays 
  delete [] J_k;
  delete [] J_lambda_a;

  delete [] J_lambda_k;
  delete [] J_lambda_k_a;
  delete [] J_lambda_k_b; 
}

// [[Rcpp::export]]
bool isTauSignificantlyChanged(double tauPrecision, NumericMatrix tau, NumericMatrix prevTau) 
{
  int numOfVertices = tau.nrow();
  int numOfClasses = tau.ncol();
  for (int i = 0; i < numOfVertices; i++)
    for (int k = 0; k < numOfClasses; k++)
      if (fabs(tau(i, k) - prevTau(i, k)) > tauPrecision)
        return true;
      return false;
}


void updateTau(NumericMatrix& new_tau, NumericMatrix& stat, const NumericMatrix& tau, 
               const NumericMatrix& logPi, NumericMatrix& temp, int numOfVertices, int numOfClasses) {
  
  // Iterates through non-zero values of stat
  
  for (int i = 0; i < numOfVertices; i++)
    for (int l = 0; l < numOfClasses; l++)
      temp(i, l) = 0;
  //Rcpp::print(temp);
  for (int i = 0; i < numOfVertices; i++) {
    for (int j = 0; j < numOfVertices; j++) {
      if (stat(i,j) != 0) {
        for (int l = 0; l < numOfClasses; l++)
          temp(i, l) += tau(j, l);
      }
    }
  }
  //Rcpp::print(temp);
  for (int i = 0; i < numOfVertices; i++) 
    for (int k = 0; k < numOfClasses; k++) 
      for (int l = 0; l < numOfClasses; l++) 
        new_tau(i, k) += logPi(k, l) * temp(i, l);
  

}

void updateTauByNegativeReflection(NumericMatrix& new_tau, NumericMatrix& stat,
                                   const NumericMatrix& tau, const NumericMatrix& logPi, NumericMatrix& temp,
                                   int numOfVertices, int numOfClasses) {
  
  NumericVector tauL(numOfClasses);
  sumDoubleMatrixByRow(tau, tauL);
  
  for (int i = 0; i < numOfVertices; i++)
    for (int l = 0; l < numOfClasses; l++)
      temp(i, l) = tauL(l);
  
  for (int i = 0; i < numOfVertices; i++) {
    for (int j = 0; j < numOfVertices; j++) {
      if (stat(i,j) != 0) {
        for (int l = 0; l < numOfClasses; l++)
          temp(i, l) -= tau(j, l);
      }
    }
  }
  
  for (int i = 0; i < numOfVertices; i++)
    for (int k = 0; k < numOfClasses; k++)
      for (int l = 0; l < numOfClasses; l++)
        new_tau(i, k) += logPi(k, l) * temp(i, l);
}

// [[Rcpp::export]]
NumericMatrix easy_E_step(int numOfVertices, int numOfClasses,
                          NumericVector alpha, NumericMatrix pi, NumericMatrix stat,
                          NumericMatrix tau) {
  
  NumericMatrix new_tau(numOfVertices, numOfClasses);
  for (int i = 0; i < numOfVertices; i++)
    for (int k = 0; k < numOfClasses; k++)
      new_tau(i, k) = log(alpha(k));
  
  NumericMatrix temp(numOfVertices, numOfClasses);
  
  
  NumericMatrix logPi(numOfClasses, numOfClasses);
  logMatrix(pi, logPi);
  
  updateTau(new_tau, stat, tau, logPi, temp, numOfVertices, numOfClasses);
  
  normalizeLogTau2Tau(new_tau, minTau);
  
  return new_tau;
}

// [[Rcpp::export]]
NumericMatrix runFixedPointEstimationEStep(int numOfVertices, int numOfClasses,
                                           NumericVector alpha, NumericMatrix pi, NumericMatrix stat00,
                                           NumericMatrix stat01, NumericMatrix stat10, NumericMatrix stat11,
                                           NumericMatrix tau) {
  
  NumericMatrix new_tau(numOfVertices, numOfClasses);
  for (int i = 0; i < numOfVertices; i++)
    for (int k = 0; k < numOfClasses; k++)
      new_tau(i, k) = log(alpha(k));
  
  NumericMatrix temp(numOfVertices, numOfClasses);
  
  
  NumericMatrix pi10(numOfClasses, numOfClasses);
  NumericMatrix pi11(numOfClasses, numOfClasses);
  NumericMatrix pi00(numOfClasses, numOfClasses);
  for (int k = 0; k < numOfClasses; k++) {
    for (int l = 0; l < numOfClasses; l++) {
      pi11(k,l) = pi(k,l)*pi(l,k);
      pi10(k,l) = pi(k,l)*(1 - pi(l,k));
      pi00(k,l) = 1 - pi11(k,l) - pi10(k,l) - pi10(l,k);
    }
  }
  
  NumericMatrix logPi10(numOfClasses, numOfClasses);
  logMatrix(pi10, logPi10);
  
  updateTau(new_tau, stat10, tau, logPi10, temp, numOfVertices, numOfClasses);
  
  NumericMatrix logPi10_t(numOfClasses, numOfClasses);
  logTransposedMatrix(pi10, logPi10_t);
  updateTau(new_tau, stat01, tau, logPi10_t, temp, numOfVertices, numOfClasses);

  NumericMatrix logPi11(numOfClasses, numOfClasses);
  logMatrix(pi11, logPi11);
  updateTau(new_tau, stat11, tau, logPi11, temp, numOfVertices, numOfClasses);

  NumericMatrix logPi00(numOfClasses, numOfClasses);
  logMatrix(pi00, logPi00);
  updateTauByNegativeReflection(new_tau, stat00, tau, logPi00, temp, numOfVertices, numOfClasses);

  normalizeLogTau2Tau(new_tau, minTau);

  return new_tau;
}

// [[Rcpp::export]]
List calculateStats(NumericMatrix& network, NumericMatrix& stat00, NumericMatrix& stat01,
                    NumericMatrix& stat10, NumericMatrix& stat11){
  int numOfVertices = network.nrow();
  for (int i = 0; i < numOfVertices; i++) {
    for (int j = 0; j < numOfVertices; j++) {
      if (network(i,j) == 1 && network(j,i) == 0) {
        stat10(i, j) =  1;
        stat00(i, j) =  1;
        stat01(j, i) =  1;
        stat00(j, i) =  1;
      }
      // (1,1): stat11 is symmetric
      if (i < j && network(i,j)== 1 && network(j,i) == 1) {
        stat11(i, j) =  1;
        stat00(i, j) =  1;
        stat11(j, i) =  1;
        stat00(j, i) =  1;
      }
    }
  }
  return List::create(stat00,stat01,stat10,stat11);
}

// [[Rcpp::export]]     
NumericMatrix runFixedPointEstimationEStepMM(int numOfVertices, int numOfClasses,
                                             NumericVector alpha, NumericMatrix pi, 
                                             NumericMatrix stat00, NumericMatrix stat01,
                                             NumericMatrix stat10, NumericMatrix stat11,
                                             NumericMatrix tau, NumericMatrix& network)
{
  //cout << "MMBinaryReciprocityModel::runFixedPointEstimationEStep()" << endl;
  
  NumericMatrix A = NumericMatrix(numOfVertices, numOfClasses);
  NumericMatrix s = NumericMatrix(numOfVertices, numOfClasses);
  
  //Calculate pi's for reciprocity model
  NumericMatrix pi10(numOfClasses, numOfClasses);
  NumericMatrix pi11(numOfClasses, numOfClasses);
  NumericMatrix pi00(numOfClasses, numOfClasses);
  for (int k = 0; k < numOfClasses; k++) {
    for (int l = 0; l < numOfClasses; l++) {
      pi11(k,l) = pi(k,l)*pi(l,k);
      pi10(k,l) = pi(k,l)*(1 - pi(l,k));
      pi00(k,l) = 1 - pi11(k,l) - pi10(k,l) - pi10(l,k);
    }
  }
  
  
  // Calculate the quadratic coefficients
  
  // Compute the norm term, i.e. \pi_kl^0
  NumericMatrix logPi00 = NumericMatrix(numOfClasses, numOfClasses);
  logMatrix(pi00, logPi00);
  NumericVector tauL = NumericVector(numOfClasses);
  sumDoubleMatrixByRow(tau, tauL);
  for (int i = 0; i < numOfVertices; i++)
    for (int k = 0; k < numOfClasses; k++) {
      A(i, k) = 0;
      for (int l = 0; l < numOfClasses; l++)
        A(i, k) += (logPi00(k, l) * (tauL(l) - tau(i, l)));
    }
    
    
    // y_ij = 1, y_ji = 0 AND y_ij = 0, y_ji = 1
    NumericMatrix logPi10 = NumericMatrix(numOfClasses, numOfClasses);
  for (int k = 0; k < numOfClasses; k++)
    for (int l = 0; l < numOfClasses; l++)
      logPi10(k, l) = log(pi10(k, l) / pi00(k, l));
  
  for (int i = 0; i < numOfVertices; i++) {
    for (int j = 0; j < numOfVertices; j++) {
      if (network(i,j) == 1 && network(j,i) == 0) {
        for (int k = 0; k < numOfClasses; k++)
          for (int l = 0; l < numOfClasses; l++) {
            A(i, k) += tau(j, l) * logPi10(k, l);
            A(j, l) += tau(i, k) * logPi10(l, k);
          }
      }
    }
  }
  
  // y_ij = 1, y_ji = 1
  NumericMatrix logPi11 = NumericMatrix(numOfClasses, numOfClasses);
  for (int k = 0; k < numOfClasses; k++)
    for (int l = 0; l < numOfClasses; l++)
      logPi11(k, l) = log(pi11(k, l) / pi00(k, l));
  
  for (int i = 0; i < numOfVertices; i++) {
    for (int j = 0; j < numOfVertices; j++) {
      if (network(i,j) == 1 && network(j,i) == 1) {
        for (int k = 0; k < numOfClasses; k++)
          for (int l = 0; l < numOfClasses; l++)
            A(i, k) += tau(j, l) * logPi11(k, l);
      }
    }
  }
  
  // Finalize by subtracting half of from 1 dividing tau_{ik}
  for (int i = 0; i < numOfVertices; i++)
    for (int k = 0; k < numOfClasses; k++) {
      // In theory, A(i, k) must be negative or 0.
      if (A(i, k) > 0) // In reality, A(i, k) can be greater than 0 because of numerical precision.
        A(i, k) = 0; // Therefore, we cut it off to 0 in this case
      A(i, k) = 1 - A(i, k) / 2;
      A(i, k) /= tau(i, k);
    }
    
    // Calculate the linear coefficients
    NumericVector logAlpha(numOfClasses);
  for (int k = 0; k < numOfClasses; k++)
    logAlpha(k) = log(alpha(k));
  for (int i = 0; i < numOfVertices; i++)
    for (int k = 0; k < numOfClasses; k++)
      s(i, k) = 1 + logAlpha(k) - log(tau(i, k));
  
  NumericMatrix new_tau = NumericMatrix(numOfVertices, numOfClasses);
  solveQP(A, s, new_tau, minTau);
  
  normalizeTau(new_tau, minTau);
  
  return new_tau;
}

 
void updatePi(NumericMatrix& pi,NumericMatrix& stat, const NumericMatrix& tau,
              const NumericMatrix& sumTaus) {
  
  int numOfClasses = pi.nrow();
  int numOfVertices = stat.nrow();
  // reset to 0
  for (int k = 0; k < numOfClasses; k++)
    for (int l = 0; l < numOfClasses; l++)
      pi(k, l) = 0;
  
  // compute numerator sum
  for (int i = 0; i < numOfVertices; i++) {
    for (int j = 0; j < numOfVertices; j++) {
      if (stat(i,j) != 0) {
        for (int k = 0; k < numOfClasses; k++)
          for (int l = 0; l < numOfClasses; l++)
            pi(k, l) += tau(i, k) * tau(j, l);
      }
    }
  }
  
  // divided by denominator sum
  for (int k = 0; k < numOfClasses; k++)
    for (int l = 0; l < numOfClasses; l++)
      pi(k, l) /= sumTaus(k, l);
}


// [[Rcpp::export]]   
List runModelEstimationMStep(int numOfVertices, int numOfClasses,
                             NumericVector alpha, NumericMatrix pi, NumericMatrix stat00,
                             NumericMatrix stat01, NumericMatrix stat10, NumericMatrix stat11,
                             NumericMatrix tau) {
  
  NumericVector tauL = NumericVector(numOfClasses);
  sumDoubleMatrixByRow(tau, tauL);
  
  sumDoubleMatrixByRow(tau, alpha);
  
    for (int k = 0; k < alpha.size(); k++)
    alpha(k) /= numOfVertices;
  normalizeVector(alpha, 1e-6);
  
  NumericMatrix sumTaus(numOfClasses, numOfClasses);
  for (int k = 0; k < numOfClasses; k++)
    for (int l = 0; l < numOfClasses; l++) {
      sumTaus(k, l) = 0;
      for (int i = 0; i < numOfVertices; i++)
        sumTaus(k, l) += tau(i, k) * (tauL(l) - tau(i, l));
    }
    
  //Calculate pi's for reciprocity model
  NumericMatrix pi10(numOfClasses, numOfClasses);
  NumericMatrix pi11(numOfClasses, numOfClasses);
  NumericMatrix pi00(numOfClasses, numOfClasses);
  for (int k = 0; k < numOfClasses; k++) {
    for (int l = 0; l < numOfClasses; l++) {
      pi11(k,l) = pi(k,l)*pi(l,k);
      pi10(k,l) = pi(k,l)*(1 - pi(l,k));
      pi00(k,l) = 1 - pi11(k,l) - pi10(k,l) - pi10(l,k);
    }
  }
  updatePi(pi10, stat10, tau, sumTaus);
  
  updatePi(pi11, stat11, tau, sumTaus);
  
  for (int k = 0; k < numOfClasses; k++)
    for (int l = 0; l <= k; l++) {
      pi00(k, l) = 1 - pi10(k, l) - pi10(l, k) - pi11(k, l);
      pi00(l, k) = pi00(k, l);
    }
    
    // check min and normalize again
    for (int k = 0; k < numOfClasses; k++)
      for (int l = 0; l < numOfClasses; l++) {
        if (pi10(k, l) < minPi)
          pi10(k, l) = minPi;
        if (pi11(k, l) < minPi)
          pi11(k, l) = minPi;
        if (pi00(k, l) < minPi)
          pi00(k, l) = minPi;
      }
      for (int k = 0; k < numOfClasses; k++)
        for (int l = 0; l <= k; l++) {
          double norms = pi00(k, l) + pi10(k, l) + pi10(l, k) + pi11(k, l);
          pi10(k, l) = pi10(k, l) / norms;
          pi10(l, k) = pi10(l, k) / norms;
          pi11(k, l) = pi11(k, l) / norms;
          pi11(l, k) = pi11(k, l);
          pi00(k, l) = pi00(k, l) / norms;
          pi00(l, k) = pi00(k, l);
        }
        
        NumericMatrix pi_new(numOfClasses, numOfClasses);

        for (int k = 0; k < numOfClasses; k++)
          for (int l = 0; l < numOfClasses; l++) 
            pi_new(k, l) = pi11(k,l) + pi10(k,l);
  return List::create(pi_new, alpha);
}

NumericMatrix find_sumTaus(int numOfVertices, int numOfClasses,
                           NumericVector alpha,
                           NumericMatrix tau) {
NumericVector tauL = NumericVector(numOfClasses);
sumDoubleMatrixByRow(tau, tauL);

for (int k = 0; k < alpha.size(); k++)
  alpha(k) /= numOfVertices;
normalizeVector(alpha, 1e-6);

NumericMatrix sumTaus(numOfClasses, numOfClasses);
for (int k = 0; k < numOfClasses; k++)
  for (int l = 0; l < numOfClasses; l++) {
    sumTaus(k, l) = 0;
    for (int i = 0; i < numOfVertices; i++)
      sumTaus(k, l) += tau(i, k) * (tauL(l) - tau(i, l));
  }
  return sumTaus;
}

// [[Rcpp::export]] 
NumericMatrix easy_M_Step(int numOfVertices, int numOfClasses,
                             NumericVector alpha, NumericMatrix pi, NumericMatrix stat,
                             NumericMatrix tau) {
  NumericMatrix sumTaus(numOfClasses, numOfClasses);
  
  sumTaus = find_sumTaus(numOfVertices,numOfClasses,alpha, tau);
  
  updatePi(pi, stat, tau, sumTaus);
  
  for (int k = 0; k < numOfClasses; k++)
    for (int l = 0; l < numOfClasses; l++) {
      if (pi(k, l) < minPi)
        pi(k, l) = minPi;
    }
  
  NumericMatrix pi_new(numOfClasses, numOfClasses);
  pi_new = pi;
  
  return pi_new;
}


  
