#include <Rcpp.h>
#include <algorithm>
#include <queue>
using namespace Rcpp;

// Functions for a single observed dataset -------------------------------------

std::vector<int> MaxOneCell2x2C(int k_int, int b_int, int c_int) {
  std::vector<int> x(2);
  x[0] = 0;
  x[1] = k_int;
  if (b_int == 0 || k_int == 0) {
    return x;
  } else {
    double k = (double) k_int, b = (double) b_int, c = (double) c_int;
    x[0] = std::max(0.0, ceil( ((b * k) - c) / (b + c) ));
    x[1] = k_int - x[0];
    return x;
  }
}

//' Maximum Likelihood Estimate under the sharp null for Compliers.
//'
//' Find the maximum likelihood estimate of the 2 by 4 contingency table
//'     assuming only Compliers and Never Takers in the population,
//'     under the sharp null for Compliers
//'     and with the multivariate hypergeometric sampling distribution.
//'     This is the C++ function for the R wrapper function
//'     \code{\link{FindMLE_CONT_H0_hypergeoR}}.
//'
//' @param n_y0x0z0 Number of individuals with Y=0, X=0, Z=0.
//' @param n_y1x0z0 Number of individuals with Y=1, X=0, Z=0.
//' @param n_y0x0z1 Number of individuals with Y=0, X=0, Z=1.
//' @param n_y1x0z1 Number of individuals with Y=1, X=0, Z=1.
//' @param n_y0x1z1 Number of individuals with Y=0, X=1, Z=1.
//' @param n_y1x1z1 Number of individuals with Y=1, X=1, Z=1.
//' @return The maximum likelihood under the sharp null for Compliers,
//'     and the corresponding (possibly non-unique)
//'     total number of Compliers of response types
//'     Never Recover and Always Recover.
// [[Rcpp::export]]
List FindMLE_CONT_H0_hypergeoC( int n_y0x0z0, int n_y1x0z0,
                                int n_y0x0z1, int n_y1x0z1,
                                int n_y0x1z1, int n_y1x1z1 ){

  std::vector<int> hatCONRz0 = MaxOneCell2x2C(n_y0x0z0, n_y0x1z1, n_y0x0z1);

  std::vector<int> hatCOARz0 = MaxOneCell2x2C(n_y1x0z0, n_y1x1z1, n_y1x0z1);

  int nz1 = n_y0x0z1 + n_y1x0z1 + n_y0x1z1 + n_y1x1z1;

  double maxpH0 =
    ::Rf_lchoose(n_y0x0z1 + hatCONRz0[1], n_y0x0z1) +
    ::Rf_lchoose(n_y0x1z1 + hatCONRz0[0], n_y0x1z1) +
    ::Rf_lchoose(n_y1x1z1 + hatCOARz0[0], n_y1x1z1) +
    ::Rf_lchoose(n_y1x0z1 + hatCOARz0[1], n_y1x0z1) -
    ::Rf_lchoose(nz1 + n_y0x0z0 + n_y1x0z0, nz1) ;

  List COmle = List::create(hatCONRz0, hatCOARz0) ;

  return List::create( exp(maxpH0), COmle ) ;
}

std::vector<int> KLargestInt_pq(std::vector<int> vin, int k) {

  std::priority_queue< std::pair<int, int> > q;

  for (int i = 0; i < vin.size(); ++i) {
    q.push( std::pair<int, int>(vin[i], i) );
  }

  std::vector<int> out(k);
  for (int i = 0; i < k; ++i) {
    out[i] = q.top().second + 1;
    q.pop();
  }
  return out;
}

std::vector<int> KLargest_pq(std::vector<double> vin, int k) {

  std::priority_queue< std::pair<double, int> > q;

  for (int i = 0; i < vin.size(); ++i) {
    q.push( std::pair<double, int>(vin[i], i) );
  }

  std::vector<int> out(k);
  for (int i = 0; i < k; ++i) {
    out[i] = q.top().second + 1;
    q.pop();
  }
  return out;
}

std::vector<int> OberhoferAlgoC_int(std::vector<int> nis_i, int Nm_i) {

  double Nm = double( Nm_i );
  int k = nis_i.size();
  std::vector<double> nis(nis_i.begin(), nis_i.end());

  std::vector<double> jcur(k);
  std::fill(jcur.begin(), jcur.end(), 1.0);
  std::vector<double> fij_cur(k);
  double nsum = std::accumulate(nis.begin(), nis.end(), 0.0);

  int counter = 0;
  for (int i = 0; i < k; ++i) {
    jcur[i] += std::max(0.0, ceil( (Nm + 1) * nis[i] / nsum ) - 1.0);
    fij_cur[i] = 1.0 + nis[i] / jcur[i] ;
    counter += jcur[i] - 1.0 ;
  }

  std::vector<int> iord = KLargest_pq(fij_cur, 2);
  int icur = iord[0] - 1;
  double inext = fij_cur[iord[1] - 1];

  while (counter < Nm) {
    jcur[icur]++;
    counter++;
    fij_cur[icur] = 1 + (nis[icur] / jcur[icur]);
    if (fij_cur[icur] < inext) {
      iord = KLargest_pq(fij_cur, 2);
      icur = iord[0] - 1;
      inext = fij_cur[iord[1] - 1];
    }
  }

  std::vector<int> result(jcur.begin(), jcur.end());
  for (int i = 0; i < k; ++i) {
    result[i]-- ;
  }

  return result;
}

//' Maximum Likelihood Estimate without assuming the sharp null for Compliers.
//'
//' Find the maximum likelihood estimate of the 2 by 4 contingency table
//'     assuming only Compliers and Never Takers in the population,
//'     with the multivariate hypergeometric sampling distribution.
//'     This is the C++ function for the R wrapper function
//'     \code{\link{FindMLE_CONT_H1_hypergeoR}}.
//'
//' @inheritParams FindMLE_CONT_H0_hypergeoC
//' @return The maximum likelihood and the corresponding
//'     (possibly non-unique) number of Compliers in each Z arm.
// [[Rcpp::export]]
List FindMLE_CONT_H1_hypergeoC( int n_y0x0z0, int n_y1x0z0,
                                int n_y0x0z1, int n_y1x0z1,
                                int n_y0x1z1, int n_y1x1z1 ) {

  double maxpH1 = 0.0;
  List COmle(4);

  for (int COHEz1=0; COHEz1 <= n_y1x1z1; ++COHEz1) {
    for (int COHUz1=0; COHUz1 <= n_y0x1z1; ++COHUz1) {

      std::vector<int> y0z1(3);
      y0z1[0] = n_y0x0z1;
      y0z1[1] = n_y0x1z1 - COHUz1;
      y0z1[2] = COHEz1;

      std::vector<int> y0z0hat = OberhoferAlgoC_int( y0z1, n_y0x0z0 ) ;

      std::vector<int>  y1z1(3);
      y1z1[0] = n_y1x0z1;
      y1z1[1] = n_y1x1z1 - COHEz1;
      y1z1[2] = COHUz1;

      std::vector<int> y1z0hat = OberhoferAlgoC_int( y1z1, n_y1x0z0 ) ;

      double maxpH1_cand = 0.0;
      for (int jj=0; jj < 3; ++jj) {
        maxpH1_cand = maxpH1_cand +
          ::Rf_lchoose(y0z0hat[jj] + y0z1[jj], y0z1[jj]) +
          ::Rf_lchoose(y1z0hat[jj] + y1z1[jj], y1z1[jj]);
      }
      if (maxpH1_cand > maxpH1) {
        maxpH1 = maxpH1_cand;
        COmle[0] = y0z0hat;
        COmle[1] = y1z0hat;
        COmle[2] = y0z1;
        COmle[3] = y1z1;
      }
    }
  }

  int nz1 = n_y0x0z1 + n_y1x0z1 + n_y0x1z1 + n_y1x1z1;

  maxpH1 -= ::Rf_lchoose(nz1 + n_y0x0z0 + n_y1x0z0, nz1) ;

  return List::create( exp(maxpH1), COmle ) ;
}


//' Finite population sample space given an observed dataset.
//'
//' Sample space of all possibly observable datasets given an observed dataset,
//'     assuming only Compliers and Never Takers in the population.
//'     This is the C++ function for the R wrapper function
//'     \code{\link{AllPossiblyObsH0_CONT}}.
//'
//' @param obs_y0x0z0 Number of observed individuals with Y=0, X=0, Z=0.
//' @param obs_y1x0z0 Number of observed individuals with Y=1, X=0, Z=0.
//' @param obs_y0x0z1 Number of observed individuals with Y=0, X=0, Z=1.
//' @param obs_y1x0z1 Number of observed individuals with Y=1, X=0, Z=1.
//' @param obs_y0x1z1 Number of observed individuals with Y=0, X=1, Z=1.
//' @param obs_y1x1z1 Number of observed individuals with Y=1, X=1, Z=1.
//' @return All possibly observable datasets in a list format.
// [[Rcpp::export]]
List AllPossiblyObsH0_CONT_C(
    int obs_y0x0z0, int obs_y1x0z0,
    int obs_y0x0z1, int obs_y1x0z1, int obs_y0x1z1, int obs_y1x1z1) {

  int ny0 = obs_y0x0z0 + obs_y0x0z1 + obs_y0x1z1;
  int ny1 = obs_y1x0z0 + obs_y1x0z1 + obs_y1x1z1;
  int nz0 = obs_y0x0z0 + obs_y1x0z0;
  int nz1 = ny0 + ny1 - nz0;

  // matrix to cache lchoose results
  int nCk_dim = std::max(ny0, ny1) + 1;
  NumericMatrix nCk( nCk_dim, nCk_dim );
  std::fill( nCk.begin(), nCk.end(), -1.0 );

  double normconst = ::Rf_lchoose(nz0 + nz1, nz1);

  std::vector<int> n_y0x0z0_out;
  std::vector<int> n_y1x0z0_out;
  std::vector<int> n_y0x0z1_out;
  std::vector<int> n_y1x0z1_out;
  std::vector<int> n_y0x1z1_out;
  std::vector<int> n_y1x1z1_out;
  std::vector<double> qH0_all;

  for (int ny0z0=std::max(0, ny0 - nz1);
       ny0z0 <= std::min(ny0, nz0);
       ++ny0z0) {

    int ny1z0 = nz0 - ny0z0;
    int ny0z1 = ny0 - ny0z0;
    int ny1z1 = nz1 - ny0z1;

    for (int n_y0x0z1 = std::max(0, obs_y0x0z1 - ny0z0);
         n_y0x0z1 <= std::min(ny0z1, obs_y0x0z0 + obs_y0x0z1);
         ++n_y0x0z1) {
      for (int n_y1x0z1=std::max(0, obs_y1x0z1 - ny1z0);
           n_y1x0z1 <= std::min(ny1z1, obs_y1x0z0 + obs_y1x0z1);
           ++ n_y1x0z1) {

        // (1) save the unique datasets -------------------------------------
        int n000 = ny0z0, n100 = ny1z0,
          n001 = n_y0x0z1, n101 = n_y1x0z1,
          n011 = ny0z1 - n_y0x0z1, n111 = ny1z1 - n_y1x0z1;

        n_y0x0z0_out.push_back(n000);
        n_y1x0z0_out.push_back(n100);
        n_y0x0z1_out.push_back(n001);
        n_y1x0z1_out.push_back(n101);
        n_y0x1z1_out.push_back(n011);
        n_y1x1z1_out.push_back(n111);
        // (1) --------------------------------------------------------------

        // (2) find the max likelihood under H0 -----------------------------
        std::vector<int> hatCONRz0 = MaxOneCell2x2C(n000, n011, n001);
        std::vector<int> hatCOARz0 = MaxOneCell2x2C(n100, n111, n101);

        std::vector<int> coltots(4);
        coltots[0] = n001 + hatCONRz0[1] ;
        coltots[1] = n011 + hatCONRz0[0] ;
        coltots[2] = n111 + hatCOARz0[0] ;
        coltots[3] = n101 + hatCOARz0[1] ;

        std::vector<int> nz1s(4);
        nz1s[0] = n001 ;
        nz1s[1] = n011 ;
        nz1s[2] = n111 ;
        nz1s[3] = n101 ;

        double maxpH0 = -normconst;
        for (int jj=0; jj < 4; ++jj) {
          int nn = coltots[jj], kk = nz1s[jj] ;
          // check cache
          if( nCk(nn, kk) < 0.0 ) {
            // calculate the lchoose value
            nCk(nn, kk) = ::Rf_lchoose( nn, kk );
            nCk(nn, nn - kk) = nCk(nn, kk);
          }
          maxpH0 += nCk(nn, kk);
        }
        qH0_all.push_back( exp( maxpH0 ) ) ;
        // (2) --------------------------------------------------------------
      }
    }
  }

  List result(7);
  result[0] = n_y0x0z0_out;
  result[1] = n_y1x0z0_out;
  result[2] = n_y0x0z1_out;
  result[3] = n_y1x0z1_out;
  result[4] = n_y0x1z1_out;
  result[5] = n_y1x1z1_out;
  result[6] = qH0_all;

  return result;
}


//' Finite population sample space given an observed dataset.
//'
//' Sample space of all possibly observable datasets given an observed dataset,
//'     assuming only Compliers and Never Takers in the population.
//'     The maximum likelihood estimate is also calculated for each dataset.
//'     This is the C++ function for the R wrapper function
//'     \code{\link{AllPossiblyObsH0_CONT}}.
//'
//' @inheritParams AllPossiblyObsH0_CONT_C
//' @return All possibly observable datasets in a list format.
// [[Rcpp::export]]
List AllPossiblyObsH0qH1_CONT_C(
    int obs_y0x0z0, int obs_y1x0z0,
    int obs_y0x0z1, int obs_y1x0z1, int obs_y0x1z1, int obs_y1x1z1) {

  int ny0 = obs_y0x0z0 + obs_y0x0z1 + obs_y0x1z1;
  int ny1 = obs_y1x0z0 + obs_y1x0z1 + obs_y1x1z1;
  int nz0 = obs_y0x0z0 + obs_y1x0z0;
  int nz1 = ny0 + ny1 - nz0;

  // matrix to cache lchoose results
  int nCk_dim = ny0 + ny1 + 1;
  NumericMatrix nCk( nCk_dim, nCk_dim );
  std::fill( nCk.begin(), nCk.end(), -1.0 );

  double normconst = ::Rf_lchoose(nz0 + nz1, nz1);

  // Single map for Oberhofer MLE
  typedef std::map< std::pair< int, std::vector<int> >,
                    std::vector<int> > OberMap;
  typedef OberMap::iterator OberIter;
  OberMap ober;

  // Observed value of the GLR
  List qH0out = FindMLE_CONT_H0_hypergeoC(
    obs_y0x0z0, obs_y1x0z0, obs_y0x0z1, obs_y1x0z1, obs_y0x1z1, obs_y1x1z1 ) ;

  List qH1out = FindMLE_CONT_H1_hypergeoC(
    obs_y0x0z0, obs_y1x0z0, obs_y0x0z1, obs_y1x0z1, obs_y0x1z1, obs_y1x1z1 ) ;

  double obsG = ( (double) qH0out[0] ) / ( (double) qH1out[0] ) *
    exp( std::numeric_limits<double>::epsilon() * 1e2 ) ;

  std::vector<int> n_y0x0z0_out;
  std::vector<int> n_y1x0z0_out;
  std::vector<int> n_y0x0z1_out;
  std::vector<int> n_y1x0z1_out;
  std::vector<int> n_y0x1z1_out;
  std::vector<int> n_y1x1z1_out;
  std::vector<double> qH0_all;
  std::vector<double> qH1_all;

  for (int ny0z0=std::max(0, ny0 - nz1);
       ny0z0 <= std::min(ny0, nz0);
       ++ny0z0) {

    int ny1z0 = nz0 - ny0z0;
    int ny0z1 = ny0 - ny0z0;
    int ny1z1 = nz1 - ny0z1;

    for (int n_y0x0z1 = std::max(0, obs_y0x0z1 - ny0z0);
         n_y0x0z1 <= std::min(ny0z1, obs_y0x0z0 + obs_y0x0z1);
         ++n_y0x0z1) {
      for (int n_y1x0z1=std::max(0, obs_y1x0z1 - ny1z0);
           n_y1x0z1 <= std::min(ny1z1, obs_y1x0z0 + obs_y1x0z1);
           ++ n_y1x0z1) {

        // (1) save the unique datasets -------------------------------------
        int n000 = ny0z0, n100 = ny1z0,
          n001 = n_y0x0z1, n101 = n_y1x0z1,
          n011 = ny0z1 - n_y0x0z1, n111 = ny1z1 - n_y1x0z1;

        n_y0x0z0_out.push_back(n000);
        n_y1x0z0_out.push_back(n100);
        n_y0x0z1_out.push_back(n001);
        n_y1x0z1_out.push_back(n101);
        n_y0x1z1_out.push_back(n011);
        n_y1x1z1_out.push_back(n111);
        // (1) --------------------------------------------------------------

        // (2) find the max likelihood under H0 -----------------------------
        std::vector<int> hatCONRz0 = MaxOneCell2x2C(n000, n011, n001);
        std::vector<int> hatCOARz0 = MaxOneCell2x2C(n100, n111, n101);

        std::vector<int> coltots(4);
        coltots[0] = n001 + hatCONRz0[1] ;
        coltots[1] = n011 + hatCONRz0[0] ;
        coltots[2] = n111 + hatCOARz0[0] ;
        coltots[3] = n101 + hatCOARz0[1] ;

        std::vector<int> nz1s(4);
        nz1s[0] = n001 ;
        nz1s[1] = n011 ;
        nz1s[2] = n111 ;
        nz1s[3] = n101 ;

        double maxpH0 = -normconst;
        for (int jj=0; jj < 4; ++jj) {
          int nn = coltots[jj], kk = nz1s[jj] ;
          // check cache
          if( nCk(nn, kk) < 0.0 ) {
            // calculate the lchoose value
            nCk(nn, kk) = ::Rf_lchoose( nn, kk );
            nCk(nn, nn - kk) = nCk(nn, kk);
          }
          maxpH0 += nCk(nn, kk);
        }
        qH0_all.push_back( exp( maxpH0 ) ) ;
        // (2) --------------------------------------------------------------

        // (3) find the max likelihood under H1 (if needed) -----------------
        double maxpH1 = 0.0; // default value for less extreme datasets
        if( maxpH0 <= log(obsG) ) {

          int a, b, c, d, e, f;

          a = n000, b = n100, c = n011, d = n111, e = n001, f = n101;

          // ---------------------- FindMLEnoDE_CiterMap ----------------------------
          // Part I: condition on z1 arm ---------------------------------------------
          // Initialise hash tables
          typedef std::set< std::vector<int> > MemoMap;
          typedef MemoMap::iterator MemoIter;
          MemoMap memo;

          for (int t1 = 0; t1 <= d; ++t1) {
            for (int u1 = 0; u1 <= c; ++u1) {

              int s1 = c - u1, w1 = d - t1;

              int vlins;
              std::vector<int> vin_o(3), inc_idx(3), vino_sort(3);
              std::vector<int> vhat_o(3), vhat(3);
              OberIter v_it;

              vlins = a;
              vin_o[0] = e, vin_o[1] = s1, vin_o[2] = t1;

              // Oberhofer (start) ----------------------------------------------//
              // search map for precalculated values
              inc_idx = KLargestInt_pq(vin_o, 3); // Sorted order
              for (int idx=0; idx < 3; ++idx) {
                vino_sort[idx] = vin_o[ inc_idx[idx] - 1 ];
              }
              v_it = ober.find( std::make_pair( vlins, vino_sort ) );
              // if it exists pull the result
              if( v_it != ober.end() ) {
                vhat_o = v_it->second ;
              } else {
                // find the mle
                vhat_o = OberhoferAlgoC_int( vino_sort, vlins );
                // save
                ober.insert( std::make_pair(
                    std::make_pair( vlins, vino_sort ), vhat_o) );
              }
              for (int idx=0; idx < 3; ++idx) {
                vhat[ inc_idx[idx] - 1 ] = vhat_o[idx]; // Unsort output
              }
              // Oberhofer (end) ------------------------------------------------//
              int s0 = vhat[1], t0 = vhat[2];

              vlins = b;
              vin_o[0] = f, vin_o[1] = u1, vin_o[2] = w1;

              // Oberhofer (start) ----------------------------------------------//
              // search map for precalculated values
              inc_idx = KLargestInt_pq(vin_o, 3); // Sorted order
              for (int idx=0; idx < 3; ++idx) {
                vino_sort[idx] = vin_o[ inc_idx[idx] - 1 ];
              }
              v_it = ober.find( std::make_pair( vlins, vino_sort ) );
              // if it exists pull the result
              if( v_it != ober.end() ) {
                vhat_o = v_it->second ;
              } else {
                // find the mle
                vhat_o = OberhoferAlgoC_int( vino_sort, vlins );
                // save
                ober.insert( std::make_pair(
                    std::make_pair( vlins, vino_sort ), vhat_o) );
              }
              for (int idx=0; idx < 3; ++idx) {
                vhat[ inc_idx[idx] - 1 ] = vhat_o[idx]; // Unsort output
              }
              // Oberhofer (end) ------------------------------------------------//
              int u0 = vhat[1], w0 = vhat[2];


              // iterate back to z1 arm now ------------------------------------------
              std::vector<int> s1u1 = MaxOneCell2x2C(c, s0, u0) ;
              std::vector<int> t1w1 = MaxOneCell2x2C(d, t0, w0) ;

              std::vector<int> t1u1(2);
              t1u1[0] = t1w1[0], t1u1[1] = s1u1[1];

              memo.insert( t1u1 );
            }
          }

          // Part II: find max hypergeometric for saved combinations ----------------
          for (MemoIter miter=memo.begin(); miter!=memo.end(); ++miter) {

            std::vector<int> t1u1 = *miter;
            int t1 = t1u1[0], u1 = t1u1[1] ;
            int s1 = c - u1, w1 = d - t1;

            int vlins;
            std::vector<int> vin_o(3), inc_idx(3), vino_sort(3);
            std::vector<int> vhat_o(3), vhat(3);
            OberIter v_it;

            vlins = a;
            vin_o[0] = e, vin_o[1] = s1, vin_o[2] = t1;
            // Oberhofer (start) ----------------------------------------------//
            inc_idx = KLargestInt_pq(vin_o, 3); // Sorted order
            for (int idx=0; idx < 3; ++idx) {
              vino_sort[idx] = vin_o[ inc_idx[idx] - 1 ];
            }
            v_it = ober.find( std::make_pair( vlins, vino_sort ) );
            vhat_o = v_it->second ;
            for (int idx=0; idx < 3; ++idx) {
              vhat[ inc_idx[idx] - 1 ] = vhat_o[idx]; // Unsort output
            }
            // Oberhofer (end) ------------------------------------------------//
            int s0 = vhat[1], t0 = vhat[2];

            vlins = b;
            vin_o[0] = f, vin_o[1] = u1, vin_o[2] = w1;
            // Oberhofer (start) ----------------------------------------------//
            inc_idx = KLargestInt_pq(vin_o, 3); // Sorted order
            for (int idx=0; idx < 3; ++idx) {
              vino_sort[idx] = vin_o[ inc_idx[idx] - 1 ];
            }
            v_it = ober.find( std::make_pair( vlins, vino_sort ) );
            vhat_o = v_it->second ;
            for (int idx=0; idx < 3; ++idx) {
              vhat[ inc_idx[idx] - 1 ] = vhat_o[idx]; // Unsort output
            }
            // Oberhofer (end) ------------------------------------------------//
            int u0 = vhat[1], w0 = vhat[2];

            std::vector<int> onetab(6);
            onetab[0]=a-s0-t0;
            onetab[1]=s0;
            onetab[2]=t0;
            onetab[3]=u0;
            onetab[4]=w0;
            onetab[5]=b-u0-w0;

            std::vector<int> onetabtot = onetab;
            onetabtot[0]+=e;
            onetabtot[1]+=s1;
            onetabtot[2]+=t1;
            onetabtot[3]+=u1;
            onetabtot[4]+=w1;
            onetabtot[5]+=f;

            double mle_cand = 0.0;

            for (int idx=0; idx < 6; ++idx) {

              int nn = onetabtot[idx], kk = onetab[idx];

              // check cache
              if( nCk(nn, kk) < 0.0 ) {
                // calculate the lchoose value
                nCk(nn, kk) = ::Rf_lchoose( nn, kk );
                nCk(nn, nn - kk) = nCk(nn, kk);
              }
              mle_cand += nCk(nn, kk);
              if (mle_cand > maxpH1) {
                maxpH1 = mle_cand;
              }
            }
          }
          // ---------------------- FindMLEnoDE_CiterMap ----------------------------
          maxpH1 -= normconst ;
        }
        qH1_all.push_back( exp( maxpH1 ) ) ;
        // (3) --------------------------------------------------------------
      }
    }
  }

  List result(8);
  result[0] = n_y0x0z0_out;
  result[1] = n_y1x0z0_out;
  result[2] = n_y0x0z1_out;
  result[3] = n_y1x0z1_out;
  result[4] = n_y0x1z1_out;
  result[5] = n_y1x1z1_out;
  result[6] = qH0_all;
  result[7] = qH1_all;

  return result;

}

// Finding probabilities under H0 ---------------------------------------------

//' Exact finite population p-values under the sharp null for Compliers.
//'
//' Find the exact population-specific p-values
//'     under the sharp null for Compliers,
//'     for each compatible population with only Compliers and Never Takers.
//'     This is the C++ function for the R wrapper function
//'     \code{\link{Get_pvalues_CONT}}.
//'
//' @param n_y0x0z0_H0 Number of individuals with Y=0, X=0, Z=0.
//' @param n_y1x0z0_H0 Number of individuals with Y=1, X=0, Z=0.
//' @param n_y0x0z1_H0 Number of individuals with Y=0, X=0, Z=1.
//' @param n_y1x0z1_H0 Number of individuals with Y=1, X=0, Z=1.
//' @param n_y0x1z1_H0 Number of individuals with Y=0, X=1, Z=1.
//' @param n_y1x1z1_H0 Number of individuals with Y=1, X=1, Z=1.
//' @param n_NTy0_H0 Total number of Never Takers with Y=0.
//' @param n_CONR_H0 Total number of Compliers of response type
//'     Never Recover.
//' @param n_COAR_H0 Total number of Compliers of response type
//'     Always Recover.
//' @param n_NTy1_H0 Total number of Never Takers with Y=1.
//' @param critical_regions Matrix denoting whether each
//'     possibly observed dataset (rows) is in a critical region
//'     (columns).
//' @return A matrix of p-values correspnding to each population
//'     (row) and the critical region (columns).
// [[Rcpp::export]]
NumericMatrix GetPvalueshypergeoC_allpsi_CONT(
    std::vector<int> n_y0x0z0_H0, std::vector<int> n_y1x0z0_H0,
    std::vector<int> n_y0x0z1_H0, std::vector<int> n_y1x0z1_H0,
    std::vector<int> n_y0x1z1_H0, std::vector<int> n_y1x1z1_H0,
    std::vector<int> n_NTy0_H0,
    std::vector<int> n_CONR_H0, std::vector<int> n_COAR_H0,
    std::vector<int> n_NTy1_H0,
    IntegerMatrix critical_regions) {

  // critical_regions: (n x c) 0/1 matrix for c different critical regions
  IntegerVector critical_check( critical_regions.nrow(), 0);
  for (int ii=0; ii < critical_regions.nrow(); ++ii) {
    // Is this dataset in any of the critical regions?
    for (int cc=0; cc < critical_regions.ncol(); ++cc) {
      critical_check[ii] += critical_regions(ii, cc);
    }
  }

  // ProbH0hypergeoC_allpsi ---------------------------------------------------
  int nz1 = n_y0x0z1_H0[0] + n_y1x0z1_H0[0] + n_y0x1z1_H0[0] + n_y1x1z1_H0[0];

  double normconst = ::Rf_lchoose(
    n_NTy0_H0[0] + n_CONR_H0[0] + n_COAR_H0[0] + n_NTy1_H0[0], nz1 );

  int allns = n_y0x0z1_H0.size();

  int psisize = n_NTy0_H0.size();

  // pvalues_all: (|Psi| x c) matrix for c different critical regions
  NumericMatrix pvalues_all( psisize, critical_regions.ncol() );
  std::fill( pvalues_all.begin(), pvalues_all.end(), 0.0) ;

  std::vector<int> maxpsi(4);
  maxpsi[0] = *max_element( n_NTy0_H0.begin(), n_NTy0_H0.end() );
  maxpsi[1] = *max_element( n_CONR_H0.begin(), n_CONR_H0.end() );
  maxpsi[2] = *max_element( n_COAR_H0.begin(), n_COAR_H0.end() );
  maxpsi[3] = *max_element( n_NTy1_H0.begin(), n_NTy1_H0.end() );
  int nCk_dim = *max_element( maxpsi.begin(), maxpsi.end() ) + 1;

  NumericMatrix nCk( nCk_dim, nCk_dim );
  std::fill( nCk.begin(), nCk.end(), -1.0 );

  for (int psi=0; psi < psisize; ++psi) {

    std::vector<int> psi_fixed(4);
    psi_fixed[0] = n_NTy0_H0[psi];
    psi_fixed[1] = n_NTy1_H0[psi];
    psi_fixed[2] = n_CONR_H0[psi];
    psi_fixed[3] = n_COAR_H0[psi];

    for (int ii=0; ii < allns; ++ii) {

      if ( critical_check[ii] > 0 ) {

        int n_CONR_zz = std::max(n_y0x0z1_H0[ii]+n_y0x0z0_H0[ii]-psi_fixed[0],
                                 n_y0x1z1_H0[ii]);
        int n_COAR_zz = std::max(n_y1x0z1_H0[ii]+n_y1x0z0_H0[ii]-psi_fixed[1],
                                 n_y1x1z1_H0[ii]);

        bool check =
          (n_CONR_zz < 0) || (n_COAR_zz < 0) ||
          (psi_fixed[0] < n_y0x0z1_H0[ii]) ||
          (psi_fixed[1] < n_y1x0z1_H0[ii]) ||
          (psi_fixed[2] < n_CONR_zz) ||
          (psi_fixed[3] < n_COAR_zz) ;

        if (check == false) {

          std::vector<int> ntilde_fixed(4);
          ntilde_fixed[0] = n_y0x0z1_H0[ii];
          ntilde_fixed[1] = n_y1x0z1_H0[ii];
          ntilde_fixed[2] = n_CONR_zz;
          ntilde_fixed[3] = n_COAR_zz;

          double pH0_fixed = 0.0;

          for (int idx=0; idx < 4; ++idx) {

            int nn = psi_fixed[idx], kk = ntilde_fixed[idx];

            // check cache
            if( nCk(nn, kk) < 0.0 ) {
              // calculate the lchoose value
              nCk(nn, kk) = ::Rf_lchoose( nn, kk );
              nCk(nn, nn - kk) = nCk(nn, kk);
            }
            pH0_fixed += nCk(nn, kk);
          }
          for (int cc=0; cc < pvalues_all.ncol(); ++cc) {
            if ( critical_regions( ii, cc ) == 1 ) {
              pvalues_all( psi, cc ) += exp( pH0_fixed - normconst );
            }
          }
        }
      }
    }
  }
  return pvalues_all;
}
