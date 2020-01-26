#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double NAP_cpp(NumericVector ascores, NumericVector bscores) 
{
  int i, j;
  double counter = 0;
  int asize = ascores.size();
  int bsize = bscores.size();
  
  for(i=0; i<asize; i++)
  {
    for(j=0; j<bsize; j++)
    {
      if(ascores[i] < bscores[j])
        counter += 1;

      if(ascores[i] == bscores[j])
        counter += 0.5;
    }
  }
  
  return counter / (asize * bsize);
}
