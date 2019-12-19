#include <math.h>
#include <string.h>
#include "mex.h"

/*
  % Deterministic binormalization for symmetric A.
  % [x inform s it] = dsbin(B,x,tol,mns)
  %   B = A.^2. B is equilibrated in the 1-norm and A in the 2-norm.
  %   x is the initial guess.
  %   tol is a small tolerance factor (default 1e-3).
  %   mns is the maximum number of sweeps allowed (default 100*size).
  %   inform: 0 - Success. 1 - Failure.
  %   s is a convergence measure. If s <= 1, then the method converged.
  %   it is the number of iterations (roughly equivalent to matrix-vector
  %     products).
  % On success, diag(x) B diag(x) is equilibrated in the 1-norm.
  
  % Jul 2010. Algorithm by Oren Livne and Gene Golub. See
  %   O.E. Livne and G.H. Golub, "Scaling by Binormalization", Numerical
  %   Algorithms, 35(1), 97-120, 2004.
  % Code by Andrew M. Bradley (ambrad@cs.stanford.edu), though borrowing
  % heavily from the  implementation in Livne and Golub's associated technical
  % report.
*/
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  /* In */
  const mxArray *m_B, *m_x0;
  double tol;
  int mns;
  /* Out */
  double *x;
  int inform;
  double s;
  int it;
  /* Internal */
  const double *pr, *prd;
  mwIndex *ir, *jc, *ird, *jcd;
  mxArray *m_Bdiag, *m_b, *m_plhs[1], *m_prhs[2];
  double *Bdiag, *b, bbar, tmp, c0, c1, c2, disc, xi, di;
  int i, j, k, n, jd;

  m_B = prhs[0];
  n = mxGetM(m_B);
  if(mxGetN(m_B) != n)
    mexErrMsgTxt("B is not square and so certainly not symmetric.");
  pr = mxGetPr(m_B);
  ir = mxGetIr(m_B);
  jc = mxGetJc(m_B);
  m_x0 = prhs[1];
  tol = (nrhs < 3) ? 1e-3 : mxGetScalar(prhs[2]);
  mns = (nrhs < 4) ? 100*n : (int)mxGetScalar(prhs[3]);

  /* x = x0 */
  plhs[0] = mxDuplicateArray(m_x0);
  x = mxGetPr(plhs[0]);
  
  /* b = B*x */
  m_prhs[0] = m_B;
  m_prhs[1] = plhs[0];
  mexCallMATLAB(1, m_plhs, 2, m_prhs, "mtimes");
  m_b = m_plhs[0];
  b = mxGetPr(m_b);

  /* bbar = dot(x,b)/n; */
  bbar = 0;
  for(i = 0; i < n; i++) bbar += x[i]*b[i];
  bbar /= (double)n;

  /* diag(B) */
  m_prhs[0] = m_B;
  mexCallMATLAB(1, m_plhs, 1, m_prhs, "diag");
  m_Bdiag = m_plhs[0];
  ird = mxGetIr(m_Bdiag);
  prd = mxGetPr(m_Bdiag);
  
  it = 0;
  inform = 0;
  for(;;) {
    /* Check for convergence: s = sqrt(sum((x.*b - bbar).^2) / n) */
    s = 0;
    for(i = 0; i < n; i++) {
      tmp = x[i]*b[i] - bbar;
      s += tmp*tmp;
    }
    s = sqrt(s / n);
    if(s <= tol*bbar) break;

    /* Sweep */
    jd = 0;
    for(i = 0; i < n; i++) {
      if(i == ird[jd]) {
	/* Nonzero B(i,i) (== prd[jd]) */
	c2 = (n - 1)*prd[jd];
	c1 = (n - 2)*(b[i] - prd[jd]*x[i]);
	c0 = -prd[jd]*x[i]*x[i] + 2*b[i]*x[i] - n*bbar;
	jd++;
      }
      else {
	/* Zero B(i,i) */
	c2 = 0;
	c1 = (n - 2)*b[i];
	c0 = 2*b[i]*x[i] - n*bbar;
      }
      if(c0 >= 0) {
	inform = 1;
	break;
      }
      disc = c1*c1 - 4*c0*c2;
      xi = -2.0*c0/(c1 + sqrt(disc));
      di = xi - x[i];
      x[i] = xi;

      /* b += B(:,i)*di */
      for(j = jc[i]; j < jc[i+1]; j++)
	b[ir[j]] += pr[j]*di;

      /* bbar = dot(x,b)/n; */
      bbar = 0;
      for(j = 0; j < n; j++) bbar += x[j]*b[j];
      bbar /= (double)n;      
    }
    
    if(inform == 1) break;

    /* Number of sweeps */
    it++;
    if(it >= mns) {
      inform = 1;
      break;
    }
  }

  s = s/(tol*bbar);
  bbar = sqrt(bbar);
  for(i = 0; i < n; i++)
    x[i] = x[i]/bbar;

  if(nlhs > 1) plhs[1] = mxCreateDoubleScalar(inform);
  if(nlhs > 2) plhs[2] = mxCreateDoubleScalar(s);
  if(nlhs > 3) plhs[3] = mxCreateDoubleScalar(it);

  mxDestroyArray(m_Bdiag);
  mxDestroyArray(m_b);
}
