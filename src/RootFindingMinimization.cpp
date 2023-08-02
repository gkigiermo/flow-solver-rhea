#include "RootFindingMinimization.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////
#define _ACTIVATE_COUT_ 0
#define _TAYLORED_LNSRCH_ 1

////////// BaseRootFindingMinimization CLASS //////////

BaseRootFindingMinimization::BaseRootFindingMinimization(vector<double> &fvec_) : fvec(fvec_) {};
        
BaseRootFindingMinimization::~BaseRootFindingMinimization() {};


////////// Brent CLASS //////////

//Brent::Brent() : BaseRootFindingMinimization() {};

Brent::Brent(vector<double> &fvec_) : BaseRootFindingMinimization(fvec_) {};

Brent::~Brent() {};                                                                           

void Brent::set_ax_bx_cx(const double &ax_, const double &bx_, const double &cx_) {

    ax = ax_;
    bx = bx_;
    cx = cx_;

};

void Brent::solve(double &fxmin, vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance) {

    /// Parabolic interpolation and Brent's method in one dimension
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    const double CGOLD = 0.3819660;
    const double ZEPS = numeric_limits<double>::epsilon();

    double a,b,d=0.0,p,q,r,tol1,tol2,u,v,w,xm,etemp,fu,fv,fw,e=0.0;
    vector<double> x( 1, 0.0 );
    vector<double> u_( 1, 0.0 );

    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x[0]=w=v=bx;
    function_vector(x,fvec);
    fw=fv=fxmin=fvec[0];
    for(iter=0; iter<max_iter; ++iter) {
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tolerance*fabs(x[0])+ZEPS);
        if(fabs(x[0]-xm)<=(tol2-0.5*(b-a))) {
            xmin[0]=x[0];
            return ;
        }
        if(fabs(e)>tol1) {
            r=(x[0]-w)*(fxmin-fv);
            q=(x[0]-v)*(fxmin-fw);
            p=(x[0]-v)*q-(x[0]-w)*r;
            q=2.0*(q-r);
            if(q>0.0) {
                p=-p;
            }
            q=fabs(q);
            etemp=e;
            e=d;
            if((fabs(p)>=fabs(0.5*q*etemp))||(p<=q*(a-x[0]))||(p>=q*(b-x[0]))) {
                d=CGOLD*(e=(x[0]>=xm ? a-x[0] : b-x[0]));
            } else {
                d=p/q;
                u=x[0]+d;
                if((u-a < tol2)||(b-u < tol2)) {
                    d=SIGN(tol1,xm-x[0]);
                }
            }
        } else {
            d=CGOLD*(e=(x[0] >= xm ? a-x[0] : b-x[0]));
        }
        u=(fabs(d) >= tol1 ? x[0]+d : x[0]+SIGN(tol1,d));
        u_[0]=u;
        function_vector(u_,fvec);
        fu=fvec[0];
        if(fu<=fxmin) {
            if(u>=x[0]) a=x[0]; else b=x[0];
            shft3(v,w,x[0],u);
            shft3(fv,fw,fxmin,fu);
        } else {
            if(u<x[0]) a=u; else b=u;
            if(fu<=fw || w==x[0]) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if(fu <= fv || v==x[0] || v==w) {
                v=u;
                fv=fu;
            }
        }
    }
#if _ACTIVATE_COUT_
    cout << "Too many iterations in Brent" <<endl;
#endif
    xmin[0]=x[0];

};


////////// NewtonRaphson CLASS //////////

//NewtonRaphson::NewtonRaphson() : BaseRootFindingMinimization() {};

NewtonRaphson::NewtonRaphson(vector<double> &fvec_) : BaseRootFindingMinimization(fvec_) {};

NewtonRaphson::~NewtonRaphson() {};                                                                           

void NewtonRaphson::solve(double &fxmin, vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance) {

    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    //const int MAXITS=200;
    const int MAXITS=max_iter;
    const double TOLF=1.0e-8, TOLX=numeric_limits<double>::epsilon(), STPMX=100.0, TOLMIN=1.0e-12;
    int i, j, n = xmin.size();
    bool check;
    double den,d,f,fold,stpmax,sum,temp,test;

    vector<int> indx(n);
    vector<double> g(n), p(n), xold(n);
    vector < vector < double > > fjac(n,vector<double> (n));

    f = fmin( xmin );
    test = 0.0;
    for(i = 0; i < n; i++) {
      if( fabs( fvec[i] ) > test ) fxmin = test = fabs( fvec[i] );
    }
    if( test < 0.01*TOLF ) {
      check = false;
      return;
    }
    sum = 0.0;
    for(i = 0; i < n; i++) sum += SQR( xmin[i] );
    stpmax = STPMX*max( sqrt( sum ), double( n ) );
    for(iter = 0; iter < MAXITS; iter++) {
      fdjac(xmin,fjac);
      for( i = 0; i < n; i++ ) {
          sum = 0.0;
	  for( j = 0; j < n; j++ ) sum += fjac[j][i]*fvec[j];
	  g[i] = sum;
      }
      for(i = 0; i < n; i++) xold[i] = xmin[i];
      fold = f;
      for(i = 0; i < n; i++) p[i] = -fvec[i];
      ludcmp(fjac,indx,d);
      lubksb(fjac,indx,p);
      lnsrch(xold,fold,g,p,xmin,f,stpmax,check,fjac,indx);
      test = 0.0;
      for(i = 0; i < n; i++) {
        if( fabs( fvec[i] ) > test ) fxmin = test = fabs( fvec[i] );
      }
      if( test < TOLF ) {
        check = false;
#if _ACTIVATE_COUT_
        cout << "Newton-Raphson's minimization TOLF error: " << test << " ( iteration: " << iter << " )" << endl;
#endif
	return;
      }
      if(check) {
          test = 0.0;
	  den = max( f, 0.5*n );
	  for( i = 0; i < n; i++ ) {
	      temp = fabs( g[i])*max( fabs( xmin[i] ), 1.0 )/den;
	      if( temp > test ) test = temp;
	  }
	  check = ( test < TOLMIN );
	  return;
      }
      test = 0.0;
      for(i = 0; i < n; i++) {
	temp = ( fabs( xmin[i] - xold[i] ) )/max( fabs(xmin[i]), 1.0 );
	if( temp > test ) fxmin = test = temp;
      }
      if( test < TOLX ) {
#if _ACTIVATE_COUT_
        cout << "Newton-Raphson's minimization TOLX error: " << test << " ( iteration: " << iter << " )" << endl;
#endif
	return;
      }
    }
#if _ACTIVATE_COUT_
    cout << "MAXITS exceeded in Newton-Raphson" << endl;
#endif

};

void NewtonRaphson::lnsrch(vector<double> &xold, double &fold, vector<double> &g, vector<double> &p, vector<double> &x,double &f, double &stpmax, bool &check,vector< vector<double> > &fjac, vector<int> &indx) {

#if _TAYLORED_LNSRCH_
    int i, n = xold.size();
    double alam = 1.0e-5;

    for(i = 0; i < n; ++i) x[i] = xold[i] + alam*p[i];
#else
    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    const double ALF = 1.0e-4, TOLX = numeric_limits<double>::epsilon();
    double a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
    double rhs1, rhs2, slope = 0.0, sum = 0.0, temp, test, tmplam;
    int i, n = xold.size();

    check = false;
    for( i = 0; i < n; i++ ) sum += p[i]*p[i];
    sum = sqrt( sum );
    if( sum > stpmax ) {
        for( i = 0; i < n; i++ ) p[i] *= stpmax/sum;
    }
    for( i = 0; i < n; i++ ) slope += g[i]*p[i];
#if _ACTIVATE_COUT_
    if(slope >= 0.0) cout << "Roundoff problem in lnsrch" << endl;
#endif
    test = 0.0;
    for( i = 0; i < n; i++ ) {
        temp = fabs( p[i] )/max( fabs( xold[i] ), 1.0 );
        if( temp > test ) test = temp;
    }
    alamin = TOLX/test;
    //alam = 1.0;
    alam = 1.0e-5;
    for(;;) {
        for( i = 0; i < n; i++ ) x[i] = xold[i] + alam*p[i];
        f = fmin( x );
        if( alam < alamin ) {
            for( i = 0; i < n; i++ ) x[i] = xold[i];
            check = true;
            return;
	} else if( f <= fold + ALF*alam*slope ) {
            return;
	} else {
            if( alam == 1.0 ) {
	        tmplam = -slope/( 2.0*( f - fold - slope ) );
	    } else {
	        rhs1 = f - fold - alam*slope;
	        rhs2 = f2 - fold - alam2*slope;
	        a = ( rhs1/( alam*alam ) - rhs2/( alam2*alam2 ) )/( alam - alam2 );
	        b = ( -alam2*rhs1/( alam*alam ) + alam*rhs2/( alam2*alam2 ) )/( alam - alam2 );
	        if( a == 0.0 ) {
                    tmplam = -slope/(2.0*b);
	        } else {
	            disc = b*b - 3.0*a*slope;
		    if( disc < 0.0 ) tmplam = 0.5*alam;
		    else if( b <= 0.0 ) tmplam = ( -b + sqrt( disc ) )/( 3.0*a );
		    else tmplam = -slope/( b + sqrt( disc ) );
	        }
	        if( tmplam > 0.5*alam ) tmplam = 0.5*alam;
	    }
        }
        alam2 = alam;
        f2 = f;
        alam = max( tmplam, 0.1*alam );
    }
#endif

};

void NewtonRaphson::lubksb(vector< vector<double> > &a, vector<int> &indx, vector<double> &b) {

    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    int i,ii=0,ip,j;
    double sum;

    int n=a.size();
    for(i=0;i<n;i++) {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if(ii != 0) {
	for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
      } else if (sum != 0.0) {
	ii=i+1;
      }
      b[i]=sum;
    }
    for(i=n-1;i>=0;i--) {
      sum=b[i];
      for(j=i+1;j<n;j++) sum -= a[i][j]*b[j];
      b[i]=sum/a[i][i];
    }

};

void NewtonRaphson::ludcmp(vector< vector<double> > &a, vector<int> &indx, double &d) {

    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    const double TINY=1.0e-20;
    int i,j,k,imax=0;
    double big,dum,sum,temp;

    int n=a.size();
    vector<double> vv(n);
    d=1.0;
    for(i=0;i<n;i++) {
      big=0.0;
      for(j=0;j<n;j++) {
	if ((temp=fabs(a[i][j])) > big) big=temp;
      }
#if _ACTIVATE_COUT_
      if(big == 0.0) cout << "Singular matrix in routine ludcmp" <<endl;
#endif
      vv[i]=1.0/big;
    }
    for(j=0;j<n;j++) {
      for(i=0;i<j;i++) {
	sum=a[i][j];
	for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
	a[i][j]=sum;
      }
      big=0.0;
      for(i=j;i<n;i++) {
	sum=a[i][j];
	for(k=0;k<j;k++) sum -= a[i][k]*a[k][j];
	a[i][j]=sum;
	if((dum=vv[i]*fabs(sum)) >= big) {
	  big=dum;
	  imax=i;
	}
      }
      if(j != imax) {
	for(k=0;k<n;k++) {
	  dum=a[imax][k];
	  a[imax][k]=a[j][k];
	  a[j][k]=dum;
	}
	d = -d;
	vv[imax]=vv[j];
      }
      indx[j]=imax;
      if(a[j][j] == 0.0) a[j][j]=TINY;
      if(j != n-1) {
	dum=1.0/(a[j][j]);
	for(i=j+1;i<n;i++) a[i][j] *= dum;
      }
    }

};

void NewtonRaphson::fdjac(vector<double> &x, vector< vector<double> > &df) {

    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    const double EPS = 1.0e-8;
    double h, temp;
    int i, j, n = x.size();
    vector<double> f(n);

    for( j = 0; j < n; j++ ) {
      temp = x[j];
      h = EPS*fabs( temp );
      if( h == 0.0 ) h = EPS;
      x[j] = temp + h;
      h = x[j] - temp;
      function_vector(x,f);
      x[j] = temp;
      for( i = 0; i < n; i++ ) {
        df[i][j] = ( f[i] - fvec[i] )/h;
      }
    }   

};

double NewtonRaphson::fmin(vector<double> &x) {

    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    int i, n = x.size();
    double sum = 0.0;

    function_vector(x,fvec);
    for( i = 0; i < n; i++ ) sum += SQR( fvec[i] );

    return( 0.5*sum );

};

