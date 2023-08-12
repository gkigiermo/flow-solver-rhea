#include "RootFindingMinimization.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////
#define _ACTIVATE_COUT_ 0

////////// BaseRootFindingMinimization CLASS //////////

BaseRootFindingMinimization::BaseRootFindingMinimization(vector<double> &fvec_) : fvec(fvec_) {};
        
BaseRootFindingMinimization::~BaseRootFindingMinimization() {};

void BaseRootFindingMinimization::lnsrch(vector<double> &xold, double &fold, vector<double> &g, vector<double> &p, vector<double> &x,double &f, double &stpmax, bool &check,vector< vector<double> > &fjac) {

    /// Globally Convergent Methods for Nonlinear Systems of Equations: linear search
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
    alam = 1.0;
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

};

void BaseRootFindingMinimization::ludcmp(vector< vector<double> > &a, vector<int> &indx, double &d) {

    /// Globally Convergent Methods for Nonlinear Systems of Equations: ludcmp method
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

void BaseRootFindingMinimization::fdjac(vector<double> &x, vector< vector<double> > &df) {

    /// Globally Convergent Methods for Nonlinear Systems of Equations: fdjac method
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

void BaseRootFindingMinimization::lubksb(vector< vector<double> > &a, vector<int> &indx, vector<double> &b) {

    /// Globally Convergent Methods for Nonlinear Systems of Equations: lubksb method
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

void BaseRootFindingMinimization::qrdcmp(vector< vector<double> > &a, vector<double> &c, vector<double> &d, bool &sing) {

    /// Globally Convergent Methods for Nonlinear Systems of Equations: qrdcmp method
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    int i, j, k, n = c.size();
    double scale, sigma, sum, tau;

    sing = false;
    for( k = 0; k < (n-1); k++ ) {
        scale = 0.0;
	for( i = k; i < n; i++ ) scale = max( scale, fabs( a[i][k] ) );
	if( scale == 0.0 ) {
	    sing = true;
	    c[k] = d[k] = 0.0;
	} else {
	    for( i = k; i < n; i++ ) a[i][k] /= scale;
	    for( sum = 0.0, i = k; i < n; i++ ) sum += SQR(a[i][k]);
	    sigma = SIGN( sqrt( sum ), a[k][k] );
	    a[k][k] += sigma;
	    c[k] = sigma*a[k][k];
	    d[k] = -scale*sigma;
	    for( j = (k+1); j < n; j++ ) {
                sum = 0.0;
	        for( i = k; i < n; i++ ) sum += a[i][k]*a[i][j];
	        tau = sum/c[k];
	        for( i = k; i < n; i++ ) a[i][j] -= tau*a[i][k];
	    }
	}
    }
    d[n-1] = a[n-1][n-1];
    if( d[n-1] == 0.0 ) sing = true;

};

void BaseRootFindingMinimization::rsolv(vector< vector<double> > &a, vector<double> &d, vector<double> &b) {

    /// Globally Convergent Methods for Nonlinear Systems of Equations: rsolv method
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    int i,j, n = d.size();
    double sum;

    b[n-1] /= d[n-1];
    for( i = (n-2); i >= 0; i-- ) {
        sum = 0.0;
        for( j = (i+1); j < n; j++ ) sum += a[i][j]*b[j];
	b[i] = ( b[i] - sum )/d[i];
    }

};

void BaseRootFindingMinimization::qrupdt(vector< vector<double> > &r, vector< vector<double> > &qt, vector<double> &u, vector<double> &v) {

    /// Globally Convergent Methods for Nonlinear Systems of Equations: qrdcmp method
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    int i, k, n = u.size();
    for( k = (n-1); k >= 0; k-- ) {
        if(u[k] != 0.0) break;
    }
    if( k < 0 ) k = 0;
    for( i = (k-1); i >= 0; i-- ) {
        rotate(r,qt,i,u[i],-u[i+1]);
	if( u[i] == 0.0 ) u[i] = fabs(u[i+1]);
	else if( fabs(u[i]) > fabs(u[i+1] ) ) u[i] = fabs( u[i] )*sqrt( 1.0 + SQR( u[i+1]/u[i] ) );
	else u[i] = fabs( u[i+1])*sqrt( 1.0+SQR(u[i]/u[i+1] ) );
    }
    for( i = 0; i < n; i++ ) r[0][i] += u[0]*v[i];
    for( i = 0; i < k; i++ ) rotate(r,qt,i,r[i][i],-r[i+1][i]);

};

void BaseRootFindingMinimization::rotate(vector< vector<double> > &r, vector< vector<double> > &qt, const int &i, const double &a, const double &b) {

    /// Globally Convergent Methods for Nonlinear Systems of Equations: qrdcmp method
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.
	
    int j, n = r.size();
    double c,fact,s,w,y;

    if( a == 0.0 ) {
        c = 0.0;
	s = (b >= 0.0 ? 1.0 : -1.0);
    } else if( fabs(a) > fabs(b) ) {
        fact = b/a;
	c = SIGN( 1.0/sqrt( 1.0 + ( fact*fact ) ), a );
	s = fact*c;
    } else {
        fact = a/b;
	s = SIGN( 1.0/sqrt( 1.0 + ( fact*fact ) ), b );
	c = fact*s;
    }
    for( j = i; j < n; j++ ) {
	y = r[i][j];
	w = r[i+1][j];
	r[i][j] = c*y - s*w;
	r[i+1][j] = s*y + c*w;
    }
    for( j = 0; j < n; j++ ) {
	y = qt[i][j];
	w = qt[i+1][j];
	qt[i][j] = c*y - s*w;
	qt[i+1][j] = s*y + c*w;
    }

};

double BaseRootFindingMinimization::fmin(vector<double> &x) {

    /// Globally Convergent Methods for Nonlinear Systems of Equations: fmin method
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    int i, n = x.size();
    double sum = 0.0;

    function_vector(x,fvec);
    for( i = 0; i < n; i++ ) sum += SQR( fvec[i] );

    return( 0.5*sum );

};

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
      lnsrch(xold,fold,g,p,xmin,f,stpmax,check,fjac);
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


////////// Broyden CLASS //////////

//Broyden::Broyden() : BaseRootFindingMinimization() {};

Broyden::Broyden(vector<double> &fvec_) : BaseRootFindingMinimization(fvec_) {};

Broyden::~Broyden() {};                                                                           

void Broyden::solve(double &fxmin, vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance) {

    /// Broyden method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    //const int MAXITS = 200;
    const int MAXITS = max_iter;
    const double EPS = numeric_limits<double>::epsilon();
    const double TOLF = 1.0e-8, TOLX = EPS, STPMX = 100.0, TOLMIN = 1.0e-12;
    bool check, restrt, sing, skip;
    int i, j, k, n = xmin.size();
    double den, f, fold, stpmax, sum, temp, test;

    vector<double> c(n), d(n), fvcold(n), g(n), p(n), s(n), t(n), w(n), xold(n);
    vector < vector < double > > fjac(n,vector<double> (n));
    vector < vector < double > > qt(n,vector<double> (n));
    
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
    restrt = true;
    for(iter = 0; iter < MAXITS; iter++) {
	if (restrt) {
            fdjac(xmin,fjac);
	    qrdcmp(fjac,c,d,sing);
#if _ACTIVATE_COUT_
	    if (sing) cout << "Singular Jacobian in Broyden" << endl;
#endif
	    for( i = 0; i < n; i++ ) {
	        for( j = 0; j < n; j++ ) qt[i][j] = 0.0;
		qt[i][i] = 1.0;
	    }
	    for( k = 0; k < (n-1); k++ ) {
	        if( c[k] != 0.0) {
		    for( j = 0; j < n; j++ ) {
		        sum = 0.0;
			for( i = k; i < n; i++ ) sum += fjac[i][k]*qt[i][j];
			sum /= c[k];
			for( i = k; i< n; i++ ) qt[i][j] -= sum*fjac[i][k];
		    }
		}
	    }
	    for( i = 0; i < n; i++ ) {
	        fjac[i][i] = d[i];
		for( j = 0; j < i; j++ ) fjac[i][j]=0.0;
	    }
	} else {
	    for( i = 0; i < n; i++ ) s[i] = xmin[i] - xold[i];
	    for( i = 0; i < n; i++) {
                sum = 0.0;
	        for( j = i; j < n; j++ ) sum += fjac[i][j]*s[j];
		t[i] = sum;
	    }
	    skip = true;
	    for( i = 0; i < n; i++ ) {
                sum = 0.0;
	        for( j = 0; j < n; j++) sum += qt[j][i]*t[j];
                w[i] = fvec[i] - fvcold[i] - sum;
		if( fabs( w[i] ) >= EPS*( fabs( fvec[i]) + fabs( fvcold[i] ) ) ) skip = false;
		else w[i] = 0.0;
	    }
	    if( !skip ) {
	        for( i = 0; i < n; i++ ) {
                    sum = 0.0;
		    for( j = 0; j < n; j++ ) sum += qt[i][j]*w[j];
		    t[i] = sum;
		}
                den = 0.0;
		for( i = 0; i < n; i++ ) den += SQR( s[i] );
		for( i = 0; i < n;i++ ) s[i] /= den;
		qrupdt(fjac,qt,t,s);
		for( i = 0; i < n; i++ ) {
#if _ACTIVATE_COUT_
		    if( fjac[i][i] == 0.0 ) cout << "Singular Jacobian in Broyden" << endl;
#endif
		    d[i] = fjac[i][i];
		}
            }
	}
	for( i = 0; i < n; i++ ) {
            sum = 0.0;
	    for( j = 0; j < n; j++ ) sum += qt[i][j]*fvec[j];
	    p[i] = -sum;
	}
	for( i = (n-1); i >= 0; i-- ) {
            sum = 0.0;
	    for( j = 0; j < i; j++ ) sum -= fjac[j][i]*p[j];
	    g[i] = sum;
	}
	for( i = 0; i < n; i++ ) {
	    xold[i] = xmin[i];
	    fvcold[i] = fvec[i];
	}
	fold = f;
	rsolv(fjac,d,p);
	lnsrch(xold,fold,g,p,xmin,f,stpmax,check,fjac);
	test = 0.0;
	for( i = 0; i < n; i++ ) {
	    if( fabs( fvec[i] ) > test ) fxmin = test = fabs(fvec[i]);
        }
	if( test < TOLF ) {
	    check = false;
	    return;
	}
	if( check ) {
	    if( restrt ) {
                return;
	    } else {
	        test = 0.0;
		den = max( f, 0.5*n );
		for( i = 0; i < n; i++ ) {
		    temp = fabs( g[i] )*max( fabs( xmin[i]), 1.0 )/den;
		    if( temp > test ) fxmin = test = temp;
		}
		if( test < TOLMIN ) {
                    return;
                } else {
                    restrt = true;
		}
            }
	} else {
	    restrt = false;
	    test = 0.0;
	    for( i = 0; i < n; i++ ) {
	        temp = ( fabs( xmin[i] - xold[i] ) )/max( fabs( xmin[i] ), 1.0 );
		if( temp > test ) fxmin = test = temp;
	    }
	    if( test < TOLX ) return;
	}
    }
#if _ACTIVATE_COUT_
    cout << "MAXITS exceeded in Broyden" << endl;
#endif

};
