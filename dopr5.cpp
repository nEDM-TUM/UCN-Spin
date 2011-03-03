#include "dopr5.h"
#include <math.h>
#include <algorithm>
#include <limits>

using namespace std;

Dopr5::Dopr5(double firsthtry, unsigned int dimension, double* initialY, double* initialDYDX, double initialT, Derivatives *DERI, double atoll, double rtoll, bool dens)
: htry(firsthtry), n(dimension), y(initialY), dydx(initialDYDX), t(initialT), derivatives(DERI), atol(atoll), rtol(rtoll), errold(1.0e-4), reject(false), dense(dens)
{
	EPS=numeric_limits<double>::epsilon();

	yout = new double[n];
	yerr = new double[n];
	
	k2 = new double[n];
	k3 = new double[n];
	k4 = new double[n];
	k5 = new double[n];
	k6 = new double[n];
	
	rcont1 = new double[n];
	rcont2 = new double[n];
	rcont3 = new double[n];
	rcont4 = new double[n];
	rcont5 = new double[n];
	
	dydxnew = new double[n];
	
	tracker = derivatives->getTracker();
}

void Dopr5::reset(double firsthtry, double* initialY, double* initialDYDX, double initialT)
{
	y = initialY;
	dydx = initialDYDX;
	t = initialT;
	htry = firsthtry;
}

Dopr5::~Dopr5()
{
	delete yout;
	delete yerr;
	delete k2;
	delete k3;
	delete k4;
	delete k5;
	delete k6;
	
	delete rcont1;
	delete rcont2;
	delete rcont3;
	delete rcont4;
	delete rcont5;
	
	delete dydxnew;
}

void Dopr5::step()
{
	double h=htry;
	tracker->makeTrack(t,htry);
	for (;;)
	{
		dy(h);
		double err=error();
		if (success(err,h))
		{
			break;
		}
		else
		{
			tracker->setScaling(h/htry);
		}
		if (fabs(h) <= fabs(t)*EPS)
		{
			throw("stepsize underflow in Dopr5");
		}
	}
	if (dense)
	{
		prepare_dense(h);
	}
	for(unsigned int i=0; i<n; i++)
	{
		dydx[i] = dydxnew[i];
		y[i] = yout[i];
	}
	told=t;
	t += (hdid=h);
	htry = hnext;
	tracker->stepDone();
}

void Dopr5::dy(const double h)
{
	static const double c2=0.2,c3=0.3,c4=0.8,c5=8.0/9.0,a21=0.2,a31=3.0/40.0,
	a32=9.0/40.0,a41=44.0/45.0,a42=-56.0/15.0,a43=32.0/9.0,a51=19372.0/6561.0,
	a52=-25360.0/2187.0,a53=64448.0/6561.0,a54=-212.0/729.0,a61=9017.0/3168.0,
	a62=-355.0/33.0,a63=46732.0/5247.0,a64=49.0/176.0,a65=-5103.0/18656.0,
	a71=35.0/384.0,a73=500.0/1113.0,a74=125.0/192.0,a75=-2187.0/6784.0,
	a76=11.0/84.0,e1=71.0/57600.0,e3=-71.0/16695.0,e4=71.0/1920.0,
	e5=-17253.0/339200.0,e6=22.0/525.0,e7=-1.0/40.0;
	double ytemp[n];
	unsigned int i;
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*a21*dydx[i];
	}
	derivatives->eval(c2,ytemp,k2);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a31*dydx[i]+a32*k2[i]);
	}
	derivatives->eval(c3,ytemp,k3);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a41*dydx[i]+a42*k2[i]+a43*k3[i]);
	}
	derivatives->eval(c4,ytemp,k4);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a51*dydx[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]);
	}
	derivatives->eval(c5,ytemp,k5);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a61*dydx[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]);
	}
	derivatives->eval(1.,ytemp,k6);
	
	for (i=0;i<n;i++)
	{
		yout[i]=y[i]+h*(a71*dydx[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
	}
	derivatives->eval(1.,yout,dydxnew);
	
	for (i=0;i<n;i++)
	{
		yerr[i]=h*(e1*dydx[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*dydxnew[i]);
	}
}

double Dopr5::error()
{
	double err=0.0,sk,fact;
	for (unsigned int i=0;i<n;i++)
	{
		sk=atol+rtol*max(fabs(y[i]),fabs(yout[i]));
		fact = yerr[i]/sk;
		err += fact*fact;
	}
	return sqrt(err/n);
}

bool Dopr5::success(const double err,double &h)
{
	static const double beta=0.0,alpha=0.2-beta*0.75,safe=0.9,minscale=0.2,maxscale=10.0;
	double scale;
	if (err <= 1.0)
	{
		if (err == 0.0)
		{
			scale=maxscale;
		}
		else
		{
			scale=safe*pow(err,-alpha)*pow(errold,beta);
			if (scale<minscale) scale=minscale;
			if (scale>maxscale) scale=maxscale;
		}
		if (reject)
			hnext=h*min(scale,1.0);
		else
			hnext=h*scale;
		errold=max(err,1.0e-4);
		reject=false;
		return true;
	}
	else
	{
		scale=max(safe*pow(err,-alpha),minscale);
		h *= scale;
		reject=true;
		return false;
	}
}

void Dopr5::prepare_dense(const double h)
{
	static const double d1=-12715105075.0/11282082432.0,d3=87487479700.0/32700410799.0, d4=-10690763975.0/1880347072.0,d5=701980252875.0/199316789632.0, d6=-1453857185.0/822651844.0,d7=69997945.0/29380423.0;
	for (unsigned int i=0;i<n;i++)
	{
		rcont1[i]=y[i];
		double ydiff=yout[i]-y[i];
		rcont2[i]=ydiff;
		double bspl=h*dydx[i]-ydiff;
		rcont3[i]=bspl;
		rcont4[i]=ydiff-h*dydxnew[i]-bspl;
		rcont5[i]=h*(d1*dydx[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*dydxnew[i]);
	}
}

double Dopr5::dense_out(const int i,const double t,const double h)
{
	double s=(t-told)/h;
	double s1=1.0-s;
	return rcont1[i]+s*(rcont2[i]+s1*(rcont3[i]+s*(rcont4[i]+s1*rcont5[i])));
}
