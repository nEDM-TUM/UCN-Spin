#include "dopr.h"
#include <limits>
#include <math.h>
#include <iostream>

using namespace std;

Dopr::Dopr(double firsthtry, int dimension, double* initialY, double* initialDYDX, double initialT, Derivatives *DERI, double atoll, double rtoll, bool dens)
: told (0.0), hdid(0.0), htry(firsthtry), hnext(firsthtry), n(dimension),
 y(initialY), yout(NULL), dydx(initialDYDX), yerr(NULL), yerr2(NULL), t(initialT),
 k2(NULL), k3(NULL), k4(NULL), k5(NULL), k6(NULL), k7(NULL), k8(NULL), k9(NULL), k10(NULL),
 rcont1(NULL), rcont2(NULL), rcont3(NULL), rcont4(NULL), rcont5(NULL), rcont6(NULL), rcont7(NULL), rcont8(NULL),
 derivatives(DERI), atol(atoll), rtol(rtoll), EPS(0.0) , errold(1.0e-4), reject(false), dense(dens), tracker(NULL),
 stepsnottaken(0), hmax(0.003),
c2(0.526001519587677318785587544488e-01),
c3(0.789002279381515978178381316732e-01),
c4(0.118350341907227396726757197510e+00),
c5(0.281649658092772603273242802490e+00),
c6(0.333333333333333333333333333333e+00),
c7(0.25e+00),
c8(0.307692307692307692307692307692e+00),
c9(0.651282051282051282051282051282e+00),
c10(0.6e+00),
c11(0.857142857142857142857142857142e+00),
c14(0.1e+00),
c15(0.2e+00),
c16(0.777777777777777777777777777778e+00),

b1(5.42937341165687622380535766363e-2),
b6(4.45031289275240888144113950566e0),
b7(1.89151789931450038304281599044e0),
b8(-5.8012039600105847814672114227e0),
b9(3.1116436695781989440891606237e-1),
b10(-1.52160949662516078556178806805e-1),
b11(2.01365400804030348374776537501e-1),
b12(4.47106157277725905176885569043e-2),

bhh1(0.244094488188976377952755905512e+00),
bhh2(0.733846688281611857341361741547e+00),
bhh3(0.220588235294117647058823529412e-01),

er1(0.1312004499419488073250102996e-01),
er6(-0.1225156446376204440720569753e+01),
er7(-0.4957589496572501915214079952e+00),
er8( 0.1664377182454986536961530415e+01),
er9(-0.3503288487499736816886487290e+00),
er10(0.3341791187130174790297318841e+00),
er11(0.8192320648511571246570742613e-01),
er12(-0.2235530786388629525884427845e-01),

a21(5.26001519587677318785587544488e-2),
a31(1.97250569845378994544595329183e-2),
a32(5.91751709536136983633785987549e-2),
a41(2.95875854768068491816892993775e-2),
a43(8.87627564304205475450678981324e-2),
a51(2.41365134159266685502369798665e-1),
a53(-8.84549479328286085344864962717e-1),
a54(9.24834003261792003115737966543e-1),
a61(3.7037037037037037037037037037e-2),
a64(1.70828608729473871279604482173e-1),
a65(1.25467687566822425016691814123e-1),
a71(3.7109375e-2),
a74(1.70252211019544039314978060272e-1),
a75(6.02165389804559606850219397283e-2),
a76(-1.7578125e-2),

a81(3.70920001185047927108779319836e-2),
a84(1.70383925712239993810214054705e-1),
a85(1.07262030446373284651809199168e-1),
a86(-1.53194377486244017527936158236e-2),
a87(8.27378916381402288758473766002e-3),
a91(6.24110958716075717114429577812e-1),
a94(-3.36089262944694129406857109825e0),
a95(-8.68219346841726006818189891453e-1),
a96(2.75920996994467083049415600797e1),
a97(2.01540675504778934086186788979e1),
a98(-4.34898841810699588477366255144e1),
a101(4.77662536438264365890433908527e-1),
a104(-2.48811461997166764192642586468e0),
a105(-5.90290826836842996371446475743e-1),
a106(2.12300514481811942347288949897e1),
a107(1.52792336328824235832596922938e1),
a108(-3.32882109689848629194453265587e1),
a109(-2.03312017085086261358222928593e-2),

a111(-9.3714243008598732571704021658e-1),
a114(5.18637242884406370830023853209e0),
a115(1.09143734899672957818500254654e0),
a116(-8.14978701074692612513997267357e0),
a117(-1.85200656599969598641566180701e1),
a118(2.27394870993505042818970056734e1),
a119(2.49360555267965238987089396762e0),
a1110(-3.0467644718982195003823669022e0),
a121(2.27331014751653820792359768449e0),
a124(-1.05344954667372501984066689879e1),
a125(-2.00087205822486249909675718444e0),
a126(-1.79589318631187989172765950534e1),
a127(2.79488845294199600508499808837e1),
a128(-2.85899827713502369474065508674e0),
a129(-8.87285693353062954433549289258e0),
a1210(1.23605671757943030647266201528e1),
a1211(6.43392746015763530355970484046e-1),

a141(5.61675022830479523392909219681e-2),
a147(2.53500210216624811088794765333e-1),
a148(-2.46239037470802489917441475441e-1),
a149(-1.24191423263816360469010140626e-1),
a1410(1.5329179827876569731206322685e-1),
a1411(8.20105229563468988491666602057e-3),
a1412(7.56789766054569976138603589584e-3),
a1413(-8.298e-3),

a151(3.18346481635021405060768473261e-2),
a156(2.83009096723667755288322961402e-2),
a157(5.35419883074385676223797384372e-2),
a158(-5.49237485713909884646569340306e-2),
a1511(-1.08347328697249322858509316994e-4),
a1512(3.82571090835658412954920192323e-4),
a1513(-3.40465008687404560802977114492e-4),
a1514( 1.41312443674632500278074618366e-1),
a161(-4.28896301583791923408573538692e-1),
a166(-4.69762141536116384314449447206e0),
a167(7.68342119606259904184240953878e0),
a168(4.06898981839711007970213554331e0),
a169(3.56727187455281109270669543021e-1),
a1613(-1.39902416515901462129418009734e-3),
a1614(2.9475147891527723389556272149e0),
a1615(-9.15095847217987001081870187138e0),

d41(-0.84289382761090128651353491142e+01),
d46(0.56671495351937776962531783590e+00),
d47(-0.30689499459498916912797304727e+01),
d48(0.23846676565120698287728149680e+01),
d49(0.21170345824450282767155149946e+01),
d410(-0.87139158377797299206789907490e+00),
d411(0.22404374302607882758541771650e+01),
d412(0.63157877876946881815570249290e+00),
d413(-0.88990336451333310820698117400e-01),
d414(0.18148505520854727256656404962e+02),
d415(-0.91946323924783554000451984436e+01),
d416(-0.44360363875948939664310572000e+01),

d51(0.10427508642579134603413151009e+02),
d56(0.24228349177525818288430175319e+03),
d57(0.16520045171727028198505394887e+03),
d58(-0.37454675472269020279518312152e+03),
d59(-0.22113666853125306036270938578e+02),
d510(0.77334326684722638389603898808e+01),
d511(-0.30674084731089398182061213626e+02),
d512(-0.93321305264302278729567221706e+01),
d513(0.15697238121770843886131091075e+02),
d514(-0.31139403219565177677282850411e+02),
d515(-0.93529243588444783865713862664e+01),
d516(0.35816841486394083752465898540e+02),

d61( 0.19985053242002433820987653617e+02),
d66(-0.38703730874935176555105901742e+03),
d67(-0.18917813819516756882830838328e+03),
d68(0.52780815920542364900561016686e+03),
d69(-0.11573902539959630126141871134e+02),
d610(0.68812326946963000169666922661e+01),
d611(-0.10006050966910838403183860980e+01),
d612( 0.77771377980534432092869265740e+00),
d613(-0.27782057523535084065932004339e+01),
d614(-0.60196695231264120758267380846e+02),
d615( 0.84320405506677161018159903784e+02),
d616( 0.11992291136182789328035130030e+02),

d71(-0.25693933462703749003312586129e+02),
d76(-0.15418974869023643374053993627e+03),
d77(-0.23152937917604549567536039109e+03),
d78( 0.35763911791061412378285349910e+03),
d79( 0.93405324183624310003907691704e+02),
d710(-0.37458323136451633156875139351e+02),
d711( 0.10409964950896230045147246184e+03),
d712( 0.29840293426660503123344363579e+02),
d713(-0.43533456590011143754432175058e+02),
d714( 0.96324553959188282948394950600e+02),
d715(-0.39177261675615439165231486172e+02),
d716(-0.14972683625798562581422125276e+03)
{
	EPS=numeric_limits<double>::epsilon();

	yout = new double[n];
	yerr = new double[n];
	yerr2 = new double[n];
	
	k2 = new double[n];
	k3 = new double[n];
	k4 = new double[n];
	k5 = new double[n];
	k6 = new double[n];
	k7 = new double[n];
	k8 = new double[n];
	k9 = new double[n];
	k10 = new double[n];
	
	rcont1 = new double[n];
	rcont2 = new double[n];
	rcont3 = new double[n];
	rcont4 = new double[n];
	rcont5 = new double[n];
	rcont6 = new double[n];
	rcont7 = new double[n];
	rcont8 = new double[n];
	tracker = derivatives->getTracker();
}

void Dopr::reset(double firsthtry, double* initialY, double* initialDYDX, double initialT)
{
	y = initialY;
	dydx = initialDYDX;
	t = initialT;
	htry = firsthtry;
	stepsnottaken = 0;
}

Dopr::~Dopr()
{
	delete[] yout;
	delete[] yerr;
	delete[] yerr2;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	delete[] k5;
	delete[] k6;
	delete[] k7;
	delete[] k8;
	delete[] k9;
	delete[] k10;
	delete[] rcont1;
	delete[] rcont2;
	delete[] rcont3;
	delete[] rcont4;
	delete[] rcont5;
	delete[] rcont6;
	delete[] rcont7;
	delete[] rcont8;
}

void Dopr::step()
{
	double dydxnew[n];
	double h = htry;
	//cout << "htry = " << htry << endl;
	tracker->makeTrack(t,htry);
	//int i=0;
	for(;;)
	{
		//i++;
		dy(h);
		double err = error(h);
		if(success(err,h))
		{
			break;
		}
		else
		{
			tracker->setScaling(h/htry);
			stepsnottaken++;
		}
		if(fabs(h) <= fabs(t)*EPS)
		{
			throw("stepsize underflow in Dopr");
		}
	}
	derivatives->eval(1.0, yout, dydxnew);
	//cout << "\t hdone (after " << i << " steps) = " << h << endl;
	if (dense)
	{
		prepare_dense(h,dydxnew);
	}
	for(int i=0; i<n; i++)
	{
		dydx[i] = dydxnew[i];
		y[i] = yout[i];
	}
	told = t;
	t += (hdid = h);
	htry = hnext;
	tracker->stepDone();
}

void Dopr::dy(const double h) {
	double ytemp[n];
	int i=0;
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
		ytemp[i]=y[i]+h*(a41*dydx[i]+a43*k3[i]);
	}
	derivatives->eval(c4,ytemp,k4);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a51*dydx[i]+a53*k3[i]+a54*k4[i]);
	}
	derivatives->eval(c5,ytemp,k5);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a61*dydx[i]+a64*k4[i]+a65*k5[i]);
	}
	derivatives->eval(c6,ytemp,k6);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a71*dydx[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
	}
	derivatives->eval(c7,ytemp,k7);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a81*dydx[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);
	}
	derivatives->eval(c8,ytemp,k8);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a91*dydx[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i]);
	}
	derivatives->eval(c9,ytemp,k9);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a101*dydx[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]+a108*k8[i]+a109*k9[i]);
	}
	derivatives->eval(c10,ytemp,k10);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a111*dydx[i]+a114*k4[i]+a115*k5[i]+a116*k6[i]+a117*k7[i]+a118*k8[i]+a119*k9[i]+a1110*k10[i]);
	}
	derivatives->eval(c11,ytemp,k2);
	
	for (i=0;i<n;i++)
	{
		ytemp[i]=y[i]+h*(a121*dydx[i]+a124*k4[i]+a125*k5[i]+a126*k6[i]+a127*k7[i]+a128*k8[i]+a129*k9[i]+a1210*k10[i]+a1211*k2[i]);
	}
	derivatives->eval(1.,ytemp,k3);
	
	for (i=0;i<n;i++)
	{
		k4[i]=b1*dydx[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k2[i]+b12*k3[i];
		yout[i]=y[i]+h*k4[i];
	}
	
	for (i=0;i<n;i++)
	{
		yerr[i]=k4[i]-bhh1*dydx[i]-bhh2*k9[i]-bhh3*k3[i];
		yerr2[i]=er1*dydx[i]+er6*k6[i]+er7*k7[i]+er8*k8[i]+er9*k9[i]+er10*k10[i]+er11*k2[i]+er12*k3[i];
	}
}

double Dopr::error(const double h)
{
	double err=0.0, err2=0.0, sk=1.0, deno=1.0, fact=0.0;
	for(int i=0; i<n; i++)
	{
		sk = atol+rtol*max(fabs(y[i]),fabs(yout[i]));
		fact = yerr[i]/sk;
		err2 += fact*fact;
		fact = yerr2[i]/sk;
		err += fact*fact;
	}
	deno = err+0.01*err2;
	if(deno <= 0.0)
	{
		deno = 1.0;
	}
	return fabs(h)*err*sqrt(1.0/(n*deno));
}

bool Dopr::success(const double err, double &h)
{
	const double beta=0.0;
	const double alpha=1.0/8.0-beta*0.2, safe=0.9, minscale=0.333, maxscale=6.0;
	double scale=1.0;
	if (err <= 1.0)
	{
		if(hnext < hmax)
		{
			if (err == 0.0)
			{
				scale=maxscale;
			}
			else
			{
				scale=safe*pow(err,-alpha)*pow(errold,beta);
				if(scale<minscale)
				{
					scale=minscale;
				}
				if (scale>maxscale)
				{
					scale=maxscale;
				}
			}
			if(reject)
			{
				hnext=h*min(scale,1.0);
			}
			else
			{
				hnext=h*scale;
			}
			hnext=min(hnext,hmax);
			
		}
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

void Dopr::prepare_dense(const double h, double *dydxnew)
{
	int i=0;
	double ydiff=0.0,bspl=0.0;
	double ytemp[n];
	for (i=0;i<n;i++)
	{
		rcont1[i] = y[i];
		ydiff = yout[i]-y[i];
		rcont2[i] = ydiff;
		bspl = h*dydx[i]-ydiff;
		rcont3[i] = bspl;
		rcont4[i] = ydiff-h*dydxnew[i]-bspl;
		rcont5[i] = d41*dydx[i]+d46*k6[i]+d47*k7[i]+d48*k8[i]+d49*k9[i]+d410*k10[i]+d411*k2[i]+d412*k3[i];
		rcont6[i] = d51*dydx[i]+d56*k6[i]+d57*k7[i]+d58*k8[i]+d59*k9[i]+d510*k10[i]+d511*k2[i]+d512*k3[i];
		rcont7[i] = d61*dydx[i]+d66*k6[i]+d67*k7[i]+d68*k8[i]+d69*k9[i]+d610*k10[i]+d611*k2[i]+d612*k3[i];
		rcont8[i] = d71*dydx[i]+d76*k6[i]+d77*k7[i]+d78*k8[i]+d79*k9[i]+d710*k10[i]+d711*k2[i]+d712*k3[i];
	}
    for (i=0;i<n;i++)
    {
		ytemp[i]=y[i]+h*(a141*dydx[i]+a147*k7[i]+a148*k8[i]+a149*k9[i]+a1410*k10[i]+a1411*k2[i]+a1412*k3[i]+a1413*dydxnew[i]);
	}
    derivatives->eval(c14,ytemp,k10);
    for (i=0;i<n;i++)
    {
		ytemp[i]=y[i]+h*(a151*dydx[i]+a156*k6[i]+a157*k7[i]+a158*k8[i]+a1511*k2[i]+a1512*k3[i]+a1513*dydxnew[i]+a1514*k10[i]);
	}
    derivatives->eval(c15,ytemp,k2);
    for (i=0;i<n;i++)
    {
		ytemp[i]=y[i]+h*(a161*dydx[i]+a166*k6[i]+a167*k7[i]+a168*k8[i]+a169*k9[i]+a1613*dydxnew[i]+a1614*k10[i]+a1615*k2[i]);
	}
    derivatives->eval(c16,ytemp,k3);
	for (i=0;i<n;i++)
	{
		rcont5[i]=h*(rcont5[i]+d413*dydxnew[i]+d414*k10[i]+d415*k2[i]+d416*k3[i]);
		rcont6[i]=h*(rcont6[i]+d513*dydxnew[i]+d514*k10[i]+d515*k2[i]+d516*k3[i]);
		rcont7[i]=h*(rcont7[i]+d613*dydxnew[i]+d614*k10[i]+d615*k2[i]+d616*k3[i]);
		rcont8[i]=h*(rcont8[i]+d713*dydxnew[i]+d714*k10[i]+d715*k2[i]+d716*k3[i]);
	}
}

double Dopr::dense_out(const int i,const double time,const double h)
{
	double s=(time-told)/h;
	double s1=1.0-s;
	return rcont1[i]+s*(rcont2[i]+s1*(rcont3[i]+s*(rcont4[i]+s1*(rcont5[i]+s*(rcont6[i]+s1*(rcont7[i]+s*rcont8[i]))))));
}
