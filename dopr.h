#include "derivatives.h"

class Dopr
{
	public:
		Dopr(double, int, double*, double*, double, Derivatives*, double, double, bool);
		~Dopr();
		void reset(double, double*, double*, double);
		void step();
		void dy(const double);
		double error(const double);
		bool success(const double, double&);
		void prepare_dense(const double, double*);
		double dense_out(const int, const double, const double);
		
		inline double getHdid(){return hdid;};
		inline double getT(){return t;};
		inline int getStepsnottaken() {return stepsnottaken;};
	
	private:	
		double told, hdid, htry, hnext;
		const int n;
		double *y, *yout, *dydx, *yerr, *yerr2;
		double t;
		double *k2,*k3,*k4,*k5,*k6,*k7,*k8,*k9,*k10;
		double *rcont1, *rcont2, *rcont3, *rcont4, *rcont5, *rcont6, *rcont7, *rcont8;
		Derivatives *derivatives;
		double atol, rtol, EPS, errold;
		bool reject;
		bool dense;
		Tracking *tracker;
		int stepsnottaken;
		const double hmax;
		
		const double c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c14,c15,c16;
		const double b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3;
		const double er1,er6,er7,er8,er9,er10,er11,er12;
		const double a21,a31,a32,a41,a43,a51,a53,a54,a61,a64,a65,a71,a74,a75,a76;
		const double a81,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,a101,a104,a105;
		const double a106,a107,a108,a109,a111,a114,a115,a116,a117,a118,a119,a1110;
		const double a121,a124,a125,a126,a127,a128,a129,a1210,a1211,a141,a147,a148;
		const double a149,a1410,a1411,a1412,a1413,a151,a156,a157,a158,a1511,a1512;
		const double a1513,a1514,a161,a166,a167,a168,a169,a1613,a1614,a1615;
		const double d41,d46,d47,d48,d49,d410,d411,d412,d413,d414,d415,d416,d51,d56;
		const double d57,d58,d59,d510,d511,d512,d513,d514,d515,d516,d61,d66,d67,d68;
		const double d69,d610,d611,d612,d613,d614,d615,d616,d71,d76,d77,d78,d79;
		const double d710,d711,d712,d713,d714,d715,d716;
};
