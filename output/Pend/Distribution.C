#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>

using namespace std;

void Distribution(string filename)
{
	const double Pi = 3.141592653589793;
	stringstream linestream;
	string line;
	vector<double> Theta,Phi;
	double x,y,z,theta,phi,r,rxy;
	double theta_stat[2], phi_stat[2];
	int n=0;
	ifstream file(filename.c_str());
	file.precision(15);
	cout.precision(15);
	if(!file.good())
	{
		cerr << "Error in opening file" << endl;
		return;
	}
	while(!file.eof())
	{
		getline(file,line);
		linestream.str(line);
		linestream >> x >> y >> z;
		rxy = sqrt(x*x+y*y);
		r = sqrt(x*x+y*y+z*z);
		linestream.clear();
		theta = acos(z/r);
		phi = acos(x/rxy);
		
		/*
		if(y > 0.0)		// Need to do something like this for flip-angles > Pi
		{
			theta = 2*Pi-theta;
		}
		*/
		/*
		if(theta > Pi)
		{
			phi = phi-Pi;
		}
		*/

		Phi.push_back(phi);
		Theta.push_back(theta);
		if(n == 0)
		{
			for(int i=0; i<2; i++)
			{
				theta_stat[i] = theta;
				phi_stat[i] = phi;
			}
		}
		else
		{
			if(theta < theta_stat[0])
			{
				theta_stat[0] = theta;
			}
			else if(theta > theta_stat[1])
			{
				theta_stat[1] = theta;
			}
			if(phi < phi_stat[0])
			{
				phi_stat[0] = phi;
			}
			else if(phi > phi_stat[1])
			{
				phi_stat[1] = phi;
			}
		}
		n++;
	}
	Theta.pop_back();
	Phi.pop_back();
	cout << "n = " << n << endl;
	if(theta_stat[0] == theta_stat[1])
	{
		theta_stat[1] += 0.01*theta_stat[0];
		theta_stat[0] -= 0.01*theta_stat[0];
	}
	if(phi_stat[0] == phi_stat[1])
	{
		phi_stat[1] += 0.01*phi_stat[0];
		phi_stat[0] -= 0.01*phi_stat[0];
	}
	TH1D *hist_theta = new TH1D("hist_theta","distribution theta",500,theta_stat[0]/Pi,theta_stat[1]/Pi);
	TH1D *hist_phi = new TH1D("hist_phi","distribution phi",500,phi_stat[0]/Pi,phi_stat[1]/Pi);
	for(int i=0; i<Theta.size(); i++)
	{
		hist_theta->Fill(Theta[i]/Pi);
		hist_phi->Fill(Phi[i]/Pi);
	}
	hist_theta->GetXaxis()->SetTitle("* #pi");
	hist_phi->GetXaxis()->SetTitle("* #pi");
	TCanvas *canvas_theta = new TCanvas("canvas_theta");
	TCanvas *canvas_phi = new TCanvas("canvas_phi");
	canvas_theta->cd();
	hist_theta->Draw();
	cout << "theta min: " << theta_stat[0] << " theta max: " << theta_stat[1] << endl;
	canvas_phi->cd();
	hist_phi->Draw();
	canvas_theta->Update();
	canvas_phi->Update();
}
