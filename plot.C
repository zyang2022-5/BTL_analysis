#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <cstdlib>

vector<vector<double>> read_file(const char *filename)
{
	ifstream myfile (filename);
	string temp;
	vector<vector<double>> sd_err_l;
	vector<double> sd_l;
	vector<double> err_l;
	double sd, err;
	
		for (int i=0; i<9; i++) {
			getline(myfile, temp);
			myfile >> sd;
			std::getline(myfile, temp, ',');
			myfile >> err;
			std::getline(myfile, temp);
			cout << sd << endl;
			cout << err << endl;
			sd_l.push_back(sd);
			err_l.push_back(err);

		}
	
	sd_err_l.push_back(sd_l);
	sd_err_l.push_back(err_l);
	return sd_err_l;
}


void plot_both(vector<int>& x, vector<double>& first, vector<double>& first_err, vector<double>& second, vector<double>& second_err, const char *title)
{
	
	TGraphErrors *gr = new TGraphErrors();

	TGraphErrors *gr1 = new TGraphErrors();
	gr1->SetMarkerColor(kBlue);
	TCanvas *c1 = new TCanvas();
	gr->SetMinimum(0.);
	int n = 0;
	for (int i = 0; i<x.size(); i++) {
		n = gr->GetN();
		gr->SetPoint(n, x[i], first[i]);
		gr->SetPointError(n,0, first_err[i]);
	}
	
	gr->Draw("A*");

	n=0;
	for (int i = 0; i<x.size(); i++) {
		n = gr1->GetN();
		gr1->SetPoint(n, x[i], second[i]);
		gr1->SetPointError(n,0, second_err[i]);
	}
	gr1->Draw("*");
	
	TF1 *fit = new TF1("fit", "pol2", 0 ,0.017);
	gr->Fit("fit");
	TF1 *fit2 = new TF1("fit2", "pol2", 0, 0.017);
	fit2->SetLineColor(kBlue);
	gr1->Fit("fit2");

	gr->SetTitle(title);
}


void plot() 
{
/*
	fstream file;
	file.open("temp.csv");

	TGraph *gr = new TGraph();

	string trash, event, yield, decay_time, sigma;

	getline(file, trash);

	for (int i=0; i<11; i++) {

//		if(file.eof()) break;
		std::getline(file, event, ',');
		std::getline(file, yield, ',');
		std::getline(file, decay_time, ',');
		std::getline(file, sigma);
		
		gr->SetPoint(gr->GetN(), 1000*(std::stod(decay_time) / std::stod(yield)), std::stod(sigma));

	}
	
	file.close();
*/

//	TGraphErrors *gr = new TGraphErrors();
	vector<int> x = {1,5,10,25,50,75,100,150};
//	vector<double> Tyvek = {228.42,226.71,280.73,449.49,787.66,1212.51,1741.79,3690.93};
//	vector<double> Tyvek_err = {8.296,8.234,10.197,16.326,28.609,44.04,63.265,134.061};

	vector<double> ESR = {51.69,67.44,82.01,118.14,163.72,208.36,250,331.334};
	vector<double> ESR_err = {1.853,2.418,2.94,4.235,5.869,7.47,8.963,11.879};
	cout << "2222222222222222222222" << endl;
	vector<vector<double>> Tyvek = read_file("testfile.txt");
	
	cout << Tyvek[0][1] << endl;
	

//	plot_both(x, Tyvek[0], Tyvek[1], ESR, ESR_err, "test:x:y");
/*
	TGraphErrors *gr1 = new TGraphErrors();
	gr1->SetMarkerColor(kBlue);
	TCanvas *c1 = new TCanvas();
	gr->SetMinimum(0.);
	int n = 0;
	for (int i = 0; i <8; i++) {
		n = gr->GetN();
		gr->SetPoint(n, x[i], Tyvek[i]);
		gr->SetPointError(n,0, Tyvek_err[i]);
	}
	
	gr->Draw("A*");

	n=0;
	for (int i = 0; i <5; i++) {
		n = gr1->GetN();
		gr1->SetPoint(n, x[i], ESR[i]);
		gr1->SetPointError(n,0, ESR_err[i]);
	}
	gr1->Draw("*");
	
	TF1 *fit = new TF1("fit", "pol2", 0 ,0.017);
	gr->Fit("fit");
	TF1 *fit2 = new TF1("fit2", "pol2", 0, 0.017);
	gr1->Fit("fit2");
	gr->GetXaxis()->SetTitle("photon count");
	gr->GetYaxis()->SetTitle("sigma(ps)");
*/

}
