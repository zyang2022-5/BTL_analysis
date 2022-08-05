#include <vector>
#include <algorithm>
#include <cstdio>
#include <tuple>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

vector<tuple<vector<double>, vector<double>, double>> store_time(const char *filename)
{
	TFile *input = new TFile(filename, "read");
	TTree *tree = (TTree*)input->Get("Hits");
	TTree *tree_3 = (TTree*)input->Get("Scoring");

	double time;
	int event;
	double Z_pos;
	double E_dep;
	double muon_Z;
	tree->SetBranchAddress("ftime", &time);
	tree->SetBranchAddress("fEvent", &event);
	tree->SetBranchAddress("fZ", &Z_pos);
	tree_3->SetBranchAddress("fEdep", &E_dep);
	tree_3->SetBranchAddress("iZ", &muon_Z);

	int entries = tree->GetEntries();
	vector<tuple<vector<double>, vector<double>, double>> time_l;
	vector<double> temp_l_L;
	vector<double> temp_l_R;
	vector<int> hits;
	int test = 0;
	int init;
	int temp;
	
	while (test < entries-1) {
		
		init = test;
		temp_l_L = {};
		temp_l_R = {};
		tree->GetEntry(init);
		temp = event;
		tree_3->GetEntry(temp-1);

		for(int i = init; i < entries; i++){
			tree->GetEntry(i);
			if (temp != event) {
				break;
			}

			if (Z_pos > 20) {
				temp_l_R.push_back(time);
			} else {
				temp_l_L.push_back(time);
			}
			test = i+1;
		}
		
//		if (E_dep > 0.5) {
		time_l.push_back(make_tuple(temp_l_L, temp_l_R, muon_Z*1000.));
//		}

	}

	return time_l;
	
}

vector<int> count_photon(const char *filename) 
{
	TFile *input = new TFile(filename, "read");
	TTree *tree = (TTree*)input->Get("Hits");
	TTree *tree_3 = (TTree*)input->Get("Scoring");

	double time;
	int event;
	double Z_pos;
	double E_dep;
	tree->SetBranchAddress("ftime", &time);
	tree->SetBranchAddress("fEvent", &event);
	tree->SetBranchAddress("fZ", &Z_pos);
	tree_3->SetBranchAddress("fEdep", &E_dep);

	int entries = tree->GetEntries();
	vector<int> hits;
	int test = 0;
	int init;
	int temp;
	int counter = 0;
	
	while (test < entries-1) {
		
		init = test;
		tree->GetEntry(init);
		temp = event;
		tree_3->GetEntry(temp);
		counter = 0;

		for(int i = init; i < entries; i++){
			tree->GetEntry(i);
			if (temp != event) {
				break;
			}

			test = i+1;
			counter++;
		}
		
		if (E_dep > 0.5) {
			hits.push_back(counter);
		}

	}

	return hits;
}

//This function takes in a filename, a Ttree name, and a branch name  as the parameters and returns a vector of doubles containing the requested information in the designated branch of the given tree in the given file for every event.
vector<double> get_stuff(const char *filename, const char *treename, const char *branch_address)
{
	
	TFile *input = new TFile(filename, "read");
	TTree *tree = (TTree*)input->Get(treename);

	double content;
	tree->SetBranchAddress(branch_address, &content);

	vector<double> content_l;
	int entries = tree->GetEntries();

	for (int i=0; i<entries; i++) {
		tree->GetEntry(i);
		content_l.push_back(content);
	}
	return content_l;
}

vector<double> find_ave_path(const char *filename, int p_count_int, bool center) 
{
	TFile *input = new TFile(filename, "read");
	TTree *tree = (TTree*)input->Get("Hits");
	TTree *tree_3 = (TTree*)input->Get("Scoring");

	double path, time;
	int event;
	double Z_pos;
	double muon_Z_pos;
	tree->SetBranchAddress("ftime", &time);
	tree->SetBranchAddress("fMeanPath", &path);
	tree->SetBranchAddress("fEvent", &event);
	tree->SetBranchAddress("fZ", &Z_pos);
	tree_3->SetBranchAddress("iZ", &muon_Z_pos);

	vector<tuple<double, double>> plength_and_pos;
	vector<double> plength;
	vector<tuple<double, double>> time_l;

	int entries = tree->GetEntries();
	vector<double> temp_l_L;
	vector<double> temp_l_R;
	vector<double> temp_l_t_L;
	vector<double> temp_l_t_R;
	int test = 0;
	int init;
	int temp;

	int min_L_ind, min_R_ind;
	vector<double> low_ten_L;
	vector<double> low_ten_R;
	
	while (test < entries-1) {
		
		init = test;
		temp_l_L = {};
		temp_l_R = {};
		temp_l_t_L = {};
		temp_l_t_R = {};
		tree->GetEntry(init);
		temp = event;
		tree_3->GetEntry(event-1);

		for(int i = init; i < entries; i++){
			tree->GetEntry(i);
			if (temp != event) {
				break;
			}

			if (Z_pos > 20) {
				temp_l_R.push_back(path);
				temp_l_t_R.push_back(time);
			} else {
				temp_l_L.push_back(path);
				temp_l_t_L.push_back(time);
			}
			test = i+1;
		}
			
		min_L_ind = 0;
		min_R_ind = 0;

		low_ten_L.clear();
		low_ten_R.clear();
		
		for(int j = 0; j < p_count_int; j++) {
			min_L_ind = min_element(temp_l_t_L.begin(), temp_l_t_L.end()) - temp_l_t_L.begin();
			min_R_ind = min_element(temp_l_t_R.begin(), temp_l_t_R.end()) - temp_l_t_R.begin();
			low_ten_L.push_back(temp_l_L[min_L_ind]);
			temp_l_t_L.erase(temp_l_t_L.begin() + min_L_ind);
			temp_l_L.erase(temp_l_L.begin() + min_L_ind);
			low_ten_R.push_back(temp_l_R[min_R_ind]);
			temp_l_t_R.erase(temp_l_t_R.begin() + min_R_ind);
			temp_l_R.erase(temp_l_R.begin() + min_R_ind);


		}
//		cout << muon_Z_pos << endl;

		if (center) {
			plength.push_back(accumulate(low_ten_R.begin(), low_ten_R.end(), 0.0) / low_ten_R.size());
			plength.push_back(accumulate(low_ten_L.begin(), low_ten_L.end(), 0.0) / low_ten_L.size());
		} else {
			if (muon_Z_pos > 0) {
				plength_and_pos.push_back(make_tuple(accumulate(low_ten_R.begin(), low_ten_R.end(), 0.0) / low_ten_R.size(), (0.0285 - muon_Z_pos)*1000));
				
				plength_and_pos.push_back(make_tuple(accumulate(low_ten_L.begin(), low_ten_L.end(), 0.0) / low_ten_L.size(), (0.0285 + muon_Z_pos)*1000));
			} else if (muon_Z_pos < 0) {
				plength_and_pos.push_back(make_tuple(accumulate(low_ten_R.begin(), low_ten_R.end(), 0.0) / low_ten_R.size(), (0.0285 + abs(muon_Z_pos))*1000));
				
				plength_and_pos.push_back(make_tuple(accumulate(low_ten_L.begin(), low_ten_L.end(), 0.0) / low_ten_L.size(), (0.0285 + muon_Z_pos)*1000));
			} else if (muon_Z_pos == 0) {
				plength_and_pos.push_back(make_tuple(accumulate(low_ten_R.begin(), low_ten_R.end(), 0.0) / low_ten_R.size(), 0.0285*1000));
				
				plength_and_pos.push_back(make_tuple(accumulate(low_ten_L.begin(), low_ten_L.end(), 0.0) / low_ten_L.size(), 0.0285*1000));
			}
		}
	}

	if (center) {
		vector<double> center_ave_path;
		center_ave_path.push_back(accumulate(plength.begin(), plength.end(), 0.0) / plength.size());
		return center_ave_path;
	} else {
//		TCanvas *c1 = new TCanvas();
		TGraph *gr = new TGraph();
		
		for (int j=0; j<plength_and_pos.size(); j++) {
			gr->SetPoint(gr->GetN(), get<1>(plength_and_pos[j]), get<0>(plength_and_pos[j]));
		}
		char add[100], title[350];
		sprintf(add, "  ||  ave path length of %d ; muon hit distance (mm);ave path length", p_count_int);
		strcpy(title, filename);
		strcat(title, add);
		gr->SetTitle(title);
//		gr->Draw("A*");

		TF1 *fit = new TF1("fit", "pol4", 0, 58);
		gr->Fit("fit");


		vector<double> fit_para{fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3), fit->GetParameter(4)};
		return fit_para;
	}
}


// This function takes in a vector containing all the timing of every optical photon with 0.5 mev or above and a integer p_count_int that represents the first how many photons to average. It returns a list of doubles containing the average timing per event which is the first p_count_int number of photons time averaged.
vector<double> calc_min_time(vector<tuple<vector<double>, vector<double>, double>>& time_l, int p_count_int, const char *filename, bool center)
{
	double p_count = (double) p_count_int;
	double r_index = 1.845;
	double min_L;
	double min_R;
	double ave_path_L;
	double ave_path_R;
	double c_n = 3.0e+11 / r_index;   // mm / s
	double sigma_L, sigma_R, w_L, w_R;
	double pos_L, pos_R;
	vector<double> fit;

	fit = find_ave_path(filename, p_count_int, center);
	
	vector<double> low_ten;
	vector<double> low_ten_L;
	vector<double> low_ten_R;
	vector<double> min_time_l;
	double ave_min;

	for (int i = 0; i < time_l.size(); i++) {
		
		low_ten_L.clear();
		low_ten_R.clear();
		sigma_L = 1.0 / get<0>(time_l[i]).size();
		sigma_R = 1.0 / get<1>(time_l[i]).size();
		w_L = 1.0 / (sigma_L * sigma_L);
		w_R = 1.0 / (sigma_R * sigma_R);

		for(int j = 0; j < p_count_int; j++) {
			min_L = *min_element(get<0>(time_l[i]).begin(), get<0>(time_l[i]).end());
			min_R = *min_element(get<1>(time_l[i]).begin(), get<1>(time_l[i]).end());
			get<0>(time_l[i]).erase(std::remove(get<0>(time_l[i]).begin(), get<0>(time_l[i]).end(), min_L), get<0>(time_l[i]).end());
			low_ten_L.push_back(min_L);

			get<1>(time_l[i]).erase(std::remove(get<1>(time_l[i]).begin(), get<1>(time_l[i]).end(), min_R), get<1>(time_l[i]).end());
			low_ten_R.push_back(min_R);

		}

		if (center) {
			ave_path_L = fit[0];
			ave_path_R = fit[0];
		} else {
			pos_L = 28.5 + get<2>(time_l[i]);
			pos_R = 28.5 - get<2>(time_l[i]);

			ave_path_L = fit[0] + fit[1]*pos_L + fit[2]*pow(pos_L,2.0) + fit[3]*pow(pos_L, 3.0) + fit[4]*pow(pos_L, 4.0);
			ave_path_R = fit[0] + fit[1]*pos_R + fit[2]*pow(pos_R,2.0) + fit[3]*pow(pos_R, 3.0) + fit[4]*pow(pos_R, 4.0);
		}
		
		min_time_l.push_back((((accumulate(low_ten_L.begin(), low_ten_L.end(), 0.0)/p_count) - (ave_path_L / c_n)*1.0e12) * w_L + ((accumulate(low_ten_R.begin(), low_ten_R.end(), 0.0)/p_count) - (ave_path_R / c_n)*1.0e12) * w_R) / (w_L + w_R));

	}
	return min_time_l;
}


// This function takes in a vector containing all the timing of every optical photon with 0.5 mev or above and a integer p_count_int that represents the first how many photons to average. It returns a list of doubles containing the average timing per event which is the first p_count_int number of photons time averaged.
vector<double> calc_min_time_old(vector<tuple<vector<double>, vector<double>, double>>& time_l, int p_count_int)
{
	double p_count = (double) p_count_int;
	double min_L;
	double min_R;
	vector<double> low_ten;
	vector<double> low_ten_L;
	vector<double> low_ten_R;
	vector<double> min_time_l;
	double ave_min;

	for (int i = 0; i < time_l.size(); i++) {
		
		low_ten_L.clear();
		low_ten_R.clear();

		for(int j = 0; j < p_count_int; j++) {
			min_L = *min_element(get<0>(time_l[i]).begin(), get<0>(time_l[i]).end());
			min_R = *min_element(get<1>(time_l[i]).begin(), get<1>(time_l[i]).end());
			get<0>(time_l[i]).erase(std::remove(get<0>(time_l[i]).begin(), get<0>(time_l[i]).end(), min_L), get<0>(time_l[i]).end());
			low_ten_L.push_back(min_L);

			get<1>(time_l[i]).erase(std::remove(get<1>(time_l[i]).begin(), get<1>(time_l[i]).end(), min_R), get<1>(time_l[i]).end());
			low_ten_R.push_back(min_R);

		}
		
		min_time_l.push_back((accumulate(low_ten_L.begin(), low_ten_L.end(), 0.0) / p_count + accumulate(low_ten_L.begin(), low_ten_L.end(), 0.0) / p_count) / 2.0);

	}
	return min_time_l;
}

vector<double>  hist_simple(const char *hist_name, vector<double>& min_time_l)
{
//	TCanvas *c1 = new TCanvas();
//	gStyle->SetOptStat(0);

	double x_low = *min_element(min_time_l.begin(), min_time_l.end()) - 5;
	double x_high = *max_element(min_time_l.begin(), min_time_l.end()) + 5;
	
	vector<double> sd_err;

	TH1F *hist = new TH1F("hist", hist_name, 150, x_low, x_high);

	for (int i = 0; i < min_time_l.size(); i++) {
		hist->Fill(min_time_l[i]);

	}
	
	hist->GetXaxis()->SetTitle("time (ps)");
	hist->GetYaxis()->SetTitle("event count");
	double sigma = hist->GetRMS();
	double err = sigma/sqrt(min_time_l.size());
//	hist->Draw();	

	sd_err.push_back(sigma);
	sd_err.push_back(err);
//	c1->SaveAs(hist_name, "pdf");
	return sd_err;

}


void twod_hist(vector<double>& x, vector<double>& y, double x_s, double x_end, double y_s, double y_end, int x_bin, int y_bin, const char *title, const char *x_title, const char *y_title)
{
	int entry_count = x.size();

	double x_low = *min_element(x.begin(), x.end());
	double y_low = *min_element(y.begin(), y.end());	      double x_high = *max_element(x.begin(), x.end());
	double y_high = *max_element(y.begin(), y.end());

	TH2F *hist = new TH2F("hist", title, x_bin, x_low, x_high, y_bin, y_low, y_high);
	
	for (int i=0; i<entry_count; i++) {
		hist->Fill(x[i],y[i]);
	}

	TCanvas *c1 = new TCanvas();
	gStyle->SetPalette(kRainBow);
	hist->GetXaxis()->SetTitle(x_title);
	hist->GetXaxis()->SetRangeUser(x_s, x_end);
	hist->GetYaxis()->SetTitle(y_title);
	hist->GetYaxis()->SetRangeUser(y_s, y_end);
	hist->Draw("colz");
}

void write_file(vector<vector<double>> sd_err_l, const char *filename)
{
	ofstream myfile (filename);
	if (myfile.is_open()){
		for (int i=0; i<sd_err_l.size(); i++){
			myfile << sd_err_l[i][0] << "," << sd_err_l[i][1] << "\n";
		}
		myfile.close();
	}
}


void plot_both(vector<int>& x, vector<double>& first, vector<double>& first_err, vector<double>& second, vector<double>& second_err, const char *title)
{

        TGraphErrors *gr = new TGraphErrors();
        TGraphErrors *gr1 = new TGraphErrors();

        gr->SetMarkerColor(kRed);
	gr->SetLineStyle(1);
	gr->SetLineColor(kRed);
	gr->SetMarkerStyle(22);

       	gr1->SetMarkerColor(kBlue);
	gr1->SetLineStyle(2);
	gr1->SetLineColor(kBlue);
	gr1->SetMarkerStyle(kStar);

        TCanvas *c1 = new TCanvas();
	c1->SetLogy();
        gr->SetMinimum(7.);
	gr->SetMaximum(200000.);
        int n = 0;
        for (int i = 0; i<x.size(); i++) {
                n = gr->GetN();
                gr->SetPoint(n, x[i], first[i]);
                gr->SetPointError(n,0, first_err[i]);
        }

        gr->Draw("AP");

        n=0;
        for (int i = 0; i<x.size(); i++) {
                n = gr1->GetN();
                gr1->SetPoint(n, x[i], second[i]);
                gr1->SetPointError(n,0, second_err[i]);
        }
        gr1->Draw("P");
	


        TF1 *fit = new TF1("fit", "pol2", 0 , 40);
	fit->SetLineColor(kRed);
	fit->SetLineStyle(1);
        gr->Fit("fit");

        TF1 *fit1 = new TF1("fit1", "pol2", 0, 40);
        fit1->SetLineColor(kBlue);
	fit1->SetLineStyle(2);
        gr1->Fit("fit1");

	gr->SetTitle(title);

	auto legend = new TLegend(0.7, 0.6, 0.85, 0.75);
	legend->AddEntry(gr, "weighted ave with path");
	legend->AddEntry(gr1, "old formula no path no weight");

	legend->Draw();
}

void plot_all(vector<int>& x, vector<double>& first, vector<double>& first_err, vector<double>& second, vector<double>& second_err, vector<double>& third, vector<double>& third_err, vector<double>& fourth, vector<double>& fourth_err, vector<double>& fifth, vector<double>& fifth_err, vector<double>& sixth, vector<double>& sixth_err, vector<double>& seventh, vector<double>& seventh_err, vector<double>& eighth, vector<double>& eighth_err, const char *title)
{

        TGraphErrors *gr = new TGraphErrors();
        TGraphErrors *gr1 = new TGraphErrors();
        TGraphErrors *gr2 = new TGraphErrors();
        TGraphErrors *gr3 = new TGraphErrors();
        TGraphErrors *gr4 = new TGraphErrors();
        TGraphErrors *gr5 = new TGraphErrors();
        TGraphErrors *gr6 = new TGraphErrors();
        TGraphErrors *gr7 = new TGraphErrors();

        gr->SetMarkerColor(kRed);
	gr->SetLineStyle(1);
	gr->SetLineColor(kRed);
	gr->SetMarkerStyle(kStar);

       	gr1->SetMarkerColor(kBlue);
	gr1->SetLineStyle(1);
	gr1->SetLineColor(kBlue);
	gr1->SetMarkerStyle(kStar);

	gr2->SetMarkerColor(kRed);
	gr2->SetLineStyle(2);
	gr2->SetLineColor(kRed);
	gr2->SetMarkerStyle(kCircle);

	gr3->SetLineStyle(2);
	gr3->SetLineColor(kBlue);
	gr3->SetMarkerColor(kBlue);
	gr3->SetMarkerStyle(kCircle);

	gr4->SetLineStyle(1);
	gr4->SetLineColor(3);
	gr4->SetMarkerColor(3);
	gr4->SetMarkerStyle(kStar);

	gr5->SetLineStyle(2);
	gr5->SetLineColor(3);
	gr5->SetMarkerColor(3);
	gr5->SetMarkerStyle(kCircle);
	
	gr6->SetLineStyle(2);
	gr6->SetLineColor(6);
	gr6->SetMarkerColor(6);
	gr6->SetMarkerStyle(22);

	gr7->SetLineStyle(2);
	gr7->SetLineColor(7);
	gr7->SetMarkerColor(7);
	gr7->SetMarkerStyle(22);


        TCanvas *c1 = new TCanvas();
	c1->SetLogy();
        gr->SetMinimum(7.);
	gr->SetMaximum(200000.);
        int n = 0;
        for (int i = 0; i<x.size(); i++) {
                n = gr->GetN();
                gr->SetPoint(n, x[i], first[i]);
                gr->SetPointError(n,0, first_err[i]);
        }

        gr->Draw("AP");

        n=0;
        for (int i = 0; i<x.size(); i++) {
                n = gr1->GetN();
                gr1->SetPoint(n, x[i], second[i]);
                gr1->SetPointError(n,0, second_err[i]);
        }
        gr1->Draw("P");
	
        for (int i = 0; i<x.size(); i++) {
                n = gr2->GetN();
                gr2->SetPoint(n, x[i], third[i]);
                gr2->SetPointError(n,0, third_err[i]);
        }
        gr2->Draw("P");

        for (int i = 0; i<x.size(); i++) {
                n = gr3->GetN();
                gr3->SetPoint(n, x[i], fourth[i]);
                gr3->SetPointError(n,0, fourth_err[i]);
        }
        gr3->Draw("P");

        for (int i = 0; i<x.size(); i++) {
                n = gr4->GetN();
                gr4->SetPoint(n, x[i], fifth[i]);
                gr4->SetPointError(n,0, fifth_err[i]);
        }
        gr4->Draw("P");


        for (int i = 0; i<x.size(); i++) {
                n = gr5->GetN();
                gr5->SetPoint(n, x[i], sixth[i]);
                gr5->SetPointError(n,0, sixth_err[i]);
        }
        gr5->Draw("P");


        for (int i = 0; i<x.size(); i++) {
                n = gr6->GetN();
                gr6->SetPoint(n, x[i], seventh[i]);
                gr6->SetPointError(n,0, seventh_err[i]);
        }
        gr6->Draw("P");

        for (int i = 0; i<x.size(); i++) {
                n = gr7->GetN();
                gr7->SetPoint(n, x[i], eighth[i]);
                gr7->SetPointError(n,0, eighth_err[i]);
        }
        gr7->Draw("P");


        TF1 *fit = new TF1("fit", "pol2", 0 , 40);
        gr->Fit("fit");
        TF1 *fit1 = new TF1("fit1", "pol2", 0, 40);
        fit1->SetLineColor(kBlue);
        gr1->Fit("fit1");
        
        TF1 *fit2 = new TF1("fit2", "pol2", 0, 40);
        fit2->SetLineColor(kRed);
	fit2->SetLineStyle(2);
        gr2->Fit("fit2");
 
        TF1 *fit3 = new TF1("fit3", "pol2", 0, 40);
        fit3->SetLineColor(kBlue);
	fit3->SetLineStyle(2);
        gr3->Fit("fit3");

        TF1 *fit4 = new TF1("fit4", "pol2", 0, 0.40);
        fit4->SetLineColor(3);
	fit4->SetLineStyle(1);
        gr4->Fit("fit4");

        TF1 *fit5 = new TF1("fit5", "pol2", 0, 40);
        fit5->SetLineColor(3);
	fit5->SetLineStyle(2);
        gr5->Fit("fit5");

        TF1 *fit6 = new TF1("fit6", "pol2", 0, 40);
        fit6->SetLineColor(6);
	fit6->SetLineStyle(2);
        gr6->Fit("fit6");

        TF1 *fit7 = new TF1("fit7", "pol2", 0, 40);
        fit7->SetLineColor(7);
	fit7->SetLineStyle(2);
        gr7->Fit("fit7");

	gr->SetTitle(title);

	auto legend = new TLegend(0.7, 0.6, 0.85, 0.75);
	legend->AddEntry(gr, "center Ideal tyvek");
	legend->AddEntry(gr1, "center ESR");
	legend->AddEntry(gr2, "rand Ideal tyvek");
	legend->AddEntry(gr3, "rand ESR");
	legend->AddEntry(gr4, "99.9 specular");
	legend->AddEntry(gr5, "rand 99.9 specular");
	legend->AddEntry(gr6, "airgap ri=1.0 lambertian");
	legend->AddEntry(gr7, "airgap ri=1.6 lambertain");

	legend->Draw();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////Main Function//////////////////////////////////////////////////////////////////////////////////////////////////////////


void extr_time()
{
//	vector<double> fit_para = find_ave_path("lamb_gap_rindex1.6_muon_rand_500entries0.root", 100);
//	printf("%f \n %f", fit_para[0], fit_para[1]);

	vector<tuple<vector<double>, vector<double>, double>> time_l_ESR = store_time("ESR_muon_center_500_entries0.root");
	
	vector<tuple<vector<double>, vector<double>, double>> time_l_Tyvek = store_time("ideal_tyvek_muon_center_500_entries0.root");

	vector<tuple<vector<double>, vector<double>, double>> time_l_ESR_rand = store_time("ESR_muon_500_entries0.root");
	
	vector<tuple<vector<double>, vector<double>, double>> time_l_Tyvek_rand = store_time("tyvek_muon_500_entries0.root");
	vector<tuple<vector<double>, vector<double>, double>> time_l_99_spec = store_time("99_dot_9_percent_specular_muon_center0.root");
	vector<tuple<vector<double>, vector<double>, double>> time_l_99_spec_rand = store_time("99_specular_muon_randZ_500entries0.root");
	
	vector<tuple<vector<double>, vector<double>, double>> time_l_airgap = store_time("lamb_airgap_rindex1_muon_rand_500entries0.root");
	vector<tuple<vector<double>, vector<double>, double>> time_l_1_6_spike = store_time("lamb_gap_rindex1.6_muon_rand_500entries0.root");
	
	
	vector<double> sd_l_ESR, err_l_ESR, sd_l_tyvek, err_l_tyvek, sd_l_spec, err_l_spec, sd_l_airgap, err_l_airgap, sd_l_1_6spike, err_l_1_6spike;

	vector<double> sd_l_ESR_rand, err_l_ESR_rand, sd_l_tyvek_rand, err_l_tyvek_rand, sd_l_spec_rand, err_l_spec_rand;

	vector<double> sd_l_ESR_old, err_l_ESR_old, sd_l_tyvek_old, err_l_tyvek_old, sd_l_spec_old, err_l_spec_old, sd_l_airgap_old, err_l_airgap_old, sd_l_1_6spike_old, err_l_1_6spike_old;

	vector<double> sd_l_ESR_rand_old, err_l_ESR_rand_old, sd_l_tyvek_rand_old, err_l_tyvek_rand_old, sd_l_spec_rand_old, err_l_spec_rand_old;

	vector<int> photon_num = {1,5,10,25,50,75,100,150};
	for (int j : photon_num) {

	vector<double>min_time_l_ESR = calc_min_time(time_l_ESR, j, "ESR_muon_center_500_entries0.root", true);
	vector<double> min_time_l_Tyvek = calc_min_time(time_l_Tyvek, j, "ideal_tyvek_muon_center_500_entries0.root", true);
	vector<double> min_time_l_ESR_rand = calc_min_time(time_l_ESR_rand, j, "ESR_muon_500_entries0.root", false);
	vector<double> min_time_l_Tyvek_rand = calc_min_time(time_l_Tyvek_rand, j, "tyvek_muon_500_entries0.root", false);
	vector<double> min_time_l_spec = calc_min_time(time_l_99_spec, j, "99_dot_9_percent_specular_muon_center0.root", true);
	vector<double> min_time_l_spec_rand = calc_min_time(time_l_99_spec_rand, j, "99_specular_muon_randZ_500entries0.root", false);
	vector<double> min_time_l_airgap = calc_min_time(time_l_airgap, j, "lamb_airgap_rindex1_muon_rand_500entries0.root", false);
	vector<double> min_time_l_1_6_spike = calc_min_time(time_l_1_6_spike, j, "lamb_gap_rindex1.6_muon_rand_500entries0.root", false);
	

	vector<double> min_time_l_ESR_old = calc_min_time_old(time_l_ESR, j);
	vector<double> min_time_l_Tyvek_old = calc_min_time_old(time_l_Tyvek, j);
	vector<double> min_time_l_ESR_rand_old = calc_min_time_old(time_l_ESR_rand, j);
	vector<double> min_time_l_Tyvek_rand_old = calc_min_time_old(time_l_Tyvek_rand, j);
	vector<double> min_time_l_spec_old = calc_min_time_old(time_l_99_spec, j);
	vector<double> min_time_l_spec_rand_old = calc_min_time_old(time_l_99_spec_rand, j);
	vector<double> min_time_l_airgap_old = calc_min_time_old(time_l_airgap, j);
	vector<double> min_time_l_1_6_spike_old = calc_min_time_old(time_l_1_6_spike, j);

	char name1[100], name2[100], name3[100], name4[100], name5[100], name6[100], name7[100], name8[100];
	sprintf(name1, "ESR average timing of first %d", j);
	sprintf(name2, "tyvek average timing of first %d", j);
	sprintf(name3, "rand ESR average timing of first %d", j);
	sprintf(name4, "rand tyvek average timing of the first %d", j);
	sprintf(name5, "99.9 reflective specular average timing of the first %d", j);
	sprintf(name6, "rand 99.9 reflective specular average timing of the first %d", j);
	sprintf(name7, "rand airgap rindex 1.0 average timing of the first %d", j);
	sprintf(name8, "rand airgap rindex 1.6 average timing of the first %d", j);


	vector<double> temp_tup_ESR = hist_simple(name1, min_time_l_ESR);
	vector<double> temp_tup_tyvek = hist_simple(name2, min_time_l_Tyvek);
	vector<double> temp_tup_tyvek_rand = hist_simple(name4, min_time_l_Tyvek_rand);
	vector<double> temp_tup_ESR_rand = hist_simple(name3, min_time_l_ESR_rand);
	vector<double> temp_tup_spec = hist_simple(name5, min_time_l_spec);
	vector<double> temp_tup_spec_rand = hist_simple(name6, min_time_l_spec_rand);
	vector<double> temp_tup_airgap_rand = hist_simple(name7, min_time_l_airgap);
	vector<double> temp_tup_1_6_spike_rand = hist_simple(name8, min_time_l_1_6_spike);
	sd_l_ESR.push_back(temp_tup_ESR[0]);
	err_l_ESR.push_back(temp_tup_ESR[1]);
	sd_l_tyvek.push_back(temp_tup_tyvek[0]);
	err_l_tyvek.push_back(temp_tup_tyvek[1]);
	sd_l_spec.push_back(temp_tup_spec[0]);
	err_l_spec.push_back(temp_tup_spec[1]);

	sd_l_ESR_rand.push_back(temp_tup_ESR_rand[0]);
	err_l_ESR_rand.push_back(temp_tup_ESR_rand[1]);
	sd_l_tyvek_rand.push_back(temp_tup_tyvek_rand[0]);
	err_l_tyvek_rand.push_back(temp_tup_tyvek_rand[1]);
	sd_l_spec_rand.push_back(temp_tup_spec_rand[0]);
	err_l_spec_rand.push_back(temp_tup_spec_rand[1]);
	sd_l_airgap.push_back(temp_tup_airgap_rand[0]);
	err_l_airgap.push_back(temp_tup_airgap_rand[1]);
	sd_l_1_6spike.push_back(temp_tup_1_6_spike_rand[0]);
	err_l_1_6spike.push_back(temp_tup_1_6_spike_rand[1]);



	vector<double> temp_tup_ESR_old = hist_simple(name1, min_time_l_ESR_old);
	vector<double> temp_tup_tyvek_old = hist_simple(name2, min_time_l_Tyvek_old);
	vector<double> temp_tup_tyvek_rand_old = hist_simple(name4, min_time_l_Tyvek_rand_old);
	vector<double> temp_tup_ESR_rand_old = hist_simple(name3, min_time_l_ESR_rand_old);
	vector<double> temp_tup_spec_old = hist_simple(name5, min_time_l_spec_old);
	vector<double> temp_tup_spec_rand_old = hist_simple(name6, min_time_l_spec_rand_old);
	vector<double> temp_tup_airgap_rand_old = hist_simple(name7, min_time_l_airgap_old);
	vector<double> temp_tup_1_6_spike_rand_old = hist_simple(name8, min_time_l_1_6_spike_old);
	sd_l_ESR_old.push_back(temp_tup_ESR_old[0]);
	err_l_ESR_old.push_back(temp_tup_ESR_old[1]);
	sd_l_tyvek_old.push_back(temp_tup_tyvek_old[0]);
	err_l_tyvek_old.push_back(temp_tup_tyvek_old[1]);
	sd_l_spec_old.push_back(temp_tup_spec_old[0]);
	err_l_spec_old.push_back(temp_tup_spec_old[1]);

	sd_l_ESR_rand_old.push_back(temp_tup_ESR_rand_old[0]);
	err_l_ESR_rand_old.push_back(temp_tup_ESR_rand_old[1]);
	sd_l_tyvek_rand_old.push_back(temp_tup_tyvek_rand_old[0]);
	err_l_tyvek_rand_old.push_back(temp_tup_tyvek_rand_old[1]);
	sd_l_spec_rand_old.push_back(temp_tup_spec_rand_old[0]);
	err_l_spec_rand_old.push_back(temp_tup_spec_rand_old[1]);
	sd_l_airgap_old.push_back(temp_tup_airgap_rand_old[0]);
	err_l_airgap_old.push_back(temp_tup_airgap_rand_old[1]);
	sd_l_1_6spike_old.push_back(temp_tup_1_6_spike_rand_old[0]);
	err_l_1_6spike_old.push_back(temp_tup_1_6_spike_rand_old[1]);
	}
	vector<int> x = {1,5,10,25,50,75,100,150};

//	plot_all(x, sd_l_tyvek, err_l_tyvek, sd_l_ESR, err_l_ESR, sd_l_tyvek_rand, err_l_tyvek_rand, sd_l_ESR_rand, err_l_ESR_rand, sd_l_spec, err_l_spec, sd_l_spec_rand, err_l_spec_rand, sd_l_airgap, err_l_airgap, sd_l_1_6spike, err_l_1_6spike, "Time Resolution vs. Threshold; Number of Early Photons Averaged; Time Resolution (ps)");
	
	plot_both(x, sd_l_tyvek, err_l_tyvek, sd_l_tyvek_old, err_l_tyvek_old, "tyvek center; Number of Early Photons Averaged; Time Resolution (ps)");
	plot_both(x, sd_l_ESR, err_l_ESR, sd_l_ESR_old, err_l_ESR_old, "ESR center; Number of Early Photons Averaged; Time Resolution (ps)");
	plot_both(x, sd_l_tyvek_rand, err_l_tyvek_rand, sd_l_tyvek_rand_old, err_l_tyvek_rand_old, "tyvek rand; Number of Early Photons Averaged; Time Resolution (ps)");
	plot_both(x, sd_l_ESR_rand, err_l_ESR_rand, sd_l_ESR_rand_old, err_l_ESR_rand_old, "ESR rand; Number of Early Photons Averaged; Time Resolution (ps)");
	plot_both(x, sd_l_spec, err_l_spec, sd_l_spec_old, err_l_spec_old, "99.9 specular center; Number of Early Photons Averaged; Time Resolution (ps)");
	plot_both(x, sd_l_spec_rand, err_l_spec_rand, sd_l_spec_rand_old, err_l_spec_rand_old, "99.9 specular rand; Number of Early Photons Averaged; Time Resolution (ps)");
	plot_both(x, sd_l_airgap, err_l_airgap, sd_l_airgap_old, err_l_airgap_old, "airgap ri=1.0 rand; Number of Early Photons Averaged; Time Resolution (ps)");
	plot_both(x, sd_l_1_6spike, err_l_1_6spike, sd_l_1_6spike_old, err_l_1_6spike_old, "airgap ri=1.6 rand; Number of Early Photons Averaged; Time Resolution (ps)");

	
	/*
	vector<double> X_pos = get_stuff("ideal_tyvek_muon_center_500_entries0.root", "ScoringKilled", "fX");
	vector<double> Y_pos = get_stuff("ideal_tyvek_muon_center_500_entries0.root", "ScoringKilled", "fY");
	vector<double> Z_pos = get_stuff("ideal_tyvek_muon_center_500_entries0.root", "ScoringKilled", "fZ");

	twod_hist(X_pos, Y_pos, -1.7, 1.7, -1.7, 1.7, 1500, 1500, "lambertian reflector photon killed position(x,y)", "X pos", "Y pos");
	twod_hist(X_pos, Z_pos, -1.7, 1.7, -31., 31., 1500, 150, "lambertian reflector photon killed position(x,z)", "X pos", "Z pos");
	twod_hist(Y_pos, Z_pos, -1.7, 1.7, -31., 31., 1500, 150, "lambertian reflector photon killed position(y,z)", "Y pos", "Z pos");
*/

/*
//	TH1F *ESR = new TH1F("ESR", "ESR hits count", 150, *min_element(EST_hits.begin(),ESR_hits.end()), *max_element(ESR_hits.begin(), ESR_hits.end()));
	
	int start = 0;
	int end = 3700;

	TH1F *ESR = new TH1F("ESR", "photon hits count", 150, start, end);
//	TH1F *Tyvek = new TH1F("Tyvek", "Tyvek hits count", 150, *min_element(Tyvek_hits.begin(),Tyvek_hits.end()), *max_element(Tyvek_hits.begin(), Tyvek_hits.end()));

	TH1F *Tyvek = new TH1F("Tyvek", "Tyvek hits count", 150, start, end);
	for (int i = 0; i < ESR_hits.size(); i++) {
		
		ESR->Fill(ESR_hits[i]);
	}

	for (int i = 0; i < Tyvek_hits.size(); i++) {
		
		Tyvek->Fill(Tyvek_hits[i]);
	}
	ESR->SetLineColor(kRed);
	Tyvek->SetLineColor(kBlack);
	ESR->GetXaxis()->SetTitle("hits count");
	ESR->GetYaxis()->SetTitle("event count");
	ESR->Draw("SAME");	
	Tyvek->Draw("SAME");
	auto legend = new TLegend(0.7, 0.6, 0.85, 0.75);
	legend->AddEntry(ESR, "ESR hits");
	legend->AddEntry(Tyvek, "Tyvek hits");
	legend->Draw();
*/
//	input->Close();
}
