#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// defaults
	int seed_min = 1;
	int seed_max = 10;

	// cross-section models
	struct CSModel
	{
		string label;
		double A_ref;
		double rho_ref;
	};

	vector<CSModel> csModels = {
		{ "exp,rho=0.10", 580., 0.10 },
		{ "exp,rho=0.14", 580., 0.14 },
		{ "non-exp,rho=0.10", 580., 0.10 },
		{ "non-exp,rho=0.14", 580., 0.14 },
	};

	// normalisation adjustments
/*
	vector<double> normAdjustments = {
		0.800,
		0.850,
		0.880,
		0.900,
		0.930,
		0.950,
		0.970,
		0.990,
		1.000,
		1.010,
		1.030,
		1.050,
		1.070,
		1.100,
		1.120,
		1.150,
		1.200,
	};
*/

	vector<double> normAdjustments = {
		0.800,
		0.900,
		0.950,
		1.000,
		1.050,
		1.100,
		1.200,
	};

	// prepare output
	// TODO
	//TFile *f_out = TFile::Open("process_fits.root", "recreate");
	TFile *f_out = TFile::Open("process_fits2.root", "recreate");

	// make histograms
	vector<TH1D *> csHistograms;

	for (auto &m : csModels)
	{
		printf("* %s\n", m.label.c_str());
		TDirectory *dir_model = f_out->mkdir(m.label.c_str());

		TGraph *g_mean_De_rho = new TGraph(); g_mean_De_rho->SetName("g_mean_De_rho"); g_mean_De_rho->SetTitle(";norm adjustment;mean_De_rho");
		TGraph *g_rms_De_rho = new TGraph(); g_rms_De_rho->SetName("g_rms_De_rho"); g_rms_De_rho->SetTitle(";norm adjustment;rms_De_rho");
		TGraph *g_mean_de_A = new TGraph(); g_mean_de_A->SetName("g_mean_de_A"); g_mean_de_A->SetTitle(";norm adjustment;mean_de_A");
		TGraph *g_rms_de_A = new TGraph(); g_rms_de_A->SetName("g_rms_de_A"); g_rms_de_A->SetTitle(";norm adjustment;rms_de_A");

		for (double normAdjustment : normAdjustments)
		{
			char buf[100];
			sprintf(buf, "%.3f", normAdjustment);
			TDirectory *dir_norm = dir_model->mkdir(buf);

			TH1D *h_De_rho = new TH1D("", ";De_rho", 50, -0.05, +0.05);
			TH1D *h_de_A = new TH1D("", ";de_A", 50, -0.25, +0.25);

			for (int seed = seed_min; seed <= seed_max; ++seed)
			{
				// TODO
				//sprintf(buf, "fits/%s/%.3f/fit_seed%i.root", m.label.c_str(), normAdjustment, seed);
				sprintf(buf, "fits2/%s/%.3f/fit_seed%i.root", m.label.c_str(), normAdjustment, seed);
				TFile *f_in = TFile::Open(buf);

				TGraph *g_results = (TGraph *) f_in->Get("g_results");

				double A, A_e;
				g_results->GetPoint(0, A, A_e);

				double rho, rho_e;
				g_results->GetPoint(3, rho, rho_e);

				double de_A = (A - m.A_ref) / m.A_ref;
				double De_rho = rho - m.rho_ref;

				//printf("seed=%2i, A=%.1f, rho=%.3f | de_A=%.3f, De_rho=%.3f\n", seed, A, rho, de_A, De_rho);

				h_de_A->Fill(de_A);
				h_De_rho->Fill(De_rho);

				delete f_in;
			}

			gDirectory = dir_norm;
			h_De_rho->Write("h_De_rho");
			h_de_A->Write("h_de_A");

			int idx = g_mean_De_rho->GetN();
			g_mean_De_rho->SetPoint(idx, normAdjustment, h_De_rho->GetMean());
			g_rms_De_rho->SetPoint(idx, normAdjustment, h_De_rho->GetRMS());
			g_mean_de_A->SetPoint(idx, normAdjustment, h_de_A->GetMean());
			g_rms_de_A->SetPoint(idx, normAdjustment, h_de_A->GetRMS());
		}

		gDirectory = dir_model;

		g_mean_De_rho->Write();
		g_rms_De_rho->Write();
		g_mean_de_A->Write();
		g_rms_de_A->Write();
	}

	// clean up
	delete f_out;

	return 0;
}
