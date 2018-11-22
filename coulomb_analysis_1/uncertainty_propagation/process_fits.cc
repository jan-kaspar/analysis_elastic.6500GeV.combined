#include "TH1D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "../../stat.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

void ProcessOneGroup(const string &base_dir,
		const vector<string> &unc_labels,
		unsigned int seed_min, unsigned int seed_max,
		const vector<string> &fit_labels,
		TFile *f_out)
{
	printf("\n>> %s\n", base_dir.c_str());
	
	for (const auto &fit_label : fit_labels)
	{
		TDirectory *d_fit_label = f_out->mkdir(fit_label.c_str());

		printf("    %s\n", fit_label.c_str());

		for (const auto &unc_label : unc_labels)
		{
			TDirectory *d_unc_mode = d_fit_label->mkdir(unc_label.c_str());

			TH1D *h_de_rho = new TH1D("h_de_rho", "", 40, -0.10, +0.10);
			TH1D *h_de_a = new TH1D("h_de_a", "", 20, -150., +150.);
			TH1D *h_de_si_tot = new TH1D("h_de_si_tot", "", 40, -20., +20.);

			TGraph *g_de_rho_vs_de_a = new TGraph(); g_de_rho_vs_de_a->SetName("g_de_rho_vs_de_a");

			Stat stat(3);

			for (unsigned int seed = seed_min; seed <= seed_max; ++seed)
			{
				// get reco data
				char file[500];
				sprintf(file, "%s/%s/seed%u/%s/fit.root", base_dir.c_str(), unc_label.c_str(), seed, fit_label.c_str());
				TFile *f = TFile::Open(file);
				if (f == NULL)
					continue;

				TGraph *g_fit_data = (TGraph *) f->Get("g_fit_data");

				const double *v = g_fit_data->GetY();
				const double rho = v[4];
				const double a = v[6];
				const double si_tot = v[14];

				// calculate simu-reco deltas
				double rho_0 = 0.09, a_0 = 0., si_tot_0 = 110.;

				if (fit_label.find("approach1") == 0) rho_0 = 8.794886E-02, si_tot_0 = 1.117773E+02;
				if (fit_label.find("approach2_step_c") == 0) rho_0 = 8.596895E-02, si_tot_0 = 1.113469E+02;
				if (fit_label.find("approach3_step_d") == 0) rho_0 = 8.478408E-02, si_tot_0 = 1.103178E+02;
				if (fit_label.find("approach3_step_f") == 0) rho_0 = 9.905548E-02, si_tot_0 = 1.093493E+02;

				const double de_rho = rho - rho_0;
				const double de_a = a - a_0;
				const double de_si_tot = si_tot - si_tot_0;

				// skip outliers
				//if (abs(de_rho) > 0.05 || abs(de_a > 150.))
				if (abs(de_rho) > 0.05 || abs(de_si_tot > 10.))
					continue;

				h_de_rho->Fill(de_rho);
				h_de_a->Fill(de_a);
				h_de_si_tot->Fill(de_si_tot);

				int idx = g_de_rho_vs_de_a->GetN();
				g_de_rho_vs_de_a->SetPoint(idx, de_a, de_rho);

				stat.Fill(de_rho, de_a, de_si_tot);

				delete f;
			}

			gDirectory = d_unc_mode;

			h_de_rho->Write();
			h_de_a->Write();
			h_de_si_tot->Write();

			printf("        %20s: RMS(rho) = %.4f, RMS(si_tot) = %.3f\n", unc_label.c_str(), h_de_rho->GetRMS(), h_de_si_tot->GetRMS());
		}
	}
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// prepare output
	TFile *f_out = TFile::Open("process_fits.root", "recreate");

	{
		// settings
		string base_dir = "data";

		vector<string> unc_labels = {
			//"stat",
			//"stat+syst",
			"stat+syst+norm",
			"norm"
		};

		unsigned int seed_min = 1, seed_max = 100;

		vector<string> fit_labels = {
			"approach1",
			"approach2_step_c",
			"approach3_step_d",
			"approach3_step_f",
		};

		gDirectory = f_out;
		ProcessOneGroup(base_dir,unc_labels, seed_min, seed_max, fit_labels, f_out);
	}

	// clean up
	delete f_out;

	return 0;
}
