#include "TH1D.h"
#include "TFile.h"
#include "TGraph.h"

#include "../../stat.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// settings
	vector<string> models = {
		"exp1,rho=0.10",
		"exp1,rho=0.14",
		"exp3,rho=0.10",
		"exp3,rho=0.14",
	};

	vector<string> unc_modes = {
		"stat",
		"stat+syst",
		"stat+syst+norm",
	};

	vector<string> norms = {
		"0.00",
	};

	vector<string> binnings = {
		"ob-2-10-0.05",
	};

	unsigned int seed_min = 1, seed_max = 100;

	vector<string> t_maxs = {
		"0.15"
	};

	// open model file
	TFile *f_in_model = TFile::Open("build_models.root");

	// prepare output
	TFile *f_out = TFile::Open("process_fits.root", "recreate");

	// process all
	for (const auto &model : models)
	{
		printf("* %s\n", model.c_str());

		TDirectory *d_model = f_out->mkdir(model.c_str());

		// get model data
		TGraph *g_model_data = (TGraph *) f_in_model->Get((model + "/g_data").c_str());
		const double a0 = g_model_data->GetY()[0];
		const double rho0 = g_model_data->GetY()[1];

		for (const auto &unc_mode : unc_modes)
		{
			TDirectory *d_unc_mode = d_model->mkdir(unc_mode.c_str());

			for (const auto &norm : norms)
			{
				TDirectory *d_norm = d_unc_mode->mkdir(("norm"+norm).c_str());

				for (const auto &binning : binnings)
				{
					TDirectory *d_binning = d_norm->mkdir(binning.c_str());

					for (const auto &t_max : t_maxs)
					{
						TDirectory *d_t_max = d_binning->mkdir(("tmax" + t_max).c_str());

						TH1D *h_de_rho = new TH1D("h_de_rho", "", 20, -0.03, +0.03);
						TH1D *h_de_a = new TH1D("h_de_a", "", 20, -60., +60.);

						Stat stat(2);

						for (unsigned int seed = seed_min; seed <= seed_max; ++seed)
						{
							string simuDir = model+"/"+unc_mode+"/norm"+norm+"/"+binning;
							string fitDir = model.substr(0, 4) + "," + unc_mode + "," + t_max;

							char file[500];
							sprintf(file, "simu-fit-study/%s/seed%u/%s/fit.root", simuDir.c_str(), seed, fitDir.c_str());
							TFile *f = TFile::Open(file);
							TGraph *g_fit_data = (TGraph *) f->Get("g_fit_data");

							const double *v = g_fit_data->GetY();
							const double rho = v[4];
							const double a = v[6];

							h_de_rho->Fill(rho - rho0);
							h_de_a->Fill(a - a0);

							stat.Fill(rho - rho0, a - a0);

							delete f;
						}

						gDirectory = d_t_max;

						h_de_rho->Write();
						h_de_a->Write();

						TGraph *g_data = new TGraph();
						g_data->SetPoint(0, stat.GetMean(0), stat.GetMeanUnc(0));
						g_data->SetPoint(1, stat.GetStdDev(0), stat.GetStdDevUnc(0));
						g_data->SetPoint(2, stat.GetMean(1), stat.GetMeanUnc(1));
						g_data->SetPoint(3, stat.GetStdDev(1), stat.GetStdDevUnc(1));
						g_data->Write("g_data");
					}
				}
			}
		}
	}

	// clean up
	delete f_in_model;
	delete f_out;

	return 0;
}
