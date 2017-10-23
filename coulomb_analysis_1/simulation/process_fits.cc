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
		const vector<string> &models, const vector<string> &unc_modes, const vector<string> norms,
		const vector<string> binnings, unsigned int seed_min, unsigned int seed_max, const vector<string> t_maxs,
		TFile *f_in_model, TFile *f_out)
{
	printf("\n>> %s\n", base_dir.c_str());

	TDirectory *d_base = f_out->mkdir(base_dir.c_str());

	for (const auto &model : models)
	{
		printf("* %s\n", model.c_str());

		TDirectory *d_model = d_base->mkdir(model.c_str());

		// get model data
		TGraph *g_model_data = (TGraph *) f_in_model->Get((model + "/g_data").c_str());
		const double a0 = g_model_data->GetY()[0];
		const double rho0 = g_model_data->GetY()[1];

		for (const auto &unc_mode : unc_modes)
		{
			TDirectory *d_unc_mode = d_model->mkdir(unc_mode.c_str());

			for (const auto &binning : binnings)
			{
				TDirectory *d_binning = d_unc_mode->mkdir(binning.c_str());

				for (const auto &t_max : t_maxs)
				{
					TDirectory *d_t_max = d_binning->mkdir(("tmax" + t_max).c_str());

					TGraphErrors *g_rho_bias_vs_norm = new TGraphErrors();
					TGraphErrors *g_rho_rms_vs_norm = new TGraphErrors();
					TGraphErrors *g_a_bias_vs_norm = new TGraphErrors();
					TGraphErrors *g_a_rms_vs_norm = new TGraphErrors();

					for (const auto &norm : norms)
					{
						TDirectory *d_norm = d_t_max->mkdir(("norm"+norm).c_str());
						gDirectory = d_norm;

						TH1D *h_de_rho = new TH1D("h_de_rho", "", 20, -0.03, +0.03);
						TH1D *h_de_a = new TH1D("h_de_a", "", 20, -60., +60.);

						Stat stat(2);

						for (unsigned int seed = seed_min; seed <= seed_max; ++seed)
						{
							string simuDir = model+"/"+unc_mode+"/norm"+norm+"/"+binning;
							string fitDir = model.substr(0, 4) + "," + unc_mode + "," + t_max;

							char file[500];
							sprintf(file, "%s/%s/seed%u/%s/fit.root", base_dir.c_str(), simuDir.c_str(), seed, fitDir.c_str());
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

						gDirectory = d_norm;

						h_de_rho->Write();
						h_de_a->Write();

						TGraph *g_data = new TGraph();
						g_data->SetPoint(0, stat.GetMean(0), stat.GetMeanUnc(0));
						g_data->SetPoint(1, stat.GetStdDev(0), stat.GetStdDevUnc(0));
						g_data->SetPoint(2, stat.GetMean(1), stat.GetMeanUnc(1));
						g_data->SetPoint(3, stat.GetStdDev(1), stat.GetStdDevUnc(1));
						g_data->Write("g_data");

						const int idx = g_rho_bias_vs_norm->GetN();
						const double de_norm = atof(norm.c_str());

						g_rho_bias_vs_norm->SetPoint(idx, de_norm, stat.GetMean(0));
						g_rho_bias_vs_norm->SetPointError(idx, 0., stat.GetMeanUnc(0));

						g_rho_rms_vs_norm->SetPoint(idx, de_norm, stat.GetStdDev(0));
						g_rho_rms_vs_norm->SetPointError(idx, 0., stat.GetStdDevUnc(0));

						g_a_bias_vs_norm->SetPoint(idx, de_norm, stat.GetMean(1));
						g_a_bias_vs_norm->SetPointError(idx, 0., stat.GetMeanUnc(1));

						g_a_rms_vs_norm->SetPoint(idx, de_norm, stat.GetStdDev(1));
						g_a_rms_vs_norm->SetPointError(idx, 0., stat.GetStdDevUnc(1));
					}

					gDirectory = d_t_max;
					g_rho_bias_vs_norm->Write("g_rho_bias_vs_norm");
					g_rho_rms_vs_norm->Write("g_rho_rms_vs_norm");
					g_a_bias_vs_norm->Write("g_a_bias_vs_norm");
					g_a_rms_vs_norm->Write("g_a_rms_vs_norm");
				}
			}
		}
	}
}

//----------------------------------------------------------------------------------------------------

int main()
{

	// open model file
	TFile *f_in_model = TFile::Open("build_models.root");

	// prepare output
	TFile *f_out = TFile::Open("process_fits.root", "recreate");

	// process fit study
	{
		// settings
		string base_dir = "simu-fit-study";

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

		gDirectory = f_out;
		ProcessOneGroup(base_dir, models, unc_modes, norms, binnings, seed_min, seed_max, t_maxs, f_in_model, f_out);
	}

	// process normalisation study
	{
		// settings
		string base_dir = "simu-norm-study";

		vector<string> models = {
			"exp1,rho=0.06",
			"exp1,rho=0.10",
			"exp1,rho=0.14",
			"exp3,rho=0.06",
			"exp3,rho=0.10",
			"exp3,rho=0.14",
		};

		vector<string> unc_modes = {
			"none",
		};

		vector<string> norms = {
			"+0.10",
			"+0.08",
			"+0.06",
			"+0.04",
			"+0.02",
			"0.00",
			"-0.02",
			"-0.04",
			"-0.06",
			"-0.08",
			"-0.10"
		};

		vector<string> binnings = {
			"ob-2-10-0.05",
		};

		unsigned int seed_min = 1, seed_max = 1;

		vector<string> t_maxs = {
			"0.15"
		};

		gDirectory = f_out;
		ProcessOneGroup(base_dir, models, unc_modes, norms, binnings, seed_min, seed_max, t_maxs, f_in_model, f_out);
	}

	// clean up
	delete f_in_model;
	delete f_out;

	return 0;
}
