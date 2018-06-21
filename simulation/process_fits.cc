#include "TH1D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "../stat.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

void ProcessOneGroup(const string &base_dir,
		const vector<string> &models, const vector<string> &unc_modes, const vector<string> &norms,
		const vector<string> &binnings, unsigned int seed_min, unsigned int seed_max,
		const vector<string> &fitMethods,
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
		const double a_0 = g_model_data->GetY()[0];
		const double rho_0 = g_model_data->GetY()[1];
		const double si_tot_0 = g_model_data->GetY()[2];

		for (const auto &unc_mode : unc_modes)
		{
			TDirectory *d_unc_mode = d_model->mkdir(unc_mode.c_str());

			for (const auto &binning : binnings)
			{
				TDirectory *d_binning = d_unc_mode->mkdir(binning.c_str());

				for (const auto &fitMethod : fitMethods)
				{
					TDirectory *d_fit_meth = d_binning->mkdir(fitMethod.c_str());

					const string &expN = model.substr(0, 4);

					string t_max = "t_max_unknown";
					if (expN == "exp1") t_max = "0.07";
					if (expN == "exp3") t_max = "0.15";

					{
						TDirectory *d_t_max = d_fit_meth->mkdir(("tmax" + t_max).c_str());

						TGraphErrors *g_rho_bias_vs_norm = new TGraphErrors();
						TGraphErrors *g_rho_rms_vs_norm = new TGraphErrors();
						TGraphErrors *g_a_bias_vs_norm = new TGraphErrors();
						TGraphErrors *g_a_rms_vs_norm = new TGraphErrors();

						for (const auto &norm : norms)
						{
							TDirectory *d_norm = d_t_max->mkdir(("norm"+norm).c_str());
							gDirectory = d_norm;

							TH1D *h_de_rho = new TH1D("h_de_rho", "", 20, -0.05, +0.05);
							TH1D *h_de_a = new TH1D("h_de_a", "", 20, -150., +150.);
							TH1D *h_de_si_tot = new TH1D("h_de_si_tot", "", 20, -10., +10.);

							TGraph *g_de_rho_vs_de_a = new TGraph(); g_de_rho_vs_de_a->SetName("g_de_rho_vs_de_a");

							Stat stat(3);

							for (unsigned int seed = seed_min; seed <= seed_max; ++seed)
							{
								// get reco data
								string simuDir = model+"/"+unc_mode+"/norm"+norm+"/"+binning;
								string fitDir = model.substr(0, 4) + "," + unc_mode + "," + t_max;

								char file[500];
								sprintf(file, "%s/%s/seed%u/%s/%s/fit.root", base_dir.c_str(), simuDir.c_str(), seed, fitMethod.c_str(), fitDir.c_str());
								TFile *f = TFile::Open(file);
								TGraph *g_fit_data = (TGraph *) f->Get("g_fit_data");

								const double *v = g_fit_data->GetY();
								const double rho = v[4];
								const double a = v[6];
								const double si_tot = v[14];

								// calculate simu-reco deltas
								const double de_rho = rho - rho_0;
								const double de_a = a - a_0;
								const double de_si_tot = si_tot - si_tot_0;

								// skip outliers
								//if (abs(de_rho) > 0.05 || abs(de_a > 150.))
								//	continue;

								h_de_rho->Fill(de_rho);
								h_de_a->Fill(de_a);
								h_de_si_tot->Fill(de_si_tot);

								int idx = g_de_rho_vs_de_a->GetN();
								g_de_rho_vs_de_a->SetPoint(idx, de_a, de_rho);

								stat.Fill(de_rho, de_a, de_si_tot);

								delete f;
							}

							gDirectory = d_norm;

							h_de_rho->Write();
							h_de_a->Write();
							h_de_si_tot->Write();
							g_de_rho_vs_de_a->Write();

							TGraph *g_data = new TGraph();
							g_data->SetPoint(0, stat.GetMean(0), stat.GetMeanUnc(0));
							g_data->SetPoint(1, stat.GetStdDev(0), stat.GetStdDevUnc(0));
							g_data->SetPoint(2, stat.GetMean(1), stat.GetMeanUnc(1));
							g_data->SetPoint(3, stat.GetStdDev(1), stat.GetStdDevUnc(1));
							g_data->SetPoint(4, stat.GetMean(2), stat.GetMeanUnc(2));
							g_data->SetPoint(5, stat.GetStdDev(2), stat.GetStdDevUnc(2));
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
			//"exp1,rho=0.14",
			//"exp2,rho=0.10",
			//"exp2,rho=0.14",
			"exp3,rho=0.10",
			//"exp3,rho=0.14",
		};

		vector<string> unc_modes = {
			//"stat",
			//"stat+syst",
			"stat+syst+norm",
			//"norm"
		};

		vector<string> norms = {
			"0.00",
		};

		vector<string> binnings = {
			"ob-2-10-0.05",
			//"ob-3-5-0.05",
		};

		unsigned int seed_min = 1, seed_max = 500;

		vector<string> fitMethods = {
			"coulomb_analysis_1",
			"coulomb_analysis_1_weighted",
			"coulomb_analysis_1_weighted_norm_lim",
			"coulomb_analysis_1_weighted_norm_unlim",
		};

		gDirectory = f_out;
		ProcessOneGroup(base_dir, models, unc_modes, norms, binnings, seed_min, seed_max, fitMethods, f_in_model, f_out);
	}

#if 0
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
#endif

	// clean up
	delete f_in_model;
	delete f_out;

	return 0;
}
