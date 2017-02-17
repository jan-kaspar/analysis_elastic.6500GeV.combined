#include "include/Elegent/Constants.h"
#include "include/Elegent/CoulombInterference.h"

#include "../HadronicFitModel.h"

#include "TFile.h"
#include "TH1D.h"

#include "TRandom3.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

TH1D* BuildGeneratingHistogram(HadronicFitModel *m, TH1D *h_in, double t_min, double t_max)
{
	model = m;
	//model->Print();

	// create new empty histogram with the same binning as h_in
	TH1D *h_out = new TH1D(*h_in);
	for (int bi = 1; bi <= h_out->GetNbinsX(); ++bi)
	{
		h_out->SetBinContent(bi, 0.);
		h_out->SetBinError(bi, 0.);
	}

	// determine bin range
	int bi_min = h_in->GetXaxis()->FindBin(t_min);
	int bi_max = h_in->GetXaxis()->FindBin(t_max);

	//printf("%i, %i\n", bi_min, bi_max);

	// determine central values and uncertainties
	for (int bi = bi_min; bi <= bi_max; ++bi)
	{
		double le = h_out->GetBinLowEdge(bi);
		double he = le + h_out->GetBinWidth(bi);

		int n_div = 10;
		double w = (he - le)/ 10;

		double S = 0.;
		for (int si = 0; si < n_div; ++si)
		{
			double c = le + w * (0.5 + si);
			S += w * cnts->sig_fac * coulomb->Amp(-c).Rho2();
		}

		double y_cen = S / (he - le);

		double y_unc_rel = h_in->GetBinError(bi) / h_in->GetBinContent(bi);

		// TODO
		// special settings
		if (bi == 5)
			y_unc_rel = 0.05;
		if (bi == 6)
			y_unc_rel = 0.02;

		double y_unc = y_cen * y_unc_rel;

		h_out->SetBinContent(bi, y_cen);
		h_out->SetBinError(bi, y_unc);

		//printf("    t=%f, rel unc=%f\n", (he+le)/2., y_unc_rel);
	}


	return h_out;
}

//----------------------------------------------------------------------------------------------------

TH1D* SimulateErrors(TH1D *h_gen, double t_min, double t_max)
{
	// create new empty histogram with the same binning as h_in
	TH1D *h_out = new TH1D(*h_gen);

	// determine bin range
	int bi_min = h_out->GetXaxis()->FindBin(t_min);
	int bi_max = h_out->GetXaxis()->FindBin(t_max);

	for (int bi = bi_min; bi <= bi_max; ++bi)
	{
		double y_cen = h_out->GetBinContent(bi);
		double y_unc = h_out->GetBinError(bi);

		double y_err = gRandom->Gaus() * y_unc;

		h_out->SetBinContent(bi, y_cen + y_err);
	}

	return h_out;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// defaults
	string dataModelFile = "data_model/merged.root";
	string dataModelObject = "ob-1-30-0.05/merged/combined/h_dsdt";

	double t_min = 8.1E-4;
	double t_max = 0.2;

	int seed_min = 1;
	int seed_max = 10;

	// init Elegent
	Constants::Init(2*6500, cnts->mPP);
    cnts->Print();

	coulomb->ffType = coulomb->ffPuckett;
	coulomb->Print();

	// cross-section models
	struct CSModel
	{
		string label;
		HadronicFitModel *hfm;
		TH1D *h_gen;
	};

	vector<CSModel> csModels;

	HadronicFitModel *hfm = new HadronicFitModel();

	coulomb->mode = CoulombInterference::mKL;
	
	hfm->t1 = 0.2;
	hfm->t2 = 0.5;

	hfm->a = 1.84E9;
	hfm->b1 = 10.2;
	hfm->b2 = 0.;
	hfm->b3 = 0.;

	hfm->phaseMode = HadronicFitModel::pmConstant;

	hfm->p0 = M_PI/2. - atan(0.10);
	csModels.push_back({"exp,rho=0.10", new HadronicFitModel(*hfm), NULL});

	hfm->p0 = M_PI/2. - atan(0.14);
	csModels.push_back({"exp,rho=0.14", new HadronicFitModel(*hfm), NULL});

	hfm->b1 = 10.2;
	hfm->b2 = 4.4;
	hfm->b3 = 10.;

	hfm->p0 = M_PI/2. - atan(0.10);
	csModels.push_back({"non-exp,rho=0.10", new HadronicFitModel(*hfm), NULL});

	hfm->p0 = M_PI/2. - atan(0.14);
	csModels.push_back({"non-exp,rho=0.14", new HadronicFitModel(*hfm), NULL});

	// load data_model
	TFile *f_in = TFile::Open(dataModelFile.c_str());
	TH1D *h_in = (TH1D *) f_in->Get(dataModelObject.c_str());

	// prepare output
	// TODO
	//TFile *f_out = TFile::Open("simulation.root", "recreate");
	TFile *f_out = TFile::Open("simulation2.root", "recreate");

	// make histograms
	vector<TH1D *> csHistograms;

	for (auto &m : csModels)
	{
		printf("* %s\n", m.label.c_str());
		gDirectory = f_out->mkdir(m.label.c_str());

		// build generating histogram
		m.h_gen = BuildGeneratingHistogram(m.hfm, h_in, t_min, t_max);
		m.h_gen->SetLineColor(2);
		m.h_gen->Write("h_gen");

		// simulate errors
		for (int seed = seed_min; seed <= seed_max; ++seed)
		{
			gRandom->SetSeed(seed);
			TH1D *h_sim = SimulateErrors(m.h_gen, t_min, t_max);

			char buf[100];
			sprintf(buf, "h_sim_seed%i", seed);
			h_sim->SetLineColor(4);
			h_sim->Write(buf);
		}
	}

	// clean up
	delete f_out;
	delete f_in;

	return 0;
}
