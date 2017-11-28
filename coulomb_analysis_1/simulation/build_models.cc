#include "include/Elegent/Constants.h"
#include "include/Elegent/CoulombInterference.h"

#include "../../HadronicFitModel.h"

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TSpline.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

TSpline* BuildSpline(TGraph *g)
{
	TSpline3 *s = new TSpline3("", g->GetX(), g->GetY(), g->GetN());
	s->SetName(g->GetName());
	return s;
}

//----------------------------------------------------------------------------------------------------

void SampleModel(Model *m, const string &label, const vector<TH1D*> &b_hists)
{
	printf("* %s\n", label.c_str());

	model = m;

	TDirectory *d_top = gDirectory;

	TDirectory *d_model = d_top->mkdir(label.c_str());
	gDirectory = d_model;

	// save model data
	TGraph *g_data = new TGraph();
	TComplex amp0 = model->Amp(0.);
	g_data->SetPoint(0, 0, cnts->sig_fac * amp0.Rho2());
	g_data->SetPoint(1, 1, amp0.Re() / amp0.Im());
	g_data->Write("g_data");

	// make graphs
	TGraph *g_dsdt_H = new TGraph();
	TGraph *g_dsdt_CH = new TGraph();

	for (double mt = 1E-4; mt < 1.1; )
	{
		coulomb->mode = CoulombInterference::mPH;
		TComplex F_H = coulomb->Amp(-mt);

		coulomb->mode = CoulombInterference::mKL;
		TComplex F_CH = coulomb->Amp(-mt);

		int idx = g_dsdt_H->GetN();
		g_dsdt_H->SetPoint(idx, mt, cnts->sig_fac * F_H.Rho2());
		g_dsdt_CH->SetPoint(idx, mt, cnts->sig_fac * F_CH.Rho2());

		double dmt = 1E-2;
		if (mt < 0.2) dmt = 5E-3;
		if (mt < 0.05) dmt = 1E-3;
		if (mt < 0.004) dmt = 1E-4;
		if (mt < 0.001) dmt = 1E-5;
		mt += dmt;
	}

	g_dsdt_H->Write("g_dsdt_H");
	g_dsdt_CH->Write("g_dsdt_CH");

	// make histograms
	TSpline *s = BuildSpline(g_dsdt_CH);

	for (const auto &h : b_hists)
	{
		gDirectory = d_model->mkdir(h->GetName());

		TH1D *h_dsdt = new TH1D(*h);

		for (int bi = 1; bi <= h_dsdt->GetNbinsX(); ++bi)
		{
			// integrate function over bin
			const double l = h_dsdt->GetBinLowEdge(bi);
			const double w = h_dsdt->GetBinWidth(bi);

			const unsigned int n_div = 100;

			double S = 0.;
			for (unsigned int id = 0; id < n_div; id++)
			{
				const double mt = l + (0.5 + id) * w / n_div;
				S += s->Eval(mt);
			}

			S /= n_div;

			h_dsdt->SetBinContent(bi, S);
			h_dsdt->SetBinError(bi, 0.);
		}

		h_dsdt->Write("h_dsdt_CH");
	}


	gDirectory = d_top;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// binnings
	vector<string> binnings;
	binnings.push_back("ob-1-20-0.05");
	binnings.push_back("ob-2-10-0.05");
	binnings.push_back("ob-3-5-0.05");

	// load binnings
	TFile *f_in = TFile::Open("/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/DS-merged/merged.root");
	vector<TH1D *> b_hists;
	for (const auto &b : binnings)
	{
		TH1D *h = (TH1D *) f_in->Get((b + "/merged/combined/h_dsdt").c_str());
		h->SetName(b.c_str());
		b_hists.push_back(h);
	}

	// init Elegent
	Constants::Init(2*6500, cnts->mPP);
    cnts->Print();

	coulomb->ffType = coulomb->ffPuckett;
	coulomb->mode = CoulombInterference::mKL;

	coulomb->Print();

	// prepare output
	TFile *f_out = TFile::Open("build_models.root", "recreate");

	// init hadronic model
	HadronicFitModel *hfm = new HadronicFitModel();

	hfm->hts = sqrt(1.06);	// takes into account recent normalisation change

	hfm->t1 = 0.2;
	hfm->t2 = 0.5;

	hfm->phaseMode = HadronicFitModel::pmConstant;

	// sample exp1 models
	hfm->a = 19.5E8;
	hfm->b1 = 10.23;
	hfm->b2 = 0.;
	hfm->b3 = 0.;

	hfm->p0 = M_PI/2. - atan(0.06);
	SampleModel(hfm, "exp1,rho=0.06", b_hists);

	hfm->p0 = M_PI/2. - atan(0.10);
	SampleModel(hfm, "exp1,rho=0.10", b_hists);

	hfm->p0 = M_PI/2. - atan(0.14);
	SampleModel(hfm, "exp1,rho=0.14", b_hists);

	// sample exp2 models
	hfm->a = 19.5E8;
	hfm->b1 = 10.51;
	hfm->b2 = 1.95;
	hfm->b3 = 0.;

	hfm->p0 = M_PI/2. - atan(0.06);
	SampleModel(hfm, "exp2,rho=0.06", b_hists);

	hfm->p0 = M_PI/2. - atan(0.10);
	SampleModel(hfm, "exp2,rho=0.10", b_hists);

	hfm->p0 = M_PI/2. - atan(0.14);
	SampleModel(hfm, "exp2,rho=0.14", b_hists);

	// sample exp3 models
	hfm->a = 19.5E8;
	hfm->b1 = 10.63;
	hfm->b2 = 3.96;
	hfm->b3 = 9.59;

	hfm->p0 = M_PI/2. - atan(0.06);
	SampleModel(hfm, "exp3,rho=0.06", b_hists);

	hfm->p0 = M_PI/2. - atan(0.10);
	SampleModel(hfm, "exp3,rho=0.10", b_hists);

	hfm->p0 = M_PI/2. - atan(0.14);
	SampleModel(hfm, "exp3,rho=0.14", b_hists);

	// clean up
	delete f_out;

	return 0;
}
