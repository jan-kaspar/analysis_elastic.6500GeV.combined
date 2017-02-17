#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"

#include "TH1D.h"
#include "TGraph.h"

#include <vector>

using namespace std;

namespace Fit
{

unsigned int B_degree;

TH1D *h_fit;

int bin_fit_min, bin_fit_max;

HadronicFitModel *hfm;

TGraph *interpolatedPsiRe, *interpolatedPsiIm;

//----------------------------------------------------------------------------------------------------

struct Results
{
	double A = 0., A_e = 0.;
	double B = 0., B_e = 0.;

	double p0 = 0., p0_e = 0.;
	double rho = 0., rho_e = 0.;
};

//----------------------------------------------------------------------------------------------------

template <class T>
void SetModelParameters(const T &par, bool all_phase_parameters = true)
{
	hfm->a = par[0] * 1E8;

	if (B_degree > 0) hfm->b1 = par[1];
	if (B_degree > 1) hfm->b2 = par[2];
	if (B_degree > 2) hfm->b3 = par[3];
	if (B_degree > 3) hfm->b4 = par[4];
	if (B_degree > 4) hfm->b5 = par[5];
	if (B_degree > 5) hfm->b6 = par[6];
	if (B_degree > 6) hfm->b7 = par[7];
	if (B_degree > 7) hfm->b8 = par[8];
	if (B_degree > 8) hfm->b9 = par[9];

	hfm->p0 = par[B_degree + 1];
	if (all_phase_parameters)
	{
		hfm->p_A = par[B_degree + 2];
		hfm->p_ka = par[B_degree + 3];
		hfm->p_tm = par[B_degree + 4];
	}
}

//----------------------------------------------------------------------------------------------------

void InterpolatePsi()
{
	printf(">> Fit::InterpolatePsi\n");

	interpolatedPsiRe->Set(0);
	interpolatedPsiIm->Set(0);

	for (int bi = bin_fit_min; bi <= bin_fit_max; ++bi)
	{
		double x = h_fit->GetBinCenter(bi);

		TComplex Psi = coulomb->Psi_KL(-x);

		int idx = interpolatedPsiRe->GetN();
		interpolatedPsiRe->SetPoint(idx, x, Psi.Re());
		interpolatedPsiIm->SetPoint(idx, x, Psi.Im());
	}
}

//----------------------------------------------------------------------------------------------------

double F_fit(double mt)
{
	// amplitude components
	TComplex F_C = coulomb->Amp_pure(-mt);
	TComplex F_H = hfm->Amp(-mt);
	
	TComplex Psi(interpolatedPsiRe->Eval(mt), interpolatedPsiIm->Eval(mt));
	TComplex F_T = F_C + F_H * TComplex::Exp(i*Psi);

	return cnts->sig_fac * F_T.Rho2();
}

//----------------------------------------------------------------------------------------------------

class S2_FCN : public ROOT::Minuit2::FCNBase
{
	public:
		S2_FCN() {}

  		double operator() (const std::vector<double> &) const;
  		double Up() const { return 1.; }
};

//----------------------------------------------------------------------------------------------------

double S2_FCN::operator() (const std::vector<double> &par) const
{
	SetModelParameters(par, false);

	double S2 = 0.;
	for (int bi = bin_fit_min; bi <= bin_fit_max; ++bi)
	{
		double x = h_fit->GetBinCenter(bi);

		double y = h_fit->GetBinContent(bi);
		double y_unc = h_fit->GetBinError(bi);

		double f = F_fit(x);

		double rd = (y - f) / y_unc;

		S2 += rd * rd;
	}

	return S2;
}

//----------------------------------------------------------------------------------------------------

void MakePlots()
{
	printf(">> MakePlots\n");

	TGraph *g_dsdt_C = new TGraph();
	TGraph *g_dsdt_H = new TGraph();
	TGraph *g_dsdt_CH = new TGraph();

	for (double mt = 0.5E-4; mt <= 1.;)
	{
		double dsdt_C = cnts->sig_fac * coulomb->Amp_pure(-mt).Rho2();
		double dsdt_H = cnts->sig_fac * model->Amp(-mt).Rho2();
		double dsdt_CH = cnts->sig_fac * coulomb->Amp(-mt).Rho2();

		int idx = g_dsdt_C->GetN();
		g_dsdt_C->SetPoint(idx, mt, dsdt_C);
		g_dsdt_H->SetPoint(idx, mt, dsdt_H);
		g_dsdt_CH->SetPoint(idx, mt, dsdt_CH);

		double dmt = 0.01;
		if (mt < 0.2)
			dmt = 0.001;
		if (mt < 0.004)
			dmt = 0.0001;
		if (mt < 0.001)
			dmt = 0.00005;

		mt += dmt;
	}

	g_dsdt_C->Write("g_dsdt_C");
	g_dsdt_H->Write("g_dsdt_H");
	g_dsdt_CH->Write("g_dsdt_CH");
}

//----------------------------------------------------------------------------------------------------

unsigned int RunFit(Results &results)
{
	printf(">> RunFit\n");

	// defaults
	double a_def = 1.84E9;
	double b1_def = 10.2;
	double p0_def = M_PI/2. - atan(0.12);

	unsigned int N_it = 3;

	// initialise hadronic model
	hfm = new HadronicFitModel();
	hfm->t1 = 0.2;
	hfm->t2 = 0.5;
	hfm->phaseMode = HadronicFitModel::pmConstant;

	model = hfm;

	// initialize storage for interpolated Psi function
	interpolatedPsiRe = new TGraph(); interpolatedPsiRe->SetName("interpolatedPsiRe");
	interpolatedPsiIm = new TGraph(); interpolatedPsiIm->SetName("interpolatedPsiIm");

	// select CNI formula
	coulomb->mode = CoulombInterference::mKL;

	// initialise fitter
	TFitterMinuit *minuit = new TFitterMinuit();
	S2_FCN fcn;
	minuit->SetMinuitFCN(&fcn);
	minuit->SetPrintLevel(10);
	minuit->CreateMinimizer();

	// initial point - modulus
	char buf[200];
	minuit->SetParameter(0, "a", a_def / 1E8, 0.7, 0., 0.);	// without the factor 1E8

	for (unsigned int i = 1; i <= B_degree; i++)
	{
		sprintf(buf, "b%i", i);

		double val = 0., unc = 2.;
		double lim_low=0., lim_high=0.;
		if (i == 1) val = b1_def;
		if (i == 2) val = 0.;
		if (i == 3) val = 0.;
		if (i == 4) val = 0.;

		minuit->SetParameter(i, buf, val, unc, lim_low, lim_high);
	}
	
	// initial point - phase
	double p0_lim_min = 0.;
	double p0_lim_max = 0.;
	minuit->SetParameter(B_degree+1, "p0", p0_def, 0.01, p0_lim_min, p0_lim_max);

	// run minimisation
	for (unsigned int it = 0; it < N_it; it++)
	{
		printf("* iteration %i\n", it);

		// transfer parameters to hfm
		unsigned int n_par = minuit->GetNumberTotalParameters();
		vector<double> pars;
		for (unsigned int i = 0; i < n_par; i++)
			pars.push_back(minuit->GetParameter(i));

		SetModelParameters(pars);

		// interpolate Psi
		InterpolatePsi();

		// start minimisation
		minuit->Minimize();
	}

	// save results
	results.p0 = minuit->GetParameter(B_degree+1);
	results.p0_e = minuit->GetParError(B_degree+1);

	results.rho = cos(results.p0) / sin(results.p0);
	results.rho_e = fabs(1. / sin(results.p0) / sin(results.p0)) * results.p0_e;

	results.A = cnts->sig_fac * minuit->GetParameter(0) * 1E8 * minuit->GetParameter(0) * 1E8;
	results.A_e = 2. * cnts->sig_fac * minuit->GetParameter(0) * 1E8 * minuit->GetParError(0) * 1E8;

	results.B = 2.*minuit->GetParameter(1);
	results.B_e = 2.*minuit->GetParError(1);

	// transfer parameters to hfm
	unsigned int n_par = minuit->GetNumberTotalParameters();
	vector<double> pars;
	for (unsigned int i = 0; i < n_par; i++)
		pars.push_back(minuit->GetParameter(i));

	SetModelParameters(pars);

	// save plots
	MakePlots();

	// clean up
	delete hfm;

	delete interpolatedPsiRe;
	delete interpolatedPsiIm;

	// TODO
	//delete minuit;

	return 0;
}

};
