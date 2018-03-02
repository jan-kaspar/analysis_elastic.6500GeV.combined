#include "NumericalIntegration.h"

#include "include/Elegent/Constants.h"

#include "TSpline.h"
#include "TFile.h"
#include "TGraph.h"

#include <string>

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

// input for FB transformation
TSpline *s_dsdt_H, *s_phase_H;

//----------------------------------------------------------------------------------------------------

TComplex Amp_H(double mt)
{
	// mt = -t

	double mod = sqrt(s_dsdt_H->Eval(mt) / cnts->sig_fac);
	double ph = s_phase_H->Eval(mt);

	//printf("   dsdt=%E | mod=%E, ph=%E\n", s_dsdt_H->Eval(mt), mod, ph);

	//if (mt < 0.001)
	//.	printf("==> %E\n", mod);
	
	return mod * TComplex(cos(ph), sin(ph));
}

//----------------------------------------------------------------------------------------------------

TComplex FBInteg(double mt, double *par, const void *)
{
	double b = par[0];

	TComplex v = TMath::BesselJ0(b * sqrt(mt)) * Amp_H(mt);

	//printf("mt=%E, b=%E --> amp=%E + i%E --> integ=%E + i%E\n", mt, b, Amp_H(mt).Re(), Amp_H(mt).Im(), v.Re(), v.Im());

	return v;
}

//----------------------------------------------------------------------------------------------------

unsigned long integ_ws_size, integ_ws_mb_size;
gsl_integration_workspace *integ_ws, *integ_ws_mb;

//----------------------------------------------------------------------------------------------------

/// b in fm
TComplex Prf(double b)
{
	double par[] = { b / cnts->hbarc };	// b in GeV^-1
	double mt_max = 3.0;
	double rel_err = 1E-2;

	return ComplexIntegrate(FBInteg, par, NULL, 0., mt_max, 0., rel_err, integ_ws_size, integ_ws, "Prf")
		/ 4. / cnts->p_cms / cnts->sqrt_s;
}

//----------------------------------------------------------------------------------------------------

enum { dtEl, dtTot, dtInel } distType;

/// b in fm
double MeanBInteg(double b, double *par, const void *)
{
	double power = par[0];

	double A_sq = 0;
	switch (distType)
	{
		case dtEl:
			A_sq = Prf(b).Rho2();
			break;

		case dtTot:
			A_sq = Prf(b).Im();
			break;

		case dtInel:
			A_sq = Prf(b).Im() - Prf(b).Rho2();
			break;
	}

	return b * pow(b, power) * A_sq;
}

//----------------------------------------------------------------------------------------------------

void DoTransformation()
{
	TGraph *g_A_re = new TGraph(); g_A_re->SetName("g_A_re"); g_A_re->SetTitle(";b   (fm);Re A(b)"); g_A_re->SetLineColor(2);
	TGraph *g_A_im = new TGraph(); g_A_im->SetName("g_A_im"); g_A_im->SetTitle(";b   (fm);Im A(b)"); g_A_im->SetLineColor(4);
	TGraph *g_A_mod2 = new TGraph(); g_A_mod2->SetName("g_A_mod2"); g_A_mod2->SetTitle(";b   (fm);|A(b)|^{2}"); g_A_mod2->SetLineColor(1);
	TGraph *g_A_mod2_b = new TGraph(); g_A_mod2_b->SetName("g_A_mod2_b"); g_A_mod2_b->SetTitle(";b   (fm);b |A(b)|^{2}"); g_A_mod2_b->SetLineColor(8);
	
	// sample graphs
	// b in fm
	for (double b = 0.; b <= 5.; b += 0.05)
	{
		TComplex A = Prf(b);

		//printf("%2.3f | %E, %E\n", b, A.Re(), A.Im());

		int idx = g_A_re->GetN();

		g_A_re->SetPoint(idx, b, A.Re());
		g_A_im->SetPoint(idx, b, A.Im());
		g_A_mod2->SetPoint(idx, b, A.Rho2());
		g_A_mod2_b->SetPoint(idx, b, A.Rho2() * b);
	}

	g_A_re->Write();
	g_A_im->Write();
	g_A_mod2->Write();
	g_A_mod2_b->Write();

	// calculate moments of b
	double par[] = { 0. };
	double b_max = 20.;
	double rel_err = 0.001;

	distType = dtEl; // elastic
	printf("* elastic\n");

	par[0] = 2.;
	double mbi2 = RealIntegrate(MeanBInteg, par, NULL, 0., b_max, 0., rel_err, integ_ws_mb_size, integ_ws_mb, "mb2");

	par[0] = 0.;
	double mbi0 = RealIntegrate(MeanBInteg, par, NULL, 0., b_max, 0., rel_err, integ_ws_mb_size, integ_ws_mb, "mb0");

	double b_el_ms = mbi2 / mbi0;
	double b_el_rms = sqrt(b_el_ms);

	// the factor 10 = sq_hbarc / hbarc / hbarc
	double si_el = 8.*M_PI * mbi0 * 10.;	// in mb

	printf("si_el: %.3f mb\n", si_el);
	printf("RMS of b_el: %.3f fm\n", b_el_rms);
	

	distType = dtInel; // inel
	printf("\n* inelastic\n");

	par[0] = 2.; mbi2 = RealIntegrate(MeanBInteg, par, NULL, 0., b_max, 0., rel_err, integ_ws_mb_size, integ_ws_mb, "mb2");
	par[0] = 0.; mbi0 = RealIntegrate(MeanBInteg, par, NULL, 0., b_max, 0., rel_err, integ_ws_mb_size, integ_ws_mb, "mb0");
	printf("mbi2=%E, mbi0=%E\n", mbi2, mbi0);

	double b_inel_ms = mbi2 / mbi0;
	double b_inel_rms = sqrt(b_inel_ms);
	printf("RMS of b_inel: %.3f fm\n", b_inel_rms);
	

	distType = dtTot; // total
	printf("\n* total\n");
	
	par[0] = 2.; mbi2 = RealIntegrate(MeanBInteg, par, NULL, 0., b_max, 0., rel_err, integ_ws_mb_size, integ_ws_mb, "mb2");
	par[0] = 0.; mbi0 = RealIntegrate(MeanBInteg, par, NULL, 0., b_max, 0., rel_err, integ_ws_mb_size, integ_ws_mb, "mb0");
	printf("mbi2=%E, mbi0=%E\n", mbi2, mbi0);

	double b_tot_ms = mbi2 / mbi0;
	double b_tot_rms = sqrt(b_tot_ms);
	printf("RMS of b_tot: %.3f fm\n", b_tot_rms);


	TGraph *g_b_mom = new TGraph();
	g_b_mom->SetName("g_b_mom");
	g_b_mom->SetPoint(0, 0, si_el);
	g_b_mom->SetPoint(1, 1, b_el_rms);
	g_b_mom->SetPoint(2, 2, b_tot_rms);
	g_b_mom->SetPoint(3, 3, b_inel_rms);
	g_b_mom->Write();
}

//----------------------------------------------------------------------------------------------------
	
TSpline* BuildSpline(TGraph *g)
{
	TSpline3 *s = new TSpline3("", g->GetX(), g->GetY(), g->GetN());
	s->SetName(g->GetName());
	return s;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: impactParameterDistributions <input file> <output file>\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// parse command line
	if (argc != 3)
	{
		printf("ERROR: wrong number of parameters.\n");
		PrintUsage();
		return 1;
	}

	string fn_in = argv[1];
	string fn_out = argv[2];

	// init constants
	Constants::Init(2*6500., cnts->mPP);
    cnts->Print();

	// get input
	TFile *f_in = TFile::Open(fn_in.c_str());
	if (f_in == NULL)
	{
		printf("ERROR: can't open file `%s'.\n", fn_in.c_str());
		return 2;
	}

	TGraph *g_dsdt_H = (TGraph *) f_in->Get("g_fit_H");
	TGraph *g_phase_H = (TGraph *) f_in->Get("g_Phase_H");

	if (!g_dsdt_H || !g_phase_H)
	{
		printf("ERROR: can't load input: %p, %p\n", g_dsdt_H, g_phase_H);
		return 3;
	}

	s_dsdt_H = BuildSpline(g_dsdt_H);
	s_phase_H = BuildSpline(g_phase_H);

	// prepare output
	TFile *f_out = TFile::Open(fn_out.c_str(), "recreate");

	g_dsdt_H->Write("g_dsdt_H");
	g_phase_H->Write("g_Phase_H");

	// do transformation
	integ_ws_size = 100;
	integ_ws = gsl_integration_workspace_alloc(integ_ws_size);

	integ_ws_mb_size = 100;
	integ_ws_mb = gsl_integration_workspace_alloc(integ_ws_mb_size);

	DoTransformation();

	delete f_out;
	delete f_in;

	return 0;
}
