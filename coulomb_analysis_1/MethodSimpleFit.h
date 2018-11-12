#include "TMatrixDSymEigen.h"
#include "TDirectory.h"

namespace MethodSimpleFit
{

bool useNormalisationFitParameter;
bool useNormalisationChiSqTerm;
bool useNormalisationLimit;
bool useNormalisationFromA;
double A_p_value_fix;

bool useEtaFixed;
double eta_value_fix;

bool useAFixed;
double a_value_fix;

bool useB1Fixed;
double b1_value_fix;

bool useB2Fixed;
double b2_value_fix;

bool useB3Fixed;
double b3_value_fix;

bool useRhoFixed;
double rho_value_fix;

bool useInterpolatedPsi;

double t_min_fit, t_max_fit;

TF1 *f_fit;

unsigned int N_ii;

bool debug_fit;
TVectorD eig_values;
TMatrixD E;
TGraphErrors *g_de;
TGraphErrors *g_cont;
TGraph *g_largest_cont_idx;

//----------------------------------------------------------------------------------------------------

double EtaFromA(double a)
{
	return A_p_value_fix / (cnts->sig_fac * a * 1E8 * a * 1E8);
}

//----------------------------------------------------------------------------------------------------

template <class T>
double F_fit(double mt, const T &par)
{
	// transfer parameters to FitModel
	SetModelParameters(par, false);

	double norm_corr = 1.;
	if (useNormalisationFitParameter)
		norm_corr = par[par_off_norm];
	if (useNormalisationFromA)
		norm_corr = EtaFromA(par[par_off_a]);

	// amplitude components
	TComplex F_C = coulomb->Amp_pure(-mt);
	TComplex F_H = hfm->Amp(-mt);
	TComplex F_T = 0.;

	// amplitude choice
	if (coulomb->mode == coulomb->mPC)
		F_T = F_C;

	if (coulomb->mode == coulomb->mPH)
		F_T = F_H;

	if (coulomb->mode == coulomb->mKL)
	{
		const TComplex Psi = (useInterpolatedPsi) ? TComplex(interpolatedPsiRe->Eval(mt), interpolatedPsiIm->Eval(mt)) : - coulomb->Phi_SWY(-mt);
		F_T = F_C + F_H * TComplex::Exp(i*Psi);
	}

	if (coulomb->mode == coulomb->mSWY)
	{
		const TComplex Psi = - coulomb->Phi_SWY(-mt);
		F_T = F_C + F_H * TComplex::Exp(i*Psi);
	}

	/*
	printf("mt=%.1E | FC: re=%+.1E, im=%+.1E, FH: re=%+.1E, im=%+.1E, Psi: re=%+.1E, im=%+.1E, FT: re=%+.1E, im=%+.1E | ds/dt = %.1E\n",
		mt, F_C.Re(), F_C.Im(), F_H.Re(), F_H.Im(), Psi.Re(), Psi.Im(), F_T.Re(), F_T.Im(), cnts->sig_fac * F_T.Rho2());
	*/

	return norm_corr * cnts->sig_fac * F_T.Rho2();
}

//----------------------------------------------------------------------------------------------------

double f_fit_imp(double x[], double par[])
{
	return F_fit(x[0], par);
}

//----------------------------------------------------------------------------------------------------

void InterpolatePsi()
{
	printf(">> MethodSimpleFit::InterpolatePsi\n");

	interpolatedPsiRe->Set(0);
	interpolatedPsiIm->Set(0);

	// make list of points
	vector<double> a_mt;
	for (unsigned int i = 0; i < data_coll.size(); ++i)
		a_mt.push_back(data_coll[i].x_repr);

	sort(a_mt.begin(), a_mt.end());

	for (double mt = a_mt.back() + 0.05; mt < 1.1; mt += 0.05)
		a_mt.push_back(mt);

	for (unsigned int i = 0; i < a_mt.size(); i++)
	{
		double mt = a_mt[i];

		TComplex Psi = coulomb->Psi_KL(-mt);

		interpolatedPsiRe->SetPoint(i, mt, Psi.Re());
		interpolatedPsiIm->SetPoint(i, mt, Psi.Im());
		//printf("\t%E| %E, %E\n", mt, Psi.Re(), Psi.Im());
	}
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
	// propagate parameters to fit function
	for (unsigned int i = 0; i < par.size(); i++)
		f_fit->SetParameter(i, par[i]);

	// make vector of differences
	int dim = data_coll_unc.GetNrows();
	TVectorD diff(dim);
	for (int i = 0; i < dim; i++)
	{
		if (data_coll[i].x_left < t_min_fit || data_coll[i].x_right > t_max_fit)
			diff(i) = 0.;
		else
			diff(i) = data_coll[i].y - f_fit->Eval(data_coll[i].x_repr);
	}

	// calculate diff^T * cov_matrix^-1 * diff
	double S2 = 0.;
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			S2 += diff(i) * data_coll_unc_inv(i, j) * diff(j);

	//printf("\tS2 = %E\n", S2);

	// debugging code
	if (debug_fit)
	{
		// transform diff vector into eigen coordinates
		TVectorD de(E * diff);

		// make graph, one point per bin - central value de(i), sigma: corresponding uncertainty 1./sqrt( eig_values(i) )
		g_de = new TGraphErrors();
		g_de->SetName("g_de");
		g_de->SetTitle(";eigenvector index");
		g_de->SetMarkerStyle(20);

		g_cont = new TGraphErrors();
		g_cont->SetName("g_cont");
		g_cont->SetTitle(";eigenvector index;(delta / sigma)^{2}");
		g_cont->SetMarkerStyle(20);

		double S2_alt = 0.;
		vector<pair<unsigned int, double>> v_idx_cont;
		for (int i = 0; i < dim; i++)
		{
			double unc2 = eig_values(i);
			double unc = (unc2 > 0.) ? 1./sqrt(unc2) : 0.;

			g_de->SetPoint(i, i, de(i));
			g_de->SetPointError(i, 0, unc);

			double cont2 = de(i) * de(i) * unc2;
			//double cont = sqrt(cont2);

			v_idx_cont.emplace_back(i, cont2);

			g_cont->SetPoint(i, i, cont2);

			S2_alt += cont2;
		}

		sort(v_idx_cont.begin(), v_idx_cont.end(), [](const pair<unsigned int, double> &left, const pair<unsigned int, double> &right) {
			return left.second > right.second;
		});

		g_largest_cont_idx = new TGraph();
		g_largest_cont_idx->SetName("g_largest_cont_idx");
		g_largest_cont_idx->SetTitle(";eigenvector index;contribution to #chi^{2}");

		for (const auto & p : v_idx_cont)
		{
			int idx = g_largest_cont_idx->GetN();
			g_largest_cont_idx->SetPoint(idx, p.first, p.second);
		}

		/*
		unsigned int n_b = par.size() - 1;
		int ndf = dim - (1 + n_b);
		double prob = TMath::Prob(S2, ndf);
		double sigma_eq = sqrt(2.) * TMath::ErfcInverse(prob);

		printf("\t\tS2 = %.3E, S2_alt = %.3E; ndf = %i, S2 / ndf = %.3f; Prob(S2, ndf) = %.2E (%.2f sigmas)\n",
			S2, S2_alt, ndf, S2 / ndf, prob, sigma_eq);
		*/
	}

	// additional S2 contribution from normalisation
	if (useNormalisationChiSqTerm)
	{
		double norm_corr = 1.;
		if (useNormalisationFitParameter)
			norm_corr = par[par_off_norm];
		if (useNormalisationFromA)
			norm_corr = EtaFromA(par[par_off_a]);

		const double de_norm = (norm_corr - 1.) / 0.055;
		S2 += de_norm * de_norm;
	}

	return S2;
}

//----------------------------------------------------------------------------------------------------

void SavePlotsBeforeFit(S2_FCN &fcn, TFitterMinuit *minuit)
{
	TDirectory *topDirectory = gDirectory;

	// eigen decomposition of inverse covariance matrix
	TMatrixDSymEigen eig_decomp(data_coll_unc_inv);

	int dim = data_coll_unc_inv.GetNrows();
	eig_values.ResizeTo(dim);
	eig_values = eig_decomp.GetEigenValues();

	E.ResizeTo(dim, dim);
	E = eig_decomp.GetEigenVectors();	// in columns
	E = E.Transpose(E);					// in rows

	// save eigenvectors as graphs
	gDirectory = topDirectory->mkdir("eigenvectors");
	for (int i = 0; i < dim; i++)
	{
		char buf[50];
		sprintf(buf, "eigenvector %i", i);
		TGraph *g = new TGraph();
		g->SetName(buf);
		g->SetTitle(";bin t_repr");
		g->SetMarkerStyle(20);

		for (int j = 0; j < dim; j++)
		{
			const double x = data_coll[j].x_repr;
			const double y = E(i, j);
			g->SetPoint(j, x, y);
		}

		g->Write();
	}

	gDirectory = topDirectory;

	// get parameters
	int n_par = minuit->GetNumberTotalParameters();
	vector<double> par(n_par);
	for (int i = 0; i < n_par; i++)
		par[i] = minuit->GetParameter(i);

	// calculate quantities in eigenspace
	debug_fit = true;
	fcn(par);
	debug_fit = false;

	// save plots
	g_de->SetName("g_de before fit"); g_de->SetTitle("before fit;eigenvector idx;#delta"); g_de->SetLineColor(2); g_de->SetMarkerColor(2);
	g_de->Write();

	g_cont->SetName("g_cont before fit"); g_cont->SetTitle("before fit;eigenvector idx;contribution to #chi^{2}"); g_cont->SetLineColor(2); g_cont->SetMarkerColor(2);
	g_cont->Write();
}

//----------------------------------------------------------------------------------------------------

void SavePlotsAfterFit(S2_FCN &fcn, TFitterMinuit *minuit)
{
	// get parameters
	int n_par = minuit->GetNumberTotalParameters();
	vector<double> par(n_par);
	for (int i = 0; i < n_par; i++)
		par[i] = minuit->GetParameter(i);

	// calculate quantities in eigenspace
	debug_fit = true;
	fcn(par);
	debug_fit = false;

	// save plots
	g_de->SetName("g_de after fit"); g_de->SetTitle("after fit;eigenvector idx;#delta"); g_de->SetLineColor(4); g_de->SetMarkerColor(4);
	g_de->Write();

	g_cont->SetName("g_cont after fit"); g_cont->SetTitle("after fit;eigenvector idx;contribution to #chi^{2}"); g_cont->SetLineColor(4); g_cont->SetMarkerColor(4);
	g_cont->Write();

	g_largest_cont_idx->SetName("g_largest_cont_idx after fit");
	g_largest_cont_idx->Write();
}

//----------------------------------------------------------------------------------------------------

void AnalyzeCompatibilityWithZero(const TVectorD &val, const TMatrixD &cov)
{
	printf(">> AnalyzeCompatibilityWithZero\n");

	int dim = val.GetNrows();

	printf("\tvalues and uncertainties\n");
	for (int i = 0; i < dim; i++)
		printf("\t\t%.3f +- %.3f (= %.3f si)\n", val(i), sqrt(cov(i, i)), fabs(val(i)) / sqrt(cov(i, i)));

	printf("\tcorrelation\n");
	for (int i = 0; i < dim; i++)
	{
		printf("\t");

		for (int j = 0; j < dim; j++)
		{
			double D = cov(i, i) * cov(j, j);
			double rho = (D > 0) ? cov(i, j) / sqrt(D) : 0.;
			printf("%+10.3f", rho);
		}

		printf("\n");
	}

	printf("\tcombined significances\n");

	TMatrixD cov_I(TMatrixD::kInverted, cov);

	double S2 = 0.;
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			S2 += val(i) * cov_I(i, j) * val(j);

	printf("\t\tS2 = %E\n", S2);

	double prob = TMath::Prob(S2, dim);
	printf("\t\tprob = %E\n", prob);
	double sig = sqrt(2.) * TMath::ErfcInverse(prob);
	printf("\t\tsig = %.1f\n", sig);
}

//----------------------------------------------------------------------------------------------------

unsigned int RunFit(const string & /*settings*/, Results &results, FILE *f_out_tex)
{
	printf(">> MethodSimpleFit::RunFit\n");

	SaveInputData();

	TDirectory *topDirectory = gDirectory;

	// set parameter offsets
	unsigned int n_fit_parameters = 0;

	if (useNormalisationFitParameter)
	{
		par_off_norm = 0;
		par_off_a = 1;
		par_off_b = 2;
		par_off_p0 = B_degree + 2;
		par_off_pAdd = B_degree + 3;
		n_fit_parameters = B_degree + 3;
	} else {
		par_off_norm = -1;
		par_off_a = 0;
		par_off_b = 1;
		par_off_p0 = B_degree + 1;
		par_off_pAdd = B_degree + 2;
		n_fit_parameters = B_degree + 2;
	}

	// initialize storage for interpolated Psi function
	interpolatedPsiRe = new TGraph(); interpolatedPsiRe->SetName("interpolatedPsiRe");
	interpolatedPsiIm = new TGraph(); interpolatedPsiIm->SetName("interpolatedPsiIm");

	// initialize fit function
	f_fit = new TF1("f_fit", f_fit_imp, 1E-5, 1., n_fit_parameters);
	f_fit->SetNpx(1000);

	t_min_fit = 0.;
	t_max_fit = 0.25;
	N_ii = 3; // default
	//N_ii = 10; // test

	// initialise fitter
	TFitterMinuit *minuit = new TFitterMinuit();
	S2_FCN fcn;
	minuit->SetMinuitFCN(&fcn);
	minuit->SetPrintLevel(10);
	minuit->CreateMinimizer();

	debug_fit = false;

	double m_chi2, m_edm, m_errdef;
	int m_n_par_var, m_n_par_tot, m_ndf;

	// initial point - normalisation
	if (useNormalisationFitParameter)
	{
		if (useNormalisationLimit)
			minuit->SetParameter(par_off_norm, "eta", 1., 0.01, 1. - 0.055, 1. + 0.055);
		else
			minuit->SetParameter(par_off_norm, "eta", 1., 0.01, 0., 0.);
	}

	if (useEtaFixed)
	{
		minuit->SetParameter(par_off_norm, "eta", eta_value_fix, 0., 0., 0.);
		minuit->FixParameter(par_off_norm);
	}

	// initial point - modulus
	char buf[200];
	minuit->SetParameter(par_off_a, "a", hfm->a / 1E8, 0.7, 0., 0.);	// without the factor 1E8

	if (useAFixed)
	{
		minuit->SetParameter(par_off_a, "a", a_value_fix, 0., 0., 0.);
		minuit->FixParameter(par_off_a);
	}

	for (unsigned int i = 1; i <= B_degree; i++)
	{
		sprintf(buf, "b%i", i);

		double val = 0., unc = 5.;
		double lim_low=0., lim_high=0.;
		if (i == 1) val = hfm->b1;
		if (i == 2) val = hfm->b2;
		if (i == 3) val = hfm->b3;
		if (i == 4) val = hfm->b4;

		if (i == 5) { val = hfm->b5; unc = 20.; lim_low = 0.; lim_high = 0.; }
		if (i == 6) { val = hfm->b6; unc = 20.; lim_low = 0.; lim_high = 0.; }
		if (i == 7) { val = hfm->b7; unc = 20.; lim_low = 0.; lim_high = 0.; }
		if (i == 8) { val = hfm->b8; unc = 20.; lim_low = 0.; lim_high = 0.; }
		/*
		if (i == 9) { val = hfm->b9; unc = 20.; lim_low = -100.; lim_high = +100.; }
		*/

		minuit->SetParameter(par_off_b + i - 1, buf, val, unc, lim_low, lim_high);
	}
	
	if (useB1Fixed)
	{
		minuit->SetParameter(par_off_b+0, "b1", b1_value_fix, 0., 0., 0.);
		minuit->FixParameter(par_off_b+0);
	}

	if (useB2Fixed)
	{
		minuit->SetParameter(par_off_b+1, "b2", b2_value_fix, 0., 0., 0.);
		minuit->FixParameter(par_off_b+1);
	}

	if (useB3Fixed)
	{
		minuit->SetParameter(par_off_b+2, "b3", b3_value_fix, 0., 0., 0.);
		minuit->FixParameter(par_off_b+2);
	}

	// initial point - phase
	minuit->SetParameter(par_off_p0, "p0", hfm->p0, 0.01, p0_lim_min, p0_lim_max);

	if (useRhoFixed)
	{
		minuit->SetParameter(par_off_p0, "p0", M_PI/2. - atan(rho_value_fix), 0., 0., 0.);
		minuit->FixParameter(par_off_p0);
	}

	// ------------------------------ F_CH fits, chosen formula

	printf("\n>> setting chosen interference formula\n");
	coulomb->mode = chosenCIMode;
	coulomb->Print();

	// TODO: uncomment ?
	/*
	bool release_p0 = (chosenCIMode != CoulombInterference::mPH);
	if (release_p0)
	{
		minuit->ReleaseParameter(par_off_p0);
		printf("* p0 released\n");
	} else {
		minuit->FixParameter(par_off_p0);
		printf("* p0 fixed\n");
	}
	*/

	useInterpolatedPsi = false;

	for (unsigned int ii = 0; ii < N_ii; ii++)
	{
		if (ii == 0)
		{
			BuildUncertaintyMatrix(NULL);
		} else {
			AdjustBinRepresentativePoints(f_fit);
			BuildUncertaintyMatrix(f_fit);

			InterpolatePsi();
			useInterpolatedPsi = true;
		}

		// update high |t| normalisation
		{
			double norm_corr = 1.;

			if (useNormalisationFitParameter)
				norm_corr = minuit->GetParameter(par_off_norm);

			if (useNormalisationFromA)
				norm_corr = EtaFromA(minuit->GetParameter(par_off_a));

			hfm->hts = 1./sqrt(norm_corr);
		}

		printf(">> setting hts = %.5f\n", hfm->hts);

		printf("\n\n>> F_CH, iteration %u\n", ii);

		char buf[50];
		sprintf(buf, "F_CH, iteration %u", ii);
		gDirectory = topDirectory->mkdir(buf);


		if (useInterpolatedPsi)
		{
			interpolatedPsiRe->Write();
			interpolatedPsiIm->Write();
		}

		SavePlotsBeforeFit(fcn, minuit);

		printf(">> state before fit\n");
		PrintMinuit2State(minuit);

		minuit->Minimize();

		printf(">> state after fit\n");
		PrintMinuit2State(minuit);

		minuit->GetStats(m_chi2, m_edm, m_errdef, m_n_par_var, m_n_par_tot);
		m_ndf = data_coll.size() - m_n_par_var;
		printf("chi^2 / ndf = %.3E / %i = %.3f\n", m_chi2, m_ndf, m_chi2 / m_ndf);

		SavePlotsAfterFit(fcn, minuit);
	}

	// ------------------------------ save last fit input

	gDirectory = topDirectory;

	TGraphErrors *g_data_coll_unc_stat = new TGraphErrors(); g_data_coll_unc_stat->SetName("g_data_coll_unc_stat");
	TGraphErrors *g_data_coll_unc_syst = new TGraphErrors(); g_data_coll_unc_syst->SetName("g_data_coll_unc_syst");
	TGraphErrors *g_data_coll_unc_full = new TGraphErrors(); g_data_coll_unc_full->SetName("g_data_coll_unc_full");

	TGraph *g_dataset_idx = new TGraph(); g_dataset_idx->SetName("g_dataset_idx");

	for (unsigned int i = 0; i < data_coll.size(); i++)
	{
		double stat_unc = data_coll[i].y_stat_unc;
		double full_unc = sqrt(data_coll_unc(i, i));
		double syst_unc = sqrt(full_unc*full_unc - stat_unc*stat_unc);

		g_data_coll_unc_stat->SetPoint(i, data_coll[i].x_repr, data_coll[i].y);
		g_data_coll_unc_stat->SetPointError(i, 0., stat_unc);

		g_data_coll_unc_syst->SetPoint(i, data_coll[i].x_repr, data_coll[i].y);
		g_data_coll_unc_syst->SetPointError(i, 0., syst_unc);

		g_data_coll_unc_full->SetPoint(i, data_coll[i].x_repr, data_coll[i].y);
		g_data_coll_unc_full->SetPointError(i, 0., full_unc);

		g_dataset_idx->SetPoint(i, data_coll[i].dataset, 0.);
	}

	g_data_coll_unc_stat->Write();
	g_data_coll_unc_syst->Write();
	g_data_coll_unc_full->Write();

	g_dataset_idx->Write();

	data_coll_unc.Write("data_coll_unc");

	// ------------------------------ print covariance matrix

	printf("\n>> covariance matrix\n");
	int n_par = minuit->GetNumberFreeParameters();
	for (int i = 0; i < n_par; i++)
	{
		printf("\t");
		for (int j = 0; j < n_par; j++)
		{
			printf("%12.3E,", minuit->GetCovarianceMatrixElement(i, j));
		}
		printf("\n");
	}

	// ------------------------------ draw fit result

	TCanvas *c, *c_rel, *c_relC;
	PlotDataExt(c, c_rel, c_relC);

	TGraph *g_fit_C = new TGraph(); g_fit_C->SetName("g_fit_C"); g_fit_C->SetLineColor(1);
	TGraph *g_fit_CH = new TGraph(); g_fit_CH->SetName("g_fit_CH"); g_fit_CH->SetLineColor(1);
	TGraph *g_fit_CH_rel = new TGraph(); g_fit_CH_rel->SetName("g_fit_CH_rel"); g_fit_CH_rel->SetLineColor(1);
	TGraph *g_fit_CH_relC = new TGraph(); g_fit_CH_relC->SetName("g_fit_CH_relC"); g_fit_CH_relC->SetLineColor(1);

	TGraph *g_fit_CH_Zv = new TGraph(); g_fit_CH_Zv->SetName("g_fit_CH_Zv"); g_fit_CH_Zv->SetLineColor(8);

	TGraph *g_fit_H = new TGraph(); g_fit_H->SetName("g_fit_H"); g_fit_H->SetLineColor(6);
	TGraph *g_fit_H_rel = new TGraph(); g_fit_H_rel->SetName("g_fit_H_rel"); g_fit_H_rel->SetLineColor(6);
	TGraph *g_fit_H_relC = new TGraph(); g_fit_H_relC->SetName("g_fit_H_relC"); g_fit_H_relC->SetLineColor(6);

	TGraph *g_Phase_H = new TGraph(); g_Phase_H->SetName("g_Phase_H"); g_Phase_H->SetLineColor(1);

	TGraph *g_ref = new TGraph(); g_ref->SetName("g_ref"); g_ref->SetLineColor(1);
	TGraph *g_refC = new TGraph(); g_refC->SetName("g_refC"); g_refC->SetLineColor(1);

	double dt = 1E-5;
	for (double t = 0. + dt; t <= 1.; t += dt)
	{
		if (t > 0.01)
			dt = 1E-4;

		if (t > 0.02)
			dt = 1E-3;

		double y_ref = A_ref * exp(-B_ref * t);
		double y_refC = y_ref + cnts->sig_fac * coulomb->Amp_pure(-t).Rho2();

		int idx = g_fit_CH->GetN();

		g_ref->SetPoint(idx, t, y_ref);
		g_refC->SetPoint(idx, t, y_refC);

		coulomb->mode = chosenCIMode;
		double si_CH = f_fit->Eval(t);
		g_fit_CH->SetPoint(idx, t, si_CH);
		g_fit_CH_rel->SetPoint(idx, t, si_CH / y_ref - 1.);
		g_fit_CH_relC->SetPoint(idx, t, si_CH / y_refC - 1.);

		coulomb->mode = coulomb->mPH;
		double si_H = f_fit->Eval(t);
		g_fit_H->SetPoint(idx, t, si_H);
		g_fit_H_rel->SetPoint(idx, t, si_H / y_ref - 1.);
		g_fit_H_relC->SetPoint(idx, t, (si_H - y_ref) / y_refC);

		coulomb->mode = coulomb->mPC;
		double si_C = f_fit->Eval(t);
		g_fit_C->SetPoint(idx, t, si_C);

		double phase = hfm->Amp(-t).Theta();
		g_Phase_H->SetPoint(idx, t, phase);

		g_fit_CH_Zv->SetPoint(idx, t, (si_CH - si_C - si_H) / (si_C + si_H));
	}

	g_ref->Write();
	g_refC->Write();

	g_fit_H->Write();
	g_fit_C->Write();
	g_fit_CH->Write();
	g_Phase_H->Write();

	gPad = c; g_fit_CH->Draw("l"); g_fit_H->Draw("l");
	gPad = c_rel; g_fit_CH_rel->Draw("l"); g_fit_H_rel->Draw("l");
	gPad = c_relC; g_fit_CH_relC->Draw("l"); g_fit_H_relC->Draw("l"); g_fit_CH_Zv->Draw("l");

	c->Write();
	c_rel->Write();
	c_relC->Write();

	coulomb->mode = chosenCIMode;
	f_fit->SetRange(1E-5, t_max_fit);
	f_fit->SetLineColor(2);
	f_fit->Write();

	// ------------------------------ calculate cross-sections and RMS of b

	// propagate parameters to model
	f_fit->Eval(0.);
	hfm->Print();

	coulomb->mode = coulomb->mPH;
	// TODO: this can be wrong, eta is included !
	double si_el = f_fit->Integral(0., 1.5);

	double a = cnts->sig_fac * minuit->GetParameter(par_off_a) * 1E8 * minuit->GetParameter(par_off_a) * 1E8;
	double p0 = minuit->GetParameter(par_off_p0);
	double rho = cos(p0) / sin(p0);

	double si_tot = sqrt( 16.*cnts->pi * cnts->sq_hbarc / (1. + rho * rho) * a );

	double si_inel = si_tot - si_el;

	TF1 *f_t_der_amp_sq = new TF1("f_t_der_amp_sq", f_t_der_amp_sq_imp, 0., 3., 0);
	/*
	f_t_der_amp_sq->SetNpx(10000);
	integrationMode = imModulus; double int_t_der_amp_sq_mod = f_t_der_amp_sq->Integral(0., 2., (double *) NULL, 1E-7); f_t_der_amp_sq->Write("f_int_mod");
	integrationMode = imPhase; double int_t_der_amp_sq_phase = f_t_der_amp_sq->Integral(0., 2., (double *) NULL, 1E-7); f_t_der_amp_sq->Write("f_int_phase");
	integrationMode = imFullSum; double int_t_der_amp_sq_fullsum = f_t_der_amp_sq->Integral(0., 2., (double *) NULL, 1E-7); f_t_der_amp_sq->Write("f_int_full_sum");
	*/
	integrationMode = imFull; double int_t_der_amp_sq_full = f_t_der_amp_sq->Integral(0., 2., (double *) NULL, 1E-7); //f_t_der_amp_sq->Write("f_int_full");

	/*
	printf("***** %E, %E, %E | %E\n", int_t_der_amp_sq_mod, int_t_der_amp_sq_phase, int_t_der_amp_sq_full,
			int_t_der_amp_sq_mod + int_t_der_amp_sq_phase);
	*/

	double int_amp_sq = si_el / cnts->sig_fac;
	double b_ms_el = 4. * int_t_der_amp_sq_full / int_amp_sq;

	double ep = 1E-5;
	double b_ms_tot = 4. * ( log(hfm->Amp(0.).Im()) - log(hfm->Amp(-ep).Im()) ) / ep;

	double b_ms_inel = si_tot / si_inel * b_ms_tot - si_el / si_inel * b_ms_el;

	double b_rms_el = sqrt(b_ms_el) * cnts->hbarc;		// conversion 1/GeV to fm
	double b_rms_tot = sqrt(b_ms_tot) * cnts->hbarc;
	double b_rms_inel = sqrt(b_ms_inel) * cnts->hbarc;

	printf("\n");
	printf(">> cross-sections: el = %.1f, inel = %.1f, tot = %.1f mb\n", si_el, si_inel, si_tot);
	printf(">> RMS of b: el = %.3f, inel = %.3f, tot = %.3f fm\n", b_rms_el, b_rms_inel, b_rms_tot);

	// ------------------------------ get fit quality

	minuit->GetStats(m_chi2, m_edm, m_errdef, m_n_par_var, m_n_par_tot);
	m_ndf = data_coll.size() - m_n_par_var;

	double prob = TMath::Prob(m_chi2, m_ndf);
	double sigma_eq = sqrt(2.) * TMath::ErfcInverse(prob);

	printf("\n>> fit quality: prob = %.3E, sig = %.3E\n\n", prob, sigma_eq);

	// ------------------------------ fill results

	results.method = 0;
	results.metric = 0;

	results.chi_sq = m_chi2;
	results.ndf = m_ndf;
	results.prob = prob;
	results.sig = sigma_eq;
	results.quality = results.chi_sq / results.ndf;

	results.p0 = minuit->GetParameter(par_off_p0);
	results.p0_e = minuit->GetParError(par_off_p0);

	results.rho = cos(results.p0) / sin(results.p0);
	results.rho_e = fabs(1. / sin(results.p0) / sin(results.p0)) * results.p0_e;

	results.a = cnts->sig_fac * minuit->GetParameter(par_off_a) * 1E8 * minuit->GetParameter(par_off_a) * 1E8;
	results.a_e = 2. * cnts->sig_fac * minuit->GetParameter(par_off_a) * 1E8 * minuit->GetParError(par_off_a) * 1E8;

	results.B = 2.*minuit->GetParameter(par_off_b);
	results.B_e = 2.*minuit->GetParError(par_off_b);

	results.si_el = si_el;
	results.si_inel = si_inel;
	results.si_tot = si_tot;

	results.b_rms_el = b_rms_el;
	results.b_rms_inel = b_rms_inel;
	results.b_rms_tot = b_rms_tot;

	// ------------------------------ error propagation

	// indices in GetCovarianceMatrixElement function refer only to non-fixed parameters
	int par_uncoff_p0 = par_off_p0;
	if (useB1Fixed)
		par_uncoff_p0--;
	if (useB2Fixed)
		par_uncoff_p0--;
	if (useB3Fixed)
		par_uncoff_p0--;

	// TODO: uncomment
	/*
	double V_pa_pa = minuit->GetCovarianceMatrixElement(par_off_a, par_off_a);
	double V_pa_p0 = minuit->GetCovarianceMatrixElement(par_off_a, par_uncoff_p0);
	double V_p0_p0 = minuit->GetCovarianceMatrixElement(par_uncoff_p0, par_uncoff_p0);

	printf("sqrt(V_pa_pa) = %.3f\n", sqrt(V_pa_pa));
	printf("sqrt(V_p0_p0) = %.3f\n", sqrt(V_p0_p0));
	printf("correlation pa and p0 = %.3f\n", V_pa_p0 / sqrt(V_pa_pa) / sqrt(V_p0_p0));

	double sc_a = 2. * cnts->sig_fac * minuit->GetParameter(par_off_a) * 1E8 * 1E8;
	double sc_rho = 1. / sin(results.p0) / sin(results.p0);

	double V_a_a = sc_a * V_pa_pa * sc_a;
	double V_a_rho = sc_a * V_pa_p0 * sc_rho;
	double V_rho_rho = sc_rho * V_p0_p0 * sc_rho;

	printf("sqrt(V_a_a) = %.3f\n", sqrt(V_a_a));
	printf("sqrt(V_rho_rho) = %.3f\n", sqrt(V_rho_rho));
	printf("correlation a and rho = %.3f\n", V_a_rho / sqrt(V_a_a) / sqrt(V_rho_rho));

	double der_a = si_tot/2. * 1./results.a;
	double der_rho = si_tot/2. * 2.*results.rho / (1. + results.rho*results.rho);

	double V_si_tot = der_a * V_a_a * der_a + 2.* der_a * V_a_rho * der_rho + der_rho * V_rho_rho * der_rho;
	double si_tot_unc = sqrt(V_si_tot);
	*/
	double si_tot_unc = 0.;

	// ------------------------------ save fit data

	TGraph *g_fit_data = new TGraph();
	g_fit_data->SetPoint(0, 0., m_chi2);
	g_fit_data->SetPoint(1, 0., m_ndf);
	g_fit_data->SetPoint(2, 0., prob);
	g_fit_data->SetPoint(3, 0., sigma_eq);

	g_fit_data->SetPoint(4, 0., results.rho);
	g_fit_data->SetPoint(5, 0., results.rho_e);

	g_fit_data->SetPoint(6, 0., results.a);
	g_fit_data->SetPoint(7, 0., results.a_e);

	g_fit_data->SetPoint(8, 0., results.B);
	g_fit_data->SetPoint(9, 0., results.B_e);

	double eta = 1., eta_unc = 0.;
	if (useNormalisationFitParameter)
	{
		eta = minuit->GetParameter(par_off_norm);
		eta_unc = minuit->GetParError(par_off_norm);
	}
	if (useNormalisationFromA)
	{
		eta = EtaFromA(minuit->GetParameter(par_off_a));
	}
	g_fit_data->SetPoint(10, 0., eta);
	g_fit_data->SetPoint(11, 0., eta_unc);

	g_fit_data->SetPoint(12, 0., data_coll.size());

	const double A_p = eta * cnts->sig_fac * minuit->GetParameter(par_off_a) * 1E8 * minuit->GetParameter(par_off_a) * 1E8;
	g_fit_data->SetPoint(13, 0., A_p);

	g_fit_data->SetPoint(14, 0., si_tot);
	g_fit_data->SetPoint(15, 0., si_tot_unc);

	g_fit_data->Write("g_fit_data");

	// ------------------------------ save fit data in TeX format

	fprintf(f_out_tex, "$n_{\\rm points} = %lu$, $\\chi^2/\\hbox{ndf} =  %.3f / (%lu - %u) = %.3f$ \\\\ \n",
		data_coll.size(), m_chi2, data_coll.size(), m_n_par_var, m_chi2 / m_ndf);
	fprintf(f_out_tex, "$\\rho = %.4f \\pm %.4f$, $\\si_{\\rm tot} = (%.2f \\pm %.2f)\\un{mb}$, $A' \\equiv \\et\\, \\d\\si^{\\rm N}/\\d t|_0 = %.2f\\un{mb/GeV^2}$, $\\eta = %.5f \\pm %.5f$ \\\\ \n",
		results.rho, results.rho_e, si_tot, si_tot_unc, A_p, eta, eta_unc);

	for (unsigned int i = 0; i < n_fit_parameters; ++i)
	{
		if (i > 0)
			fprintf(f_out_tex, ", ");

		fprintf(f_out_tex, "%u/%s = %.5E $\\pm$ %.5E\n", i, minuit->GetParName(i), minuit->GetParameter(i), minuit->GetParError(i));
	}

	return 0;
}


};	// namespace
