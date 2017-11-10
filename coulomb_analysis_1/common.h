#include <vector>
#include <string>

#include "TCanvas.h"
#include "TSpline.h"
#include "TH2D.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

struct DataSetInfo
{
	string tag;

	string file_hist, obj_hist;
	TH1D *hist;

	string file_unc, obj_unc;
	TMatrixD *unc;

	double t_max;
	double norm_unc_add;
};

//----------------------------------------------------------------------------------------------------

struct InputData
{
	double norm_unc_global;
	vector<DataSetInfo> dataSetInfo;
};

//----------------------------------------------------------------------------------------------------

struct BinData
{
	unsigned int dataset;
	int bin;
	double x_left, x_right, x_repr;
	double y, y_stat_unc;
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

// recognized datasets
vector<DataSetInfo> stdDataSets;

void InitStdDataSets()
{
	vector<string> binnings;
	binnings.push_back("ob-1-20-0.05");
	binnings.push_back("ob-2-10-0.05");
	binnings.push_back("ob-3-5-0.05");

	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		DataSetInfo ds;
		ds.tag = "2500-2rp-" + binnings[bi];

		ds.file_hist = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/DS-merged/merged.root";
		ds.obj_hist = binnings[bi] + "/merged/combined/h_dsdt";

		ds.file_unc = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/systematics/matrix.root";
		ds.obj_unc = "matrices/all-but-norm/" + binnings[bi] + "/cov_mat";

		ds.t_max = 0.25;
		ds.norm_unc_add = 0.;

		stdDataSets.push_back(ds);
	}

	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		DataSetInfo ds;
		ds.tag = "simu-" + binnings[bi];

		ds.file_hist = "simu.root";
		ds.obj_hist = "h_dsdt";

		ds.file_unc = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/systematics/matrix.root";
		ds.obj_unc = "matrices/all-but-norm/" + binnings[bi] + "/cov_mat";

		ds.t_max = 0.25;
		ds.norm_unc_add = 0.;

		stdDataSets.push_back(ds);
	}
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

InputData inputData; 

vector<BinData> data_coll;

bool use_stat_unc, use_syst_unc;
TMatrixDSym data_coll_rel_unc;
TMatrixDSym data_coll_unc;
TMatrixDSym data_coll_unc_inv;

unsigned int B_degree;

CoulombInterference::CIMode chosenCIMode;

HadronicFitModel::PhaseMode phaseMode;

TGraph *interpolatedPsiRe = NULL, *interpolatedPsiIm = NULL;
TSpline3 *splinePsiRe = NULL, *splinePsiIm = NULL;

HadronicFitModel *hfm = NULL;

// fit limits
double p0_lim_min, p0_lim_max;
double p_A_lim_min, p_A_lim_max;
double p_ka_lim_min, p_ka_lim_max;
double p_tm_lim_min, p_tm_lim_max;

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

void AdjustBinRepresentativePoints(TF1 *f_fit)
{
	printf(">> AdjustBinRepresentativePoints(%p)\n", f_fit);

	for (unsigned int i = 0; i < data_coll.size(); i++)
	{
		double l = data_coll[i].x_left;
		double r = data_coll[i].x_right;

		double w = r - l;
		double V = f_fit->Integral(l, r) / w;

		double xr;
		while (r - l > w/10000)
		{
			xr = (r+l)/2.;
			if (f_fit->Eval(xr) < V)
			{
				r -= (r-l)/2.;
			} else {
				l += (r-l)/2.;
			}
		}
		xr = (r+l)/2.;

		//printf("\t\ti=%2u, bit=%2u : c=%.5f, x=%.5f\n", i, bi, c, xr);

		data_coll[i].x_repr = xr;
	}
}

//----------------------------------------------------------------------------------------------------

void BuildUncertaintyMatrix(TF1 *f, const vector<BinData> data, const TMatrixDSym &rel_unc, TMatrixDSym &unc, TMatrixDSym &unc_inv)
{
	printf(">> BuildUncertaintyMatrix(%p)\n", f);

	unc.ResizeTo(rel_unc);

	for (int i = 0; i < unc.GetNrows(); i++)
	{
		for (int j = 0; j < unc.GetNcols(); j++)
		{
			double y_ref_i = 0., y_ref_j = 0.;
			if (f == NULL)
			{
				// no smoothing
				y_ref_i = data[i].y;
				y_ref_j = data[j].y;
			} else {
				y_ref_i = f->Eval(data[i].x_repr);
				y_ref_j = f->Eval(data[j].x_repr);
			}
			
			double el_syst = y_ref_i * rel_unc(i, j) * y_ref_j;

			double el_stat = (i == j) ? pow(data[j].y_stat_unc, 2.) : 0.;

			double sum = 0.;
			if (use_stat_unc)
				sum += el_stat;
			if (use_syst_unc)
				sum += el_syst;

			unc(i, j) = sum;
		}
	}

	/*
	// discard off-diagonal terms
	printf(">> BuildUncertaintyMatrix: discarding off-diagonal terms in unc\n");
	for (int i = 0; i < unc.GetNrows(); i++)
	{
		for (int j = 0; j < unc.GetNcols(); j++)
		{
			if (i != j)
				unc(i, j) = 0.;
		}
	}
	*/

	unc_inv.ResizeTo(unc);
	unc_inv = unc;
	unc_inv.Invert();
}

//----------------------------------------------------------------------------------------------------

void BuildUncertaintyMatrix(TF1 *f)
{
	BuildUncertaintyMatrix(f, data_coll, data_coll_rel_unc, data_coll_unc, data_coll_unc_inv);
}

//----------------------------------------------------------------------------------------------------
// Psi-interpolation tools
//----------------------------------------------------------------------------------------------------

void BuildPsiSplines(bool wide_range = false)
{
	// code below interpolation corresponds to KL formula
	if (coulomb->mode != CoulombInterference::mKL)
		return;

	interpolatedPsiRe->Set(0);
	interpolatedPsiIm->Set(0);

	double dmt = 0.001;
	double mt_max = (wide_range) ? 1.5 : 0.23;
	for (double mt = 4E-4; mt < mt_max; mt += dmt)
	{
		TComplex Psi = coulomb->Psi_KL(-mt);

		int idx = interpolatedPsiRe->GetN();
		interpolatedPsiRe->SetPoint(idx, mt, Psi.Re());
		interpolatedPsiIm->SetPoint(idx, mt, Psi.Im());

		if (mt > 0.002)
			dmt = 0.002;

		if (mt > 0.01)
			dmt = 0.005;

		if (mt > 0.03)
			dmt = 0.01;

		if (mt > 0.07)
			dmt = 0.02;

		if (mt > 0.185)
			dmt = 0.01;
	}

	splinePsiRe = new TSpline3("splinePsiRe", interpolatedPsiRe->GetX(), interpolatedPsiRe->GetY(), interpolatedPsiRe->GetN());
	splinePsiIm = new TSpline3("splinePsiIm", interpolatedPsiIm->GetX(), interpolatedPsiIm->GetY(), interpolatedPsiIm->GetN());
}

//----------------------------------------------------------------------------------------------------

void ReleasePsiSplines()
{
	// code below interpolation corresponds to KL formula
	if (coulomb->mode != CoulombInterference::mKL)
		return;

	delete splinePsiRe; splinePsiRe = NULL;
	delete splinePsiIm; splinePsiIm = NULL;
}

//----------------------------------------------------------------------------------------------------
// Plotting tools
//----------------------------------------------------------------------------------------------------

double A_ref = 519.5;
double B_ref = 19.38;

//----------------------------------------------------------------------------------------------------

TGraphErrors* DataToGraph(bool relative = false, const char *name="g_data", int color = 1)
{
	TGraphErrors *g_data = new TGraphErrors();
	g_data->SetName(name);
	g_data->SetLineColor(color);
	g_data->SetMarkerColor(color);
	g_data->SetMarkerStyle(20);
	g_data->SetMarkerSize(0.8);

	for (unsigned int i = 0; i < data_coll.size(); i++)
	{
		int idx = g_data->GetN();

		double x = data_coll[i].x_repr;
		double y = data_coll[i].y;
		double y_unc = data_coll[i].y_stat_unc;

		double y_ref = A_ref * exp(-x * B_ref);

		if (relative)
		{
			y = y / y_ref - 1.;
			y_unc = y_unc / y_ref;
		}

		g_data->SetPoint(idx, x, y);
		g_data->SetPointError(idx, 0., y_unc);
	}

	return g_data;
}

//----------------------------------------------------------------------------------------------------

TGraphErrors* FitToGraph(TF1 *ff, bool relative = false, const char *name="g_fit", int color = 2)
{
	TGraphErrors *g_fit = new TGraphErrors();
	g_fit->SetName(name);
	g_fit->SetLineColor(color);
	g_fit->SetMarkerColor(color);

	double dt = 1E-5;
	for (double t = 0. + dt; t <= 4.0; t += dt)
	{
		if (t > 0.01)
			dt = 1E-4;

		if (t > 0.02)
			dt = 1E-3;
		
		if (t > 0.2)
			dt = 0.002;
		
		if (t > 1.0)
			dt = 0.01;

		double y = ff->Eval(t);
		
		double y_ref = A_ref * exp(-B_ref * t);

		if (relative)
		{
			y = y / y_ref - 1.;
		}

		int idx = g_fit->GetN();
		g_fit->SetPoint(idx, t, y);
		g_fit->SetPointError(idx, 0., 0.);
	}

	return g_fit;
}

//----------------------------------------------------------------------------------------------------

void SaveInputData()
{
	unsigned int prev_ds = 0;
	TGraphErrors *g_data = NULL;
	for (unsigned int i = 0; i < data_coll.size(); i++)
	{
		unsigned int ds = data_coll[i].dataset;
		if (g_data == NULL || ds != prev_ds)
		{
			if (g_data != NULL)
				g_data->Write();

			g_data = new TGraphErrors();
			char buf[50];
			sprintf(buf, "g_input_data%i", ds);
			g_data->SetName(buf);
			int colors[] = { 2, 4, 8, 6, 10 };
			g_data->SetMarkerColor(colors[ds]);
			g_data->SetLineColor(colors[ds]);
			g_data->SetMarkerStyle(20);
			g_data->SetMarkerSize(0.6);
			g_data->SetTitle(inputData.dataSetInfo[ds].tag.c_str());

			prev_ds = ds;
		}

		int idx = g_data->GetN();

		g_data->SetPoint(idx, data_coll[i].x_repr, data_coll[i].y);
		g_data->SetPointError(idx, 0., data_coll[i].y_stat_unc);
	}

	if (g_data != NULL)
		g_data->Write();
}

//----------------------------------------------------------------------------------------------------

void PlotData(TCanvas* &c, TCanvas* &c_rel)
{
	c = new TCanvas("fit canvas");
	c_rel = new TCanvas("fit canvas, rel");

	unsigned int prev_ds = 0;
	TGraphErrors *g_data = NULL, *g_data_rel = NULL;
	for (unsigned int i = 0; i < data_coll.size(); i++)
	{
		unsigned int ds = data_coll[i].dataset;
		if (g_data == NULL || ds != prev_ds)
		{
			if (g_data != NULL)
			{
				gPad = c; g_data->Draw((prev_ds == 0) ? "ap" : "p");
				gPad = c_rel; g_data_rel->Draw((prev_ds == 0) ? "ap" : "p");
			}

			g_data = new TGraphErrors();
			char buf[50];
			sprintf(buf, "g_data%i", ds);
			int colors[] = { 2, 4, 8, 6, 10 };
			g_data->SetMarkerColor(colors[ds]);
			g_data->SetLineColor(colors[ds]);
			g_data->SetMarkerStyle(20);
			g_data->SetMarkerSize(0.6);
			g_data->SetName(buf);
			g_data->SetTitle(inputData.dataSetInfo[ds].tag.c_str());
			
			g_data_rel = new TGraphErrors(*g_data);
			sprintf(buf, "g_data_rel%i", ds);
			g_data_rel->SetName(buf);

			prev_ds = ds;
		}

		int idx = g_data->GetN();

		g_data->SetPoint(idx, data_coll[i].x_repr, data_coll[i].y);
		g_data->SetPointError(idx, 0., data_coll[i].y_stat_unc);

		double y_ref = A_ref * exp(-B_ref * data_coll[i].x_repr);

		g_data_rel->SetPoint(idx, data_coll[i].x_repr, data_coll[i].y / y_ref - 1.);
		g_data_rel->SetPointError(idx, 0., data_coll[i].y_stat_unc / y_ref);
	}

	if (g_data)
	{
		gPad = c; g_data->Draw((prev_ds == 0) ? "ap" : "p");
		gPad = c_rel; g_data_rel->Draw((prev_ds == 0) ? "ap" : "p");
	}
}

//----------------------------------------------------------------------------------------------------

TH2D* BuildLimitFrame(const vector<TGraphErrors *> &coll, double margin = 0.1)
{
	double x_min = +1E100, x_max = -1E100;
	double y_min = +1E100, y_max = -1E100;

	for (const auto &g : coll)
	{
		for (int i = 0; i < g->GetN(); i++)
		{
			double x, y;
			g->GetPoint(i, x, y);

			double y_unc = g->GetErrorY(i);

			x_min = min(x_min, x);
			x_max = max(x_max, x);

			y_min = min(y_min, y-y_unc);
			y_max = max(y_max, y+y_unc);
		}
	}

	double x_width = x_max - x_min;
	x_min -= margin * x_width;
	x_max += margin * x_width;

	double y_width = y_max - y_min;
	y_min -= margin * y_width;
	y_max += margin * y_width;

	return new TH2D("", "", 50, x_min, x_max, 50, y_min, y_max);
}

//----------------------------------------------------------------------------------------------------

void PlotDataExt(TCanvas* &c, TCanvas* &c_rel, TCanvas* &c_relC)
{
	// build graphs
	unsigned int prev_ds = 0;

	vector<TGraphErrors *> coll_g_data, coll_g_data_rel, coll_g_data_relC;
	TGraphErrors *g_data=NULL, *g_data_rel=NULL, *g_data_relC=NULL;

	for (unsigned int i = 0; i < data_coll.size(); i++)
	{
		unsigned int ds = data_coll[i].dataset;
		if (coll_g_data.empty() || ds != prev_ds)
		{
			g_data = new TGraphErrors();
			char buf[50];
			sprintf(buf, "g_data%i", ds);
			int colors[] = { 2, 4, 8, 6, 10 };
			g_data->SetMarkerColor(colors[ds]);
			g_data->SetLineColor(colors[ds]);
			g_data->SetMarkerStyle(20);
			g_data->SetMarkerSize(0.6);
			g_data->SetName(buf);
			g_data->SetTitle(inputData.dataSetInfo[ds].tag.c_str());
			coll_g_data.push_back(g_data);
			
			g_data_rel = new TGraphErrors(*g_data);
			sprintf(buf, "g_data_rel%i", ds);
			g_data_rel->SetName(buf);
			coll_g_data_rel.push_back(g_data_rel);

			g_data_relC = new TGraphErrors(*g_data);
			sprintf(buf, "g_data_relC%i", ds);
			g_data_relC->SetName(buf);
			coll_g_data_relC.push_back(g_data_relC);

			prev_ds = ds;
		}

		int idx = g_data->GetN();

		g_data->SetPoint(idx, data_coll[i].x_repr, data_coll[i].y);
		g_data->SetPointError(idx, 0., data_coll[i].y_stat_unc);

		double y_ref = A_ref * exp(-B_ref * data_coll[i].x_repr);

		g_data_rel->SetPoint(idx, data_coll[i].x_repr, data_coll[i].y / y_ref - 1.);
		g_data_rel->SetPointError(idx, 0., data_coll[i].y_stat_unc / y_ref);

		// reference including Coulomb

		y_ref = A_ref * exp(-B_ref * data_coll[i].x_repr);
		y_ref += cnts->sig_fac * coulomb->Amp_pure(-data_coll[i].x_repr).Rho2();

		g_data_relC->SetPoint(idx, data_coll[i].x_repr, data_coll[i].y / y_ref - 1.);
		g_data_relC->SetPointError(idx, 0., data_coll[i].y_stat_unc / y_ref);
	}
	
	// build canvases
	c = new TCanvas("fit canvas");
	c->SetGridx(1); c->SetGridy(1);
	BuildLimitFrame(coll_g_data)->Draw();

	c_rel = new TCanvas("fit canvas, rel");
	c_rel->SetGridx(1); c_rel->SetGridy(1);
	BuildLimitFrame(coll_g_data_rel)->Draw();

	c_relC = new TCanvas("fit canvas, relC");
	c_relC->SetGridx(1); c_relC->SetGridy(1);
	BuildLimitFrame(coll_g_data_relC)->Draw();

	for (unsigned int idx = 0; idx < coll_g_data.size(); ++idx)
	{
		gPad = c; coll_g_data[idx]->Draw("p");
		gPad = c_rel; coll_g_data_rel[idx]->Draw("p");
		gPad = c_relC; coll_g_data_relC[idx]->Draw("p");
	}
}

//----------------------------------------------------------------------------------------------------
// parameter-comparison fit tools
//----------------------------------------------------------------------------------------------------

vector<TF1 *> feature_fcn;

vector<double> t_points;
TMatrixD dot_product_matrix;

vector<TVectorD> vec_set_fcn, vec_set_on;

TMatrixD GS_mat_T;

TVectorD data_rel;
TMatrixD data_unc_rel;
TVectorD data_fcn_coef;

//----------------------------------------------------------------------------------------------------

double Dot(const TVectorD &a, const TVectorD &b)
{
	double S = 0.;

	for (unsigned int i = 0; i < t_points.size(); i++)
		for (unsigned int j = 0; j < t_points.size(); j++)
			S += a(i) * dot_product_matrix(i, j) * b(j);

	return S;
}

//----------------------------------------------------------------------------------------------------

void GramSchmidt(const vector<TVectorD> &v_in, vector<TVectorD> &v_out, TMatrixD &mat)
{
	mat.ResizeTo(v_in.size(), v_in.size());

	for (unsigned int vi = 0; vi < v_in.size(); vi++)
	{
		//printf("\tvi = %u\n", vi);

		TVectorD v = v_in[vi];

		mat(vi, vi) = 1.;

		// remove projections
		for (unsigned int ri = 0; ri < vi; ri++)
		{
			//printf("\t\tremoving %u\n", ri);
			double proj = Dot(v_out[ri], v);
			v -= proj * v_out[ri];
	
			for (unsigned int si = 0; si <= ri; si++)
				mat(vi, si) -= proj * mat(ri, si);
		}

		// scale to unit size
		double norm = Dot(v, v);
		v *= 1./sqrt(norm);

		for (unsigned int si = 0; si <= vi; si++)
			mat(vi, si) *= 1./sqrt(norm);

		v_out.push_back(v);
	}
}

//----------------------------------------------------------------------------------------------------

double Eval(const vector<TF1 *> &functions, const TVectorD &data_fcn_coef, double x)
{
	double S = 0.;

	for (unsigned int fi = 0; fi < functions.size(); fi++)
	{
		S += functions[fi]->Eval(x) * data_fcn_coef(fi);
	}

	return S;
}

//----------------------------------------------------------------------------------------------------

void BuildVectors()
{
	// sampling points
	t_points.clear();
	for (unsigned int i = 0; i < data_coll.size(); i++)
		t_points.push_back(data_coll[i].x_repr);

	// function (feature) vectors
	vec_set_fcn.clear();
	for (unsigned int fi = 0; fi < feature_fcn.size(); fi++)
	{
		TVectorD v(t_points.size());
		
		for (unsigned int pi = 0; pi < t_points.size(); pi++)
			v(pi) = feature_fcn[fi]->Eval(t_points[pi]);

		vec_set_fcn.push_back(v);
	}

	// orthonormalise the set
	TMatrixD GS_mat;
	
	//printf("\n* GramSchmidt\n");
	vec_set_on.clear();
	GramSchmidt(vec_set_fcn, vec_set_on, GS_mat);
	GS_mat_T.ResizeTo(GS_mat);
	GS_mat_T.Transpose(GS_mat);
}

//----------------------------------------------------------------------------------------------------

void InitFeatures(const string &function_list)
{
	size_t prev = 0;
	while (true)
	{
		size_t pos = function_list.find_first_of(";", prev);
		string substr = function_list.substr(prev, (pos == string::npos) ? string::npos : pos-prev);

		feature_fcn.push_back(new TF1(substr.c_str(), substr.c_str()));

		prev = pos+1;
		if (pos == string::npos)
			break;
	}

	printf("\n* feature_fcn\n");
	for (unsigned int i = 0; i < feature_fcn.size(); i++)
		printf("\t%2u: %s\n", i, feature_fcn[i]->GetName());
}

//----------------------------------------------------------------------------------------------------
// tools for calculating RMS's of b
//----------------------------------------------------------------------------------------------------

enum IntegrationMode { imModulus, imPhase, imFull, imFullSum } integrationMode;

double f_t_der_amp_sq_imp(double x[], double [])
{
	double t = -x[0];

	double ep = 1E-7;

	double r = 0.;

	if (integrationMode == imModulus)
	{
		double mod_der = ( hfm->Amp(t).Rho() - hfm->Amp(t - ep).Rho() ) / ep;
		r = fabs(t) * mod_der * mod_der;
	}

	if (integrationMode == imPhase)
	{
		double phase_der = ( hfm->Amp(t).Theta() - hfm->Amp(t - ep).Theta() ) / ep;
		r = fabs(t) * hfm->Amp(t).Rho2() * phase_der * phase_der;
	}

	if (integrationMode == imFull)
	{
		TComplex amp_der = ( hfm->Amp(t) - hfm->Amp(t - ep) ) / ep;
		r = fabs(t) * amp_der.Rho2();
	}

	if (integrationMode == imFullSum)
	{
		double mod_der = ( hfm->Amp(t).Rho() - hfm->Amp(t - ep).Rho() ) / ep;
		double phase_der = ( hfm->Amp(t).Theta() - hfm->Amp(t - ep).Theta() ) / ep;
		r = fabs(t) * (mod_der * mod_der  +  hfm->Amp(t).Rho2() * phase_der * phase_der);
	}

	return r;
}

//----------------------------------------------------------------------------------------------------
// Minuit1 tools
//----------------------------------------------------------------------------------------------------

void PrintMinuit1State(TF1 *ff)
{
	for (int i = 0; i < ff->GetNpar(); i++)
	{
		double min, max;
		ff->GetParLimits(i, min, max);

		printf("%8s: value = %+8.3f, unc = %+8.3f, min = %+8.3f, max = %+8.3f", ff->GetParName(i),
			ff->GetParameter(i), ff->GetParError(i), min, max);

		if (min == max && min != 0.)
			printf(", fixed");

		if (max > min)
			printf(", limited");

		printf("\n");
	}

	return;
}

//----------------------------------------------------------------------------------------------------
// Minuit2 tools
//----------------------------------------------------------------------------------------------------

void PrintMinuit2State(TFitterMinuit *minuit)
{
	for (int pi = 0; pi < minuit->GetNumberTotalParameters(); pi++)
	{
		char name[20];
		double value, verr, vlow, vhigh;
		minuit->GetParameter(pi, name, value, verr, vlow, vhigh);
		printf("\t%i, %6s, value=%+7.3f, unc=%6.3f, min=%+7.3f, max=%+7.3f", pi, minuit->GetParName(pi),
			value, verr, vlow, vhigh);

		if (vhigh > vlow)
			printf(", limited");

		if (minuit->IsFixed(pi))
			printf(", fixed");

		printf("\n");
	}
}
