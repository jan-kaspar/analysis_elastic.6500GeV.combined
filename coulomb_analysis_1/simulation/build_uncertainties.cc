#include "include/Elegent/Constants.h"
#include "include/Elegent/CoulombInterference.h"

#include "../../HadronicFitModel.h"

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TSpline.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

TH1D* MakeRelUncHist(TH1D *h_in)
{
	TH1D *h_out = new TH1D(*h_in);

	for (int bi = 1; bi <= h_out->GetNbinsX(); ++bi)
	{
		double v = h_out->GetBinContent(bi);
		double u = h_out->GetBinError(bi);

		double u_r = 0.;
		if (v > 0.)
			u_r = u / v;

		h_out->SetBinContent(bi, u_r);
		h_out->SetBinError(bi, 0.);
	}

	return h_out;
}
		
//----------------------------------------------------------------------------------------------------

void MakeSystGenerator(const TMatrixDSym &cov_mat, const TH1D *h_stddev_ref)
{
	printf("%i, %i\n", cov_mat.GetNrows(), h_stddev_ref->GetNbinsX());

	TMatrixDSymEigen eig_decomp(cov_mat);
	TVectorD eig_values(eig_decomp.GetEigenValues());
	TMatrixDSym S(cov_mat.GetNrows());
	for (int i = 0; i < cov_mat.GetNrows(); i++)
		S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
	TMatrixD opt_per_gen = eig_decomp.GetEigenVectors() * S;

	opt_per_gen.Write("m_syst_gen");

	// validation: regenerate h_stddev
	TH1D *h_stddev_regen = new TH1D(*h_stddev_ref);
	for (int bi = 1; bi <= h_stddev_regen->GetNbinsX(); bi++)
	{
		const int i = bi - 1;

		double S = 0;
		for (int j = 0; j < h_stddev_regen->GetNbinsX(); j++)
		{
			const double e = opt_per_gen(i, j);
			S += e*e;
		}

		h_stddev_regen->SetBinContent(bi, sqrt(S));
	}

	h_stddev_regen->Write("h_stddev_regen");
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// defaults
	string hist_file_name = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/DS-merged/merged.root";

	// TODO: update !!
	string syst_file_name = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/systematics-old-ni/DS-fill5313/matrix_numerical_integration.root";

	/*
		ds.file_hist = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/DS-merged/merged.root";
		ds.obj_hist = binnings[bi] + "/merged/combined/h_dsdt";

		ds.file_unc = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/systematics-old-ni/DS-fill5313/matrix_numerical_integration.root";
		ds.obj_unc = "matrices/all-but-norm/combined/" + binnings[bi] + "/cov_mat";
	*/

	// binnings
	vector<string> binnings;
	binnings.push_back("ob-1-20-0.05");
	binnings.push_back("ob-2-10-0.05");
	binnings.push_back("ob-3-5-0.05");

	// prepare input
	TFile *f_in_hist = TFile::Open(hist_file_name.c_str());
	TFile *f_in_syst = TFile::Open(syst_file_name.c_str());

	// prepare output
	TFile *f_out = TFile::Open("build_uncertainties.root", "recreate");

	for (const auto &binning : binnings)
	{
		gDirectory = f_out->mkdir(binning.c_str());

		// statistical uncertainties
		TH1D *h_data = (TH1D *) f_in_hist->Get((binning + "/merged/combined/h_dsdt").c_str());
		MakeRelUncHist(h_data)->Write("h_rel_stat_unc");

		// systematic uncertainties
		string path_syst = "matrices/all-but-norm/combined";
		TH1D *h_syst_stddev_orig = (TH1D *) f_in_syst->Get((path_syst + "/" + binning + "/h_stddev").c_str());
		h_syst_stddev_orig->Write("h_syst_stddev_orig");

		TMatrixDSym *m_syst_cov = (TMatrixDSym *) f_in_syst->Get((path_syst + "/" + binning + "/cov_mat").c_str());
		MakeSystGenerator(*m_syst_cov, h_syst_stddev_orig);
	}

	// clean up
	delete f_out;

	return 0;
}
