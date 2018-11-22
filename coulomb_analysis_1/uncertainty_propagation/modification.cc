#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TSpline.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "../../command_line_tools.h"

using namespace std;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: simulation [option] [option] ...\n");
	printf("OPTIONS:\n");
	printf("    -unc-file <file>            \n");
	printf("    -binning <string>           \n");
	printf("    -apply-stat-err <bool>      \n");
	printf("    -apply-syst-err <bool>      \n");
	printf("    -apply-norm-err <bool>      \n");
	printf("    -seed <int>                 \n");
	printf("    -norm-unc-sigma <double>    \n");
	printf("    -norm-unc-bias <double>     \n");
	printf("    -output <file>              \n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string fn_uncertainties = "";

	string binning = "";

	bool apply_stat_err = false;
	bool apply_syst_err = false;
	bool apply_norm_err = false;

	unsigned int seed = 1;

	double norm_unc_sigma = 0.;
	double norm_unc_bias = 0.;

	string fn_output = "";

	// parse command line
	for (int argi = 1; (argi < argc) && !cl_error; argi++)
	{
		if (TestStringParameter(argc, argv, argi, "-unc-file", fn_uncertainties)) continue;

		if (TestStringParameter(argc, argv, argi, "-binning", binning)) continue;
		
		if (TestBoolParameter(argc, argv, argi, "-apply-stat-err", apply_stat_err)) continue;
		if (TestBoolParameter(argc, argv, argi, "-apply-syst-err", apply_syst_err)) continue;
		if (TestBoolParameter(argc, argv, argi, "-apply-norm-err", apply_norm_err)) continue;

		if (TestUIntParameter(argc, argv, argi, "-seed", seed)) continue;

		if (TestDoubleParameter(argc, argv, argi, "-norm-unc-sigma", norm_unc_sigma)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-norm-unc-bias", norm_unc_bias)) continue;

		if (TestStringParameter(argc, argv, argi, "-output", fn_output)) continue;

		printf("ERROR: '%s' argument not understood.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 2;
	}

	// print configuration
	printf("* configuration\n");
	printf("    fn_uncertainties = %s\n", fn_uncertainties.c_str());
	printf("    binning = %s\n", binning.c_str());

	printf("    apply_stat_err = %u\n", apply_stat_err);
	printf("    apply_syst_err = %u\n", apply_syst_err);
	printf("    apply_norm_err = %u\n", apply_norm_err);

	printf("    seed = %u\n", seed);

	printf("    norm_unc_sigma = %f\n", norm_unc_sigma);
	printf("    norm_unc_bias = %f\n", norm_unc_bias);

	printf("    fn_output = %s\n", fn_output.c_str());

	// input configuration validation
	if (binning.empty() || fn_output.empty())
	{
		printf("ERROR: some input not specified.\n");
		PrintUsage();
		return 3;
	}

	// get start histogram
	TFile *f_in_model = TFile::Open("/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/DS-merged/merged.root");
	TH1D *h_dsdt_start = (TH1D *) f_in_model->Get((binning + "/merged/combined/h_dsdt").c_str());

	gDirectory = NULL;
	TH1D *h_dsdt = new TH1D(*h_dsdt_start);
	h_dsdt->Sumw2();

	// set random seed
	gRandom->SetSeed(seed);

	// set statistical uncertainties
	{
		TFile *f_in_unc = TFile::Open(fn_uncertainties.c_str());
		TH1D *h_rel_stat_unc = (TH1D *) f_in_unc->Get((binning + "/h_rel_stat_unc").c_str());

		for (int bi = 1; bi <= h_dsdt->GetNbinsX(); bi++)
		{
			const double u_rel = h_rel_stat_unc->GetBinContent(bi);
			const double v_true = h_dsdt_start->GetBinContent(bi);

			h_dsdt->SetBinError(bi, u_rel * v_true);
		}

		delete f_in_unc;
	}

	// apply statistical errors
	if (apply_stat_err)
	{
		TFile *f_in_unc = TFile::Open(fn_uncertainties.c_str());
		TH1D *h_rel_stat_unc = (TH1D *) f_in_unc->Get((binning + "/h_rel_stat_unc").c_str());

		for (int bi = 1; bi <= h_dsdt->GetNbinsX(); bi++)
		{
			const double u_rel = h_rel_stat_unc->GetBinContent(bi);
			const double v = h_dsdt->GetBinContent(bi);
			const double v_true = h_dsdt_start->GetBinContent(bi);

			const double R = gRandom->Gaus();

			h_dsdt->SetBinContent(bi, v + R * u_rel * v_true);
			h_dsdt->SetBinError(bi, u_rel * v_true);
		}

		delete f_in_unc;
	}

	// apply systematic errors
	if (apply_syst_err)
	{
		TFile *f_in_unc = TFile::Open(fn_uncertainties.c_str());
		TMatrixD *m_syst_gen = (TMatrixD *) f_in_unc->Get((binning + "/m_syst_gen").c_str());

		const unsigned int dim = m_syst_gen->GetNrows();
		TVectorD v_R(dim);
		for (unsigned int i = 0; i < dim; i++)
		{
			v_R(i) = gRandom->Gaus();
		}
		TVectorD v_rel_err = (*m_syst_gen) * v_R;

		for (int bi = 1; bi <= h_dsdt->GetNbinsX(); bi++)
		{
			const double v = h_dsdt->GetBinContent(bi);
			const double v_true = h_dsdt_start->GetBinContent(bi);

			const unsigned int i = bi - 1;

			h_dsdt->SetBinContent(bi, v + v_true * v_rel_err(i));
		}

		delete f_in_unc;
	}

	// apply normalisation errors
	double scale = 1. + norm_unc_bias;
	if (apply_norm_err)
		scale += gRandom->Gaus() * norm_unc_sigma;
	h_dsdt->Scale(scale);

	// prepare output
	TFile *f_out = TFile::Open(fn_output.c_str(), "recreate");

	h_dsdt->Write("h_dsdt");

	// clean up
	delete f_in_model;
	delete f_out;

	return 0;
}
