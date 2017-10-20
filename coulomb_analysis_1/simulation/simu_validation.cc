#include "TFile.h"
#include "TH1D.h"

#include "../../stat.h"

// TODO: clean
/*
#include "TGraph.h"
#include "TSpline.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "../command_line_tools.h"
*/


using namespace std;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: simu_validation <output> <ref simu> <simu1> <simu2> <...>\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	if (argc < 4)
	{
		printf("ERROR: need at least 4 arguments\n");
		PrintUsage();
		return 1;
	}

	// defaults
	string fn_output = argv[1];
	string fn_ref = argv[2];

	// load reference
	TFile *f_in_ref = TFile::Open(fn_ref.c_str());
	TH1D *h_ref = (TH1D *) f_in_ref->Get("h_dsdt");

	// prepare data
	Stat st(h_ref->GetNbinsX());

	for (int ai = 3; ai < argc; ai++)
	{
		TFile *f_in = TFile::Open(argv[ai]);
		TH1D *h = (TH1D *) f_in->Get("h_dsdt");

		vector<double> diff(h_ref->GetNbinsX());
		for (int bi = 1; bi <= h_ref->GetNbinsX(); bi++)
		{
			const unsigned int i = bi - 1;

			const double v_ref = h_ref->GetBinContent(bi);
			const double v = h->GetBinContent(bi);

			diff[i] = v - v_ref;
		}

		st.Fill(diff);

		delete f_in;
	}

	// save output
	TFile *f_out = TFile::Open(fn_output.c_str(), "recreate");

	TH1D *h_bias = new TH1D(*h_ref);
	TH1D *h_unc = new TH1D(*h_ref);
	for (int bi = 1; bi <= h_ref->GetNbinsX(); bi++)
	{
		const unsigned int i = bi - 1;

		h_bias->SetBinContent(bi, st.GetMean(i));
		h_bias->SetBinError(bi, st.GetMeanUnc(i));

		h_unc->SetBinContent(bi, st.GetStdDev(i));
		h_unc->SetBinError(bi, st.GetStdDevUnc(i));
	}

	h_bias->Write("h_bias");
	h_unc->Write("h_unc");

	TH1D *h_bias_rel = new TH1D(*h_bias);
	h_bias_rel->Divide(h_ref);
	h_bias_rel->Write("h_bias_rel");

	TH1D *h_unc_rel = new TH1D(*h_unc);
	h_unc_rel->Divide(h_ref);
	h_unc_rel->Write("h_unc_rel");

	// clean up
	delete f_in_ref;
	delete f_out;

	return 0;
}
