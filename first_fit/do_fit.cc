#include "include/Elegent/Constants.h"
#include "include/Elegent/CoulombInterference.h"

#include "../HadronicFitModel.h"

#include "fit.h"

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"

#include "TRandom3.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

void RunFit(TH1D *h, double t_min, double t_max, Fit::Results &results)
{
	// settings
	Fit::B_degree = 3;
	Fit::h_fit = h;

	Fit::bin_fit_min = h->GetXaxis()->FindBin(t_min);
	Fit::bin_fit_max = h->GetXaxis()->FindBin(t_max);

	// run fit
	Fit::RunFit(results);
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: normalize_and_fit <input file> <histogram path> <normAdjustment> <output file>\n");
}

//----------------------------------------------------------------------------------------------------

int main(/*int argc, char **argv*/)
{
	// defaults
	double t_min = 8.1E-4;
	double t_max = 0.2;

	string inputFileName = "../../beta2500/2rp/DS-merged/merged.root";
	string histPath = "ob-2-10-0.05/merged/combined/h_dsdt";

	string outputFileName = "do_fit.root";

	// init Elegent
	Constants::Init(2*6500, cnts->mPP);
    cnts->Print();

	coulomb->ffType = coulomb->ffPuckett;
	coulomb->Print();

	// get input
	TFile *f_in = TFile::Open(inputFileName.c_str());
	if (f_in == NULL)
	{
		printf("ERROR: can't open input file '%s'.\n", inputFileName.c_str());
		return 2;
	}
	
	TH1D *h_in = (TH1D *) f_in->Get(histPath.c_str());

	if (h_in == NULL)
	{
		printf("ERROR: can't load input histogram '%s'.\n", histPath.c_str());
		return 3;
	}

	// prepare output
	TFile *f_out = TFile::Open(outputFileName.c_str(), "recreate");

	h_in->Write("h_in");
	
	// run fit
	Fit::Results fr;
	RunFit(h_in, t_min, t_max, fr);

	// save results
	TGraph *g_results = new TGraph();
	g_results->SetName("g_results");

	g_results->SetPoint(0, fr.A, fr.A_e);
	g_results->SetPoint(1, fr.B, fr.B_e);
	g_results->SetPoint(2, fr.p0, fr.p0_e);
	g_results->SetPoint(3, fr.rho, fr.rho_e);

	g_results->Write();

	// clean up
	delete f_out;
	delete f_in;

	return 0;
}
