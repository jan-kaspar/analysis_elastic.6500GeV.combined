#include <vector>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintParameters(TF1 *f)
{
	for (signed int i = 0; i < f->GetNpar(); i++)
	{
		printf("\t\tconst double P%i = %+E;\n", i, f->GetParameter(i));
	}
}

//----------------------------------------------------------------------------------------------------

void SaveFitGraph(TF1 *f)
{
	TGraph *g = new TGraph();
	g->SetName("g_fit");

	for (double t = 0.; t <= 1.0; t += 1E-3)
	{
		g->SetPoint(g->GetN(), t, f->Eval(t));
	}

	g->Write();
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input
	TFile *f_in = TFile::Open("/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/DS-merged/merged.root");
	TH1D *h_in = (TH1D *) f_in->Get("ob-2-10-0.05/merged/combined/h_dsdt");

	// prepare output
	TFile *f_out = TFile::Open("fit.root", "recreate");

	// settings
	const double t_fit_min = 0.01;
	const double t_fit_max = 0.95;

	// fit 1
	{
		printf("\n\n-------------------- fit 1 --------------------\n");
		gDirectory = f_out->mkdir("fit 1");

		TF1 *ff = new TF1("ff1", "[0]*exp([1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x) + [5]*exp([6]*x + [7]*x*x)");
		ff->SetParameter(0, 6.11442e+02);
		ff->SetParameter(1, -2.07544e+01);
		ff->SetParameter(2, 1.01559e+00);
		ff->SetParameter(3, 2.23444e+01);
		ff->SetParameter(4, -9.65895e+01);
		ff->SetParameter(5, 2.89226e-04);
		ff->SetParameter(6, 1.44707e+01);
		ff->SetParameter(7, -1.09700e+01);

		h_in->Fit(ff, "", "", t_fit_min, t_fit_max);

		PrintParameters(ff);
		SaveFitGraph(ff);

		h_in->Write();
	}

	// fit 2
	{
		printf("\n\n-------------------- fit 2 --------------------\n");
		gDirectory = f_out->mkdir("fit 2");

		TF1 *ff = new TF1("ff2", "([0] + [1]*x) * exp([2]*x + [3]*x*x + [4]*x*x*x) + ([5] + [6]*x) * exp([7]*x)");
		ff->SetParameter(0, 6.24949e+02);
		ff->SetParameter(1, -2.56314e+02);
		ff->SetParameter(2, -2.04532e+01);
		ff->SetParameter(3, 8.49336e+00);
		ff->SetParameter(4, -1.60850e+01);
		ff->SetParameter(5, -1.11034e+01);
		ff->SetParameter(6, 2.25886e+01);
		ff->SetParameter(7, -7.02090e+00);

		h_in->Fit(ff, "", "", t_fit_min, t_fit_max);

		PrintParameters(ff);
		SaveFitGraph(ff);

		h_in->Write();
	}

	// fit 3
	{
		printf("\n\n-------------------- fit 3 --------------------\n");
		gDirectory = f_out->mkdir("fit 3");

		TF1 *ff = new TF1("ff3", "([0] + [1]*x) * exp([2]*x + [3]*x*x + [4]*x*x*x) + ([5] + [6]*x + [7]*x*x) * exp([8]*x + [9]*x*x)");
		ff->SetParameter(0, 7.16305e+02);
		ff->SetParameter(1, -2.37871e+02);
		ff->SetParameter(2, -1.96623e+01);
		ff->SetParameter(3, 9.34281e+00);
		ff->SetParameter(4, -1.50302e+01);
		ff->SetParameter(5, -1.02707e+02);
		ff->SetParameter(6, 8.08324e+01);
		ff->SetParameter(7, 2.20613e+02);
		ff->SetParameter(8, -1.29148e+01);
		ff->SetParameter(9, 3.09810e+00);

		h_in->Fit(ff, "", "", t_fit_min, t_fit_max);

		PrintParameters(ff);
		SaveFitGraph(ff);

		h_in->Write();
	}

	// fit 4
	{
		printf("\n\n-------------------- fit 4 --------------------\n");
		gDirectory = f_out->mkdir("fit 4");

		TF1 *ff = new TF1("ff4", "[1]*[1]*exp(2*[2]*x + 2*[3]*x*x + 2*[4]*x*x*x) + 2 * cos([0]) * [1]*exp([2]*x + [3]*x*x + [4]*x*x*x) * [5]*exp([6]*x) + [5]*[5]*exp(2*[6]*x)");
		ff->SetParameter(0, 2.65511e+00);
		ff->SetParameter(1, 2.55649e+01);
		ff->SetParameter(2, -1.02703e+01);
		ff->SetParameter(3, 4.42715e+00);
		ff->SetParameter(4, -6.83600e+00);
		ff->SetParameter(5, 9.00437e-01);
		ff->SetParameter(6, -2.16005e+00);

		h_in->Fit(ff, "", "", t_fit_min, t_fit_max);

		PrintParameters(ff);
		SaveFitGraph(ff);

		h_in->Write();
	}	

	// clean up
	delete f_out;

	return 0;
}
