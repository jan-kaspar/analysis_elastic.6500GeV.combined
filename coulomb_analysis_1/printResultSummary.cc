#include <vector>
#include <string>

#include "../command_line_tools.h"
#include "classes.h"

#include "TFile.h"
#include "TGraph.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintOne(const Results &r, TGraph *g_data)
{
	//const double eta = g_data->GetY()[10];
	//const double eta_unc = g_data->GetY()[11];
	const double n_points = g_data->GetY()[12];
	const double A_p = g_data->GetY()[13];
	const double si_tot = g_data->GetY()[14];
	const double si_tot_unc = g_data->GetY()[15];

	printf("rho=%.4f+-%.4f, si_tot=%.2f+-%.2f mb, n_points=%3.0f, chisq/ndf=%.2f, dsigma_H/dt|0=%.1f mb/GeV^2",
		r.rho, r.rho_e,
		si_tot, si_tot_unc,
		n_points, r.chi_sq_norm,
		A_p
	);
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: printResultTableShort [option] ...\n");
	printf("OPTIONS:\n");
	printf("    -input <directory>     select input directory\n");
	printf("    -add-fit <string>      add fit to output\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// default parameters
	string inputDir = "";

	vector<string> fits;

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestStringParameter(argc, argv, argi, "-input", inputDir)) continue;

		if (strcmp(argv[argi], "-add-fit") == 0)
		{
			fits.push_back(argv[++argi]);
			continue;
		}
		
		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// print fit data
	for (const auto &f : fits)
	{
		string dir = inputDir + "/" + f;

		string file = dir + "/fit.out";
		Results r;
		r.Load(file);

		file = dir + "/fit.root";
		TFile *f_in = TFile::Open(file.c_str());
		if (f_in == NULL)
		{
			printf("ERROR\n");
			continue;
		}

		TGraph *g_fit_data = (TGraph *) f_in->Get("g_fit_data");

		printf("%s: ", f.c_str());
		PrintOne(r, g_fit_data);
		printf("\n");

		delete f_in;
	}
	return 0;
}
