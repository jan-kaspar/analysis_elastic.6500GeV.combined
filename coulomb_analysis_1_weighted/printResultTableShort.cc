#include <vector>
#include <string>

#include "classes.h"

#include "../command_line_tools.h"

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
	printf("    -input <directory>    select input directory\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// default parameters
	string inputDir = "fits";

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestStringParameter(argc, argv, argi, "-input", inputDir)) continue;
		
		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// define combinations
	struct Combination
	{
		string expN;
		string binning;
		string tRange;
	};

	vector<Combination> combinations = {
		{ "exp1", "2500-2rp-ob-2-10-0.05", "t_max=0.07" },
		{ "exp3", "2500-2rp-ob-2-10-0.05", "t_max=0.15" },
	};

	printf("input directory: %s\n", inputDir.c_str());

	for (const auto &c : combinations)
	{
		string dir = inputDir + "/" + c.binning + "/" + c.expN + "," + c.tRange;

		string file = dir + "/fit.out";
		Results r;
		r.Load(file);

		file = dir + "/fit.root";
		TFile *f_in = TFile::Open(file.c_str());
		TGraph *g_fit_data = (TGraph *) f_in->Get("g_fit_data");

		printf("%s, %s, %s: ", c.expN.c_str(), c.binning.c_str(), c.tRange.c_str());
		PrintOne(r, g_fit_data);
		printf("\n");

		delete f_in;
	}
	return 0;
}
