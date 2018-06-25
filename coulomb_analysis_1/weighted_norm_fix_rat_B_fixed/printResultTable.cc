#include <vector>
#include <string>

#include "../classes.h"
#include "../command_line_tools.h"

#include "TFile.h"
#include "TGraph.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintOne(const Results &r, TGraph *g_data, char *buf)
{
	const double eta = g_data->GetY()[10];
	const double eta_unc = g_data->GetY()[11];
	const double n_points = g_data->GetY()[12];
	const double A_p = g_data->GetY()[13];
	const double si_tot = g_data->GetY()[14];

	sprintf(buf, "%.0f, %.2f, %.5f+-%.5f, %.3f, %.1f, %.1f", n_points, r.chi_sq_norm, r.rho, r.rho_e, eta, A_p, si_tot);
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	// TODO
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

	// define matrix
	vector< vector<string> > rows = {
		//{ "exp1", "2500-2rp-ob-1-20-0.05" },
		{ "exp1", "2500-2rp-ob-2-10-0.05" },
	};

	vector< vector<string> > cols = {
		{ "t_max=0.005" },
	};

	unsigned int width_row_id = 30;
	unsigned int width_col = 55;

	//--------------------

	char fs_row_id[10];
	sprintf(fs_row_id, "%%%us", width_row_id);

	char fs_col[10];
	sprintf(fs_col, "%%%us", width_col);

	// print col ID
	printf(fs_row_id, "");
	for (unsigned int ci = 0; ci < cols.size(); ci++)
	{
		string col_id;
		for (unsigned int cci = 0; cci < cols[ci].size(); cci++)
		{
			if (!col_id.empty())
				col_id += ","; 
			col_id += cols[ci][cci];
		}
		
		printf(fs_col, col_id.c_str());
	}
	printf("\n");

	for (unsigned int ri = 0; ri < rows.size(); ri++)
	{
		if (rows[ri].size() == 0)
		{
			printf("\n");
			continue;
		}

		// print row ID
		string row_id;
		for (unsigned int rci = 0; rci < rows[ri].size(); rci++)
		{
			if (!row_id.empty())
				row_id += ",";
			row_id += rows[ri][rci];
		}
		printf(fs_row_id, row_id.c_str());


		for (unsigned int ci = 0; ci < cols.size(); ci++)
		{
			string dir = inputDir + "/";

			dir += rows[ri][1] + "/" + rows[ri][0] + "," + cols[ci][0];
			string file = dir + "/fit.out";

			Results r;
			r.Load(file);

			file = dir + "/fit.root";
			TFile *f_out = TFile::Open(file.c_str());

			TGraph *g_fit_data = (TGraph *) f_out->Get("g_fit_data");

			char buf[100];
			PrintOne(r, g_fit_data, buf);
			printf(fs_col, buf);

			delete f_out;
		}

		printf("\n");
	}

	return 0;
}
