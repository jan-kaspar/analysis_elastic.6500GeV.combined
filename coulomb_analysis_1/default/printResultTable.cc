#include <vector>
#include <string>

#include "../classes.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintOne(const Results &r, char *buf)
{
	sprintf(buf, "%.2f, %.3f+-%.3f, %.2f", r.chi_sq_norm, r.rho, r.rho_e, r.si_tot);
}

//----------------------------------------------------------------------------------------------------

int main()
{
	vector< vector<string> > rows = {
		{ "exp1", "2500-2rp-ob-1-20-0.05" },
		{ "exp1", "2500-2rp-ob-2-10-0.05" },
		{ "exp1", "2500-2rp-ob-3-5-0.05" },
		{},
		{ "exp2", "2500-2rp-ob-1-20-0.05" },
		{ "exp2", "2500-2rp-ob-2-10-0.05" },
		{ "exp2", "2500-2rp-ob-3-5-0.05" },
		{},
		{ "exp3", "2500-2rp-ob-1-20-0.05" },
		{ "exp3", "2500-2rp-ob-2-10-0.05" },
		{ "exp3", "2500-2rp-ob-3-5-0.05" },
	};

	vector< vector<string> > cols = {
		{ "t_max=0.07" },
		//{ "t_max=0.13" },
		{ "t_max=0.15" },
		//{ "t_max=0.17" },
	};

	unsigned int width_row_id = 35;
	unsigned int width_col = 30;

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
			string dir = "fits/";

			dir += rows[ri][1] + "/" + rows[ri][0] + "," + cols[ci][0];
			string file = dir + "/fit.out";

			Results r;
			r.Load(file);

			char buf[100];
			PrintOne(r, buf);
			printf(fs_col, buf);
		}

		printf("\n");
	}

	return 0;
}
