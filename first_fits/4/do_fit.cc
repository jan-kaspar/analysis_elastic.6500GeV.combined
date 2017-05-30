#include "include/Elegent/Constants.h"
#include "include/Elegent/CoulombInterference.h"

#include "../HadronicFitModel.h"

#include "classes.h"
#include "common.h"

#include "MethodSimpleFit.h"

#include <vector>
#include <string>

#include "TFile.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	// TODO
	//printf("USAGE: normalize_and_fit <input file> <histogram path> <normAdjustment> <output file>\n");
}

//----------------------------------------------------------------------------------------------------

int main(/*int argc, char **argv*/)
{
	// settings
	string inputDataSpec = "2500-2rp-ob-2-10-0.05";
	B_degree = 3;
	chosenCIMode = CoulombInterference::mKL;
	use_stat_unc = true;
	use_syst_unc = true;

	string rootOutputFileName = "do_fit.root";
	string resultOutputFileName = "do_fit.out";

	// prepare list of standard datasets
	InitStdDataSets();

	// compile input data list
	inputData.norm_unc_global = 0.;

	char buf[300];
	strcpy(buf, inputDataSpec.c_str());
	char *pch;
	pch = strtok(buf, ",");
	printf(">> input datasets:\n");
	while (pch != NULL)
	{
		printf ("\t%s\n", pch);

		vector<DataSetInfo>::iterator it = stdDataSets.begin();
		for (; it != stdDataSets.end(); ++it)
		{
			if (it->tag.compare(pch) == 0)
				break;
		}

		if (it == stdDataSets.end())
		{
			printf("ERROR: dataset `%s' not recognized.\n", pch);
			return 10;
		}

		inputData.dataSetInfo.push_back(*it);

		pch = strtok(NULL, ",");
	}

	// init hadronic model
	hfm = new HadronicFitModel();
	hfm->t1 = 0.2;
	hfm->t2 = 0.5;

	hfm->phaseMode = HadronicFitModel::pmConstant;

	hfm->modulusHighTVariant = 2;

	hfm->a = 1.84E9;
	hfm->b1 = 10.2;
	hfm->b2 = 0.;
	hfm->b3 = 0.;

	double rho_init = 0.12; // default value of rho
	hfm->p0 = M_PI/2. - atan(rho_init);

	model = hfm;

	// print settings
	printf(">> settings\n");
	printf("\tN_b = %u\n", B_degree);
	printf("\tchosenCIMode = %u\n", chosenCIMode);
	printf("\tphaseMode = %u\n", phaseMode);
	printf("\tuse_stat_unc = %i, use_syst_unc = %i\n", use_stat_unc, use_syst_unc);

	// select t range
	double t_min_data_coll, t_max_data_coll;
	t_min_data_coll = 1E-5;
	t_max_data_coll = 0.20;
	printf(">> t_min_data_coll = %.2E, t_max_data_coll = %.2E\n", t_min_data_coll, t_max_data_coll);

	// load data
	printf(">> loading data\n");

	for (unsigned int dsi = 0; dsi < inputData.dataSetInfo.size(); dsi++)
	{
		printf("\tdataset %u (%s)\n", dsi, inputData.dataSetInfo[dsi].tag.c_str());
		// data histogram
		TFile *f_in = new TFile(inputData.dataSetInfo[dsi].file_hist.c_str());
		if (f_in->IsZombie())
		{
			printf("ERROR: can't open file `%s'.\n", inputData.dataSetInfo[dsi].file_hist.c_str());
			return 1;
		}

		printf("\t\tdata file: %s\n", inputData.dataSetInfo[dsi].file_hist.c_str());

		TH1D *h = (TH1D *) f_in->Get(inputData.dataSetInfo[dsi].obj_hist.c_str());
		if (!h)
		{
			printf("ERROR: can't load object `%s' from file `%s'.\n", inputData.dataSetInfo[dsi].obj_hist.c_str(), inputData.dataSetInfo[dsi].file_hist.c_str());
			return 2;
		}

		printf("\t\tdata object: %s\n", inputData.dataSetInfo[dsi].obj_hist.c_str());
		
		inputData.dataSetInfo[dsi].hist = h;

		// systematic uncertainty matrix
		f_in = new TFile(inputData.dataSetInfo[dsi].file_unc.c_str());
		if (f_in->IsZombie())
		{
			printf("ERROR: can't open file `%s'.\n", inputData.dataSetInfo[dsi].file_unc.c_str());
			return 3;
		}

		printf("\t\tuncertainty file: %s\n", inputData.dataSetInfo[dsi].file_unc.c_str());

		TMatrixD *mat = (TMatrixD *) f_in->Get(inputData.dataSetInfo[dsi].obj_unc.c_str());
		if (!mat)
		{
			printf("ERROR: can't load object `%s' from file `%s'.\n", inputData.dataSetInfo[dsi].obj_unc.c_str(), inputData.dataSetInfo[dsi].file_unc.c_str());
			return 4;
		}

		printf("\t\tuncertainty object: %s\n", inputData.dataSetInfo[dsi].obj_unc.c_str());

		inputData.dataSetInfo[dsi].unc = mat;
	}

	printf(">> input data loaded succesfully\n");
	// ------------------------------ build measurement collection

	printf(">> building measurement collection\n");

	for (unsigned int dsi = 0; dsi < inputData.dataSetInfo.size(); dsi++)
	{
		unsigned int dataset_points = 0;
		TH1D *h = inputData.dataSetInfo[dsi].hist;
		for (int bi = 1; bi <= h->GetNbinsX(); bi++)
		{
			double c = h->GetBinCenter(bi);
			double l = h->GetBinLowEdge(bi);
			double r = l + h->GetBinWidth(bi);
			double v = h->GetBinContent(bi);
			double v_u = h->GetBinError(bi);

			if (v == 0.)
				continue;

			if (c < t_min_data_coll)
				continue;

			if (c > t_max_data_coll || c > inputData.dataSetInfo[dsi].t_max)
				continue;

			BinData bd;
			bd.dataset = dsi;
			bd.bin = bi;
			bd.x_left = l;
			bd.x_right = r;
			bd.x_repr = (l+r)/2.;
			bd.y = v;
			bd.y_stat_unc = v_u;

			data_coll.push_back(bd);
			dataset_points++;
		}

		printf("\tpoints from dataset %u : %u\n", dsi, dataset_points);
	}
	
	printf("\ttotal points: %lu points\n", data_coll.size());
	
	// ------------------------------ build systematic covariance matrix related to the measurement collection

	data_coll_rel_unc.ResizeTo(data_coll.size(), data_coll.size());

	for (unsigned int i = 0; i < data_coll.size(); i++)
	{
		for (unsigned int j = 0; j < data_coll.size(); j++)
		{
			const BinData &bi = data_coll[i];
			const BinData &bj = data_coll[j];

			double el_sum = 0.;

			// global normalisation
			double el_norm_g = inputData.norm_unc_global * inputData.norm_unc_global;
			el_sum += el_norm_g;

			if (bi.dataset == bj.dataset)
			{
				const DataSetInfo &dsi = inputData.dataSetInfo[bi.dataset];

				// dataset normalisation
				double el_norm_add = dsi.norm_unc_add * dsi.norm_unc_add;
				el_sum += el_norm_add;
			
				// dataset systematics (matrix elements indexed from 0, histogram bins from 1)
				double el_anal = (*dsi.unc)(bi.bin-1, bj.bin-1);
				el_sum += el_anal;

				if (i == j)
				{
					printf("idx = %3u | ds = %u, bin = %3i | x_repr = %.5f, y = %.3E | rel unc: stat = %2.1f%%, norm_g = %.2f%%, norm_add = %.2f%%, anal = %.2f%%\n",
						i, bi.dataset, bi.bin, bi.x_repr, bi.y,
						bi.y_stat_unc / bi.y * 100., sqrt(el_norm_g)*100., sqrt(el_norm_add)*100., sqrt(el_anal)*100.);
					/*
					printf("idx = %3u | ds = %u, bin = %3i | t_left = %.4f, t_repr = %.4f, t_right = %.4f\n",
						i, bi.dataset, bi.bin, bi.x_left, bi.x_repr, bi.x_right);
					*/
				}
			}

			data_coll_rel_unc(i, j) = el_sum;
		}
	}

	// initialize calculation engine
	Constants::Init(2*6500., cnts->mPP);
    cnts->Print();

	coulomb->mode = chosenCIMode;
	coulomb->ffType = coulomb->ffPuckett;
	coulomb->precision = 1E-3;
	coulomb->Print();

	// prepare output
	TFile *f_out = new TFile(rootOutputFileName.c_str(), "recreate");
	gDirectory = f_out;

	// enforce synchronisation of output buffers
	cout.sync_with_stdio(true);

	// print initial model parameters
	hfm->Print();

	// set fit parameter limits
	p0_lim_min = 1.32; p0_lim_max = 1.67;
	p_A_lim_min = 0.1; p_A_lim_max = 8.;
	p_ka_lim_min = 1.1; p_ka_lim_max = 50.;
	p_tm_lim_min = -0.8; p_tm_lim_max = -0.10;

	printf(">> default phase limits\n");
	printf("\tp0_lim_min=%.3f, p0_lim_max=%.3f\n", p0_lim_min, p0_lim_max);
	printf("\tp_A_lim_min=%.3f, p_A_lim_max=%.3f\n", p_A_lim_min, p_A_lim_max);
	printf("\tp_ka_lim_min=%.3f, p_ka_lim_max=%.3f\n", p_ka_lim_min, p_ka_lim_max);
	printf("\tp_tm_lim_min=%.3f, p_tm_lim_max=%.3f\n", p_tm_lim_min, p_tm_lim_max);

	// execute fit
	Results r;
	unsigned int ec = 100;
	
	string method_settings = "";
	ec = MethodSimpleFit::RunFit(method_settings, r);

	if (ec == 0)
	{
		printf("\n>> final results:\n");
		r.Print();
		r.Save(resultOutputFileName.c_str());
	} else {
		printf("\n>> ERROR in RunFit: code = %u\n", ec);
	}

	delete f_out;
	return 0;
}
