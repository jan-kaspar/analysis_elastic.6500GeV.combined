#include "include/Elegent/Constants.h"
#include "include/Elegent/CoulombInterference.h"

#include "../HadronicFitModel.h"
#include "../command_line_tools.h"
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
	printf("USAGE: do_fit [option] [option] ...\n");
	printf("OPTIONS:\n");
	printf("    -input <string>\n");

	printf("    -t-min <float>\n");
	printf("    -t-max <float>\n");
	printf("    -binc <integer>\n");
	printf("    -bins <string>							comma separated list of bin indeces (counted from 1)\n");

	printf("    -b-degree <int>\n");
	printf("    -htv <int>\n");
	printf("    -use-stat-unc <bool>\n");
	printf("    -use-syst-unc <bool>\n");
	printf("    -use-norm-unc <bool>\n");

	printf("    -cni-formula <string>                    choice of CNI formula (KL or SWY)\n");

	printf("    -reweight-low-t-points <bool>            whether low-|t| points should be reweighted (default)\n");
	printf("    -reweight-low-t-points-meth1 <bool>      whether low-|t| points should be reweighted (method 1)\n");
	printf("    -reweight-low-t-points-meth2 <bool>      whether low-|t| points should be reweighted (method 2)\n");

	printf("    -use-normalisation-fit-parameter\n");
	printf("    -use-normalisation-chisq-term\n");
	printf("    -use-normalisation-limit\n");

	printf("    -use-normalisation-from-a\n");
	printf("    -Ap-value\n");

	printf("    -use-eta-fixed\n");
	printf("    -eta-value\n");

	printf("    -use-a-fixed\n");
	printf("    -a-value\n");

	printf("    -use-b1-fixed\n");
	printf("    -b1-value\n");

	printf("    -use-b2-fixed\n");
	printf("    -b2-value\n");

	printf("    -use-b3-fixed\n");
	printf("    -b3-value\n");

	printf("    -use-rho-fixed\n");
	printf("    -rho-value\n");

	printf("    -output <file>\n");
	printf("    -results <file>\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string inputDataSpec = "2500-2rp-ob-3-5-0.05";
	B_degree = 3;
	MethodSimpleFit::cniFormula = "KL";
	use_stat_unc = true;
	use_syst_unc = true;
	bool use_norm_unc = false;

	MethodSimpleFit::t_min_data_coll = 8E-4;
	MethodSimpleFit::t_max_data_coll = 0.17;
	
	string binsSpec = "";
	unsigned int bi_increment = 1;

	double t_tr1 = 0.2;
	double t_tr2 = 0.4;

	unsigned int modulusHighTVariant = 2;

	bool reweight_low_t_points = false;
	bool reweight_low_t_points_meth1 = false;
	bool reweight_low_t_points_meth2 = false;

	MethodSimpleFit::useNormalisationFitParameter = false;
	MethodSimpleFit::useNormalisationChiSqTerm = false;
	MethodSimpleFit::useNormalisationLimit = false;
	MethodSimpleFit::useNormalisationFromA = false;
	MethodSimpleFit::A_p_value_fix = 0.;

	MethodSimpleFit::useEtaFixed = false;
	MethodSimpleFit::eta_value_fix = 1.;

	MethodSimpleFit::useAFixed = false;
	MethodSimpleFit::a_value_fix = 18.4;

	MethodSimpleFit::useB1Fixed = false;
	MethodSimpleFit::b1_value_fix = 10.2;

	MethodSimpleFit::useB2Fixed = false;
	MethodSimpleFit::b2_value_fix = 0.;

	MethodSimpleFit::useB3Fixed = false;
	MethodSimpleFit::b3_value_fix = 0.;

	MethodSimpleFit::useRhoFixed = false;
	MethodSimpleFit::rho_value_fix = 0.12;

	string rootOutputFileName = "do_fit.root";
	string resultOutputFileName = "do_fit.out";
	string texOutputFileName = "do_fit.tex";

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestStringParameter(argc, argv, argi, "-input", inputDataSpec)) continue;
		
		if (TestDoubleParameter(argc, argv, argi, "-t-min", MethodSimpleFit::t_min_data_coll)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-t-max", MethodSimpleFit::t_max_data_coll)) continue;

		if (TestUIntParameter(argc, argv, argi, "-binc", bi_increment)) continue;
		if (TestStringParameter(argc, argv, argi, "-bins", binsSpec)) continue;

		if (TestUIntParameter(argc, argv, argi, "-b-degree", B_degree)) continue;

		if (TestUIntParameter(argc, argv, argi, "-htv", modulusHighTVariant)) continue;

		if (TestDoubleParameter(argc, argv, argi, "-t-tr1", t_tr1)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-t-tr2", t_tr2)) continue;

		if (TestBoolParameter(argc, argv, argi, "-use-stat-unc", use_stat_unc)) continue;
		if (TestBoolParameter(argc, argv, argi, "-use-syst-unc", use_syst_unc)) continue;
		if (TestBoolParameter(argc, argv, argi, "-use-norm-unc", use_norm_unc)) continue;

		if (TestStringParameter(argc, argv, argi, "-cni-formula", MethodSimpleFit::cniFormula)) continue;

		if (TestBoolParameter(argc, argv, argi, "-reweight-low-t-points", reweight_low_t_points)) continue;
		if (TestBoolParameter(argc, argv, argi, "-reweight-low-t-points-meth1", reweight_low_t_points_meth1)) continue;
		if (TestBoolParameter(argc, argv, argi, "-reweight-low-t-points-meth2", reweight_low_t_points_meth2)) continue;

		if (TestBoolParameter(argc, argv, argi, "-use-normalisation-fit-parameter", MethodSimpleFit::useNormalisationFitParameter)) continue;
		if (TestBoolParameter(argc, argv, argi, "-use-normalisation-chisq-term", MethodSimpleFit::useNormalisationChiSqTerm)) continue;
		if (TestBoolParameter(argc, argv, argi, "-use-normalisation-limit", MethodSimpleFit::useNormalisationLimit)) continue;
		if (TestBoolParameter(argc, argv, argi, "-use-normalisation-from-a", MethodSimpleFit::useNormalisationFromA)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-Ap-value", MethodSimpleFit::A_p_value_fix)) continue;
		
		if (TestBoolParameter(argc, argv, argi, "-use-eta-fixed", MethodSimpleFit::useEtaFixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-eta-value", MethodSimpleFit::eta_value_fix)) continue;
		
		if (TestBoolParameter(argc, argv, argi, "-use-a-fixed", MethodSimpleFit::useAFixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-a-value", MethodSimpleFit::a_value_fix)) continue;
		
		if (TestBoolParameter(argc, argv, argi, "-use-b1-fixed", MethodSimpleFit::useB1Fixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-b1-value", MethodSimpleFit::b1_value_fix)) continue;
		
		if (TestBoolParameter(argc, argv, argi, "-use-b2-fixed", MethodSimpleFit::useB2Fixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-b2-value", MethodSimpleFit::b2_value_fix)) continue;
		
		if (TestBoolParameter(argc, argv, argi, "-use-b3-fixed", MethodSimpleFit::useB3Fixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-b3-value", MethodSimpleFit::b3_value_fix)) continue;
		
		if (TestBoolParameter(argc, argv, argi, "-use-rho-fixed", MethodSimpleFit::useRhoFixed)) continue;
		if (TestDoubleParameter(argc, argv, argi, "-rho-value", MethodSimpleFit::rho_value_fix)) continue;

		if (TestStringParameter(argc, argv, argi, "-output", rootOutputFileName)) continue;
		if (TestStringParameter(argc, argv, argi, "-results", resultOutputFileName)) continue;
		if (TestStringParameter(argc, argv, argi, "-tex", texOutputFileName)) continue;

		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// input validation
	if (MethodSimpleFit::useNormalisationFitParameter && MethodSimpleFit::useNormalisationFromA)
	{
		printf("ERROR: useNormalisationFitParameter and useNormalisationFromA cannot be both true.\n");
		return 2;
	}

	// prepare list of standard datasets
	InitStdDataSets();

	// compile input data list
	inputData.norm_unc_global = (use_norm_unc) ? 0.055 : 0.;

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
	hfm->t1 = t_tr1;
	hfm->t2 = t_tr2;

	hfm->phaseMode = HadronicFitModel::pmConstant; // default
	//hfm->phaseMode = HadronicFitModel::pmBailly; hfm->p_td = -0.53; // test
	//hfm->phaseMode = HadronicFitModel::pmStandard; hfm->p_t0 = -0.50; hfm->p_tau = 0.1; // test

	hfm->modulusHighTVariant = modulusHighTVariant;

	hfm->hts = sqrt(1.);
	//hfm->hts *= 1.05;	// simulates the uncertainty of the high-|t| part

	string init_point_desc = "default";
	hfm->a = MethodSimpleFit::a_value_fix * 1E8;
	hfm->b1 = MethodSimpleFit::b1_value_fix;
	hfm->b2 = MethodSimpleFit::b2_value_fix;
	hfm->b3 = MethodSimpleFit::b3_value_fix;
	hfm->p0 = M_PI/2. - atan(MethodSimpleFit::rho_value_fix);

	//init_point_desc = "centre"; hfm->a = 1.84E9; hfm->b1 = 10.2; hfm->b2 = 0.; hfm->b3 = 0.; hfm->p0 = M_PI/2. - atan(0.12);
	//init_point_desc = "test 1"; hfm->a = 1.84E9; hfm->b1 = 10.2; hfm->b2 = 0.; hfm->b3 = 0.; hfm->p0 = M_PI/2. - atan(0.06);
	//init_point_desc = "test 2"; hfm->a = 1.84E9; hfm->b1 = 9.9; hfm->b2 = 0.; hfm->b3 = 0.; hfm->p0 = M_PI/2. - atan(0.12);
	//init_point_desc = "test 3"; hfm->a = 1.70E9; hfm->b1 = 10.2; hfm->b2 = 0.; hfm->b3 = 0.; hfm->p0 = M_PI/2. - atan(0.12);
	printf(">> initial point: %s\n", init_point_desc.c_str());

	model = hfm;

	// validate CNI mode
	if (MethodSimpleFit::cniFormula == "KL")
		chosenCIMode = CoulombInterference::mKL;
	else if (MethodSimpleFit::cniFormula == "SWY")
		chosenCIMode = CoulombInterference::mSWY;
	else if (MethodSimpleFit::cniFormula == "PH")
		chosenCIMode = CoulombInterference::mPH;
	else {
		printf("ERROR: unknown CNI formula '%s'.\n", MethodSimpleFit::cniFormula.c_str());
		return 1;
	}

	if (chosenCIMode == CoulombInterference::mPH && !MethodSimpleFit::useRhoFixed)
	{
		printf("ERROR: fits with PH formula should have fixed rho.\n");
		return 2;
	}

	// print settings
	printf(">> settings\n");
	printf("    N_b = %u\n", B_degree);
	printf("    chosenCIMode = %u\n", chosenCIMode);
	printf("    phaseMode = %u\n", phaseMode);
	printf("    use_stat_unc = %i, use_syst_unc = %i\n", use_stat_unc, use_syst_unc);
	printf("    t_min_data_coll = %.2E, t_max_data_coll = %.2E\n", MethodSimpleFit::t_min_data_coll, MethodSimpleFit::t_max_data_coll);
	printf("    bi_increment = %u\n", bi_increment);

	// prepare output
	TFile *f_out = new TFile(rootOutputFileName.c_str(), "recreate");

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

		gDirectory = f_out;
		char buf[100];
		sprintf(buf, "h_input_dataset_%u", dsi);
		h->Write(buf);

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

		// TODO
		/*
		vector<int> binIndeces;
		if (binsSpec == "")
		{
		} else {
		}
		*/

		int bi_start_inc = 8;

		for (int bi = 1; bi <= h->GetNbinsX(); bi += (bi >= bi_start_inc) ? bi_increment : 1)
		{
			double c = h->GetBinCenter(bi);
			double l = h->GetBinLowEdge(bi);
			double r = l + h->GetBinWidth(bi);

			const double scale_correction = 1.;
			//const double scale_correction = 31.1 / 31.0;

			double v = h->GetBinContent(bi) * scale_correction;
			double v_u = h->GetBinError(bi) * scale_correction;

			if (v == 0.)
				continue;

			if (c < MethodSimpleFit::t_min_data_coll)
				continue;

			if (c > MethodSimpleFit::t_max_data_coll || c > inputData.dataSetInfo[dsi].t_max)
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

				/*
				if (i == j)
				{
					const double unc_rel_tot = sqrt(pow(bi.y_stat_unc / bi.y, 2.) + el_sum);

					printf("idx = %3u | ds = %u, bin = %3i | x_repr = %.5f, y = %.3E | rel unc: stat = %4.1f%%, norm_g = %.2f%%, norm_add = %.2f%%, anal = %.2f%%, tot = %.2f%%\n",
						i, bi.dataset, bi.bin, bi.x_repr, bi.y,
						bi.y_stat_unc / bi.y * 100., sqrt(el_norm_g)*100., sqrt(el_norm_add)*100., sqrt(el_anal)*100., unc_rel_tot*100.);
				}
				*/
			}

			double sc = 1.;

			if (reweight_low_t_points_meth2)
			{
				if (bi.bin >= 5 && bi.bin <= 6) sc *= 0.10;
				if (bj.bin >= 5 && bj.bin <= 6) sc *= 0.10;

				if (bi.bin >= 7 && bi.bin <= 8) sc *= 0.10;
				if (bj.bin >= 7 && bj.bin <= 8) sc *= 0.10;
			}

			data_coll_rel_unc(i, j) = el_sum * sc;
		}
	}

	// ------------------------------ apply re-weighting

	if (reweight_low_t_points)
	{
		for (unsigned int i = 0; i < data_coll.size(); i++)
		{
			//BinData &dp = data_coll[i];

			// default: first two points, stat_unc reduced to 10%
			/*
			if (dp.x_repr < 1.1E-3)
				dp.y_stat_unc *= 0.10;
			*/

			/*
			if (dp.bin >= 5 && dp.bin <= 8)
			{
				const double unc_syst_rel = sqrt(data_coll_rel_unc(i, i));
				const double unc_tot_rel = 0.02;
				const double unc_stat_rel = sqrt(unc_tot_rel*unc_tot_rel - unc_syst_rel*unc_syst_rel);
				dp.y_stat_unc = unc_stat_rel * dp.y;
			}
			*/

			/*
			if (dp.bin >= 5 && dp.bin <= 8)
			{
				const double unc_syst = sqrt(data_coll_rel_unc(i, i)) * dp.y;
				const double unc_tot = 12.; // mb/GeV^2
				const double unc_stat = sqrt(unc_tot*unc_tot - unc_syst*unc_syst);
				dp.y_stat_unc = unc_stat;
			}
			*/

			/*
			if (dp.bin >= 5 && dp.bin <= 6)
				dp.y_stat_unc *= 0.10;

			if (dp.bin >= 7 && dp.bin <= 8)
				dp.y_stat_unc *= 0.10;
			*/

			/*
			if (dp.bin >= 7)
				dp.y_stat_unc *= 100.;
			*/

			/*
			if (dp.bin >= 5 && dp.bin <= 8)
			{
				double sen = 0.;

				// green: A_C / A_H
				//if (dp.bin == 5) sen = 0.724;
				//if (dp.bin == 6) sen = 0.602;
				//if (dp.bin == 7) sen = 0.515;
				//if (dp.bin == 8) sen = 0.444;
				//const double scale = 13.;

				// red: si_C / si_H
				//if (dp.bin == 5) sen = 0.524;
				//if (dp.bin == 6) sen = 0.362;
				//if (dp.bin == 7) sen = 0.265;
				//if (dp.bin == 8) sen = 0.197;
				//const double scale = 10.;

				// red squared: (si_C / si_H)^2
				if (dp.bin == 5) sen = 0.2746; // 0.524*0.524
				if (dp.bin == 6) sen = 0.1310; // 0.362*0.362
				if (dp.bin == 7) sen = 0.0702; // 0.265*0.265
				if (dp.bin == 8) sen = 0.0388; // 0.197*0.197
				const double scale = 7.;

				const double unc_syst = sqrt(data_coll_rel_unc(i, i)) * dp.y;
				const double unc_tot = scale / sqrt(sen);
				const double unc_stat = sqrt(unc_tot*unc_tot - unc_syst*unc_syst);

				dp.y_stat_unc = unc_stat;

				const double unc_dc_tot = sqrt(unc_syst*unc_syst + unc_stat*unc_stat);
				printf("* %u, %.3f, %.3f\n", dp.bin, unc_dc_tot, 1./unc_dc_tot/unc_dc_tot*1e3 / 0.796 * 0.039);
			}
			*/

			/*
			if (dp.bin >= 5 && dp.bin <= 8)
			{
				double sen = 0.;

				// red squared: (si_C / si_H)^2
				if (dp.bin == 5) sen = 0.2746; // 0.524*0.524
				if (dp.bin == 6) sen = 0.1310; // 0.362*0.362
				if (dp.bin == 7) sen = 0.0702; // 0.265*0.265
				if (dp.bin == 8) sen = 0.0388; // 0.197*0.197
				const double scale = 3.3;

				const double unc_syst = sqrt(data_coll_rel_unc(i, i)) * dp.y;
				const double unc_tot = sqrt(unc_syst*unc_syst + dp.y_stat_unc*dp.y_stat_unc);

				const double unc_w = scale / sqrt(sen);

				const double unc_comb = sqrt(unc_tot*unc_tot + unc_w*unc_w);

				const double unc_stat = sqrt(unc_comb*unc_comb - unc_syst*unc_syst);

				dp.y_stat_unc = unc_stat;

				printf("* bin = %u, unc_tot = %.3f, unc_w = %.3f, unc_comb = %.3f \n", dp.bin, unc_tot, unc_w, unc_comb);
			}
			*/

			/*
			if (dp.x_repr < 1.6E-3)	// first ?? points
			if (dp.x_repr < 1.8E-3)	// first ?? points
				//dp.y_stat_unc *= 0.10;
				dp.y_stat_unc = 1.;
			*/

			/*
			// introduce bias
			if (dp.bin >= 5 && dp.bin <= 8)
			{
				const double unc_syst = 10. * sqrt(data_coll_rel_unc(i, i)) * dp.y;
				dp.y += unc_syst;
			}
			*/
		}
	}

	if (reweight_low_t_points_meth1)
	{
		for (unsigned int i = 0; i < data_coll.size(); i++)
		{
			BinData &dp = data_coll[i];

			if (dp.bin >= 5 && dp.bin <= 8)
			{
				double sen = 0.;

				// red squared: (si_C / si_H)^2
				if (dp.bin == 5) sen = 0.2746; // 0.524*0.524
				if (dp.bin == 6) sen = 0.1310; // 0.362*0.362
				if (dp.bin == 7) sen = 0.0702; // 0.265*0.265
				if (dp.bin == 8) sen = 0.0388; // 0.197*0.197
				const double scale = 7.;

				const double unc_syst = sqrt(data_coll_rel_unc(i, i)) * dp.y;
				const double unc_tot = scale / sqrt(sen);
				const double unc_stat = sqrt(unc_tot*unc_tot - unc_syst*unc_syst);

				dp.y_stat_unc = unc_stat;

				const double unc_dc_tot = sqrt(unc_syst*unc_syst + unc_stat*unc_stat);
				printf("* %u, %.3f, %.3f\n", dp.bin, unc_dc_tot, 1./unc_dc_tot/unc_dc_tot*1e3 / 0.796 * 0.039);
			}
		}
	}

	if (reweight_low_t_points_meth2)
	{
		for (unsigned int i = 0; i < data_coll.size(); i++)
		{
			BinData &dp = data_coll[i];

			if (dp.bin >= 5 && dp.bin <= 6)
				dp.y_stat_unc *= 0.10;

			if (dp.bin >= 7 && dp.bin <= 8)
				dp.y_stat_unc *= 0.10;
		}
	}


	// ------------------------------ print dataset summary

	for (unsigned int i = 0; i < data_coll.size(); i++)
	{
		const BinData &bi = data_coll[i];

		const double unc_rel_stat = bi.y_stat_unc / bi.y;
		const double unc_rel_syst = sqrt(data_coll_rel_unc(i, i));
		const double unc_rel_tot = sqrt(unc_rel_stat*unc_rel_stat + unc_rel_syst*unc_rel_syst);

		printf("idx = %3u | ds = %u, bin = %3i | x_repr = %.5f, y = %.3E | rel unc: stat = %5.2f%%, syst = %5.2f%%, tot = %5.2f%%\n",
			i, bi.dataset, bi.bin, bi.x_repr, bi.y,
			unc_rel_stat*100., unc_rel_syst*100., unc_rel_tot*100.
		);
	}

	// initialize calculation engine
	Constants::Init(2*6500., cnts->mPP);
    cnts->Print();

	coulomb->mode = chosenCIMode;
	coulomb->ffType = coulomb->ffPuckett;
	//coulomb->ffType = coulomb->ffArrington;
	//coulomb->ffType = coulomb->ffBorkowski;
	coulomb->precision = 1E-3;
	coulomb->Print();

	// prepare output
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

	// prepare ouput tex file
	FILE *f_out_tex = fopen(texOutputFileName.c_str(), "w");

	// execute fit
	Results r;
	unsigned int ec = 100;
	
	string method_settings = "";
	ec = MethodSimpleFit::RunFit(method_settings, r, f_out_tex);

	if (ec == 0)
	{
		printf("\n>> final results:\n");
		r.Print();
		r.Save(resultOutputFileName.c_str());
	} else {
		printf("\n>> ERROR in RunFit: code = %u\n", ec);
	}

	// clean up
	fclose(f_out_tex);
	delete f_out;

	return 0;
}
