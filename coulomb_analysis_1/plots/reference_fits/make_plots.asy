import root;
import pad_layout;

string topDir = "../../";

string it_dir = "F_C+H, iteration 2";

string unc_file = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/systematics/matrix.root";
string unc_types[], unc_labels[];
pen unc_pens[];
unc_types.push("all"); unc_pens.push(yellow+opacity(0.5)); unc_labels.push("full syst.~unc.~band");
unc_types.push("all-but-norm"); unc_pens.push(heavygreen); unc_labels.push("syst.~unc.~without normalisation");

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

void PlotRelDiff(RootObject data, RootObject fit, pen p, marker m=nomarker)
{
	if (data.InheritsFrom("TGraph"))
	{
		int n_points = data.iExec("GetN");
		for (int i = 0; i < n_points; ++i)
		{
			real ax[] = {0.};
			real ay[] = {0.};
			data.vExec("GetPoint", i, ax, ay);
	
			real x = ax[0];
			real y = ay[0];
			real y_unc = data.rExec("GetErrorY", i);
	
			real y_fit = fit.rExec("Eval", x);
	
			real y_rel = (y - y_fit) / y_fit;
			real y_rel_unc = y_unc / y_fit;
	
			draw((x, y_rel), m);
			draw((x, y_rel-y_rel_unc)--(x, y_rel+y_rel_unc), p+squarecap);
		}
	}

	if (data.InheritsFrom("TH1"))
	{
		for (int bi = 1; bi <= data.iExec("GetNbinsX"); ++bi)
		{
			real x = data.rExec("GetBinCenter", bi);

			if (x < TH1_x_min || x > TH1_x_max)
				continue;

			real w = data.rExec("GetBinWidth", bi);
			real c = data.rExec("GetBinContent", bi);
			real u = data.rExec("GetBinError", bi);

			real y_fit = fit.rExec("Eval", x);
	
			real y_rel = (c - y_fit) / y_fit;
			real y_rel_unc = u / y_fit;
	
			//draw((x, y_rel), mCi+2pt+p);
			draw((x-w/2, y_rel)--(x+w/2, y_rel), p+squarecap);
			draw((x, y_rel-y_rel_unc)--(x, y_rel+y_rel_unc), p+squarecap);
		}
	}
}

//----------------------------------------------------------------------------------------------------

void DrawRelUncBand(RootObject bc, RootObject relUnc, RootObject ref, pen p)
{
	if (relUnc.InheritsFrom("TH1"))
	{
		guide g_u, g_b;

		for (int bi = 1; bi < relUnc.iExec("GetNbinsX"); ++bi)
		{
			real c = relUnc.rExec("GetBinCenter", bi);
			real w = relUnc.rExec("GetBinWidth", bi);
			real rel_unc = relUnc.rExec("GetBinContent", bi);

			real band_cen = bc.rExec("Eval", c);
			real y_ref = ref.rExec("Eval", c);

			real y_up = band_cen * (1. + rel_unc) / y_ref - 1.;
			real y_dw = band_cen * (1. - rel_unc) / y_ref - 1.;

			g_u = g_u -- Scale((c-w/2, y_up)) -- Scale((c+w/2, y_up));
			g_b = g_b -- Scale((c-w/2, y_dw)) -- Scale((c+w/2, y_dw));
		}

		g_b = reverse(g_b);
		filldraw(g_u--g_b--cycle, p, nullpen);
	}
}

//----------------------------------------------------------------------------------------------------

void MakeFitPlots(string f)
{
	xSizeDef = 10cm;
	ySizeDef = 8cm;

	RootObject data_ih = RootGetObject(f, "h_input_dataset_0");
	RootObject data_unc_stat = RootGetObject(f, "g_data_coll_unc_stat");
	RootObject data_unc_full = RootGetObject(f, "g_data_coll_unc_full");
	RootObject fit = RootGetObject(f, "fit canvas|g_fit_CH");

	//----------

	NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
	currentpad.xTicks = LeftTicks(0.005, 0.001);

	draw(data_ih, "eb", gray);
	draw(data_unc_full, "p", heavygreen+squarecap+3pt);
	draw(data_unc_stat, "p", blue+squarecap, mCi+blue+2pt);
	draw(fit, "l", red+1pt);

	limits((0, 400), (0.02, 900), Crop);

	AttachLegend(NE, NE);

	NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
	currentpad.xTicks = LeftTicks(0.05, 0.01);
	scale(Linear, Log);

	draw(data_ih, "eb", gray);
	draw(data_unc_full, "p", heavygreen+squarecap+3pt);
	draw(data_unc_stat, "p", blue+squarecap, mCi+blue+1pt);
	draw(fit, "l", red+1pt);

	limits((0, 10), (0.2, 900), Crop);

	AttachLegend(NE, NE);

	NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
	scale(Linear, Log);
	//currentpad.xTicks = LeftTicks(0.005, 0.001);

	draw(data_ih, "eb", gray);
	draw(data_unc_full, "p", heavygreen+squarecap+3pt);
	draw(data_unc_stat, "p", blue+squarecap, mCi+blue+1pt);
	draw(fit, "l", red+1pt);

	limits((0, 1e-3), (1., 1e3), Crop);

	AttachLegend(NE, NE);

	//----------

	NewPad(false);
	AddToLegend("points included in fit", mCi+3pt+blue);
	AddToLegend("points excluded from fit", mCi+3pt+gray);
	AddToLegend("statistical uncertainties", scale(0.05, 1)*mSq+7pt);
	AddToLegend("stat~ + syst.~uncertainties", scale(0.4, 1)*mSq+7pt+heavygreen);
	AddToLegend("fit", red);

	AttachLegend();

	//----------

	TH1_x_min = 0.005;

	ySizeDef = 4cm;

	NewRow();

	NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
	currentpad.xTicks = LeftTicks(0.005, 0.001);
	currentpad.yTicks = RightTicks(0.01, 0.005);
	PlotRelDiff(data_ih, fit, gray);
	PlotRelDiff(data_unc_full, fit, heavygreen+3pt);
	PlotRelDiff(data_unc_stat, fit, blue, mCi+2pt+blue);
	limits((0, -0.04), (0.02, +0.04), Crop);

	NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
	currentpad.xTicks = LeftTicks(0.05, 0.01);
	currentpad.yTicks = RightTicks(0.01, 0.005);
	PlotRelDiff(data_ih, fit, gray);
	PlotRelDiff(data_unc_full, fit, heavygreen+2pt);
	PlotRelDiff(data_unc_stat, fit, blue, mCi+1pt+blue);
	limits((0, -0.04), (0.2, +0.04), Crop);

	NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
	PlotRelDiff(data_ih, fit, gray);
	PlotRelDiff(data_unc_full, fit, heavygreen+2pt);
	PlotRelDiff(data_unc_stat, fit, blue, mCi+1pt+blue);
	limits((0, -0.5), (1, +0.5), Crop);
}

//----------------------------------------------------------------------------------------------------

void MakeRelativePlots(string f, string binning, real t_max = 0.25)
{
	xSizeDef = 10cm;
	ySizeDef = 8cm;

	//--------------------

	string ref_label = "\hbox{ref} = 633 \e^{-20.4\,|t|}";

	NewPad("$|t|\ung{GeV^2}$", "${\d\sigma/\d t - \hbox{ref}\over\hbox{ref}}\ ,\quad " + ref_label + "$");
	currentpad.xTicks = LeftTicks(0.05, 0.01);
	currentpad.yTicks = RightTicks(0.05, 0.01);

	RootObject fit = RootGetObject(f, "g_fit_CH");
	RootObject ref = RootGetObject(f, "g_ref");

	for (int ui : unc_types.keys)
	{
		RootObject relUnc = RootGetObject(unc_file, "matrices/" + unc_types[ui] + "/" + binning + "/h_stddev");

		DrawRelUncBand(fit, relUnc, ref, unc_pens[ui]);
		AddToLegend(unc_labels[ui], mSq+6pt+unc_pens[ui]);
	}

	draw(RootGetObject(f, "fit canvas, rel|g_fit_H_rel"), "l", blue+1pt, "fit, hadronic component only");
	draw(RootGetObject(f, "fit canvas, rel|g_fit_CH_rel"), "l", red+1pt, "fit, all components");

	draw(RootGetObject(f, "fit canvas, rel|g_data_rel0"), "p", black, mCi+black+1pt, "data with stat.~unc.");

	limits((0, -0.05), (t_max, 0.15), Crop);

	AttachLegend(shift(0, 0)*BuildLegend(NE, vSkip=-1mm), NE);

	//--------------------

	string ref_label = "\hbox{ref} = 633 \e^{-20.4\,|t|} + {\d\sigma^{\rm C}\over\d t}";

	NewPad("$|t|\ung{GeV^2}$", "${\d\sigma/\d t - \hbox{ref}\over\hbox{ref}}\ ,\quad " + ref_label + "$");
	currentpad.xTicks = LeftTicks(0.05, 0.01);
	currentpad.yTicks = RightTicks(0.05, 0.01);

	RootObject fit = RootGetObject(f, "g_fit_CH");
	RootObject ref = RootGetObject(f, "g_refC");

	for (int ui : unc_types.keys)
	{
		RootObject relUnc = RootGetObject(unc_file, "matrices/" + unc_types[ui] + "/" + binning + "/h_stddev");

		DrawRelUncBand(fit, relUnc, ref, unc_pens[ui]);
		AddToLegend(unc_labels[ui], mSq+6pt+unc_pens[ui]);
	}

	draw(RootGetObject(f, "fit canvas, relC|g_fit_H_relC"), "l", blue+1pt, "fit, hadronic component only");
	draw(RootGetObject(f, "fit canvas, relC|g_fit_CH_relC"), "l", red+1pt, "fit, all components");

	draw(RootGetObject(f, "fit canvas, relC|g_data_relC0"), "p", black, mCi+black+1pt, "data with stat.~unc.");

	limits((0, -0.15), (t_max, 0.05), Crop);

	AttachLegend(shift(0, 10)*BuildLegend(SE, vSkip=-1mm), SE);
}

//----------------------------------------------------------------------------------------------------

void MakeComponentPlots(string f)
{
	xSizeDef = 10cm;
	ySizeDef = 8cm;

	NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
	currentpad.xTicks = LeftTicks(0.005, 0.001);

	AddToLegend("<{\it data}:");
	AddToLegend("$\be^*=2500\un{m}$");

	AddToLegend("<{\it fit components}:");

	draw(RootGetObject(f, "g_fit_C"), "l", blue+1pt, "Coulomb only");
	draw(RootGetObject(f, "g_fit_H"), "l", red+1pt, "hadronic only");
	draw(RootGetObject(f, "g_fit_CH"), "l", heavygreen+1pt, "Coulomb $\oplus$ hadronic");

	draw(RootGetObject(f, "g_input_data0"), "p", black+1pt, mCi+1pt);

	limits((0, 400), (0.02, 900), Crop);

	AttachLegend();
}

//----------------------------------------------------------------------------------------------------

void MakeContributionPlots(string f)
{
	xSizeDef = 10cm;
	ySizeDef = 8cm;

	NewPad("eigenvector index", "$\de$");
	draw(RootGetObject(f, it_dir + "/g_de after fit"), "p", blue, mCi+1pt+blue);

	NewPad("eigenvector index", "contribution to $\ch^2$");
	draw(RootGetObject(f, it_dir + "/g_cont after fit"), "p", red, mCi+1pt+red);

	//----------

	NewRow();

	NewPad("$|t|\ung{GeV^2}$", "eigenvector composition");

	RootObject glc = RootGetObject(f, it_dir + "/g_largest_cont_idx after fit");
	for (int evn = 0; evn < 4; ++evn)
	{
		real ax[] = {0};
		real ay[] = {0};
		glc.vExec("GetPoint", evn, ax, ay);

		real idx = ax[0];
		string name = format("eigenvector %.0f", idx);
		pen p = StdPen(evn+1);
		draw(RootGetObject(f, it_dir + "/eigenvectors/" + name), "l,p", p, mCi+1pt+p, name);
	}

	AttachLegend("eigenvectors with largest contrib.~to $\ch^2$");

	/*
	frame f = BuildLegend("eigenvectors with largest contrib.~to $\ch^2$");

	NewPad(false);
	attach(shift(100, 0) * f);
	*/
}

//----------------------------------------------------------------------------------------------------

string expns[] = { "exp1", "exp2", "exp3" };

string binnings[] = { "ob-1-20-0.05", "ob-2-10-0.05", "ob-3-5-0.05" };

string tmaxs[] = { "0.07", "0.15" };

for (string expn : expns)
{
	for (string binning : binnings)
	{
		for (string tmax : tmaxs)
		{
			string f_in = topDir + "fits/2500-2rp-"+binning+"/"+expn+",t_max="+tmax+"/fit.root";
			string f_out = binning+"_"+expn+"_"+tmax+".pdf";

			write("* " + f_in);

			MakeFitPlots(f_in);
			NewPage();
			MakeRelativePlots(f_in, binning);
			NewPage();
			MakeComponentPlots(f_in);
			NewPage();
			MakeContributionPlots(f_in);

			GShipout(f_out);
		}
	}
}
