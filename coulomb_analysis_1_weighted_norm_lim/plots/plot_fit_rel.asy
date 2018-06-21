import root;
import pad_layout;

string topDir = "../";

string f = topDir + "fits/2500-2rp-ob-2-10-0.05/exp3,t_max=0.15/fit.root";
//string f = topDir + "fits/2500-2rp-ob-2-10-0.05/exp1,t_max=0.07/fit.root";

string it_dir = "F_C+H, iteration 2";

string unc_file = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV/beta2500/2rp/systematics/matrix.root";
string unc_types[], unc_labels[];
pen unc_pens[];
unc_types.push("all"); unc_pens.push(yellow+opacity(0.5)); unc_labels.push("full syst.~unc.~band");
unc_types.push("all-but-norm"); unc_pens.push(heavygreen); unc_labels.push("syst.~unc.~without normalisation");

string binning = "ob-2-10-0.05";

TGraph_errorBar = None;

xSizeDef = 15cm;
ySizeDef = 8cm;

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

void MakeRelativePlot(string f, string binning, real t_max)
{
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
}

//----------------------------------------------------------------------------------------------------

void MakeRelativePlot_C(string f, string binning, real t_max)
{
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

MakeRelativePlot(f, binning, 0.25);

NewRow();

MakeRelativePlot_C(f, binning, 0.25);
