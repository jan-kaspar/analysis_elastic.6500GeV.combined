import root;
import pad_layout;

string f = "do_fit.root";

xSizeDef = 10cm;
ySizeDef = 8cm;

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

void PlotRelDiff(RootObject hist, RootObject fit, pen p)
{
	int n_bins = hist.iExec("GetNbinsX");
	for (int bi = 1; bi <= n_bins; ++bi)
	{
		real x = hist.rExec("GetBinCenter", bi);
		real w = hist.rExec("GetBinWidth", bi);
		real y = hist.rExec("GetBinContent", bi);
		real y_unc = hist.rExec("GetBinError", bi);

		real y_fit = fit.rExec("Eval", x);

		real y_rel = (y - y_fit) / y_fit;
		real y_rel_unc = y_unc / y_fit;

		draw((x-w/2, y_rel)--(x+w/2, y_rel), p);
		draw((x, y_rel-y_rel_unc)--(x, y_rel+y_rel_unc), p);
	}
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
currentpad.xTicks = LeftTicks(0.005, 0.001);

AddToLegend("<{\bf elastic scattering}, $\sqrt{s} = 13\un{TeV}$");

AddToLegend("<{\it data}:");
AddToLegend("TOTEM preliminary", mCi+2pt+(black+1pt));
AddToLegend("$\be^*=2500\un{m}$");

AddToLegend("<{\it fit components}:");

draw(RootGetObject(f, "g_fit_C"), "l", blue+1pt, "Coulomb only");
draw(RootGetObject(f, "g_fit_H"), "l", red+1pt, "hadronic only");
draw(RootGetObject(f, "g_fit_CH"), "l", heavygreen+1pt, "Coulomb $\oplus$ hadronic");

draw(RootGetObject(f, "g_input_data0"), "p", black+1pt, mCi+1pt);

limits((0, 400), (0.02, 900), Crop);

AttachLegend(NE, NE);
GShipout(hSkip=3mm, vSkip=0mm);
