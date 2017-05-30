import root;
import pad_layout;

string f = "do_fit_1.root";

string variant = "variant 2";

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

draw(RootGetObject(f, variant + "/h_in"), "eb", black+1pt);
//draw(RootGetObject(f, variant + "/g_fit"), "p", heavygreen+1pt);
draw(RootGetObject(f, variant + "/g_dsdt_CH"), "l", red+1pt);

limits((0, 400), (0.02, 900), Crop);

AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
scale(Linear, Log);

draw(RootGetObject(f, variant + "/h_in"), "eb", black+1pt);
draw(RootGetObject(f, variant + "/g_dsdt_CH"), "l", red+1pt);

limits((0, 10), (0.2, 900), Crop);

AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
scale(Linear, Log);
//currentpad.xTicks = LeftTicks(0.005, 0.001);

draw(RootGetObject(f, variant + "/h_in"), "eb", black+1pt);
draw(RootGetObject(f, variant + "/g_dsdt_CH"), "l", red+1pt);

limits((0, 1e-3), (1., 1e3), Crop);

AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

ySizeDef = 4cm;

RootObject hist = RootGetObject(f, variant + "/h_in");
RootObject fit = RootGetObject(f, variant + "/g_dsdt_CH");

NewRow();

NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
currentpad.xTicks = LeftTicks(0.005, 0.001);
PlotRelDiff(hist, fit, blue);
limits((0, -0.04), (0.02, +0.04), Crop);

NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
PlotRelDiff(hist, fit, blue);
limits((0, -0.04), (0.2, +0.04), Crop);

NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
PlotRelDiff(hist, fit, blue);
limits((0, -0.5), (1, +0.5), Crop);


GShipout(hSkip=3mm, vSkip=0mm);
