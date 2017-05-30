import root;
import pad_layout;

string f = "do_fit.root";

xSizeDef = 10cm;
ySizeDef = 8cm;

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

void PlotRelDiff(RootObject data, RootObject fit, pen p)
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

		draw((x, y_rel), mCi+2pt+p);
		draw((x, y_rel-y_rel_unc)--(x, y_rel+y_rel_unc), p);
	}
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
currentpad.xTicks = LeftTicks(0.005, 0.001);

draw(RootGetObject(f, "fit canvas|g_data0"), "p", black, mCi+2pt);
draw(RootGetObject(f, "fit canvas|g_fit_CH"), "l", red+1pt);

limits((0, 400), (0.02, 900), Crop);

AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
scale(Linear, Log);

draw(RootGetObject(f, "fit canvas|g_data0"), "p", black, mCi+1pt);
draw(RootGetObject(f, "fit canvas|g_fit_CH"), "l", red+1pt);

limits((0, 10), (0.2, 900), Crop);

AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------

/*
NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
scale(Linear, Log);
//currentpad.xTicks = LeftTicks(0.005, 0.001);

draw(RootGetObject(f, "h_in"), "eb", black+1pt);
draw(RootGetObject(f, "g_dsdt_CH"), "l", red+1pt);

limits((0, 1e-3), (1., 1e3), Crop);

AttachLegend(NE, NE);
*/

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

ySizeDef = 4cm;

RootObject data = RootGetObject(f, "fit canvas|g_data0");
RootObject fit = RootGetObject(f, "fit canvas|g_fit_CH");

NewRow();

NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
currentpad.xTicks = LeftTicks(0.005, 0.001);
PlotRelDiff(data, fit, blue);
limits((0, -0.04), (0.02, +0.04), Crop);

NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
PlotRelDiff(data, fit, blue);
limits((0, -0.04), (0.2, +0.04), Crop);

/*
NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
PlotRelDiff(hist, fit, blue);
limits((0, -0.5), (1, +0.5), Crop);

*/


GShipout(hSkip=3mm, vSkip=0mm);
