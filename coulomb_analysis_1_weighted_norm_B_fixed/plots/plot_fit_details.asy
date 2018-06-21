import root;
import pad_layout;

string topDir = "../";

string f = "fits/2500-2rp-ob-2-10-0.05/exp1,t_max=0.005/fit.root";


xSizeDef = 10cm;
ySizeDef = 8cm;

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
//----------------------------------------------------------------------------------------------------

RootObject data_ih = RootGetObject(topDir + f, "h_input_dataset_0");
RootObject data_unc_stat = RootGetObject(topDir + f, "g_data_coll_unc_stat");
RootObject data_unc_full = RootGetObject(topDir + f, "g_data_coll_unc_full");
RootObject fit = RootGetObject(topDir + f, "fit canvas|g_fit_CH");

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
currentpad.xTicks = LeftTicks(0.005, 0.001);

draw(data_ih, "eb", gray);
draw(data_unc_full, "p", heavygreen+squarecap+3pt);
draw(data_unc_stat, "p", blue+squarecap, mCi+blue+2pt);
draw(fit, "l", red+1pt);

limits((0, 400), (0.02, 900), Crop);

AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
currentpad.xTicks = LeftTicks(0.05, 0.01);
scale(Linear, Log);

draw(data_ih, "eb", gray);
draw(data_unc_full, "p", heavygreen+squarecap+3pt);
draw(data_unc_stat, "p", blue+squarecap, mCi+blue+1pt);
draw(fit, "l", red+1pt);

limits((0, 10), (0.2, 900), Crop);

AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
scale(Linear, Log);
//currentpad.xTicks = LeftTicks(0.005, 0.001);

draw(data_ih, "eb", gray);
draw(data_unc_full, "p", heavygreen+squarecap+3pt);
draw(data_unc_stat, "p", blue+squarecap, mCi+blue+1pt);
draw(fit, "l", red+1pt);

limits((0, 1e-3), (1., 1e3), Crop);

AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------

NewPad(false);
AddToLegend("points included in fit", mCi+3pt+blue);
AddToLegend("points excluded from fit", mCi+3pt+gray);
AddToLegend("statistical uncertainties", scale(0.05, 1)*mSq+7pt);
AddToLegend("stat~ + syst.~uncertainties", scale(0.4, 1)*mSq+7pt+heavygreen);
AddToLegend("fit", red);

AttachLegend();

//----------------------------------------------------------------------------------------------------

TH1_x_min = 0.0015;

ySizeDef = 4cm;

NewRow();

NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
currentpad.xTicks = LeftTicks(0.005, 0.001);
currentpad.yTicks = RightTicks(0.01, 0.005);
PlotRelDiff(data_ih, fit, gray);
PlotRelDiff(data_unc_full, fit, heavygreen+3pt);
PlotRelDiff(data_unc_stat, fit, blue, mCi+2pt+blue);
limits((0, -0.05), (0.02, +0.05), Crop);

NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
currentpad.xTicks = LeftTicks(0.05, 0.01);
currentpad.yTicks = RightTicks(0.01, 0.005);
PlotRelDiff(data_ih, fit, gray);
PlotRelDiff(data_unc_full, fit, heavygreen+2pt);
PlotRelDiff(data_unc_stat, fit, blue, mCi+1pt+blue);
limits((0, -0.05), (0.2, +0.05), Crop);

NewPad("$|t|\ung{GeV^2}$", "$(\hbox{data} - \hbox{fit}) / \hbox{fit}$");
PlotRelDiff(data_ih, fit, gray);
PlotRelDiff(data_unc_full, fit, heavygreen+2pt);
PlotRelDiff(data_unc_stat, fit, blue, mCi+1pt+blue);
limits((0, -0.5), (1, +0.5), Crop);


GShipout(hSkip=3mm, vSkip=0mm);
