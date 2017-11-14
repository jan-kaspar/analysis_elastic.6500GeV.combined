import root;
import pad_layout;

string f = "../fits/2500-2rp-ob-2-10-0.05/exp3,t_max=0.15/fit.root";

string dir = "F_C+H, iteration 2";

TGraph_errorBar = None;

xSizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

NewPad("eigenvector index", "$\de$");
draw(RootGetObject(f, dir + "/g_de after fit"), "p", blue, mCi+1pt+blue);

NewPad("eigenvector index", "contribution to $\ch^2$");
draw(RootGetObject(f, dir + "/g_cont after fit"), "p", red, mCi+1pt+red);

//----------------------------------------------------------------------------------------------------

NewRow();

NewPad("$|t|\ung{GeV^2}$", "eigenvector composition");

RootObject glc = RootGetObject(f, dir + "/g_largest_cont_idx after fit");
for (int evn = 0; evn < 4; ++evn)
{
	real ax[] = {0};
	real ay[] = {0};
	glc.vExec("GetPoint", evn, ax, ay);

	real idx = ax[0];
	string name = format("eigenvector %.0f", idx);
	pen p = StdPen(evn+1);
	draw(RootGetObject(f, dir + "/eigenvectors/" + name), "l,p", p, mCi+1pt+p, name);
}

frame f = BuildLegend("eigenvectors with largest contrib.~to $\ch^2$");

NewPad(false);
attach(shift(60, 0) * f);

GShipout(hSkip=3mm, vSkip=0mm);
