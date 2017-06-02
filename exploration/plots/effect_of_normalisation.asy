import root;
import pad_layout;

string topDir = "../";

string f = topDir + "exploration.root";

string dirs[];
real d_etas[];
pen d_pens[];
string d_labels[];

/*
dirs.push("a=1.84E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.45,KL"); d_etas.push(1.2); d_pens.push(blue);
dirs.push("a=2.21E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.45,KL"); d_etas.push(1.0); d_pens.push(red+dashed);
dirs.push("a=2.21E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.41,KL"); d_etas.push(1.0); d_pens.push(heavygreen+longdashed);
*/

dirs.push("a=1.84E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.45,KL"); d_etas.push(1.04); d_pens.push(blue); d_labels.push("$1.04\ {\cal A}^{\rm C+H}(a_0, b_0, \rh = 0.12)$");
dirs.push("a=1.91E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.45,KL"); d_etas.push(1.0); d_pens.push(red+dashed); d_labels.push("${\cal A}^{\rm C+H}(1.04\ a_0, b_0, \rh = 0.12)$");
dirs.push("a=1.88E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.43,KL"); d_etas.push(1.02); d_pens.push(heavygreen+longdashed); d_labels.push("$1.02\ {\cal A}^{C+H}(1.02\ a_0, b_0, \rh = 0.14)$");

xSizeDef = 10cm;
xTicksDef = LeftTicks(0.005, 0.001);

//----------------------------------------------------------------------------------------------------

RootObject ref_obj;
real ref_eta;

void DrawWithRef(RootObject obj, real eta, pen p)
{
	int N = obj.iExec("GetN");

	guide g;

	for (int idx = 0; idx < N; ++idx)
	{
		real ax[] = { 0. };
		real ay[] = { 0. };
		obj.vExec("GetPoint", idx, ax, ay);

		real y = ay[0] * eta * eta;

		real y_ref = ref_obj.rExec("Eval", ax[0]) * ref_eta * ref_eta;

		real y_rel = y / y_ref - 1.;

		g = g--(ax[0], y_rel);
	}

	draw(g, p);
}

//----------------------------------------------------------------------------------------------------

RootObject ref_obj2;

void DrawWithRefSum(RootObject obj, real eta, pen p, string label="")
{
	int N = obj.iExec("GetN");

	guide g;

	for (int idx = 0; idx < N; ++idx)
	{
		real ax[] = { 0. };
		real ay[] = { 0. };
		obj.vExec("GetPoint", idx, ax, ay);

		real y = ay[0] * eta * eta;

		real y_ref = (ref_obj.rExec("Eval", ax[0]) + ref_obj2.rExec("Eval", ax[0])) * ref_eta * ref_eta;

		real y_rel = y / y_ref - 1.;

		g = g--(ax[0], y_rel);
	}

	draw(g, p, label);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

real t_max = 0.02;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm H}/\d t\ung{mb/GeV^2}$");
scale(Linear, Log);

for (int diri : dirs.keys)
{
	real eta = d_etas[diri];
	draw(shift(0, log10(eta*eta)), RootGetObject(f, dirs[diri] + "/g_dsdt_H"), d_pens[diri]);
}

limits((0, 3e2), (t_max, 2e3), Crop);
yaxis(XEquals(8e-4, false), black);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H}/\d t\ung{mb/GeV^2}$");
scale(Linear, Log);

for (int diri : dirs.keys)
{
	real eta = d_etas[diri];
	draw(shift(0, log10(eta*eta)), RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), d_pens[diri]);
}

limits((0, 3e2), (t_max, 2e3), Crop);
yaxis(XEquals(8e-4, false), black);

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

NewRow();

ref_obj = RootGetObject(f, dirs[0] + "/g_dsdt_H");
ref_eta = d_etas[0];

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H} / \hbox{ref} - 1$, ref = H");
currentpad.yTicks = RightTicks(0.05, 0.01);

for (int diri : dirs.keys)
{
	real eta = d_etas[diri];
	DrawWithRef(RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), eta, d_pens[diri]);
}

limits((0, -0.05), (t_max, 0.05), Crop);
yaxis(XEquals(8e-4, false), black);

//----------------------------------------------------------------------------------------------------

ref_obj = RootGetObject(f, dirs[0] + "/g_dsdt_H");
ref_obj2 = RootGetObject(f, dirs[0] + "/g_dsdt_C");
ref_eta = d_etas[0];

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H} / \hbox{ref} - 1$, ref = H + C");
currentpad.yTicks = RightTicks(0.05, 0.01);

for (int diri : dirs.keys)
{
	real eta = d_etas[diri];
	DrawWithRefSum(RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), eta, d_pens[diri], d_labels[diri]);
}

limits((0, -0.20), (t_max, 0.01), Crop);
yaxis(XEquals(8e-4, false), black);

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

frame f_legend = BuildLegend();

NewPad(false);
add(f_legend);

GShipout(vSkip=1mm);
