import root;
import pad_layout;

string topDir = "../";

string f = topDir + "exploration.root";

string dirs[];
pen d_pens[];

dirs.push("a=1.88E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.41,KL"); d_pens.push(black);
dirs.push("a=1.88E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.43,KL"); d_pens.push(red);
dirs.push("a=1.88E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.45,KL"); d_pens.push(blue);
dirs.push("a=1.88E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.47,KL"); d_pens.push(heavygreen);
dirs.push("a=1.88E+09,b1=10.20,b2=0.00,b3=0.00,constant,p0=+1.49,KL"); d_pens.push(magenta);

xSizeDef = 10cm;
xTicksDef = LeftTicks(0.005, 0.001);

real t_min = 8e-4;

//----------------------------------------------------------------------------------------------------

RootObject ref_obj;
real ref_eta;

void DrawWithRef(RootObject obj, pen p)
{
	int N = obj.iExec("GetN");

	guide g;

	for (int idx = 0; idx < N; ++idx)
	{
		real ax[] = { 0. };
		real ay[] = { 0. };
		obj.vExec("GetPoint", idx, ax, ay);

		real y = ay[0];

		real y_ref = ref_obj.rExec("Eval", ax[0]);

		real y_rel = y / y_ref - 1.;

		g = g--(ax[0], y_rel);
	}

	draw(g, p);
}

//----------------------------------------------------------------------------------------------------

RootObject ref_obj2;

void DrawWithRefSum(RootObject obj, pen p, string label="")
{
	int N = obj.iExec("GetN");

	guide g;

	for (int idx = 0; idx < N; ++idx)
	{
		real ax[] = { 0. };
		real ay[] = { 0. };
		obj.vExec("GetPoint", idx, ax, ay);

		real y = ay[0];

		real y_ref = (ref_obj.rExec("Eval", ax[0]) + ref_obj2.rExec("Eval", ax[0]));

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
	draw(RootGetObject(f, dirs[diri] + "/g_dsdt_H"), d_pens[diri]);
}

limits((0, 3e2), (t_max, 2e3), Crop);
yaxis(XEquals(t_min, false), dashed);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H}/\d t\ung{mb/GeV^2}$");
scale(Linear, Log);

for (int diri : dirs.keys)
{
	draw(RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), d_pens[diri]);
}

limits((0, 3e2), (t_max, 2e3), Crop);
yaxis(XEquals(t_min, false), dashed);

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

NewRow();

ref_obj = RootGetObject(f, dirs[0] + "/g_dsdt_H");

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H} / \hbox{ref} - 1$, ref = H");
currentpad.yTicks = RightTicks(0.05, 0.01);

for (int diri : dirs.keys)
{
	DrawWithRef(RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), d_pens[diri]);
}

limits((0, -0.05), (t_max, 0.05), Crop);
yaxis(XEquals(t_min, false), dashed);

//----------------------------------------------------------------------------------------------------

ref_obj = RootGetObject(f, dirs[0] + "/g_dsdt_H");
ref_obj2 = RootGetObject(f, dirs[0] + "/g_dsdt_C");

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H} / \hbox{ref} - 1$, ref = H + C");
currentpad.yTicks = RightTicks(0.05, 0.01);

for (int diri : dirs.keys)
{
	DrawWithRefSum(RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), d_pens[diri], dirs[diri]);
}

limits((0, -0.20), (t_max, 0.01), Crop);
yaxis(XEquals(t_min, false), dashed);

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

frame f_legend = BuildLegend();

NewPad(false);
add(f_legend);

GShipout(vSkip=1mm);
