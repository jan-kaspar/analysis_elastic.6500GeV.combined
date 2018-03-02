import root;
import pad_layout;

string topDir = "../";

string f = topDir + "exploration3.root";

string dirs[], d_labels[];
pen d_pens[];

dirs.push("a=1.84E+09,b1=10.20,b2=4.40,b3=10.00,p0=+1.47,p_A=+0.00,p_ka=+2.00,p_tm=-0.30"); d_pens.push(red); d_labels.push("const");
dirs.push("a=1.84E+09,b1=10.20,b2=4.40,b3=10.00,p0=+1.47,p_A=+3.00,p_ka=+2.00,p_tm=-0.30"); d_pens.push(blue); d_labels.push("peripheral 1");
dirs.push("a=1.84E+09,b1=10.20,b2=4.40,b3=10.00,p0=+1.47,p_A=+3.00,p_ka=+2.00,p_tm=-0.20"); d_pens.push(heavygreen); d_labels.push("peripheral 2");
dirs.push("a=1.84E+09,b1=10.20,b2=4.40,b3=10.00,p0=+1.47,p_A=+3.00,p_ka=+1.50,p_tm=-0.30"); d_pens.push(magenta); d_labels.push("peripheral 3");

xSizeDef = 8cm;
//xTicksDef = LeftTicks(0.005, 0.001);

real t_min = 8e-4;
real t_max = 0.3;

//----------------------------------------------------------------------------------------------------

RootObject ref_obj;
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

NewPad("$|t|\ung{GeV^2}$", "hadr.~amp.~phase");
currentpad.yTicks = RightTicks(1., 0.5);

for (int diri : dirs.keys)
{
	draw(RootGetObject(f, dirs[diri] + "/g_FH_Theta"), d_pens[diri], d_labels[diri]);
}

//limits((0, -0.15), (t_max, 0.05), Crop);
yaxis(XEquals(t_min, false), dashed);
//AttachLegend(E, E);

//----------------------------------------------------------------------------------------------------

ref_obj = RootGetObject(f, dirs[0] + "/g_dsdt_H");
ref_obj2 = RootGetObject(f, dirs[0] + "/g_dsdt_C");

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H} / \hbox{ref} - 1$, ref = H + C");
currentpad.yTicks = RightTicks(0.05, 0.01);

for (int diri : dirs.keys)
{
	DrawWithRefSum(RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), d_pens[diri], d_labels[diri]);
}

limits((0, -0.15), (t_max, 0.05), Crop);
//yaxis(XEquals(t_min, false), dashed);
AttachLegend(E, E);
