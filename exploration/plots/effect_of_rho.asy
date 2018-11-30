import root;
import pad_layout;

string topDir = "../";

string f = topDir + "exploration.root";

string dirs[], d_labels[];
pen d_pens[];

//dirs.push("a=1.94E+09,b1=10.39,b2=0.00,b3=0.00,constant,p0=+1.43,KL"); d_pens.push(red); d_labels.push("$\rh = 0.14$");
dirs.push("a=1.94E+09,b1=10.39,b2=0.00,b3=0.00,constant,p0=+1.47,KL"); d_pens.push(blue); d_labels.push("$\rh = 0.10$");
dirs.push("a=1.94E+09,b1=10.39,b2=0.00,b3=0.00,constant,p0=+1.48,KL"); d_pens.push(red); d_labels.push("$\rh = 0.09$");
dirs.push("a=1.94E+09,b1=10.39,b2=0.00,b3=0.00,constant,p0=+1.49,KL"); d_pens.push(heavygreen); d_labels.push("$\rh = 0.08$");
//dirs.push("a=1.94E+09,b1=10.39,b2=0.00,b3=0.00,constant,p0=+1.51,KL"); d_pens.push(heavygreen); d_labels.push("$\rh = 0.06$");

xSizeDef = 12cm;
xTicksDef = LeftTicks(0.005, 0.001);

real t_min = 8e-4;
real t_max = 0.02;

//----------------------------------------------------------------------------------------------------

RootObject ref_obj;
RootObject ref_obj2;

//----------------------------------------------------------------------------------------------------

void DrawWithRef(RootObject obj, pen p, string label="")
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

		//if (abs(ax[0] - 0.015) < 5e-4)
		//	write(y_rel);

		g = g--(ax[0], y_rel);
	}

	draw(g, p, label);
}


//----------------------------------------------------------------------------------------------------

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

		real y_ref = ref_obj.rExec("Eval", ax[0]) + ref_obj2.rExec("Eval", ax[0]);

		real y_rel = y / y_ref - 1.;

		//if (abs(ax[0] - 0.015) < 5e-4)
		//	write(y_rel);

		g = g--(ax[0], y_rel);
	}

	draw(g, p, label);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

ref_obj = RootGetObject(f, dirs[0] + "/g_dsdt_H");
ref_obj2 = RootGetObject(f, dirs[0] + "/g_dsdt_C");

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H} / \hbox{ref} - 1$, ref = H + C");
currentpad.yTicks = RightTicks(0.05, 0.01);

for (int diri : dirs.keys)
{
	DrawWithRefSum(RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), d_pens[diri], d_labels[diri]);
}

limits((0, -0.17), (t_max, 0.02), Crop);
yaxis(XEquals(t_min, false), dashed);
AttachLegend(E, E);

//----------------------------------------------------------------------------------------------------

NewRow();

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H} / \hbox{ref} - 1$, ref = H");
currentpad.yTicks = RightTicks(0.05, 0.01);

for (int diri : dirs.keys)
{
	DrawWithRef(RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), d_pens[diri], d_labels[diri]);
}

limits((0, -0.02), (t_max, 0.05), Crop);
yaxis(XEquals(t_min, false), dashed);
AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------

NewRow();

NewPad("$|t|\ung{GeV^2}$", "$\d\si^{\rm C+H} / \hbox{ref} - 1$, ref = H");
currentpad.yTicks = RightTicks(0.005, 0.001);

for (int diri : dirs.keys)
{
	DrawWithRef(RootGetObject(f, dirs[diri] + "/g_dsdt_CH"), d_pens[diri], d_labels[diri]);
}

limits((0, -0.02), (0.05, 0.0), Crop);
yaxis(XEquals(t_min, false), dashed);
AttachLegend(SE, SE);
