import root;
import pad_layout;

string f = "fit.root";

string dirs[];
pen d_pens[];
dirs.push("fit 1"); d_pens.push(red);
dirs.push("fit 2"); d_pens.push(blue);
dirs.push("fit 3"); d_pens.push(heavygreen);
dirs.push("fit 4"); d_pens.push(cyan);

//----------------------------------------------------------------------------------------------------

real a_ref = 633.;
real b_ref = 20.4;

void DrawRel(RootObject obj, pen p)
{
	if (obj.InheritsFrom("TGraph"))
	{
		guide g;

		for (int i = 0; i < obj.iExec("GetN"); ++i)
		{
			real ax[] = {0.};
			real ay[] = {0.};
			obj.vExec("GetPoint", i, ax, ay);

			real y_ref = a_ref * exp(-b_ref * ax[0]);

			real y_rel = (ay[0] - y_ref) / y_ref;
			
			g = g--Scale((ax[0], y_rel));
		}

		draw(g, p);
	}

	if (obj.InheritsFrom("TH1"))
	{
		for (int bi = 1; bi <= obj.iExec("GetNbinsX"); ++bi)
		{
			real c = obj.rExec("GetBinCenter", bi);
			real w = obj.rExec("GetBinWidth", bi);
			real v = obj.rExec("GetBinContent", bi);
			real u = obj.rExec("GetBinError", bi);

			real y_ref = a_ref * exp(-b_ref * c);
			real v_rel = (v - y_ref) / y_ref;
			real u_rel = u / y_ref;

			draw(Scale((c-w/2, v_rel))--Scale((c+w/2, v_rel)), p);
			draw(Scale((c, v_rel-u_rel))--Scale((c, v_rel+u_rel)), p);
		}
	}
}

//----------------------------------------------------------------------------------------------------

void DrawAll()
{
	draw(RootGetObject(f, dirs[0]+"/h_dsdt"), "eb", black);

	for (int di : dirs.keys)
		draw(RootGetObject(f, dirs[di]+"/g_fit"), "def", d_pens[di]);
}

//----------------------------------------------------------------------------------------------------

void DrawAllRel()
{
	DrawRel(RootGetObject(f, dirs[0]+"/h_dsdt"), black);

	for (int di : dirs.keys)
		DrawRel(RootGetObject(f, dirs[di]+"/g_fit"), d_pens[di]);
}

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
scale(Linear, Log);
DrawAll();
limits((0, 1e-2), (1, 1e3), Crop);

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
scale(Linear, Log);
DrawAll();
limits((0, 1e0), (0.3, 1e3), Crop);

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
scale(Linear, Log);
DrawAll();
limits((0, 1e2), (0.02, 1e3), Crop);

//----------------------------------------------------------------------------------------------------

NewRow();

NewPad(false);

NewPad("$|t|\ung{GeV^2}$", "$(\d\si/\d t - ref) / ref$");
DrawAllRel();
limits((0, -0.1), (0.3, 0.1), Crop);
