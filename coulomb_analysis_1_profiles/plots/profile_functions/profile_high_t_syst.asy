import root;
import pad_layout;

string topDir = "../../";

string dirs[];
pen d_pens[];
dirs.push("fit_for_profile/def"); d_pens.push(blue);
dirs.push("fit_for_profile/corr1"); d_pens.push(red+dashed);
dirs.push("fit_for_profile/corr2"); d_pens.push(heavygreen+longdashed);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
scale(Linear, Log);

for (int di : dirs.keys)
{
	string f = topDir + dirs[di] + "/impactParameterDistributions.root";
	pen p = d_pens[di];
	draw(RootGetObject(f, "g_dsdt_H"), p);
}

limits((0, 1e-8), (3, 1e3), Crop);
//limits((0.3, 1e-2), (0.7, 1e0), Crop);

//----------------------------------------------------------------------------------------------------

NewPad("$b\ung{fm}$", "$|P(b)|^2$");

for (int di : dirs.keys)
{
	string f = topDir + dirs[di] + "/impactParameterDistributions.root";
	pen p = d_pens[di];
	draw(RootGetObject(f, "g_A_mod2"), p);
}

limits((0, 0), (3, 0.35), Crop);
