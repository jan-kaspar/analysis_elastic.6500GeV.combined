import root;
import pad_layout;

string topDir = "../../simulation/";

string types[] = {
	"stat",
	"syst",
	"norm",
	//"norm-bias",
};

string f_unc = topDir + "/build_uncertainties.root";
string binning = "ob-2-10-0.05";

//----------------------------------------------------------------------------------------------------

for (int ti : types.keys)
{
	NewRow();

	NewPad(false);
	label("{\SetFontSizesXX " + types[ti] + "}");

	string f = topDir + "simu-validation/" + types[ti] + "/validation.root";

	// bias
	for (int i = 0; i < 2; ++i)
	{
		NewPad("$|t|\ung{GeV^2}$", "relative bias");
		draw(RootGetObject(f, "h_bias_rel"), "eb", heavygreen);

		if (types[ti] == "norm-bias")
			draw((0, 0.10)--(1., 0.10), red+1pt);

		if (i == 1)
			limits((0, -0.002), (0.2, 0.002), Crop);
	}

	// fluctuation
	for (int i = 0; i < 2; ++i)
	{
		NewPad("$|t|\ung{GeV^2}$", "relative fluctuation");
		draw(RootGetObject(f, "h_unc_rel"), "eb", blue);

		if (types[ti] == "stat")
			draw(RootGetObject(f_unc, binning + "/h_rel_stat_unc"), "vl", red+1pt);

		if (types[ti] == "syst")
			draw(RootGetObject(f_unc, binning + "/h_syst_stddev_regen"), "vl", red+1pt);

		if (types[ti] == "norm")
			draw((0, 0.05)--(1., 0.05), red+1pt);

		if (i == 1)
			limits((0, 0), (0.2, 0.05), Crop);
	}
}

GShipout(vSkip=1mm);
