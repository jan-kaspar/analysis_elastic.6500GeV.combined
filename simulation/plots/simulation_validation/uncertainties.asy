import root;
import pad_layout;

string topDir = "../../";

string f = topDir + "build_uncertainties.root";

string binnings[] = {
	"ob-1-20-0.05",
	"ob-2-10-0.05",
	"ob-3-5-0.05",
};

//----------------------------------------------------------------------------------------------------

for (int bi : binnings.keys)
{
	NewRow();

	NewPad(false);
	label(binnings[bi]);

	NewPad("$|t|\ung{GeV^2}$", "rel.~stat.~uncertainty");
	draw(RootGetObject(f, binnings[bi]+"/h_rel_stat_unc"), "d0,eb", red+1pt);

	NewPad("$|t|\ung{GeV^2}$", "rel.~stat.~uncertainty");
	draw(RootGetObject(f, binnings[bi]+"/h_rel_stat_unc"), "d0,eb,vl", red+1pt);
	limits((0, 0.), (0.2, 0.02), Crop);

	NewPad("$|t|\ung{GeV^2}$", "rel.~syst.~uncertainty");
	draw(RootGetObject(f, binnings[bi]+"/h_syst_stddev_orig"), "d0,vl", heavygreen, "orig");
	draw(RootGetObject(f, binnings[bi]+"/h_syst_stddev_regen"), "d0,vl", blue+dashed, "regen");

	NewPad("$|t|\ung{GeV^2}$", "rel.~syst.~uncertainty");
	draw(RootGetObject(f, binnings[bi]+"/h_syst_stddev_orig"), "d0,vl", heavygreen+1pt, "orig");
	draw(RootGetObject(f, binnings[bi]+"/h_syst_stddev_regen"), "d0,vl", blue+dashed+1pt, "regen");
	limits((0, 0.), (0.2, 0.03), Crop);
	AttachLegend();
}

GShipout(hSkip=1mm);
