import root;
import pad_layout;

string topDir = "../../simulation/";

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
	draw(RootGetObject(f, binnings[bi]+"/h_rel_stat_unc"), "d0,eb", red);

	NewPad("$|t|\ung{GeV^2}$", "rel.~syst.~uncertainty");
	draw(RootGetObject(f, binnings[bi]+"/h_syst_stddev_orig"), "d0", red, "orig");
	draw(RootGetObject(f, binnings[bi]+"/h_syst_stddev_regen"), "d0", blue+dashed, "regen");
	AttachLegend();
}
