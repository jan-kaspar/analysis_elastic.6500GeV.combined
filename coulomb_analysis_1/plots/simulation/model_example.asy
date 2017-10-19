import root;
import pad_layout;

string topDir = "../../simulation/";

string f = topDir + "build_models.root";

string binnings[] = {
	"ob-1-20-0.05",
	"ob-2-10-0.05",
	"ob-3-5-0.05",
};

//----------------------------------------------------------------------------------------------------

void DrawOne(string model)
{
	NewRow();

	NewPad(false);
	label(model);

	NewPad("$|t|\un{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
	scale(Linear, Log);
	draw(RootGetObject(f, model + "/g_dsdt_H"), blue, "hadronic");
	draw(RootGetObject(f, model + "/g_dsdt_CH"), red, "Coulomb + hadronic");
	AttachLegend();


	NewPad("$|t|\un{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
	scale(Linear, Log);
	AddToLegend("<Coulomb + hadronic");
	for (int bi : binnings.keys)
	{
		draw(RootGetObject(f, model + "/" + binnings[bi] + "/h_dsdt_CH"), StdPen(bi+1), binnings[bi]);
	}
	AttachLegend();
}

//----------------------------------------------------------------------------------------------------

DrawOne("exp1,rho=0.10");

DrawOne("exp3,rho=0.14");
