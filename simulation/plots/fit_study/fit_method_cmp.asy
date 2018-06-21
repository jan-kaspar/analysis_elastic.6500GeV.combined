import root;
import pad_layout;

string topDir = "../../";

string f = topDir + "process_fits.root";

string base_dir = "simu-fit-study/<model>/stat+syst+norm/ob-2-10-0.05/";

string models[];
pen m_pens[];
models.push("exp1,rho=0.10"); m_pens.push(blue);
models.push("exp3,rho=0.10"); m_pens.push(red);

string fit_methods[] = {
	"coulomb_analysis_1",
	"coulomb_analysis_1_weighted",
	"coulomb_analysis_1_weighted_norm_lim",
	"coulomb_analysis_1_weighted_norm_unlim",
};

string quantities[], q_labels[], q_formats[];
int q_indeces[];
quantities.push("de_rho"); q_labels.push("\De\rh"); q_formats.push("%#.3f"); q_indeces.push(1);
quantities.push("de_a"); q_labels.push("\De a\ung{mb/GeV^2}"); q_formats.push("%#.1f"); q_indeces.push(3);
quantities.push("de_si_tot"); q_labels.push("\De\si_{\rm tot}\ung{mb}"); q_formats.push("%#.2f"); q_indeces.push(5);

xSizeDef = 10cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

string tMax(string model)
{
	string expN = substr(model, 0, 4);

	string t_max = "t_max_unknown";
	if (expN == "exp1") t_max = "0.07";
	if (expN == "exp3") t_max = "0.15";

	return t_max;
}

//----------------------------------------------------------------------------------------------------

for (int fmi : fit_methods.keys)
{
	NewRow();

	NewPad(false);
	label("{\SetFontSizesXX " + replace(fit_methods[fmi], "_", "\_") + "}");

	for (int qi : quantities.keys)
	{
		NewPad("$" + q_labels[qi] + "$");
		for (int mi : models.keys)
		{
			string dir = replace(base_dir, "<model>", models[mi]) + fit_methods[fmi] + "/tmax" + tMax(models[mi]) + "/norm0.00/";
			RootObject obj = RootGetObject(f, dir+"h_" + quantities[qi]);
			draw(obj, "vl,eb", m_pens[mi], models[mi]);

			real si = obj.rExec("GetRMS"), n = obj.rExec("GetEntries"), si_unc = si / sqrt(2. * n);
			string fmt = q_formats[qi];
			//AddToLegend(format("RMS = $" + fmt + " \pm", si) + format(" " + fmt + "$", si_unc));

			RootObject obj = RootGetObject(f, dir+"g_data");
			real ax[] = {0.};
			real ay[] = {0.};
			obj.vExec("GetPoint", q_indeces[qi], ax, ay);
			si = ax[0];
			si_unc = ay[0];

			AddToLegend(format("RMS = $" + fmt + " \pm", si) + format(" " + fmt + "$", si_unc));
		}
		AttachLegend();
	}
}

GShipout(vSkip=0mm);
