import root;
import pad_layout;

string topDir = "../";

string f = topDir + "process_fits.root";

string fit_labels[];
fit_labels.push("approach1");
fit_labels.push("approach2_step_c");
fit_labels.push("approach3_step_d");
fit_labels.push("approach3_step_f");

string unc_labels[];
pen unc_pens[];
unc_labels.push("norm"); unc_pens.push(red);
unc_labels.push("stat+syst+norm"); unc_pens.push(blue+1pt);

string quantities[], q_labels[];
quantities.push("rho"); q_labels.push("\De\rh");
quantities.push("si_tot"); q_labels.push("\De\si_{\rm tot}\ung{mb}");

xSizeDef = 10cm;
ySizeDef = 10cm;

//----------------------------------------------------------------------------------------------------

NewPad(false);
for (int fi : fit_labels.keys)
	NewPadLabel(replace(fit_labels[fi], "_", "\_"));

for (int qi : quantities.keys)
{
	NewRow();

	NewPadLabel("$" + q_labels[qi] + "$");

	for (int fi : fit_labels.keys)
	{
		NewPad("$" + q_labels[qi] + "$");

		for (int ui : unc_labels.keys)
		{
			RootObject obj = RootGetObject(f, fit_labels[fi] + "/" + unc_labels[ui] + "/h_de_" + quantities[qi]);
			draw(obj, "vl", unc_pens[ui]);
			AddToLegend(unc_labels[ui], unc_pens[ui]);
			AddToLegend(format("RMS = %#.4f", obj.rExec("GetRMS")));
		}
		
		AttachLegend();
	}
}
