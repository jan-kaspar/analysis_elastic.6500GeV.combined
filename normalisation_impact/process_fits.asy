import root;
import pad_layout;

string models[] = {
	"exp,rho=0.10",
	"exp,rho=0.14",
	"non-exp,rho=0.10",
	"non-exp,rho=0.14",
};

string quantities[];
string q_label[];
real q_min[], q_max[];

quantities.push("mean_De_rho"); q_label.push("mean of $\De\rh$"); q_min.push(-0.05); q_max.push(+0.05);
quantities.push("rms_De_rho"); q_label.push("RMS of $\De\rh$"); q_min.push(0.002); q_max.push(0.003);
quantities.push("mean_de_A"); q_label.push("mean of $\de A$"); q_min.push(-0.2); q_max.push(+0.2);
quantities.push("rms_de_A"); q_label.push("RMS of $\de A$"); q_min.push(8e-4); q_max.push(13e-4);

// TODO
string f = "process_fits.root";
//string f = "process_fits2.root";

xSizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

for (int qi : quantities.keys)
{
	if (qi == 2)
		NewRow();

	NewPad("mult.~normalisation error", q_label[qi]);

	for (int mi : models.keys)
	{
		pen p = StdPen(mi);

		RootObject obj = RootGetObject(f, models[mi] + "/g_" + quantities[qi]);
		draw(obj, "l,p", p, mCi+2pt+p, models[mi]);

		real diff = obj.rExec("Eval", 1.04) - obj.rExec("Eval", 1.);
		write(diff);
	}

	ylimits(q_min[qi], q_max[qi], Crop);

	if (quantities[qi] == "mean_De_rho")
	{
		currentpad.yTicks = RightTicks(0.01, 0.005);
	}
}

AddToLegend("<$\De\rho = \rho_{\rm fit} - \rho_{\rm sim}$");
AddToLegend("<$\de A = (A_{\rm fit} - A_{\rm sim}) / A_{\rm sim}$");

frame f_legend = BuildLegend();

NewPad(false);
add(f_legend);

GShipout(vSkip=1mm);
