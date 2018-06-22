import root;
import pad_layout;

string topDir = "../../simulation/";

string f = topDir + "process_fits.root";


//----------------------------------------------------------------------------------------------------

void DrawOne(string configs[], pen c_pens[])
{
	string baseDir = "simu-norm-study/";

	NewRow();

	NewPad(false);
	for (int ci : configs.keys)
	{
		AddToLegend(configs[ci], c_pens[ci]);
	}
	AttachLegend();
	
	NewPad("mult.~norm.~error", "mean $\De\rh$");
	for (int ci : configs.keys)
	{
		RootObject graph = RootGetObject(f, baseDir + configs[ci]+"/g_rho_bias_vs_norm");

		draw(graph, "l", c_pens[ci]);
	}
	
	NewPad("mult.~norm.~error", "mean $\De a\ung{mb}$");
	for (int ci : configs.keys)
	{
		RootObject graph = RootGetObject(f, baseDir + configs[ci]+"/g_a_bias_vs_norm");

		draw(graph, "l", c_pens[ci]);
	}
	
}

//----------------------------------------------------------------------------------------------------

string configs[];
pen c_pens[];
configs.push("exp1,rho=0.06/none/ob-2-10-0.05/tmax0.15"); c_pens.push(black+dashed);
configs.push("exp1,rho=0.10/none/ob-2-10-0.05/tmax0.15"); c_pens.push(blue+dashed);
configs.push("exp1,rho=0.14/none/ob-2-10-0.05/tmax0.15"); c_pens.push(red+dashed);

configs.push("exp3,rho=0.06/none/ob-2-10-0.05/tmax0.15"); c_pens.push(black);
configs.push("exp3,rho=0.10/none/ob-2-10-0.05/tmax0.15"); c_pens.push(blue);
configs.push("exp3,rho=0.14/none/ob-2-10-0.05/tmax0.15"); c_pens.push(red);

DrawOne(configs, c_pens);

