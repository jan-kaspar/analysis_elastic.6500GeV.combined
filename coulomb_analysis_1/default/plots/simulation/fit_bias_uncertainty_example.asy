import root;
import pad_layout;

string topDir = "../../simulation/";

string f = topDir + "process_fits.root";


//----------------------------------------------------------------------------------------------------

void DrawOne(string caption, string configs[], string c_labels[])
{
	string baseDir = "simu-fit-study/";

	NewRow();

	NewPad(false);
	for (int ci : configs.keys)
	{
		pen p = StdPen(ci);
		AddToLegend(c_labels[ci], p);
	}
	AttachLegend(caption);
	
	NewPad("$\De\rh$");
	for (int ci : configs.keys)
	{
		RootObject hist = RootGetObject(f, baseDir + configs[ci]+"/h_de_rho");
		RootObject data = RootGetObject(f, baseDir + configs[ci]+"/g_data");

		real ax[] = {0.};
		real ay[] = {0.};

		data.vExec("GetPoint", 0, ax, ay); real mean = ax[0], mean_unc = ay[0];
		data.vExec("GetPoint", 1, ax, ay); real stddev = ax[0], stddev_unc = ay[0];
		
		string label = format("mean = $%#+.4f$", mean) + format("$\pm %#.4f$", mean_unc) +
			format(", RMS = $%#.4f$", stddev) + format("$\pm %#.4f$", stddev_unc);

		pen p = StdPen(ci);
		draw(hist, "vl", p, label);
	}
	AttachLegend(S, N);
	
	NewPad("$\De a\ung{mb}$");
	for (int ci : configs.keys)
	{
		RootObject hist = RootGetObject(f, baseDir + configs[ci]+"/h_de_a");
		RootObject data = RootGetObject(f, baseDir + configs[ci]+"/g_data");

		real ax[] = {0.};
		real ay[] = {0.};

		data.vExec("GetPoint", 2, ax, ay); real mean = ax[0], mean_unc = ay[0];
		data.vExec("GetPoint", 3, ax, ay); real stddev = ax[0], stddev_unc = ay[0];
		
		string label = format("mean = $%#+.2f$", mean) + format("$\pm %#.2f$", mean_unc) +
			format(", RMS = $%#.2f$", stddev) + format("$\pm %#.2f$", stddev_unc);

		pen p = StdPen(ci);
		draw(hist, "vl", p, label);
	}
	AttachLegend(S, N);
}

//----------------------------------------------------------------------------------------------------

string configs[], c_labels[];
configs.push("exp3,rho=0.10/stat/ob-2-10-0.05/tmax0.15/norm0.00"); c_labels.push("stat");
configs.push("exp3,rho=0.10/stat+syst/ob-2-10-0.05/tmax0.15/norm0.00"); c_labels.push("stat+syst");
configs.push("exp3,rho=0.10/stat+syst+norm/ob-2-10-0.05/tmax0.15/norm0.00"); c_labels.push("stat+syst+norm");

DrawOne("exp3, $\rh=0.10$, ob-2-10-0.05, $|t|_{\rm max} = 0.15\un{GeV^2}$",  configs, c_labels);
