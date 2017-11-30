import root;
import pad_layout;

string topDir = "../../simulation/";

string f = topDir + "process_fits.root";


//----------------------------------------------------------------------------------------------------

void DrawOne(string configs[])
{
	string baseDir = "simu-fit-study/";

	NewRow();

	NewPad(false);
	for (int ci : configs.keys)
	{
		pen p = StdPen(ci);
		AddToLegend(configs[ci], p);
	}
	AttachLegend();
	
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

string configs[] = {
	"exp1,rho=0.10/stat/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp1,rho=0.10/stat+syst/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp1,rho=0.10/stat+syst+norm/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp1,rho=0.10/norm/ob-3-5-0.05/tmax0.15/norm0.00",
};
DrawOne(configs);

string configs[] = {
	"exp1,rho=0.14/stat/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp1,rho=0.14/stat+syst/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp1,rho=0.14/stat+syst+norm/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp1,rho=0.14/norm/ob-3-5-0.05/tmax0.15/norm0.00",
};
DrawOne(configs);

string configs[] = {
	"exp2,rho=0.10/stat/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp2,rho=0.10/stat+syst/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp2,rho=0.10/stat+syst+norm/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp2,rho=0.10/norm/ob-3-5-0.05/tmax0.15/norm0.00",
};
DrawOne(configs);

string configs[] = {
	"exp2,rho=0.14/stat/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp2,rho=0.14/stat+syst/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp2,rho=0.14/stat+syst+norm/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp2,rho=0.14/norm/ob-3-5-0.05/tmax0.15/norm0.00",
};
DrawOne(configs);


string configs[] = {
	"exp3,rho=0.10/stat/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp3,rho=0.10/stat+syst/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp3,rho=0.10/stat+syst+norm/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp3,rho=0.10/norm/ob-3-5-0.05/tmax0.15/norm0.00",
};
DrawOne(configs);

string configs[] = {
	"exp3,rho=0.14/stat/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp3,rho=0.14/stat+syst/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp3,rho=0.14/stat+syst+norm/ob-3-5-0.05/tmax0.15/norm0.00",
	"exp3,rho=0.14/norm/ob-3-5-0.05/tmax0.15/norm0.00",
};
DrawOne(configs);
