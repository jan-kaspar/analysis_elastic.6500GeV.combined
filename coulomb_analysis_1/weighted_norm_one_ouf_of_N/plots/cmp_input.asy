import root;
import pad_layout;

string topDir = "../";

string dirs[] = {
	"fits_binc_1",
	"fits_binc_5",
	"fits_binc_10",
};

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

NewPad();
scale(Linear, Log);

for (int diri : dirs.keys)
{
	string f = topDir+dirs[diri]+"/2500-2rp-ob-2-10-0.05/exp3,t_max=0.15/fit.root";
	draw(RootGetObject(f, "g_data_coll_unc_full"), "p", mCi+1pt+StdPen(diri));
}
