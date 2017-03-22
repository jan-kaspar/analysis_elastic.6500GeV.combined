import root;
import pad_layout;

string f = "do_fit.root";

string variants[] = {
	"variant 1",
	"variant 2",
	"variant 4",
};

xSizeDef = 12cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");

currentpad.xTicks = LeftTicks(0.005, 0.001);

draw(RootGetObject(f, variants[0] + "/h_in"), "eb", black+1pt);

for (int vi : variants.keys)
{
	draw(RootGetObject(f, variants[vi] + "/g_dsdt_CH"), "l", StdPen(vi), variants[vi]);
}

limits((0, 400), (0.02, 900), Crop);

AttachLegend(NE, NE);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");
scale(Linear, Log);
//currentpad.xTicks = LeftTicks(0.005, 0.001);

draw(RootGetObject(f, variants[0] + "/h_in"), "eb", black+1pt);

for (int vi : variants.keys)
{
	draw(RootGetObject(f, variants[vi] + "/g_dsdt_CH"), "l", StdPen(vi), variants[vi]);
}

//limits((0, 400), (0.02, 900), Crop);

AttachLegend(NE, NE);
