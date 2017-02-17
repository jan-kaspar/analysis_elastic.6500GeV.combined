import root;
import pad_layout;

string f = "do_fit.root";

xSizeDef = 12cm;
ySizeDef = 8cm;

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t \ung{mb/GeV^2}$");

currentpad.xTicks = LeftTicks(0.005, 0.001);

AddToLegend("<{\bf elastic scattering}, $\sqrt{s} = 13\un{TeV}$");

AddToLegend("<{\it data}:");
AddToLegend("TOTEM preliminary", mPl+3pt+(black+1pt));
AddToLegend("$\be^*=2500\un{m}$");
AddToLegend("all fills merged");
AddToLegend("diagonals combined");
AddToLegend("statistical uncertainties only");

AddToLegend("<{\it example decomposition}:");

draw(RootGetObject(f, "g_dsdt_C"), "l", blue, "Coulomb only");
draw(RootGetObject(f, "g_dsdt_H"), "l", red, "hadronic only");
draw(RootGetObject(f, "g_dsdt_CH"), "l", heavygreen+1pt, "Coulomb $\oplus$ hadronic");

draw(RootGetObject(f, "h_in"), "eb", black+1pt);

limits((0, 400), (0.02, 900), Crop);

AttachLegend(NE, NE);
