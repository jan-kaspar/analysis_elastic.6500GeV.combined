import graph;
import pad_layout;

real rho_0, si_el_0;

real rho_90(real si_el)
{
	real K = si_el_0 * (1. + rho_0 * rho_0);

	return sqrt(K/si_el - 1);
}

real rho_2500(real si_el)
{
	real S = -0.2;

	return rho_0 + S * (si_el / si_el_0 - 1.);
}

string MakeLabel()
{
	return format("$\rh^0 = %#.4f$, ", rho_0) + format("$\si_{\rm el}^0 = %#.3f\un{mb}$", si_el_0);
}

//----------------------------------------------------------------------------------------------------

NewPad("$\si_{\rm el}\ung{mb}$", "$\rh$", xSize=8cm, ySize=7cm);
currentpad.xTicks = LeftTicks(0.1, 0.02);

real si_el_min = 30.7;
real si_el_max = 31.3;

AddToLegend("<$\be^* = 2500\un{m}$:");

rho_0 = 0.0979;
si_el_0 = 31.0;
draw(graph(rho_2500, si_el_min, si_el_max), red, MakeLabel());
draw((si_el_0, rho_0), mCi+red+2pt);

/*
rho_0 = 0.0952;
si_el_0 = 31.5;
draw(graph(rho_2500, si_el_min, si_el_max), magenta+dashed, MakeLabel());
draw((si_el_0, rho_0), mCi+magenta+2pt);
*/

AddToLegend("($|t|_{\rm max} = 0.15\un{GeV^2}$, med.~binning, exp3, const.~phase)");

AddToLegend("<$\be^* = 90\un{m}$:");

rho_0 = 0.10;
si_el_0 = 31.0359;
draw(graph(rho_90, si_el_min, si_el_max), blue, MakeLabel());
draw((si_el_0, rho_0), mCi+blue+2pt);

rho_0 = 0.14;
si_el_0 = 30.7437;
draw(graph(rho_90, si_el_min, si_el_max), heavygreen+dashed, MakeLabel());
draw((si_el_0, rho_0), mCi+heavygreen+2pt);

AddToLegend("<consistent combination:");

real rho_sol = 0.09758;
real si_el_sol = 31.0505;

real scale_factor = si_el_sol / 31.0359;
real si_tot_sol = 110.5753 * scale_factor;
real si_inel_sol = 79.5394 * scale_factor;

string label_sol = format("$\rh = %#.4f$, ", rho_sol) + format("$\si_{\rm el} = %#.3f\un{mb}$", si_el_sol);
draw((si_el_sol, rho_sol), nullpen, label_sol, mSt+3pt);

AddToLegend(format("$\si_{\rm tot} = %#.3f\un{mb}$", si_tot_sol));
AddToLegend(format("$\si_{\rm inel} = %#.3f\un{mb}$", si_inel_sol));

limits((si_el_min, 0.08), (si_el_max, 0.145), Crop);
//limits((31.0, 0.09), (31.1, 0.10), Crop);

AttachLegend(NW, NE);
