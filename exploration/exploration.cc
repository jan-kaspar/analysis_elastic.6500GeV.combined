#include "include/Elegent/Constants.h"
#include "include/Elegent/CoulombInterference.h"

#include "../HadronicFitModel.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TF1.h"
#include "TRandom2.h"

#include <vector>

using namespace std;

//----------------------------------------------------------------------------------------------------

void RunExploration(const string& file,
	const vector<double>& a_a,
	const vector<double>& a_b1, const vector<double>& a_b2, const vector<double>&a_b3,
	const vector<HadronicFitModel::PhaseMode> &a_pm, const vector<double>& a_p0,
	const vector<CoulombInterference::CIMode> &a_cm)
{
	HadronicFitModel *fm = new HadronicFitModel();
	model = fm;
	
	// default central parameters
	fm->p_td = -0.53;
	
	// default peripheral parameters
	fm->p_A = 5.53;
	fm->p_ka = 4.01;
	fm->p_tm = -0.310;

	// default blending parameters
	fm->t1 = 0.2;
	fm->t2 = 0.5;

	TFile *outF = new TFile(file.c_str(), "recreate");

	for (unsigned int i_a = 0; i_a < a_a.size(); i_a++)
	for (unsigned int i_b1 = 0; i_b1 < a_b1.size(); i_b1++)
	for (unsigned int i_b2 = 0; i_b2 < a_b2.size(); i_b2++)
	for (unsigned int i_b3 = 0; i_b3 < a_b3.size(); i_b3++)
		for (unsigned int i_pm = 0; i_pm < a_pm.size(); i_pm++)
		for (unsigned int i_p0 = 0; i_p0 < a_p0.size(); i_p0++)
		for (unsigned int i_cm = 0; i_cm < a_cm.size(); i_cm++)
		{
			fm->a = a_a[i_a];

			fm->b1 = a_b1[i_b1];
			fm->b2 = a_b2[i_b2];
			fm->b3 = a_b3[i_b3];

			fm->phaseMode = a_pm[i_pm];
			fm->p0 = a_p0[i_p0];
			
			//fm_exp->b1 = a_b1[i_b1];

			coulomb->mode = a_cm[i_cm];
			
			char buf[100];
			sprintf(buf, "a=%.2E,b1=%.2f,b2=%.2f,b3=%.2f,%s,p0=%+.2f,%s",
					fm->a, fm->b1, fm->b2, fm->b3, fm->GetPhaseModeString().c_str(), fm->p0,
					coulomb->GetModeString().c_str());
			gDirectory = outF->mkdir(buf);

			printf("* %s\n", buf);

			int color = 1;
	
			TGraph *g_FH_Rho2 = new TGraph(); g_FH_Rho2->SetName("g_FH_Rho2"); g_FH_Rho2->SetTitle(";|t|;|F^{H}|^{2}"); g_FH_Rho2->SetLineColor(color);
			TGraph *g_FH_Theta = new TGraph(); g_FH_Theta->SetName("g_FH_Theta"); g_FH_Theta->SetTitle(";|t|;Arg F^{H}"); g_FH_Theta->SetLineColor(color);
			TGraph *g_FH_B = new TGraph(); g_FH_B->SetName("g_FH_B"); g_FH_B->SetTitle(";|t|;B^{H}"); g_FH_B->SetLineColor(color);
			TGraph *g_dsdt_H = new TGraph(); g_dsdt_H->SetName("g_dsdt_H"); g_dsdt_H->SetTitle(";|t|;d#sigma^{H}/dt"); g_dsdt_H->SetLineColor(color);

			TGraph *g_Psi_R = new TGraph(); g_Psi_R->SetName("g_Psi_R"); g_Psi_R->SetTitle(";|t|;Re #Psi"); g_Psi_R->SetLineColor(color);
			TGraph *g_Psi_I = new TGraph(); g_Psi_I->SetName("g_Psi_I"); g_Psi_I->SetTitle(";|t|;Im #Psi"); g_Psi_I->SetLineColor(color);
			
			TGraph *g_FC_Rho2 = new TGraph(); g_FC_Rho2->SetName("g_FC_Rho2"); g_FC_Rho2->SetTitle(";|t|;|F^{C}|^{2}"); g_FC_Rho2->SetLineColor(color);
			TGraph *g_dsdt_C = new TGraph(); g_dsdt_C->SetName("g_dsdt_C"); g_dsdt_C->SetTitle(";|t|;d#sigma^{C}/dt"); g_dsdt_C->SetLineColor(color);
			
			TGraph *g_FCH_Rho2 = new TGraph(); g_FCH_Rho2->SetName("g_FCH_Rho2"); g_FCH_Rho2->SetTitle(";|t|;|F^{C+H}|^{2}"); g_FCH_Rho2->SetLineColor(color);
			TGraph *g_dsdt_CH = new TGraph(); g_dsdt_CH->SetName("g_dsdt_CH"); g_dsdt_CH->SetTitle(";|t|;d#sigma^{CH}/dt"); g_dsdt_CH->SetLineColor(color);
			
			TGraph *g_C = new TGraph(); g_C->SetName("g_C"); g_C->SetTitle(";|t|;C"); g_C->SetLineColor(color);
			TGraph *g_Zv = new TGraph(); g_Zv->SetName("g_Zv"); g_Zv->SetTitle(";|t|;C"); g_Zv->SetLineColor(color);
			
			TGraph *g_SV = new TGraph(); g_SV->SetName("g_SV"); g_SV->SetTitle(";|t|;effect of B slope variation"); g_SV->SetLineColor(color);

			for (double mt = 1e-4; mt <= 1.;)
			{
				//printf("\t%E\n", mt);
				TComplex FC = coulomb->Amp_pure(-mt);

				TComplex FH = fm->Amp(-mt);

				double ep = 1E-5;
				double FH_B = log(fm->Amp(-mt+ep).Rho2() / FH.Rho2()) / ep;

				/*
				TComplex vPsi = coulomb->Psi_KL(-mt);
				g_Psi_R->SetPoint(g_Psi_R->GetN(), mt, vPsi.Re());
				g_Psi_I->SetPoint(g_Psi_I->GetN(), mt, vPsi.Im());

				TComplex FCH = coulomb->Amp_pure(-mt) + FH * TComplex::Exp(i * vPsi);
				*/
				TComplex FCH = coulomb->Amp(-mt);
				double C = (FCH.Rho2() - FH.Rho2()) / FH.Rho2();
				double Zv = (FCH.Rho2() - FC.Rho2() - FH.Rho2()) / (FC.Rho2() + FH.Rho2());

				int idx = g_FH_Rho2->GetN();

				g_FH_Rho2->SetPoint(idx, mt, FH.Rho2());
				g_FH_Theta->SetPoint(idx, mt, FH.Theta());
				g_FH_B->SetPoint(idx, mt, FH_B);
				g_dsdt_H->SetPoint(idx, mt, cnts->sig_fac * FH.Rho2());

				g_FC_Rho2->SetPoint(idx, mt, FC.Rho2());
				g_dsdt_C->SetPoint(idx, mt, cnts->sig_fac * FC.Rho2());

				g_FCH_Rho2->SetPoint(idx, mt, FCH.Rho2());
				g_dsdt_CH->SetPoint(idx, mt, cnts->sig_fac * FCH.Rho2());

				g_C->SetPoint(idx, mt, C);
				g_Zv->SetPoint(idx, mt, Zv);

				/*
				double sv_corr = FH.Rho() / fm_exp->Amp(-mt).Rho();
				g_SV->SetPoint(g_SV->GetN(), mt, sv_corr);
				*/

				if (mt < 0.1)
					mt *= 1.1;
				else
					mt += 0.01;
			}
			
			g_FH_Rho2->Write();
			g_FH_Theta->Write();
			g_FH_B->Write();
			g_dsdt_H->Write();

			g_FC_Rho2->Write();
			g_dsdt_C->Write();
			
			//g_Psi_R->Write();
			//g_Psi_I->Write();
			g_FCH_Rho2->Write();
			g_dsdt_CH->Write();

			g_C->Write();
			g_Zv->Write();
			
			//g_SV->Write();

			/*
			c->cd(1);
			gPad->SetLogx(1);
			g_Psi_R->Draw((idx == 0) ? "al" : "l");
			
			c->cd(2);
			gPad->SetLogx(1);
			g_Psi_I->Draw((idx == 0) ? "al" : "l");
			
			c->cd(3);
			gPad->SetLogx(1);
			g_C->Draw((idx == 0) ? "al" : "l");

			idx++;
			*/
		}

	//gDirectory = outF;
	//c->Write();

	delete outF;
	delete fm;
}

//----------------------------------------------------------------------------------------------------

void RunPeripheralExploration(const string& file,
	const vector<double>& a_a,
	const vector<double>& a_b1, const vector<double>& a_b2, const vector<double>&a_b3,
	const vector<double>& a_p0, const vector<double>& a_p_A, const vector<double>& a_p_ka, const vector<double>& a_p_tm
	)
{
	HadronicFitModel *fm = new HadronicFitModel();
	model = fm;

	// default blending parameters
	fm->t1 = 0.3;
	fm->t2 = 0.5;
			
	fm->phaseMode = HadronicFitModel::pmPeripheral;

	coulomb->mode = CoulombInterference::mKL;

	TFile *outF = new TFile(file.c_str(), "recreate");

	for (unsigned int i_a = 0; i_a < a_a.size(); i_a++)
	for (unsigned int i_b1 = 0; i_b1 < a_b1.size(); i_b1++)
	for (unsigned int i_b2 = 0; i_b2 < a_b2.size(); i_b2++)
	for (unsigned int i_b3 = 0; i_b3 < a_b3.size(); i_b3++)
		for (unsigned int i_p0 = 0; i_p0 < a_p0.size(); i_p0++)
		for (unsigned int i_p_A = 0; i_p_A < a_p_A.size(); i_p_A++)
		for (unsigned int i_p_ka = 0; i_p_ka < a_p_ka.size(); i_p_ka++)
		for (unsigned int i_p_tm = 0; i_p_tm < a_p_tm.size(); i_p_tm++)
		{
			fm->a = a_a[i_a];

			fm->b1 = a_b1[i_b1];
			fm->b2 = a_b2[i_b2];
			fm->b3 = a_b3[i_b3];

			fm->p0 = a_p0[i_p0];
			fm->p_A = a_p_A[i_p_A];
			fm->p_ka = a_p_ka[i_p_ka];
			fm->p_tm = a_p_tm[i_p_tm];
			
			char buf[100];
			sprintf(buf, "a=%.2E,b1=%.2f,b2=%.2f,b3=%.2f,p0=%+.2f,p_A=%+.2f,p_ka=%+.2f,p_tm=%+.2f",
					fm->a, fm->b1, fm->b2, fm->b3, fm->p0, fm->p_A, fm->p_ka, fm->p_tm);
			gDirectory = outF->mkdir(buf);

			printf("* %s\n", buf);

			/*
			// TODO: test - remove
			TGraph *g_B_int_re = new TGraph(); g_B_int_re->SetName("g_B_int_re");
			TGraph *g_B_int_im = new TGraph(); g_B_int_im->SetName("g_B_int_im");
			TGraph *g_I_int = new TGraph(); g_I_int->SetName("g_I_int");
			double mt0 = 1E-4;
			for (double mt = 0.5001*mt0; mt <= 1.5*mt0; mt += 0.0001*mt0)
			{
				TComplex A = model->Amp(-mt0);
				double par[] = { -mt0, A.Re(), A.Im() };
				TComplex v = CoulombInterference::B_integrand(-mt, par, coulomb);

				int idx = g_B_int_re->GetN();
				g_B_int_re->SetPoint(idx, mt, v.Re());
				g_B_int_im->SetPoint(idx, mt, v.Im());

				double I_int = coulomb->I_integral(-mt, -mt0);

				g_I_int->SetPoint(idx, mt, I_int);
			}

			g_B_int_re->Write();
			g_B_int_im->Write();
			g_I_int->Write();

			continue;
			*/

			int color = 1;
			
			TGraph *g_FC_Rho2 = new TGraph(); g_FC_Rho2->SetName("g_FC_Rho2"); g_FC_Rho2->SetTitle(";|t|;|F^{C}|^{2}"); g_FC_Rho2->SetLineColor(color);
	
			TGraph *g_FH_Re = new TGraph(); g_FH_Re->SetName("g_FH_Re"); g_FH_Re->SetTitle(";|t|;Re F^{H}"); g_FH_Re->SetLineColor(color);
			TGraph *g_FH_Im = new TGraph(); g_FH_Im->SetName("g_FH_Im"); g_FH_Im->SetTitle(";|t|;Im F^{H}"); g_FH_Im->SetLineColor(color);
			TGraph *g_FH_Rho2 = new TGraph(); g_FH_Rho2->SetName("g_FH_Rho2"); g_FH_Rho2->SetTitle(";|t|;|F^{H}|^{2}"); g_FH_Rho2->SetLineColor(color);
			TGraph *g_FH_Theta = new TGraph(); g_FH_Theta->SetName("g_FH_Theta"); g_FH_Theta->SetTitle(";|t|;Arg F^{H}"); g_FH_Theta->SetLineColor(color);
			TGraph *g_FH_B = new TGraph(); g_FH_B->SetName("g_FH_B"); g_FH_B->SetTitle(";|t|;B^{H}"); g_FH_B->SetLineColor(color);
			
			/*
			TGraph *g_Bterm_R = new TGraph(); g_Bterm_R->SetName("g_Bterm_R"); g_Bterm_R->SetTitle(";|t|;Re Bterm"); g_Bterm_R->SetLineColor(color);
			TGraph *g_Bterm_I = new TGraph(); g_Bterm_I->SetName("g_Bterm_I"); g_Bterm_I->SetTitle(";|t|;Im Bterm"); g_Bterm_I->SetLineColor(color);

			TGraph *g_Psi_R = new TGraph(); g_Psi_R->SetName("g_Psi_R"); g_Psi_R->SetTitle(";|t|;Re #Psi"); g_Psi_R->SetLineColor(color);
			TGraph *g_Psi_I = new TGraph(); g_Psi_I->SetName("g_Psi_I"); g_Psi_I->SetTitle(";|t|;Im #Psi"); g_Psi_I->SetLineColor(color);
			*/
			
			TGraph *g_FCH_Rho2 = new TGraph(); g_FCH_Rho2->SetName("g_FCH_Rho2"); g_FCH_Rho2->SetTitle(";|t|;|F^{C+H}|^{2}"); g_FCH_Rho2->SetLineColor(color);
			
			TGraph *g_C = new TGraph(); g_C->SetName("g_C"); g_C->SetTitle(";|t|;C"); g_C->SetLineColor(color);
			
			TGraph *g_SV = new TGraph(); g_SV->SetName("g_SV"); g_SV->SetTitle(";|t|;effect of B slope variation"); g_SV->SetLineColor(color);

			for (double mt = 1e-4; mt <= 1.5;)
			{
				int idx = g_FC_Rho2->GetN();

				TComplex FC = coulomb->Amp_pure(-mt);
				g_FC_Rho2->SetPoint(idx, mt, FC.Rho2());

				TComplex FH = model->Amp(-mt);
				g_FH_Re->SetPoint(idx, mt, FH.Re());
				g_FH_Im->SetPoint(idx, mt, FH.Im());
				g_FH_Rho2->SetPoint(idx, mt, FH.Rho2());
				g_FH_Theta->SetPoint(idx, mt, FH.Theta());
				double ep = 1E-5;
				g_FH_B->SetPoint(g_FH_B->GetN(), mt, log(model->Amp(-mt+ep).Rho2() / FH.Rho2()) / ep);
				
				//printf("\tt=%E | F_H: re=%E, im=%E\n", mt, FH.Re(), FH.Im());


				/*
				TComplex Bterm = coulomb->B_term(-mt);
				g_Bterm_R->SetPoint(idx, mt, Bterm.Re());
				g_Bterm_I->SetPoint(idx, mt, Bterm.Im());
				TComplex vPsi = coulomb->Psi_KL(-mt);
				g_Psi_R->SetPoint(g_Psi_R->GetN(), mt, vPsi.Re());
				g_Psi_I->SetPoint(g_Psi_I->GetN(), mt, vPsi.Im());

				TComplex FCH = coulomb->Amp_pure(-mt) + FH * TComplex::Exp(i * vPsi);
				*/
				TComplex FCH = coulomb->Amp(-mt);
				double C = (FCH.Rho2() - FH.Rho2()) / FH.Rho2();

				g_FCH_Rho2->SetPoint(g_FCH_Rho2->GetN(), mt, FCH.Rho2());
				g_C->SetPoint(g_C->GetN(), mt, C);

				/*
				double sv_corr = FH.Rho() / fm_exp->Amp(-mt).Rho();
				g_SV->SetPoint(g_SV->GetN(), mt, sv_corr);
				*/

				double w = mt * 0.1;
				if (mt > 0.05) w = 0.005;
				if (mt > 0.2) w = 0.025;
				if (mt > 0.6) w = 0.05;
				
				mt += w;
			}
			
			g_FC_Rho2->Write();

			g_FH_Re->Write();
			g_FH_Im->Write();
			g_FH_Rho2->Write();
			g_FH_Theta->Write();
			g_FH_B->Write();

			//g_Bterm_R->Write();
			//g_Bterm_I->Write();
			
			//g_Psi_R->Write();
			//g_Psi_I->Write();
			g_FCH_Rho2->Write();

			g_C->Write();
		}

	delete outF;
	delete fm;
}

//----------------------------------------------------------------------------------------------------

void RunPolynomialExploration(const string& file,
	const vector<double>& a_a,
	const vector<double>& a_b1, const vector<double>& a_b2, const vector<double>&a_b3,
	const vector<double>& a_p0, const vector<double>& a_p1, const vector<double>& a_p2, const vector<double>& a_p3
	)
{
	HadronicFitModel *fm = new HadronicFitModel();
	model = fm;

	// default blending parameters
	fm->t1 = 0.2;
	fm->t2 = 0.5;
			
	fm->phaseMode = HadronicFitModel::pmPolynomial;

	coulomb->mode = CoulombInterference::mKL;

	TFile *outF = new TFile(file.c_str(), "recreate");

	for (unsigned int i_a = 0; i_a < a_a.size(); i_a++)
	for (unsigned int i_b1 = 0; i_b1 < a_b1.size(); i_b1++)
	for (unsigned int i_b2 = 0; i_b2 < a_b2.size(); i_b2++)
	for (unsigned int i_b3 = 0; i_b3 < a_b3.size(); i_b3++)
		for (unsigned int i_p0 = 0; i_p0 < a_p0.size(); i_p0++)
		for (unsigned int i_p1 = 0; i_p1 < a_p1.size(); i_p1++)
		for (unsigned int i_p2 = 0; i_p2 < a_p2.size(); i_p2++)
		for (unsigned int i_p3 = 0; i_p3 < a_p3.size(); i_p3++)
		{
			fm->a = a_a[i_a];

			fm->b1 = a_b1[i_b1];
			fm->b2 = a_b2[i_b2];
			fm->b3 = a_b3[i_b3];

			fm->p0 = a_p0[i_p0];
			fm->p1 = a_p1[i_p1];
			fm->p2 = a_p2[i_p2];
			fm->p3 = a_p3[i_p3];
			
			char buf[100];
			sprintf(buf, "a=%.2E,b1=%.2f,b2=%.2f,b3=%.2f,p0=%+.2f,p1=%+.2f,p2=%+.2f,p3=%+.2f",
					fm->a, fm->b1, fm->b2, fm->b3, fm->p0, fm->p1, fm->p2, fm->p3);
			gDirectory = outF->mkdir(buf);

			printf("* %s\n", buf);

			int color = 1;
			
			TGraph *g_FC_Rho2 = new TGraph(); g_FC_Rho2->SetName("g_FC_Rho2"); g_FC_Rho2->SetTitle(";|t|;|F^{C}|^{2}"); g_FC_Rho2->SetLineColor(color);
	
			TGraph *g_FH_Re = new TGraph(); g_FH_Re->SetName("g_FH_Re"); g_FH_Re->SetTitle(";|t|;Re F^{H}"); g_FH_Re->SetLineColor(color);
			TGraph *g_FH_Im = new TGraph(); g_FH_Im->SetName("g_FH_Im"); g_FH_Im->SetTitle(";|t|;Im F^{H}"); g_FH_Im->SetLineColor(color);
			TGraph *g_FH_Rho2 = new TGraph(); g_FH_Rho2->SetName("g_FH_Rho2"); g_FH_Rho2->SetTitle(";|t|;|F^{H}|^{2}"); g_FH_Rho2->SetLineColor(color);
			TGraph *g_FH_Theta = new TGraph(); g_FH_Theta->SetName("g_FH_Theta"); g_FH_Theta->SetTitle(";|t|;Arg F^{H}"); g_FH_Theta->SetLineColor(color);
			TGraph *g_FH_B = new TGraph(); g_FH_B->SetName("g_FH_B"); g_FH_B->SetTitle(";|t|;B^{H}"); g_FH_B->SetLineColor(color);
			
			/*
			TGraph *g_Bterm_R = new TGraph(); g_Bterm_R->SetName("g_Bterm_R"); g_Bterm_R->SetTitle(";|t|;Re Bterm"); g_Bterm_R->SetLineColor(color);
			TGraph *g_Bterm_I = new TGraph(); g_Bterm_I->SetName("g_Bterm_I"); g_Bterm_I->SetTitle(";|t|;Im Bterm"); g_Bterm_I->SetLineColor(color);

			TGraph *g_Psi_R = new TGraph(); g_Psi_R->SetName("g_Psi_R"); g_Psi_R->SetTitle(";|t|;Re #Psi"); g_Psi_R->SetLineColor(color);
			TGraph *g_Psi_I = new TGraph(); g_Psi_I->SetName("g_Psi_I"); g_Psi_I->SetTitle(";|t|;Im #Psi"); g_Psi_I->SetLineColor(color);
			*/
			
			TGraph *g_FCH_Rho2 = new TGraph(); g_FCH_Rho2->SetName("g_FCH_Rho2"); g_FCH_Rho2->SetTitle(";|t|;|F^{C+H}|^{2}"); g_FCH_Rho2->SetLineColor(color);
			
			TGraph *g_C = new TGraph(); g_C->SetName("g_C"); g_C->SetTitle(";|t|;C"); g_C->SetLineColor(color);
			
			TGraph *g_SV = new TGraph(); g_SV->SetName("g_SV"); g_SV->SetTitle(";|t|;effect of B slope variation"); g_SV->SetLineColor(color);

			for (double mt = 1e-4; mt <= 1.4;)
			{
				int idx = g_FC_Rho2->GetN();

				TComplex FC = coulomb->Amp_pure(-mt);
				g_FC_Rho2->SetPoint(idx, mt, FC.Rho2());

				TComplex FH = model->Amp(-mt);
				g_FH_Re->SetPoint(idx, mt, FH.Re());
				g_FH_Im->SetPoint(idx, mt, FH.Im());
				g_FH_Rho2->SetPoint(idx, mt, FH.Rho2());
				g_FH_Theta->SetPoint(idx, mt, FH.Theta());
				double ep = 1E-5;
				g_FH_B->SetPoint(g_FH_B->GetN(), mt, log(model->Amp(-mt+ep).Rho2() / FH.Rho2()) / ep);
				

				TComplex FCH = coulomb->Amp(-mt);
				double C = (FCH.Rho2() - FH.Rho2()) / FH.Rho2();

				g_FCH_Rho2->SetPoint(g_FCH_Rho2->GetN(), mt, FCH.Rho2());
				g_C->SetPoint(g_C->GetN(), mt, C);

				double w = mt * 0.1;
				if (mt > 0.05) w = 0.005;
				if (mt > 0.2) w = 0.01;
				if (mt > 0.6) w = 0.02;
				
				mt += w;
			}
			
			g_FC_Rho2->Write();

			g_FH_Re->Write();
			g_FH_Im->Write();
			g_FH_Rho2->Write();
			g_FH_Theta->Write();
			g_FH_B->Write();

			g_FCH_Rho2->Write();

			g_C->Write();
		}

	delete outF;
	delete fm;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	Constants::Init(2*6500, cnts->mPP);
    cnts->Print();

	coulomb->ffType = coulomb->ffPuckett;
	coulomb->Print();

	vector<double> a_a, a_b1, a_b2, a_b3, a_p0, a_p1, a_p2, a_p3, a_p_A, a_p_ka, a_p_tm;
	vector<HadronicFitModel::PhaseMode> a_pm;
	vector<CoulombInterference::CIMode> a_cm;

	//----------------------------------------------------------------------------------------------------

	a_a.clear();
	a_a.push_back(1.84E9);
	a_a.push_back(1.84E9 * 1.02);
	a_a.push_back(1.84E9 * 1.04);
	a_a.push_back(1.84E9 * 1.20);

	a_b1.clear();
	a_b1.push_back(10.2);

	a_b2.clear();
	a_b2.push_back(0.);
	
	a_b3.clear();
	a_b3.push_back(0.);

	a_p0.clear();
	//a_p0.push_back(M_PI/2. - atan(0.06));
	a_p0.push_back(M_PI/2. - atan(0.08));
	a_p0.push_back(M_PI/2. - atan(0.10));
	a_p0.push_back(M_PI/2. - atan(0.12));
	a_p0.push_back(M_PI/2. - atan(0.14));
	a_p0.push_back(M_PI/2. - atan(0.16));

	a_pm.clear();
	a_pm.push_back(HadronicFitModel::pmConstant);

	a_cm.clear();
	a_cm.push_back(CoulombInterference::mKL);
	
	RunExploration("exploration.root", a_a, a_b1, a_b2, a_b3, a_pm, a_p0, a_cm);
	
	//----------------------------------------------------------------------------------------------------

	a_a.clear();
	a_a.push_back(1.84E9);

	a_b1.clear();
	a_b1.push_back(10.2);

	a_b2.clear();
	a_b2.push_back(4.4);
	
	a_b3.clear();
	a_b3.push_back(10);

	a_p0.clear();
	a_p0.push_back(M_PI/2. - atan(0.08));
	a_p0.push_back(M_PI/2. - atan(0.10));
	a_p0.push_back(M_PI/2. - atan(0.12));
	a_p0.push_back(M_PI/2. - atan(0.14));
	a_p0.push_back(M_PI/2. - atan(0.16));

	a_pm.clear();
	a_pm.push_back(HadronicFitModel::pmConstant);

	a_cm.clear();
	a_cm.push_back(CoulombInterference::mKL);
	
	RunExploration("exploration2.root", a_a, a_b1, a_b2, a_b3, a_pm, a_p0, a_cm);

	//----------------------------------------------------------------------------------------------------

	return 0;
}
