#ifndef _HadronicFitModel_h_
#define _HadronicFitModel_h_

#include "include/Elegent/Model.h"

#include <string>

using namespace Elegent;

class HadronicFitModel : public Model
{
  public:
	/// modulus parameters
	enum ModulusMode { mmUnknown, mmExp, mmVojtech, mmVojtechHighTFixed } modulusMode;
	double a, b1, b2, b3, b4, b5, b6, b7, b8, b9;
	double a2, c1, c2, d1;
	
	/// parameters for blending the low-|t| (variable) and high-|t| (fixed) modulus
	/// the interval (t1, t2) corresponds to (-3 sigma, +3 sigma)
	double t1, t2;

	/// scale factor for the high-|t| part
	double hts;

	/// phase parameters
	enum PhaseMode { pmUnknown, pmConstant, pmBailly, pmStandard, pmPeripheral, pmPolynomial } phaseMode;
	double p0, p1, p2, p3, p_td, p_t0, p_tau, p_A, p_ka, p_tm;

    HadronicFitModel() :
		modulusMode(mmExp),
		a(0.), b1(0.), b2(0.), b3(0.), b4(0.), b5(0.), b6(0.), b7(0.), b8(0.), b9(0.),
		a2(0.), c1(0.), c2(0.), d1(0.),
		t1(0.2), t2(0.5),
		hts(1.),
		phaseMode(pmUnknown),
		p0(0.), p1(0.), p2(0.), p3(0.), p_td(0.), p_t0(0.), p_tau(0.), p_A(0.), p_ka(0.), p_tm(0.),
		integ_workspace_initialized(false)
	{
	}

    virtual ~HadronicFitModel() {}

	virtual void Init()
	{
		fullLabel.name = "hadronic fit model";
		shortLabel.name = "hfm";

		precision_t = 1E-4;
		upper_bound_t = -50.;

		if (!integ_workspace_initialized)
		{
			integ_workspace_size_t = 100;
			integ_workspace_t = gsl_integration_workspace_alloc(integ_workspace_size_t);
			integ_workspace_initialized = true;
		}
	}

    virtual void Print() const
	{
		printf(">> HadronicFitModel\n");

		printf("\tmodulus: %s\n", GetModulusModeString().c_str());
		if (modulusMode == mmExp)
			printf("\t\ta=%.3E, b1=%.3E, b2=%.3E, b3=%.3E, b4=%.3E, b5=%.3E, b6=%.3E, b7=%.3E, b8=%.3E, b9=%.3E\n",
				a, b1, b2, b3, b4, b5, b6, b7, b8, b9);
		if (modulusMode == mmVojtech || modulusMode == mmVojtechHighTFixed)
			printf("\t\ta=%.3E, a2=%.3E, b1=%.3E, b2=%.3E, b3=%.3E, c1=%.3E, c2=%.3E, d1=%.3E\n", a, a2, b1, b2, b3, c1, c2, d1);
		
		printf("\tphase: %s\n", GetPhaseModeString().c_str());
		if (phaseMode == pmConstant)
			printf("\t\tp0=%.3E\n", p0);
		if (phaseMode == pmBailly)
			printf("\t\tp0=%.3E, p_td=%.3E\n", p0, p_td);
		if (phaseMode == pmStandard)
			printf("\t\tp0=%.3E, p_t0=%.3E, p_tau=%.3E\n", p0, p_t0, p_tau);
		if (phaseMode == pmPeripheral)
			printf("\t\tp0=%.3E, p_A=%.3E, p_ka=%.3E, p_tm=%.3E\n", p0, p_A, p_ka, p_tm);
		if (phaseMode == pmPolynomial)
			printf("\t\tp0=%.3E, p1=%.3E, p2=%.3E, p3=%.3E\n", p0, p1, p2, p3);
	}

	void PrintCode(const char *name = "hfm") const
	{
		if (modulusMode == mmExp)
			printf("%s->modulusMode = HadronicFitModel::mmExp;\n%s->a=%.5E;\n%s->b1=%.5E;\n%s->b2=%.5E;\n%s->b3=%.5E;\n%s->b4=%.5E;\n%s->b5=%.5E;\n%s->b6=%.5E;\n%s->b7=%.5E;\n%s->b8=%.5E;\n%s->b9=%.5E;\n",
				name, name, a, name, b1, name, b2, name, b3, name, b4, name, b5, name, b6, name, b7, name, b8, name, b9);
		if (modulusMode == mmVojtech)
			printf("%s->modulusMode = HadronicFitModel::mmVojtech;\n%s->a=%.5E;\n%s->a2=%.5E;\n%s->b1=%.5E;\n%s->b2=%.5E;\n%s->b3=%.5E;\n%s->c1=%.5E;\n%s->c2=%.5E;\n%s->d1=%.5E;\n",
				name, name, a, name, a2, name, b1, name, b2, name, b3, name, c1, name, c2, name, d1);
		if (modulusMode == mmVojtechHighTFixed)
			printf("%s->modulusMode = HadronicFitModel::mmVojtechHighTFixed;\n%s->a=%.5E;\n%s->a2=%.5E;\n%s->b1=%.5E;\n%s->b2=%.5E;\n%s->b3=%.5E;\n%s->c1=%.5E;\n%s->c2=%.5E;\n%s->d1=%.5E;\n",
				name, name, a, name, a2, name, b1, name, b2, name, b3, name, c1, name, c2, name, d1);

		if (modulusMode == mmExp || modulusMode == mmVojtechHighTFixed)
			printf("%s->t1=%.5E;\n%s->t2=%.5E;\n", name, t1, name, t2);

		printf("\n");

		if (phaseMode == pmConstant)
			printf("%s->phaseMode = HadronicFitModel::pmConstant;\n%s->p0=%.5E;\n", name, name, p0);
		if (phaseMode == pmBailly)
			printf("%s->phaseMode = HadronicFitModel::pmBailly;\n%s->p0=%.5E;\n%s->p_td=%.5E;\n", name, name, p0, name, p_td);
		if (phaseMode == pmStandard)
			printf("%s->phaseMode = HadronicFitModel::pmStandard;\n%s->p0=%.5E;\n%s->p_t0=%.5E;\n%s->p_tau=%.5E;\n",
				name, name, p0, name, p_t0, name, p_tau);
		if (phaseMode == pmPeripheral)
			printf("%s->phaseMode = HadronicFitModel::pmPeripheral;\n%s->p0=%.5E;\n%s->p_A=%.5E;\n%s->p_ka=%.5E;\n%s->p_tm=%.5E;\n",
				name, name, p0, name, p_A, name, p_ka, name, p_tm);
		if (phaseMode == pmPolynomial)
			printf("%s->phaseMode = HadronicFitModel::pmPolynomial;\n%s->p0=%.5E;\n%s->p1=%.5E;\n%s->p2=%.5E;\n%s->p3=%.5E;\n",
				name, name, p0, name, p1, name, p2, name, p3);
	}

	std::string GetStateString() const
	{
		char buf_m[200];
		if (modulusMode == mmExp)
			sprintf(buf_m, "modulus: exp, a=%.3E, b1=%.3E, b2=%.3E, b3=%.3E, b4=%.3E, b5=%.3E, b6=%.3E, b7=%.3E, b8=%.3E, b9=%.3E", a, b1, b2, b3, b4, b5, b6, b7, b8, b9);
		if (modulusMode == mmVojtech || modulusMode == mmVojtechHighTFixed)
			sprintf(buf_m, "mudulus: Vojtech, a=%.3E, a2=%.3E, b1=%.3E, b2=%.3E, b3=%.3E, c1=%.3E, c2=%.3E, d1=%.3E", a, a2, b1, b2, b3, c1, c2, d1);

		char buf_p[200];
		if (phaseMode == pmConstant)
			sprintf(buf_p, "phase: constant, p0=%.3E", p0);
		if (phaseMode == pmBailly)
			sprintf(buf_p, "phase: Bailly, p0=%.3E, p_td=%.3E", p0, p_td);
		if (phaseMode == pmStandard)
			sprintf(buf_p, "phase: standard, p0=%.3E, p_t0=%.3E, p_tau=%.3E", p0, p_t0, p_tau);
		if (phaseMode == pmPeripheral)
			sprintf(buf_p, "phase: peripheral, p0=%.3E, p_A=%.3E, p_ka=%.3E, p_tm=%.3E", p0, p_A, p_ka, p_tm);
		if (phaseMode == pmPolynomial)
			sprintf(buf_p, "phase: polynomial, p0=%.3E, p1=%.3E, p2=%.3E, p3=%.3E", p0, p1, p2, p3);

		return std::string(buf_m) + "; " + buf_p;
	}

	std::string GetStateStringShort() const
	{
		char buf_m[200];
		if (modulusMode == mmExp)
			sprintf(buf_m, "a=%.3E, b1=%.3E, b2=%.3E, b3=%.3E, b4=%.3E, b5=%.3E, b6=%.3E, b7=%.3E, b8=%.3E, b9=%.3E | ", a, b1, b2, b3, b4, b5, b6, b7, b8, b9);
		if (modulusMode == mmVojtech || modulusMode == mmVojtechHighTFixed)
			sprintf(buf_m, "a=%.3E, a2=%.3E, b1=%.3E, b2=%.3E, b3=%.3E, c1=%.3E, c2=%.3E, d1=%.3E | ", a, a2, b1, b2, b3, c1, c2, d1);

		char buf_p[200];
		if (phaseMode == pmConstant)
			sprintf(buf_p, "p0=%.3E", p0);
		if (phaseMode == pmBailly)
			sprintf(buf_p, "p0=%.3E, p_td=%.3E", p0, p_td);
		if (phaseMode == pmStandard)
			sprintf(buf_p, "p0=%.3E, p_t0=%.3E, p_tau=%.3E", p0, p_t0, p_tau);
		if (phaseMode == pmPeripheral)
			sprintf(buf_p, "p0=%.3E, p_A=%.3E, p_ka=%.3E, p_tm=%.3E", p0, p_A, p_ka, p_tm);
		if (phaseMode == pmPolynomial)
			sprintf(buf_p, "p0=%.3E, p1=%.3E, p2=%.3E, p3=%.3E", p0, p1, p2, p3);

		return std::string(buf_m) + buf_p;
	}


    virtual std::string GetModulusModeString() const
	{
		if (modulusMode == mmExp) return "exp";
		if (modulusMode == mmVojtech) return "Vojtech";
		if (modulusMode == mmVojtechHighTFixed) return "VojtechHighTFixed";

		return "unknown";
	}

    virtual std::string GetPhaseModeString() const
	{
		if (phaseMode == pmConstant) return "constant";
		if (phaseMode == pmBailly) return "Bailly";
		if (phaseMode == pmStandard) return "standard";
		if (phaseMode == pmPeripheral) return "peripheral";
		if (phaseMode == pmPolynomial) return "polynomial";

		return "unknown";
	}

	unsigned int GetNumberOfPhaseParameters() const
	{
		if (phaseMode == pmConstant) return 1;
		if (phaseMode == pmBailly) return 1;
		if (phaseMode == pmStandard) return 1;
		//if (phaseMode == pmPeripheral) return 4;
		if (phaseMode == pmPeripheral) return 1;
		if (phaseMode == pmPolynomial) return 1;
		return 0;
	}

    virtual TComplex Amp(double t) const;

	double upper_bound_t, precision_t;

	bool integ_workspace_initialized;
	unsigned long integ_workspace_size_t;
	gsl_integration_workspace *integ_workspace_t;

	static TComplex Amp_J0(double t, double *par, const void *vobj);
    virtual TComplex Prf(double b_fm) const;
};

//----------------------------------------------------------------------------------------------------

TComplex HadronicFitModel::Amp(double t) const
{
	// ---------- modulus ----------
	
	double m = 0.;

	if (modulusMode == mmExp)
	{
		// blending region
		double t_avg = (t2 + t1)/2., t_si = (t2 - t_avg) / 3.;
		double t_opt_m1 = -t_avg + 5. * t_si;
		double t_opt_m2 = -t_avg - 5. * t_si;
	
		// main modulus part (for low |t|)
		double bPol = 0., tPow = t;
		bPol += b1 * tPow; tPow *= t;
		bPol += b2 * tPow; tPow *= t;
		bPol += b3 * tPow; tPow *= t;
		bPol += b4 * tPow; tPow *= t;
		bPol += b5 * tPow; tPow *= t;
		bPol += b6 * tPow; tPow *= t;
		bPol += b7 * tPow; tPow *= t;
		bPol += b8 * tPow; tPow *= t;
		bPol += b9 * tPow; tPow *= t;
		double m1 = a * exp(bPol);
		
		if (t > t_opt_m1)
		{
			// optimisation for very low |t|: only m1 component
			m = m1;
		} else {	
			// fixed part for higher |t|
			double x = -t;
			//double m2 = sqrt( ( 380.*exp(-18.25*x) + TMath::Erf( (x - 0.546)/145.3 ) * 2223. * exp(-6.79 * x) ) / cnts->sig_fac );
			//double m2 = ( 21.1*exp(-8.50*x - 5.01*x*x) + (TMath::Erf((x - 0.55) / 0.0811) + 1.) * 0.263 * exp(-1.71*x) ) / sqrt(cnts->sig_fac);
			double m2 = ( 21.1*exp(-8.50*x - 5.01*x*x) + (TMath::Erf((x - 0.55) / 0.0811) + 1.) * 0.469 * exp(-2.288*pow(x, 0.7)) ) / sqrt(cnts->sig_fac);
			m2 *= 0.8;

			m2 *= hts;

			if (t < t_opt_m2)
			{
				// optimisation for very high |t|: only m2 component
				m = m2;
			} else {
				// full modulus (low and high-|t| parts blended)
				double p = (TMath::Erf( (x - t_avg) / t_si / sqrt(2.) ) + 1.) / 2.;
				m = (t < -1.) ? m2 : m1*(1.-p) + m2*p;
			}
		}
	}

	if (modulusMode == mmVojtech || modulusMode == mmVojtechHighTFixed)
	{
		// variable part
		double m1 = (a + a2*t) * exp(b1*t + b2*t*t + b3*t*t*t);
		double m2 = (c1 + c2*t) * exp(d1*t);
		double m_v = m1 + m2;

		if (modulusMode == mmVojtech)
			m = m_v;
		else {
			// fixed part
			double f_a = 6.610058960000e+08;
			double f_a2 = 1.633227930000e+09;
			double f_b1 = 8.287342810000e+00;
			double f_b2 = 9.302256060000e+00;
			double f_b3 = 1.466743830000e+01;
			double f_c1 = 1.709871650000e+07;
			double f_c2 = -2.731341890000e+07;
			double f_d1 = 2.706162130000e+00;
			
			m1 = (f_a + f_a2*t) * exp(f_b1*t + f_b2*t*t + f_b3*t*t*t);
			m2 = (f_c1 + f_c2*t) * exp(f_d1*t);
			
			double m_f = m1 + m2;
			
			m_f *= hts;
	
			// blending
			double x = -t;
			double t_avg = (t2 + t1)/2., t_si = (t2 - t_avg) / 3.;
			double p = (TMath::Erf( (x - t_avg) / t_si / sqrt(2.) ) + 1.) / 2.;
			m = (t < -1.) ? m_f : m_v*(1.-p) + m_f*p;
		}
	}

	// ---------- phase ----------

	double ph = 0.;
	
	if (phaseMode == pmUnknown)
		printf("ERROR in HadronicFitModel::Amp > phase mode is `unknown'.\n");

	if (phaseMode == pmConstant)
		ph = p0;

	if (phaseMode == pmBailly)
	{
		double rho = cos(p0)/sin(p0);

		double ze = 0.;
		if (rho == 0.)
			ze = 0.;
		else {
			double r = t / p_td;
			double ep = 1E-4;
			if (r < 1. - ep) ze = atan(rho / (1. - r));
			if (r > 1. - ep && r < 1. + ep) ze = M_PI/2. + (r-1.)/rho;
			if (r > 1. + ep) ze = atan(rho / (1. - r)) + M_PI;
		}

		ph = M_PI/2. - ze;
	}

	if (phaseMode == pmStandard)
	{
		ph = p0 + atan( (fabs(t) - fabs(p_t0)) / p_tau ) - atan(- fabs(p_t0) / p_tau);
	}

	if (phaseMode == pmPeripheral)
	{
		double r = fabs(t / p_tm);
		ph = p0 - p_A * exp(p_ka * (log(r) - r + 1.));
	}

	if (phaseMode == pmPolynomial)
	{
		double ph_pol = p1*t + p2*t*t + p3*t*t*t;

		double mt = -t;
		double mu = 1., si = 0.2;
		double w = (1. - TMath::Erf( (mt - mu)/(si * 1.414) )) / 2.;

		ph = p0 + ph_pol * w;
	}

	return m * TComplex::Exp(i * ph);
}

//----------------------------------------------------------------------------------------------------

TComplex HadronicFitModel::Amp_J0(double t, double *par, const void *vobj)
{
	const HadronicFitModel *obj = (HadronicFitModel *) vobj;
	const double &b = par[0];	// impact parameter in GeV^-1

	return obj->Amp(t) * TMath::BesselJ0(b * sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex HadronicFitModel::Prf(double b_fm) const
{
	double b = b_fm / cnts->hbarc;	// b in GeV^-1
	double par[] = { b };

	TComplex I = ComplexIntegrate(Amp_J0, par, this, upper_bound_t, 0., 0., precision_t,
		integ_workspace_size_t, integ_workspace_t, "HadronicFitModel::Prf");

	return I / 4. / cnts->p_cms / cnts->sqrt_s;
}

#endif
