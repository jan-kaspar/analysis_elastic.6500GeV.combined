#include "TComplex.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#ifndef _numerical_integration_h_
#define _numerical_integration_h_

//----------------------------------------------------------------------------------------------------

typedef double (* RealFunction)(double x, double *par, const void *obj);

typedef TComplex (* ComplexFunction)(double x, double *par, const void *obj);

//----------------------------------------------------------------------------------------------------

struct RealIntegPar
{
	RealFunction fcn;
	double *parameters;
	const void *object;

	RealIntegPar(RealFunction _f, double *_p, const void *_o) : fcn(_f), parameters(_p), object(_o) {}
};

//----------------------------------------------------------------------------------------------------

double RealIntegFcn(double x, void *vpar)
{
	RealIntegPar *par = (RealIntegPar *) vpar;
	return par->fcn(x, par->parameters, par->object);
}

//----------------------------------------------------------------------------------------------------

double RealIntegrate(RealFunction fcn, double *par, const void *object,
	double from, double to,
	double abs_err, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space, const char *errorLabel="")
{
	// prepare structures
	RealIntegPar ocip(fcn, par, object);

	gsl_function F;
  	F.function = RealIntegFcn;
  	F.params = &ocip;

	// real part
	double result, unc;
	int status = gsl_integration_qag(&F, from, to, abs_err, rel_err, work_space_size, GSL_INTEG_GAUSS15, work_space, &result, &unc);
	if (status != 0) 
	{
		printf("WARNING in %s > Integration failed: %s.\n", errorLabel, gsl_strerror(status));
	}

	return result;
}

//----------------------------------------------------------------------------------------------------

struct OneCompIntegPar
{
	ComplexFunction fcn;
	double *parameters;
	const void *object;
	enum Part { pUndefined, pReal, pImaginary } part;

	OneCompIntegPar(ComplexFunction _f, double *_p, const void *_o, Part _pa = pUndefined) : fcn(_f), parameters(_p), object(_o), part(_pa) {}
};

//----------------------------------------------------------------------------------------------------

double OneCompIntegFcn(double x, void *par)
{
	OneCompIntegPar *ocip = (OneCompIntegPar *) par;
	TComplex c = ocip->fcn(x, ocip->parameters, ocip->object);

	if (ocip->part == OneCompIntegPar::pReal)
		return c.Re();
	if (ocip->part == OneCompIntegPar::pImaginary)
		return c.Im();

	return 0.;
}

//----------------------------------------------------------------------------------------------------

TComplex ComplexIntegrate(ComplexFunction fcn, double *par, const void *object,
	double from, double to,
	double abs_err, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space, const char* errorLabel="")
{
	// prepare structures
	OneCompIntegPar ocip(fcn, par, object);

	gsl_function F;
  	F.function = OneCompIntegFcn;
  	F.params = &ocip;

	// real part
	ocip.part = OneCompIntegPar::pReal;
	double result_re, unc_re;
	int status_re = gsl_integration_qag(&F, from, to, abs_err, rel_err, work_space_size, GSL_INTEG_GAUSS61, work_space, &result_re, &unc_re);
	if (status_re != 0) 
	{
		printf("WARNING in %s > Real integration failed: %s.\n", errorLabel, gsl_strerror(status_re));
	}

	// imaginary part
	ocip.part = OneCompIntegPar::pImaginary;
	double result_im, unc_im;
	int status_im = gsl_integration_qag(&F, from, to, abs_err, rel_err, work_space_size, GSL_INTEG_GAUSS61, work_space, &result_im, &unc_im);
	if (status_im != 0) 
	{
		printf("WARNING in %s > Imaginary integration failed: %s.\n", errorLabel, gsl_strerror(status_im));
	}

	// result
	return TComplex(result_re, result_im);
}

#endif
