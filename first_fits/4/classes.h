#ifndef _classes_h_
#define _classes_h_

#include <cstring>
#include <cstdio>
#include <cstdlib>

//----------------------------------------------------------------------------------------------------

struct Results
{
	unsigned int method, metric;

	// last fit
	double chi_sq, chi_sq_norm, prob, sig;
	unsigned int ndf;
	double quality;

	// preceeding fit (if applicable)
	double chi_sq_p, chi_sq_norm_p, prob_p, sig_p;
	unsigned int ndf_p;

	std::vector<double> decisivePoints;

	double p0, rho, a, B;
	double p0_e, rho_e, a_e, B_e;

	double si_el, si_inel, si_tot;				// mb
	double b_rms_el, b_rms_inel, b_rms_tot;		// fm

	Results() : method(0), metric(0),
		chi_sq(0.), chi_sq_norm(0.), prob(0.), sig(0.), ndf(0),
		quality(0.),
		chi_sq_p(0.), chi_sq_norm_p(0.), prob_p(0.), sig_p(0.),	ndf_p(0),
		p0(0.), rho(0.), a(0.), B(0.),
		p0_e(0.), rho_e(0.), a_e(0.), B_e(0.),
		si_el(0.), si_inel(0.), si_tot(0.),
		b_rms_el(0.), b_rms_inel(0.), b_rms_tot(0.)
		{}

	unsigned int Load(const std::string &fn, bool verbose=false)
	{
		FILE *f = fopen(fn.c_str(), "r");
		if (!f)
		{
			if (verbose)
				printf("ERROR in Results::Load > Can't open file `%s' for reading.\n", fn.c_str());
			return 1;
		}
	
		// determine version	
		unsigned int version = 0;

		char buf[200];
		fgets(buf, 199, f);
		if (strstr(buf, "VERSION") != NULL)
		{
			version = atoi(buf + 8);

			// pre-read next line
			fgets(buf, 199, f);
		}

		sscanf(buf, "%u%u\n", &method, &metric);
		
		if (version < 3)
		{
			fgets(buf, 199, f); sscanf(buf, "%lf%lf%lf%lf\n", &chi_sq_norm, &prob, &sig, &quality);
			chi_sq = 0.; ndf = 0;
		} else {
			fgets(buf, 199, f); sscanf(buf, "%lf%u%lf%lf%lf\n", &chi_sq, &ndf, &prob, &sig, &quality);
			chi_sq_norm = chi_sq / ndf;
			fgets(buf, 199, f); sscanf(buf, "%lf%u%lf%lf\n", &chi_sq_p, &ndf_p, &prob_p, &sig_p);
			chi_sq_norm_p = chi_sq_p / ndf_p;
		}

		fgets(buf, 199, f);
		// TODO: parse decisive points

		fscanf(f, "%lf%lf\n", &p0, &p0_e);
		fscanf(f, "%lf%lf\n", &rho, &rho_e);
		fscanf(f, "%lf%lf\n", &a, &a_e);
		fscanf(f, "%lf%lf\n", &B, &B_e);
		fscanf(f, "%lf%lf%lf\n", &si_el, &si_inel, &si_tot);
		fscanf(f, "%lf%lf%lf\n", &b_rms_el, &b_rms_inel, &b_rms_tot);
		fclose(f);

		return 0;
	}

	unsigned int Save(const std::string &fn) const
	{
		//printf(">> Results::Save\n");
		FILE *f = fopen(fn.c_str(), "w");

		if (!f)
			return 1;

		fprintf(f, "VERSION 3\n");
		fprintf(f, "%u\t%u\n", method, metric);
		fprintf(f, "%E\t%u\t%E\t%E\t%E\n", chi_sq, ndf, prob, sig, quality);
		fprintf(f, "%E\t%u\t%E\t%E\n", chi_sq_p, ndf_p, prob_p, sig_p);

		for (unsigned int i = 0; i < decisivePoints.size(); ++i)
			fprintf(f, "%E\t", decisivePoints[i]);
		fprintf(f, "\n");

		fprintf(f, "%E\t%E\n", p0, p0_e);
		fprintf(f, "%E\t%E\n", rho, rho_e);
		fprintf(f, "%E\t%E\n", a, a_e);
		fprintf(f, "%E\t%E\n", B, B_e);
		fprintf(f, "%E\t%E\t%E\n", si_el, si_inel, si_tot);
		fprintf(f, "%E\t%E\t%E\n", b_rms_el, b_rms_inel, b_rms_tot);
		fclose(f);

		return 0;
	}

	void Print() const
	{
		printf("method = %u, metric = %u\n", method, metric);
		printf("last fit     : chi^2/ndf = %5.2f/%2u = %.3f, probability = %.2E, significance = %.3f, quality = %f\n",
			chi_sq, ndf, chi_sq/ndf, prob, sig, quality);
		printf("previous fit : chi^2/ndf = %5.2f/%2u = %.3f, probability = %.2E, significance = %.3f\n",
			chi_sq_p, ndf_p, chi_sq_p/ndf_p, prob_p, sig_p);
		
		printf("decisive point |t|'s: ");
		for (unsigned int i = 0; i < decisivePoints.size(); ++i)
			printf("%E, ", decisivePoints[i]);
		printf("\n");
		
		printf("p0 = %E +- %E\n", p0, p0_e);
		printf("rho = %E +- %E\n", rho, rho_e);
		printf("a = %E +- %E\n", a, a_e);
		printf("B = %E +- %E\n", B, B_e);
		
		printf("si_el = %.2f, si_inel = %.2f, si_tot = %.2f mb\n", si_el, si_inel, si_tot);
		printf("b_rms_el = %.3f, b_rms_inel = %.3f, b_rms_tot = %.3f fm\n", b_rms_el, b_rms_inel, b_rms_tot);
	}
};

#endif
