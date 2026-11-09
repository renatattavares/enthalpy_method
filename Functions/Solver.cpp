#include "..\Headers\Include.h"

void SORt(int o, int* pp, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti) {


	int i, j, it = 0;
	double R, T_old;
	double res = 1.0;
	double kappa = 1000;
	double resmax = 1E-6;

	while (res > resmax)
	{
		res = 0.0;
		#pragma omp parallel for shared (res) private (j, T_old, R)
		for (i = 0; i < o; i++)
		{
			if (pp[i] != 0)
			{
				j = pp[i];
				T_old = T[j];
				R = b[j] + aw[j] * T[ww[j]] + ae[j] * T[ee[j]] + an[j] * T[nn[j]] + as[j] * T[ss[j]];
				T[j] = (T_old + kappa * (R / ap[j])) / (1 + kappa);
				//printf("CV = %i\t w = %5.1E\t e = %5.1E\t n = %5.1E\t s = %5.1E\n", i, aw[i] * T[ww[i]], ae[i] * T[ee[i]], an[i] * T[nn[i]], as[i] * T[ss[i]]);
				res = res + fabs(T[j] - T_old);
			}
		}
		res = res / o;
		//it++;
		//printf("Iteracao = %i\t Residuo = %13.5E\n", it, res);			  
	}
}


void SORf (int o, int* pp, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti, int* R, double* f, double* fi, double* rho, double* L, double dt, double Tmelt, double* resf) {

	int i, j;
	int N = 1;
	int can_melt = 1;
	double res = 0.0;
	double total_f = 0.0;
	double kappa = 0.23;				// Under relaxation factor
	double Y, C;
	double f_old;

	for (i = 0; i < o; i++)
	{
		if (pp[i] != 0)
		{
			j = pp[i];
			
			if (R[j] == 0 || R[j] == 5)
			{			
				f_old = f[j];
				Y = aw[j] * T[ww[j]] + ae[j] * T[ee[j]] + an[j] * T[nn[j]] + as[j] * T[ss[j]];
				C = b[j] + ((rho[j] * L[j]) / dt) * f_old;
				f[j] = (f_old + (kappa * (dt / (rho[j] * L[j])) * (Y + C - ap[j] * Tmelt))) / (1 + kappa);
			}
			else
			{
				f_old = 0.0;
				f[j] = 0.0;
				continue;
			}

			// Impose physical limit
			if (f[j] > 1.0)
			{
				f[j] = 1.0;
			}
			else if (f[j] < 0.0)
			{
				f[j] = 0.0;
			}

			// Residue
			res = res + fabs(f[j] - f_old);
			if (fabs(f[j] - f_old) > 1.0E-05)
			{
				N++;
			}
			// Behavior of F
			total_f = total_f + f[j];
			can_melt++;
		}
	}
	*resf = res / N;
	total_f = total_f / can_melt;
	//printf("Total f = %5.1E\n", total_f);
}

void SIP(int o, int N, int* pp, int* nn, int* ss, int* ee, int* ww, double* ap, double* ae, double* aw, double* an, double* as, double* b, double* T, double Nx, double Ny, double dx, double dy, double time, string folder)
{
	int i, l, m, L;
	int it = 0;
	int maxit = 1;
	double alpha = 0.95;									// TODO: Testar com valores mais altos = 0.95
	double res, resm, T_old;
	double resmax = 1E-4;
	double* c = (double*)malloc((N) * sizeof(double));
	double* d = (double*)malloc((N) * sizeof(double));
	double* e = (double*)malloc((N) * sizeof(double));
	double* r = (double*)malloc((N) * sizeof(double));
	double* s = (double*)malloc((N) * sizeof(double));
	double* u = (double*)malloc((N) * sizeof(double));
	double* res_plot = (double*)malloc((N) * sizeof(double));

	for (m = 0; m < N; m++)
	{
		c[m] = (double)0.0;
		d[m] = (double)0.0;
		e[m] = (double)0.0;
		r[m] = (double)0.0;
		s[m] = (double)0.0;
		u[m] = (double)0.0;
		res_plot[m] = (double)0.0;
	}

	//Alocação dos coeficientes
	for (m = 0; m < o; m++)
	{
		if (pp[m] != 0)
		{ 
			l = pp[m];
			c[l] = as[l] / (1.0 + alpha * s[ww[l]]);
			d[l] = aw[l] / (1.0 + alpha * u[ss[l]]);
			e[l] = 1.0 / (alpha * (c[l] * s[ww[l]] + d[l] * u[ss[l]]) - (d[l] * s[ss[l]] + c[l] * u[ww[l]]) - ap[l]);
			s[l] = (ae[l] - alpha * c[l] * s[ww[l]]) * e[l];
			u[l] = (an[l] - alpha * c[l] * s[ww[l]]) * e[l];
			//printf("CV = %i\t c = %5.1E\t d = %5.1E\t e = %5.1E\t s = %5.5E\t u = %5.1E\n", l, c[l], d[l], e[l], s[l], u[l]);
		}
	}

	//Processo Iterativo
	for (L = 1; L <= maxit; L++)
	{
		resm = 0.0;

		//Forwards Step
		for (m = 0; m < o; m++)
		{
			if (pp[m] != 0)
			{
				l = pp[m];
				res = -b[l] - (aw[l] * T[ww[l]] + ae[l] * T[ee[l]] + as[l] * T[ss[l]] + an[l] * T[nn[l]] - ap[l] * T[l]);
				r[l] = e[l] * (res - (c[l] * r[ww[l]] + d[l] * r[ss[l]]));
				//printf("pp = %i\t r = %5.1E\t res = %5.1E\t c = %5.1E\t d = %5.1E\t e = %5.1E\t r[ww[l]] = %5.1E\t r[ss[l]] = %5.1E\n", l, r[l], res, c[l], d[l], e[l], r[ww[l]], r[ss[l]]);
			}
		}

		//Backwards Step
		for (m = o - 1; m >= 0; m--)
		{
			if (pp[m] != 0)
			{
				l = pp[m];
				r[l] = r[l] - u[l] * r[ee[l]] - s[l] * r[nn[l]];
				T_old = T[l];
				T[l] = T[l] + r[l];

				resm = resm + fabs(T[l] - T_old);
				res_plot[l] = fabs(T[l] - T_old);
			}
		}

		resm = resm / o;
		//plot_res(o, Nx, Ny, dx, dy, pp, res_plot, it, time, folder);
		//printf("Iteracao = %i\t Residuo = %13.5E\n", it, resm);
		it++;

		if (resm < resmax) 
		{
			break;
		}
	}
	free(c);
	free(d);
	free(e);
	free(r);
	free(s);
	free(u);
	return;
}

void average_temperature(int o, int* pp, double* T)
{
	int l;
	double total = 0.0;

	for (int i = 0; i < o; i++)
	{
		l = pp[i];
		total = total + T[l];
	}
	printf("Temperature average = %f", total / o);
}