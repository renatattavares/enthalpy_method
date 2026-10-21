#include "..\Headers\Include.h"

void map_mesh(int type, int o, int Nx, int Ny, double dx, double dy, double D, double t, double l, double t_fin, double t_pcm, double kx_bat, double ky_bat, double k_pcm, double k_alu, double k_cpcm, double kx_gra, double ky_gra, double rho_bat, double rho_pcm, double rho_alu,	double rho_cpcm, double rho_gra, double cp_bat, double cp_pcm, double cp_alu, double cp_cpcm, double cp_gra, double L_pcm, double L_cpcm, int* pp, int* R, double* kx, double* ky, double* rho, double* cp, double* L, int lines_batt_mod, int cols_batt_mod, int* batt_pos)
{
	// Simulation timer
	clock_t start, end;
	start = clock();		// Start timer!

	if (type == 1) // Cylindrical cell 
	{
		int a, i, l, m;
		double x_center, y_center, x, y, radius;
		double cell_radius = D / 2;
		double x_cell = Nx * dx / 2;
		double y_cell = Ny * dy / 2;

		for (a = 0; a < o; a++) 
		{
			m = pp[a];
			i = m / (Ny + 2);
			l = m - i * (Ny + 2);

			// Mesh Regions and Properties
			x_center = (dx * 0.5) + (i-1) * dx;
			y_center = (dy * 0.5) + (l-1) * dy;
			x = fabs(x_center - x_cell);
			y = fabs(y_center - y_cell);
			radius = sqrt(pow(x, 2) + pow(y, 2));
			
			if (radius <= cell_radius)
			{
				R[m] = 1; // Battery Region
				kx[m] = kx_bat;
				ky[m] = ky_bat;
				rho[m] = rho_bat;
				cp[m] = cp_bat;
				L[m] = 0.0;
			}
			else 
			{
				R[m] = 0; // PCM Region
				kx[m] = k_pcm;
				ky[m] = k_pcm;
				rho[m] = rho_pcm;
				cp[m] = cp_pcm;
				L[m] = L_pcm;
			}
			//printf("CV = %4i\t x_center = %4.5f\t y_center = %4.5f\t radius = %4.5f\t R = %4i\n", m, x_center, y_center, radius, R[m]);
		}
	}
	else if(type == 2) // Pouch cell - Baseline design 
	{
		int i, l, a, m, b;
		double x_center;
		double x_batt[16];
		int batts = 8;
		int Nx_ghosts = Nx + 2;
		int	Ny_ghosts = Ny + 2;

		// Battery regions
		for (int i = 0; i < batts; i++)
		{
			x_batt[2 * i] = (t_fin) + i * (t_fin + t);
			x_batt[(2 * i) + 1] = (t_fin + t) + i * (t_fin + t);
		}

		for (b = 0; b < o; b++) 
		{
			m = pp[b];
			i = (int)m / Ny_ghosts;
			l = (int)m - i * Ny_ghosts;

			// Mesh Regions and Properties
			x_center = (dx * 0.5) + (i-1) * dx;

			for (a = 0; a < batts; a++)
			{
				if ((x_center > x_batt[2 * a]) && (x_center < x_batt[(2 * a) + 1]))
				{
					R[m] = 1; // Battery Region
					kx[m] = kx_bat;
					ky[m] = ky_bat;
					rho[m] = rho_bat;
					cp[m] = cp_bat;
					L[m] = 0.0;
					break;
				}
				else {
					R[m] = 3; // Aluminum Region
					kx[m] = k_alu;
					ky[m] = k_alu;
					rho[m] = rho_alu;
					cp[m] = cp_alu;
					L[m] = 0.0;
				}
			}
		}
	}
	else if (type == 3 || type == 4) // Battery module
	{
		if (type == 3)
		{
			cols_batt_mod = 8;
			lines_batt_mod = 1;
		}

		int i, j, a, b, c, m, z;
		double x, y;
		double x_center, y_center;
		double offsetx, offsety;
		int Nx_ghosts = Nx + 2;
		int	Ny_ghosts = Ny + 2;

		double block_width = t_fin + 2.0 * t_pcm + t;
		double totalX = cols_batt_mod * block_width + t_fin;
		double halfH = 0.5 * (static_cast<double>(lines_batt_mod) * l);
		double a_fin = t_fin;
		double a_pouch1 = a_fin + t_pcm;
		double a_cell = a_pouch1 + t;
		double a_pouch2 = a_cell + t_pcm;

		for (z = 0; z < o; z++)
		{
			m = pp[z];
			i = (int)m / Ny_ghosts;
			j = (int)m - i * Ny_ghosts;

			x_center = (dx * 0.5) + (i - 1) * dx;
			y_center = (dy * 0.5) + (j - 1) * dy;

			batt_pos[m] = min(int(floor(x_center / block_width) + 1), cols_batt_mod);
			offsetx = floor(x_center / block_width) * block_width;
			offsety = floor(y_center / l) * l;

			x = x_center - offsetx; // X center
			y = y_center - offsety; // Y center

			// decidir subsegmento por compara��es simples
			if (x > 0 and x < a_fin)
			{
				R[m] = 6; // Graphite Region
				kx[m] = kx_gra;
				ky[m] = ky_gra;
				rho[m] = rho_gra;
				cp[m] = cp_gra;
			}
			else if (x > a_fin and x < a_pouch1)
			{
				// pouch entre fin e cell -> PCM/CPCM por Y
				if (y <= l/2)
				{
					R[m] = 0; // PCM Region
					kx[m] = k_pcm;
					ky[m] = k_pcm;
					rho[m] = rho_pcm;
					cp[m] = cp_pcm;
					L[m] = L_pcm;
				}
				else
				{
					R[m] = 5; // CPCM Region
					kx[m] = k_cpcm;
					ky[m] = k_cpcm;
					rho[m] = rho_cpcm;
					cp[m] = cp_cpcm;
					L[m] = L_cpcm;
				}
			}
			else if (x > a_pouch1 and x < a_cell)
			{
				R[m] = 1; // Battery Region
				kx[m] = kx_bat;
				ky[m] = ky_bat;
				rho[m] = rho_bat;
				cp[m] = cp_bat;
			}
			else if (x > a_cell and x < a_pouch2) 
			{
				if (y <= l/2)
				{
					R[m] = 0; // PCM Region
					kx[m] = k_pcm;
					ky[m] = k_pcm;
					rho[m] = rho_pcm;
					cp[m] = cp_pcm;
					L[m] = L_pcm;
				}
				else
				{
					R[m] = 5; // CPCM Region
					kx[m] = k_cpcm;
					ky[m] = k_cpcm;
					rho[m] = rho_cpcm;
					cp[m] = cp_cpcm;
					L[m] = L_cpcm;
				}
			}
			else if (x > a_pouch2)
			{
				R[m] = 6; // Graphite Region
				kx[m] = kx_gra;
				ky[m] = ky_gra;
				rho[m] = rho_gra;
				cp[m] = cp_gra;
			}
		}
	}
	else if (type == 5)
	{
		int i, j, a, b, c, m, z;
		double x;
		double x_center;
		double offsetx;
		int Nx_ghosts = Nx + 2;
		int	Ny_ghosts = Ny + 2;

		double block_width = t_fin + t;
		double totalX = cols_batt_mod * block_width + t_fin;
		double a_fin = t_fin;
		double a_cell = a_fin + t;

		for (z = 0; z < o; z++)
		{
			m = pp[z];
			i = (int)m / Ny_ghosts;
			j = (int)m - i * Ny_ghosts;

			x_center = (dx * 0.5) + (i - 1) * dx;

			batt_pos[m] = min(int(floor(x_center / block_width) + 1), cols_batt_mod);
			offsetx = floor(x_center / block_width) * block_width;

			x = x_center - offsetx; // X center

			// decidir subsegmento por compara��es simples
			if (x > 0 and x < a_fin)
			{
				R[m] = 6; // Graphite Region
				kx[m] = kx_gra;
				ky[m] = ky_gra;
				rho[m] = rho_gra;
				cp[m] = cp_gra;
			}
			else if (x > a_fin and x < a_cell)
			{
				R[m] = 1; // Battery Region
				kx[m] = kx_bat;
				ky[m] = ky_bat;
				rho[m] = rho_bat;
				cp[m] = cp_bat;
			}
			else if (x > a_cell)
			{
				R[m] = 6; // Graphite Region
				kx[m] = kx_gra;
				ky[m] = ky_gra;
				rho[m] = rho_gra;
				cp[m] = cp_gra;
			}
		}
	}

	end = clock(); // Stop timer!
	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	cout << "\nMap time: " << time_taken << " s " << endl;
}

void assembly(int o, int* pp, int type, int Nx, int Ny, double dx, double dy, int* ww, int* ee, int* nn, int* ss, double* kx, double* ky, double* ap, double* aw, double* ae, double* an, double* as, double* su, double* sp, double* b, double Tw, double Te, double Tn, double Ts, double qw, double qe, double qn, double qs, double* T, double* Ti, double* rho, double* cp, double* L, double* f, double* fi, int* R, double Q, double dt, double h_cp, double T_cp, double h_air, double T_air, double w, bool active, double* T_fluid) 
{

	int i, j, m, l;
		
	#pragma omp parallel for shared (ap, b) private (l, i, j)
	for (m = 0; m < o; m++) 
	{

		l = pp[m];
		i = (int)l / (Ny + 2);
		j = (int)l - i * (Ny + 2);

		su[l] = 0.0;
		sp[l] = 0.0;

		aw[l] = (2 * kx[ww[l]] * kx[l]) / (dx * (kx[ww[l]] * dx + kx[l] * dx));
		ae[l] = (2 * kx[ee[l]] * kx[l]) / (dx * (kx[ee[l]] * dx + kx[l] * dx));
		an[l] = (2 * ky[nn[l]] * ky[l]) / (dy * (ky[nn[l]] * dy + ky[l] * dy));
		as[l] = (2 * ky[ss[l]] * ky[l]) / (dy * (ky[ss[l]] * dy + ky[l] * dy));

		// Boundary conditions - West
		if (i == 1)
		{
			aw[l] = 0;
			su[l] = su[l] + (h_air * T_air) / dx;
			sp[l] = sp[l] - (h_air) / dx;
		}
		else if (i == Nx) // Boundary conditions - East
		{
			ae[l] = 0;
			su[l] = su[l] + (h_air * T_air) / dx;
			sp[l] = sp[l] - (h_air) / dx;
		}

		// Boundary conditions - North
		if (j == Ny)
		{
			an[l] = 0;
			su[l] = su[l] + (h_air * T_air) / dy;
			sp[l] = sp[l] - (h_air) / dy;

		}
		else if (j == 1) // Boundary conditions - South
		{
			as[l] = 0;
			if ((type == 3 && active == true) || (type == 2))
			{
				su[l] = su[l] + (h_cp * 0.12 * T_fluid[l]) / dy; //su[l] = su[l] + (h_cp * T_cp) / dy;
				sp[l] = sp[l] - (h_cp * 0.12) / dy;
			}
		}

		//Boundary conditions - Source term
		if (R[l] == 1)
		{
			su[l] = su[l] + Q;
		}

		//Calculo do ap e do b
		ap[l] = aw[l] + ae[l] + an[l] + as[l] + ((rho[l] * cp[l]) / dt) - sp[l];
		b[l] = su[l] + ((rho[l] * cp[l] * Ti[l]) / dt) - ((rho[l] * L[l] * (f[l] - fi[l])) / dt);
		//printf("pp = %i\t aw = %5.2E\t ae = %5.2E\t an = %5.2E\t as = %5.2E\t ap = %5.2E\t su = %5.2E\t sp = %5.2E\n", l, aw[l], ae[l], an[l], as[l], ap[l], su[l], sp[l]);
	}
}


void set_mesh_problem(int type, int& N, int& Nx, int& Ny, double& Lx, double& Ly,
	double dx, double dy, double t, double l, double t_fin, double t_pcm) {

	if (type == 1)
	{
		Lx = 0.05;			// Mesh size in x direcition [m] 
		Ly = 0.05;			// Mesh size in y direcition [m]
	}
	else if (type == 2) 
	{	
		int batts = 8;
		int fins = 9;
		Lx = batts * t + fins * t_fin;			// Mesh size in x direcition [m] 
		Ly = l;									// Mesh size in y direcition [m]
	}
	else if (type == 3)
	{
		int batts = 8;
		int fins = 9;
		int pcms = 16;
		Lx = batts * t + fins * t_fin + pcms * t_pcm;	// Mesh size in x direcition [m] 
		Ly = l;											// Mesh size in y direcition [m]
	}
	Nx = (int)(Lx / dx);				// Volumes in x direction
	Ny = (int)(Ly / dy);				// Volumes in y direction 
	N = (Nx + 2) * (Ny + 2);			// Total number of CVs. Include ghost elements in boundaries.

	cout << "Mesh length in x direction [m]: " << Lx << endl;
	cout << "Mesh length in y direction [m]: " << Ly << endl;
	cout << "CV dimension in x direction [m]: " << dx << endl;
	cout << "CV dimension in y direction [m]: " << dy << endl;
	cout << "Number of cells in x direction: " << Nx << endl;
	cout << "Number of cells in y direction: " << Ny << endl;
	cout << "Total number of true cells: " << Nx * Ny << endl;
	cout << "Total number of cells: " << N << endl;
}


int path(int Nx, int Ny, int* pp, int* ee, int* ww, int* nn, int* ss) {

	int i, j, l, o, m;

	o = 0;

	for (i = 1; i <= Nx; i++)
	{
		for (j = 1; j <= Ny; j++)
		{
			pp[o] = i * (Ny + 2) + j;
			//printf("o = %i\t pp = %i\t i = %i\t j = %i\n", o, pp[o], i, j);
			o++;
		}
	}

	for (m = 0; m < o; m++)
	{
		l = pp[m];
		i = (int)l / (Ny + 2);
		j = (int)l - i * (Ny + 2);
		
		ee[l] = (i + 1) * (Ny + 2) + j;
		ww[l] = (i - 1) * (Ny + 2) + j;
		nn[l] = i * (Ny + 2) + (j + 1);
		ss[l] = i * (Ny + 2) + (j - 1);
		//printf("pp = %i\t ww = %5.1i\t ee = %5.1i\t nn = %5.1i\t ss = %5.1i\n", l, ww[l], ee[l], nn[l], ss[l]);
	}
	return o;
}

void nonlinear_cond(int o, int* pp, int* R, double* f, double* kx, double* ky, double* resk)
{
	int j;
	double res = 0.0, ka;
	double a = 1.0;

	for (int i = 0; i < o; i++)
	{
		if (pp[i] != 0 && R[pp[i]] == 0) // PCM
		{
			j = pp[i];
			ka = kx[j];

			kx[j] = 0.425 * (1 - f[j]) + 0.152 * f[j]; // Conductivity based on liquid fraction
			ky[j] = kx[j];

			res = res + fabs(kx[j] - ka);
			if (fabs(kx[j] - ka) > 1.0E-05)
			{
				a++;
			}
		}
		else if (pp[i] != 0 && R[pp[i]] == 5) // PCM
		{
			j = pp[i];
			ka = kx[j];

			kx[j] = 31.39 * (1 - f[j]) + 31.15 * f[j]; // Conductivity based on liquid fraction
			ky[j] = kx[j];

			res = res + fabs(kx[j] - ka);
			if (fabs(kx[j] - ka) > 1.0E-05)
			{
				a++;
			}
		}
	}

	*resk = res / a;
}

void update_thermal_conductivity(int* R, int o, int* pp, double* f, double* kx, double* ky, double ks_pcm, double kl_pcm, double k_cop, double porosity, double& k_pcm, double& k_cpcm)
{
	int i, l;

	for (i = 0; i < o; i++)
	{
		if (pp[i] != 0)
		{
			l = pp[i];
			if (R[l] == 0) // PCM Region
			{
				k_pcm = ks_pcm * (1 - f[l]) + kl_pcm * f[l];
				kx[l] = k_pcm;
				ky[l] = k_pcm;
			}
			else if (R[l] == 5) // CPCM Region
			{
				k_pcm = ks_pcm * (1 - f[l]) + kl_pcm * f[l];
				k_cpcm = porosity * k_pcm + (1 - porosity) * k_cop;
				kx[l] = k_cpcm;
				ky[l] = k_cpcm;
			}
		}
	}
}

void fluid_temp(int o, int* pp, int N, int Ny, double* T, int* batt_pos, double* T_fluid, double T_cp, double h_cp, double m_dot, double t, double w, double cp_liq, int cols_batt_mod)
{
	int i, j, m, l, pos;
	vector<double> Tave(cols_batt_mod);
	vector<int> num_cells(cols_batt_mod);
	vector<double> Tcorner(cols_batt_mod+1);

	// Simulation timer
	clock_t start, end;
	start = clock();		// Start timer!

	for (i = 0; i < cols_batt_mod; i++)
	{
		Tave[i] = 0.0;
		num_cells[i] = 0;
	}

	#pragma omp parallel for shared (Tave, num_cells) private (l, i, j)
	for (m = 0; m < o; m++)
	{
		l = pp[m];
		i = (int)l / (Ny + 2);
		j = (int)l - i * (Ny + 2);

		if (j == 1) // Boundary conditions - South
		{
			Tave[batt_pos[l]-1] += T[l];
			num_cells[batt_pos[l] - 1] += 1;
			//printf("pp = %i \t Batt Pos = %i\t Tave = %5.E\t num_cells = %i\n", l, batt_pos[l], Tave[batt_pos[l] - 1], num_cells[batt_pos[l] - 1]);
		}
	}

	for (i = 0; i < cols_batt_mod; i++)
	{
		Tave[i] = Tave[i] / num_cells[i];
		//printf("Batt cell = %i \t Tave = %5.3f\n", i, Tave[i]);
	}

	for (i = 0; i < (cols_batt_mod+1); i++)
	{
		if (i == 0)
		{
			Tcorner[i] = T_cp;
		}
		else
		{
			Tcorner[i] = Tcorner[i - 1] + (((h_cp * 0.12 * t * w) / (m_dot * cp_liq)) * (Tave[i-1] - Tcorner[i - 1]));
			//Tcorner[i] = Tcorner[i - 1] + (((h_cp * t * w) / (m_dot * cp_liq)) * (Tave[i-1] - Tcorner[i - 1]));
		}
		//printf("i = %i\t Tcorner = %5.3f\n", i, Tcorner[i]);
	}

	#pragma omp parallel for shared (T_fluid) private (l, i, j)
	for (m = 0; m < o; m++)
	{
		l = pp[m];
		i = (int)l / (Ny + 2);
		j = (int)l - i * (Ny + 2);

		if (j == 1) // Boundary conditions - South
		{
			pos = batt_pos[l];
			T_fluid[l] = (Tcorner[pos - 1] + Tcorner[pos]) / 2;
			//printf("i = %i\t Tfluid = %5.3f\n", i, Tcorner[i]);
		}
	}
	end = clock(); // Stop timer!
	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	//cout << "\nFluid cauculus time: " << time_taken << " s " << endl;
}
