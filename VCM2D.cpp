 #include "Headers/Include.h"

int main() {

	// ---------------------- SIMULATION SETUP ------------------------- //
	// Problem type
	const int type = 3;			// Cylindrical = 1 / Pouch Baseline = 2 / Pouch proposed = 3
	const int pouch = 1;		// DOUBLE = 1 / PCM  = 2 / CPCM = 3

	// Time Setup
	double TotalTime = 500;		// Total simulation time [s]
	double dt =	0.001;			// Time step [s]

	// -------------------------- MESH SETUP --------------------------- //
	// Control Volumes
	double dx = 0.001;			// CV dimension in x direction[m]
	double dy = 0.004;			// CV dimension in y direction[m]
	int N, Nx, Ny;
	double Lx, Ly;

	// Battery Cell
	// Cylindrical cell (Top view) - Type 1
	double D = 0.018;				// Cell diameter [m]
	double h = 0.065;				// Cell height [m]
	 
	// Pouch cell - Baseline Design (Side view) - Type 2
	double t = 0.008;				// Cell thickness [m]
	double l = 0.188;				// Cell length [m]
	double w = 0.100;				// Cell width [m]
	double t_fin = 0.001;			// Aluminum fin thickness [m]

	// Pouch cell - Proposed Design (Side view) - Type 3
	double t_pcm = 0.001;			// PCM/CPCM thickness [m]

	// ----------------- BOUNDARY AND INITIAL CONDITIONS --------------- //
	// Initial Condition
	double Tinitial = 25;			// Temperature [�C] in time = 0s
	double Q = 200000;				// W/m�

	// Temperatures for boundary conditions
	double Tn = 0.0;				// North
	double Ts = 0.0;				// South 
	double Te = 0.0;				// East
	double Tw = 0.0;				// West

	// Heat flux for boundary conditions
	double qn = 0.0;				// North
	double qs = 0.0;				// South 
	double qe = 0.0;				// East
	double qw = 0.0;				// West

	// ------------------------ PHYSICAL PROPERTIES --------------------- //
	// Option 1: RT 18 HC [T_melt = 18�C, Latent heat of fusion = 250 kJ/kg, Cp = 2 kJ/kg. K, k = 0.2 W/kg.K, rho = 880 <-> 770 kg/m� (solid/liquid)]
	// Option 2: RT 21 HC [T_melt = 21�C, Latent heat of fusion = 250 kJ/kg, Cp = 2 kJ/kg. K, k = 0.2 W/kg.K, rho = 880 <-> 770 kg/m� (solid/liquid)]

	// Battery
	double rho_bat = 2250.0;
	double kx_bat = 0.97;
	double ky_bat = 26.57;
	double cp_bat = 1230.0;
	
	// PCM
	double rho_pcm = 776.0;
	double ks_pcm = 0.425;
	double kl_pcm = 0.152;
	double k_pcm = ks_pcm;					// Initial conductivity is solid phase
	double cp_pcm = 2150.0;
	double L_pcm = 247050.0;		
	double Tmelt = 36.1;			

	// Aluminum
	double rho_alu = 2719.0;	
	double k_alu = 202.4;
	double cp_alu = 871.0;

	// Graphite
	double rho_gra = 2235.0;
	double kx_gra = 6.5;
	double ky_gra = 1659.3;
	double cp_gra = 762.5;

	// Copper
	double rho_cop = 8978.0;
	double k_cop = 387.6;
	double cp_cop = 381.0;

	// CPCM - Copper Foam
	double porosity = 0.92;			// Copper foam
	double rho_cpcm = porosity * rho_pcm + (1 - porosity) * rho_cop;
	double k_cpcm = porosity * k_pcm + (1 - porosity) * k_cop;
	double cp_cpcm = (porosity * rho_pcm * cp_pcm + (1 - porosity) * rho_cop * cp_cop) / rho_cpcm;
	double L_cpcm = (porosity * rho_pcm / rho_cpcm) * L_pcm;
	
    if (pouch == 1) // Double
    {
        rho_cpcm = porosity * rho_pcm + (1 - porosity) * rho_cop;
        k_cpcm = porosity * k_pcm + (1 - porosity) * k_cop;
        cp_cpcm = (porosity * rho_pcm * cp_pcm + (1 - porosity) * rho_cop * cp_cop) / rho_cpcm;
        L_cpcm = (porosity * rho_pcm / rho_cpcm) * L_pcm;
    }
    else if (pouch == 2) // Full PCM
    {
        rho_cpcm = rho_pcm;
        k_cpcm = k_pcm;
        cp_cpcm = cp_pcm;
        L_cpcm = L_pcm;
    }
    else if (pouch == 3) // Full CPCM
    {
        rho_cpcm = porosity * rho_pcm + (1 - porosity) * rho_cop;
        k_cpcm = porosity * k_pcm + (1 - porosity) * k_cop;
        cp_cpcm = (porosity * rho_pcm * cp_pcm + (1 - porosity) * rho_cop * cp_cop) / rho_cpcm;
        L_cpcm = (porosity * rho_pcm / rho_cpcm) * L_pcm;
        rho_pcm = rho_cpcm;
        k_pcm = k_cpcm;
        cp_pcm = cp_cpcm;
        L_pcm = L_cpcm;
    }

	// Cold Plate
	double h_cp = 2000;		// Heat transfer coefficient [W/m².K]
	double T_cp = 15;
	double m_dot = 0.008;		// Mass flow rate [kg/s]
	double cp_liq = 3285;		// Specific heat of liquid [J/kg.K]
 
	// Natural Convection
	double h_air = 2;
	double T_air = 25;

	// --------------------------- INIT SIMULATION DATA ------------------------- //
	set_mesh_problem(type, N, Nx, Ny, Lx, Ly, dx, dy, t, l, t_fin, t_pcm);
	string results = create_results_folder();

	// -------------- MEMORY ALLOCATION AND INITIALIZATION -------------- //
	int* R = (int*)malloc((N) * sizeof(int));
	int* pp = (int*)malloc((N) * sizeof(int));
	int* ee = (int*)malloc((N) * sizeof(int));
	int* ww = (int*)malloc((N) * sizeof(int));
	int* nn = (int*)malloc((N) * sizeof(int));
	int* ss = (int*)malloc((N) * sizeof(int));
	int* batt_pos = (int*)malloc((N) * sizeof(int));
	double* rho = (double*)malloc((N) * sizeof(double));
	double* cp = (double*)malloc((N) * sizeof(double));
	double* ap = (double*)malloc((N) * sizeof(double));
	double* aw = (double*)malloc((N) * sizeof(double));
	double* ae = (double*)malloc((N) * sizeof(double));
	double* an = (double*)malloc((N) * sizeof(double));
	double* as = (double*)malloc((N) * sizeof(double));
	double* Ti = (double*)malloc((N) * sizeof(double));
	double* fi = (double*)malloc((N) * sizeof(double));
	double* sp = (double*)malloc((N) * sizeof(double));
	double* su = (double*)malloc((N) * sizeof(double));
	double* kx = (double*)malloc((N) * sizeof(double));
	double* ky = (double*)malloc((N) * sizeof(double));
	double* f = (double*)malloc((N) * sizeof(double));
	double* b = (double*)malloc((N) * sizeof(double));
	double* T = (double*)malloc((N) * sizeof(double));
	double* L = (double*)malloc((N) * sizeof(double));
	double* T_fluid = (double*)malloc((N) * sizeof(double));;

	for (int i = 0; i < N; i++) 
	{
		R[i] = 0;
		pp[i] = 0;
		ee[i] = 0;
		ww[i] = 0;
		nn[i] = 0;
		ss[i] = 0;
		rho[i] = 0.0;
		cp[i] = 0.0;
		ap[i] = 0.0;
		aw[i] = 0.0;
		ae[i] = 0.0;
		an[i] = 0.0;
		as[i] = 0.0;
		Ti[i] = 0.0;
		fi[i] = 0.0;
		sp[i] = 0.0;
		su[i] = 0.0;
		kx[i] = 0.0;
		ky[i] = 0.0;
		f[i] = 0.0;
		b[i] = 0.0;
		L[i] = 0.0;
		T[i] = Tinitial;
		batt_pos[i] = 0;
		T_fluid[i] = T_cp;
	}	

	// -------------------------- INITIALIZATION ----------------------- //
	int it = 0;						// Iterations counter
	double time = 0.0;				// Current time step [s] 
	bool print_now = true;
	bool cooling_active = false;	// Liquid cooling active
	int previous_time = 0;

	// Simulation timer
	clock_t start, end;

	// Map mesh
	int o = path(Nx, Ny, pp, ee, ww, nn, ss);
	int lines_batt_mod = 1;
	int cols_batt_mod = 8;
	map_mesh(type, o, Nx, Ny, dx, dy, D, t, l, t_fin, t_pcm, 
		kx_bat, ky_bat, k_pcm, k_alu, k_cpcm, kx_gra, ky_gra, 
		rho_bat, rho_pcm, rho_alu, rho_cpcm, rho_gra, cp_bat, 
		cp_pcm, cp_alu, cp_cpcm, cp_gra, L_pcm, L_cpcm, pp, R, 
		kx, ky, rho, cp, L, lines_batt_mod, cols_batt_mod, batt_pos);

	// Store initial mesh info
	log(time, o, pp, R, T, f, results);
	plot_mesh(o, Nx, Ny,time, dx, dy, pp, ww, ee, nn, ss, R, kx, ky, rho, cp, L, batt_pos, results);
	plot_sim(o, Nx, Ny, dx, dy, pp, T, f, time, results);
	//plot_coef(o, Nx, Ny, dx, dy, pp, aw, ae, an, as, ap, b, time, results);

	// Res�duos (SOR, liquid fraction and conductivity)
	double resmax = 1.0E-4;
	double resf = 1.0;
	double resk = 0.0;

	// Start simulation!
	start = clock();		// Start timer!

	while (time < TotalTime) 
	{

		if (print_now)
		{
			printf("\n// -------------------- Time Step = %5.3fs -------------------- // \n", time);
		}

		for (int i = 0; i < N; i++)
		{
			// Copy previous time step data and check liquid cooling condition
			Ti[i] = T[i];
			fi[i] = f[i];

			if (fi[i] > 0 && cooling_active == false && type == 3)
			{
				cooling_active = true;
				printf("Liquid cooling activated\n", time);
				break;
			}
		}

		// Compute nonlinear conducticity
		//nonlinear_cond(o, pp, R, f, kx, ky, &resk);

		if (cooling_active == true)
		{
			fluid_temp(o, pp, N, Ny, T, batt_pos, T_fluid, T_cp, h_cp, m_dot, t, w, cp_liq, cols_batt_mod);
		}
		
		// Init residue and number of iterations
		resf = 1.0;
		resk = 0.0;
		it = 0;

		while (resf > resmax || resk > resmax)
		{
			// Coefficients matrix setup
			assembly(o, pp, type, Nx, Ny, dx, dy, ww, ee, nn, ss, kx, ky, ap, aw, ae, an, as, su, sp, b, Tw, Te, Tn, Ts, qw, qe, qn, qs, T, Ti, rho, cp, L, f, fi, R, Q, dt, h_cp, T_cp, h_air, T_air, w, cooling_active, T_fluid);

			// Linear system solver
			SORt(o, pp, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti);
			//SIP(o, N, pp, nn, ss, ee, ww, ap, ae, aw, an, as, b, T, Nx, Ny, dx, dy, time, results);	

			// Compute liquid fraction
			SORf(o, pp, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti, R, f, fi, rho, L, dt, Tmelt, &resf);

			// Compute nonlinear conducticity
			//nonlinear_cond(o, pp, R, f, kx, ky, &resk);

			it++;
			if (print_now)
			{
				printf("Iteration = %i\t resf = %5.1E\t resk = %5.1E\n", it, resf, resk);
			}
		}
		
		time = time + dt;
		update_thermal_conductivity(R, o, pp, f, kx, ky, ks_pcm, kl_pcm, k_cop, porosity, k_pcm, k_cpcm);
		
		// Print time step info
		//log_plot(time, o, pp, R, T, f);
		log(time, o, pp, R, T, f, results);
		if ((previous_time != int(time)) && (int(time) % 5 == 0))
		{
			plot_sim(o, Nx, Ny, dx, dy, pp, T, f, time, results);
			plot_mesh(o, Nx, Ny, time, dx, dy, pp, ww, ee, nn, ss, R, kx, ky, rho, cp, L, batt_pos, results);
		}

		if ((previous_time != int(time)) && (int(time) % 5 == 0))
		{
			print_now = true;
			previous_time = int(time);
		}
		else
		{
			previous_time = int(time);
			print_now = false;
		}
	}

	end = clock(); // Stop timer!

	// -------------------------------------------------------------------//
	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	cout << "\nTime taken by program is : " << time_taken << " s " << endl;
	// -------------------------------------------------------------------//

	average_temperature(o, pp, T);

	free(rho);
	free(pp);
	free(ee);
	free(ww);
	free(nn);
	free(ss);
	free(cp);
	free(ap);
	free(aw);
	free(ae);
	free(an);
	free(as);
	free(Ti);
	free(fi);
	free(sp);
	free(su);
	free(kx);
	free(ky);
	free(R);
	free(b);
	free(T);
	free(f);

	return 0;
}