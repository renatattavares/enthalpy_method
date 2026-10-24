 #include "Headers/Include.h"

int main() 
{
	int type, pouch;
	int N, Nx, Ny;
	int lines_batt_mod, cols_batt_mod;
	int tr_cell;
	bool variable_resistance, charge;
	double tr_q_dot;
	double Lx, Ly;
	double TotalTime, dt, SimTotalTime;
	double dx, dy;
	double D, h, t, l, w, t_fin, t_pcm;
	double Tinitial, Q;
	double Tn, Ts, Te, Tw;
	double qn, qs, qe, qw;
	double rho_bat, kx_bat, ky_bat, cp_bat;
	double rho_pcm, k_pcm, ks_pcm, kl_pcm,cp_pcm, L_pcm, Tmelt;
	double rho_cpcm, k_cpcm, cp_cpcm, L_cpcm;
	double rho_alu, k_alu, cp_alu;
	double rho_gra, kx_gra, ky_gra, cp_gra;
	double rho_cop, k_cop, cp_cop;
	double porosity;
	double h_cp, T_cp, h_air, T_air;
	double m_dot, cp_liq;
	double tr_time, tr_duration;
	bool tr_active, cooling_mission;
	double SOCinit, cell_capacity, charge_rate, max_charge_temp;
	vector<double> times;
	vector<double> q_dots;
	vector<double> current;
	
	read_data("Input.yaml", type, pouch, TotalTime, dt, dx, dy, D, h, t, l, w, t_fin, t_pcm, Tinitial, Q, Tn, Ts, Te, Tw, qn,  qs, qe,  qw, rho_bat, kx_bat, ky_bat, cp_bat, rho_pcm, k_pcm, ks_pcm, kl_pcm, cp_pcm, L_pcm, Tmelt, rho_alu, k_alu, cp_alu, rho_gra,  kx_gra,  ky_gra, cp_gra, rho_cop,  k_cop,  cp_cop, porosity, h_cp, T_cp, h_air, T_air, rho_cpcm,  k_cpcm,  cp_cpcm, L_cpcm, lines_batt_mod, cols_batt_mod, m_dot, cp_liq, times, q_dots, tr_active, tr_cell, tr_q_dot, tr_time, tr_duration, SOCinit, current, cell_capacity, variable_resistance, charge, charge_rate, max_charge_temp, cooling_mission);
	
	// --------------------------- INIT SIMULATION DATA ------------------------- //
	set_mesh_problem(type, N, Nx, Ny, Lx, Ly, dx, dy, t, l, t_fin, t_pcm, lines_batt_mod, cols_batt_mod);
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
		batt_pos[i] = 0;
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
		T_fluid[i] = T_cp;
	}
	
	// -------------------------- INITIALIZATION ----------------------- //
	double SOCi = 0.0;
	double SOC = SOCinit;
	double time = 0.0;				// Current time step [s] 
	double q_dot = 0.0;
	int it = 0;						// Iterations counter
	int k = 0;						// Mesh element index
	int cell_pos = 0;				// Battery cell index
	int previous_time = 0;
	int print_every = 10;			// Print every X seconds
	bool print_now = true;
	bool cooling_active = false;	// Liquid cooling active
	bool charging_active = false;	// Charging has started
	bool tr_started = false;		// Thermal runaway event has started
	int const total_cells = lines_batt_mod * cols_batt_mod;
	std::vector<double> T_ave_cell(total_cells);
	std::vector<double> sum_T_dV(total_cells);
	std::vector<double> sum_dV(total_cells);
	std::vector<double> cell_heat_dissipation(total_cells);

	// Simulation timer
	clock_t start, end;

	// Map mesh
	int o = path(Nx, Ny, pp, ee, ww, nn, ss);
	map_mesh(type, o, Nx, Ny, dx, dy, D, t, l, t_fin, t_pcm, kx_bat, ky_bat, k_pcm, k_alu, k_cpcm, kx_gra, ky_gra, rho_bat, rho_pcm, rho_alu, rho_cpcm, rho_gra, cp_bat, cp_pcm, cp_alu, cp_cpcm, cp_gra, L_pcm, L_cpcm, pp, R, kx, ky, rho, cp, L, lines_batt_mod, cols_batt_mod, batt_pos);

	// Store initial mesh info
	plot_sim(o, Nx, Ny, dx, dy, pp, T, f, R, time, results);
	plot_log(time, o, pp, R, T, f, SOC, results);
	//plot_mesh(o, Nx, Ny, dx, dy, pp, ww, ee, nn, ss, R, kx, ky, rho, cp, L, batt_pos, results);
	//plot_coef(o, Nx, Ny, dx, dy, pp, aw, ae, an, as, ap, b, time, results);
	
	// Resíduos (SOR, liquid fraction and conductivity)
	double resmax = 1.0E-4;
	double resf = 1.0;

	set_mission(SimTotalTime, TotalTime, type, charge, charge_rate, times, current, SOCinit, cell_capacity);

	// Start simulation!
	start = clock();		// Start timer!

	double max_batt_temp = Tinitial;	
	while (time < TotalTime || SOC < SOCinit || max_batt_temp > Tinitial)
	{
		if (print_now)
		{
			printf("\n// -------------------- Time Step = %5.3fs -------------------- // \n", time);
			if (SOC < SOCinit && time > TotalTime)
			{
				printf("SOC is %5.3f %%, less than %5.3f %%\n", SOC, SOCinit);
			}
			if (max_batt_temp > Tinitial && time > TotalTime)
			{
				printf("Max battery temperature is %5.3f degC, higher than %5.3f degC\n", max_batt_temp, Tinitial);
			}
		}

		// ----------- Copy previous time step data  ----------- //
		// SOC
		SOCi = SOC;
		// Temperature
		for (int i = 0; i < o; i++)
		{
			k = pp[i];
			if (k != 0)
			{
				Ti[k] = T[k];
				fi[k] = f[k];
			}
			/*
				if ((type == 4 || type == 5) && (variable_resistance == true) && (R[pp[i]] == 1))
				{
					k = pp[i];
					cell_pos = batt_pos[k] - 1;
					sum_T_dV[cell_pos] += T[k] * (dx * dy * w);
					sum_dV[cell_pos] += dx * dy * w;
				}
			*/
		}
			/*
			if (cooling_active == false && cooling_mission == true && fi[i] > 0 && (type == 3 || type == 4))
			{
				cooling_active = true;
				printf("Liquid cooling activated at %5.3fs\n", time);
			}
			*/
		// ---------- Check liquid cooling condition ----------- //
		if (cooling_active == false && cooling_mission == true && (type == 4 || type == 5))
		{
			cooling_active = true;
			printf("Liquid cooling activated at %5.3fs\n", time);
		}
		else if (cooling_active == false && (time > TotalTime) && (type == 4 || type == 5))
		{
			cooling_active = true;
    		h_cp = 1358;
    		T_cp = 15;    
    		m_dot = 1.07;
    		cp_liq = 3285; // Considerando 15°C
			printf("GSE liquid cooling activated at %5.3fs\n", time);
		}
		/*
		if (variable_resistance == true && (type == 4 || type == 5))
		{
			for (int i = 0; i < total_cells; i++)
			{
				T_ave_cell[i] = sum_T_dV[i] / sum_dV[i];
			}
		} 
		*/
		// ----------------------------------------------------- //

		// ------------- Check coolant temperature ------------- //
		if (cooling_active == true && (type == 4 || type == 5))
		{
			fluid_temp(o, pp, N, Ny, T, batt_pos, T_fluid, T_cp, h_cp, m_dot, t, w, cp_liq, cols_batt_mod);
		}
		// ----------------------------------------------------- //

		// ---------- Update PCM thermal conductivity ---------- //		
		update_thermal_conductivity(R, o, pp, f, kx, ky, ks_pcm, kl_pcm, k_cop, porosity, k_pcm, k_cpcm);
		// ----------------------------------------------------- //
		
		// --------- Init residue and number of iterations ------//
		resf = 1.0;
		it = 0;
		// ----------------------------------------------------- //

		while (resf > resmax)
		{
			// Get Heat generation for module simulation
			if (type == 4 || type == 5)
			{
				get_q_dot(charging_active, Tinitial, SOCinit, time, o, pp, T, R, TotalTime, SimTotalTime, times, q_dots, current, T_ave_cell, SOC, SOCi, cell_capacity, dt, t, l, w, variable_resistance, total_cells, cell_heat_dissipation, charge, charge_rate, max_charge_temp);
				/*
				if (time >= tr_time && time <= tr_time + tr_duration && tr_active == true && tr_started == false)
				{
					tr_started = true;
					printf("Thermal runaway started at %5.3fs\n", time);
				}
				else if((time < tr_time || time > (tr_time + tr_duration)) && tr_started == true)
				{
					tr_started = false;
					printf("Thermal runaway stoped at %5.3fs\n", time);
				}
				*/
			}
			
			// Coefficients matrix setup
			assembly(o, pp, type, Nx, Ny, dx, dy, ww, ee, nn, ss, kx, ky, ap, aw, ae, an, as, su, sp, b, Tw, Te, Tn, Ts, qw, qe, qn, qs, T, Ti, rho, cp, L, f, fi, R, Q, dt, h_cp, T_cp, h_air, T_air, w, cooling_active, T_fluid, q_dot, batt_pos, tr_started, tr_cell, tr_q_dot, tr_time, tr_duration, cell_heat_dissipation);

			// Linear system solver
			SORt(o, pp, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti);
			
			// Compute liquid fraction
			SORf(o, pp, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti, R, f, fi, rho, L, dt, Tmelt, &resf);

			it++;
			if (print_now)
			{
				printf("Iteration = %i\t resf = %5.5f\t SOC = %5.5f\n", it, resf, SOC);
			}
		}

		time = time + dt;
		max_batt_temp = get_max_batt_temp(o, pp, T, R);
		/*
		double max_batt_temp = get_max_batt_temp(o, pp, T, R);
		if (charge == true && ((charging_active == false && cooling_active == true) || (max_batt_temp > Tinitial && SOC >= SOCinit)))
		{
			SimTotalTime += dt;  
		}

		// ------------- Check dispatch conditions ------------- //
		
		if (charge == true && ((charging_active == false && cooling_active == true) || (max_batt_temp > Tinitial && SOC >= SOCinit)))
		{
			SimTotalTime += dt;  
		}
		*/
		// Print time step info
		plot_log(time, o, pp, R, T, f, SOC, results);
		//plot_qdot(time, SOC, it, cell_heat_dissipation, T_ave_cell, results);
		/*		
 		if ((previous_time != int(time)) && (int(time) % print_every == 0))
		{
			plot_sim(o, Nx, Ny, dx, dy, pp, T, f, R, time, results);
		} 
		*/
		if ((previous_time != int(time)) && (int(time) % print_every == 0))
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