#include "..\Headers\Include.h"

void get_q_dot(bool& charging_active, double Tinitial, double SOCinit, double time, int o, int* pp, double* T, int* R, double TotalTime, double SimTotalTime, vector<double> times, vector<double> q_dots, vector<double> current, vector<double> T_ave_cell, double& SOC, double& SOCi, double cell_capacity, double dt, double t, double l, double w, bool variable_resistance, int total_cells, vector<double>& cell_heat_dissipation, bool charge, double charge_rate, double max_charge_temp)
{
	double max_batt_temp = 0.0;
	int p;

	if (charge == true)
	{
		if (time >= TotalTime)
		{
			double max_batt_temp = get_max_batt_temp(o, pp, T, R);
			if ((max_batt_temp > max_charge_temp) || (max_batt_temp > Tinitial && SOC >= SOCinit))
			{
				for (int j = 0; j < total_cells; j++)
				{
					cell_heat_dissipation[j] = 0.0;
				}
			}
			else
			{
				if (charging_active == false)
				{
					charging_active = true;			
				}
				double i = charge_rate * (cell_capacity / 1000);                       // Current in A
				
				// Calculation of cell State of Charge 
				SOC = (SOCi/100 + ((i * (dt / 3600)) / (cell_capacity/1000))) * 100; // SOC in percentage
				
				if (variable_resistance == true)
				{
					// Dados da tabela (SOC e Rref,dis)
					vector<double> soc_tab = {20, 30, 40, 50, 60, 70, 80, 90};
					vector<double> rref_dis_tab = {12.17, 12.11, 12.09, 11.94, 12.00, 12.50, 12.65, 12.87};					

					for (int j = 0; j < total_cells; j++)
					{
						double rref = lookup_linear_clipped(SOC, soc_tab, rref_dis_tab);
						cell_heat_dissipation[j] = i * i * ((rref)/1000) / (t*l*w);
					}
				}
				/*
				if (variable_resistance == true)
				{
					// Dados da tabela (SOC e Rref,dis)
					vector<double> soc_tab = {10, 20, 30, 40, 50, 60, 70, 80};
					vector<double> rref_charge_tab = {15.72, 14.45, 14.50, 14.42, 13.95, 14.00, 14.58, 15.05};

					//vector<double> temp_tab = {26, 30, 35, 40};
					//vector<double> Tfactor_tab = {1.0, 0.881, 0.732, 0.583};
					
					for (int j = 0; j < total_cells; j++)
					{
						double rref = lookup_linear_clipped(SOC, soc_tab, rref_charge_tab);
						cell_heat_dissipation[j] = i * i * ((rref)/1000) / (t*l*w);
					}			
				}
				*/
				else
				{
					double rref = 4; // Constant resistance in mOhm
					for (int j = 0; j < total_cells; j++)
					{
						cell_heat_dissipation[j] = i * i * ((rref)/1000) / (t*l*w);
					}
				}
			}
		} 
		else
		{
			double i = select_data_in_time(time, times, current);

			// Calculation of cell State of Charge 
			SOC = (SOCi/100 - ((i * (dt / 3600)) / (cell_capacity/1000))) * 100; // SOC in percentage

			if (variable_resistance == true)
			{
				// Dados da tabela (SOC e Rref,dis)
				vector<double> soc_tab = {20, 30, 40, 50, 60, 70, 80, 90};
				vector<double> rref_dis_tab = {12.17, 12.11, 12.09, 11.94, 12.00, 12.50, 12.65, 12.87};
					
				//vector<double> temp_tab = {26, 30, 35, 40};
				//vector<double> Tfactor_tab = {1.0, 0.881, 0.732, 0.583};

				for (int j = 0; j < total_cells; j++)
				{
					double rref = lookup_linear_clipped(SOC, soc_tab, rref_dis_tab);
					cell_heat_dissipation[j] = i * i * ((rref)/1000) / (t*l*w);
				}
			}
			else
			{
				double rref = 4; 
				for (int j = 0; j < total_cells; j++)
				{
					cell_heat_dissipation[j] = i * i * ((rref)/1000) / (t*l*w);
				}
			}		
		}
	}
	else
	{
		for (int j = 0; j < total_cells; j++)
		{
			cell_heat_dissipation[j] = select_data_in_time(time, times, q_dots);
		}
	}
}

double select_data_in_time(double time, vector<double> times, vector<double> data)
{
	double selected_data = 0.0;
	int size_times = static_cast<int>(times.size());

	if (time > times[size_times - 1])
	{
		selected_data = data[size_times];
	}
	else
	{
		for (int i = 0; i < size_times; i++)
		{
			if (time <= times[i])
			{
				selected_data = data[i];
				break;
			}
		}
	}

	return selected_data;
}

double lookup_linear_clipped(double breakpoint, vector<double>& indexes, vector<double>& table_data)
{
	if (indexes.size() != table_data.size() || table_data.empty()) {
		throw invalid_argument("Tabelas inválidas");
	}

	// clipping nas bordas
	if (breakpoint <= indexes.front()) {
		return table_data.front();
	}
	if (breakpoint >= indexes.back()) {
		return table_data.back();
	}

	// encontrar intervalo [i, i+1] tal que indexes[i] <= breakpoint <= indexes[i+1]
	auto it = lower_bound(indexes.begin(), indexes.end(), breakpoint);
	size_t i1 = static_cast<size_t>(it - indexes.begin());
	size_t i0 = i1 - 1;

	double x0 = indexes[i0];
	double x1 = indexes[i1];
	double y0 = table_data[i0];
	double y1 = table_data[i1];

	// interpola??o linear
	double t = (breakpoint - x0) / (x1 - x0);
	return y0 + t * (y1 - y0);
}

void set_mission(double& SimTotalTime, double& TotalTime, int type, bool charge, double charge_rate, vector<double>& times, vector<double>& current, double& SOCinit, double cell_capacity)
{
    if (charge == true && (type == 4 || type == 5))
    {
        double soc_end = 0.0; 
        double charge_current = 0.0;
        double charge_time = 0.0;

        // Discharge
        for (int i = 0; i <= times.size(); i++)
        {
            if (i == 0)
            {
                soc_end = SOCinit - ((current[i] * (times[i] / 3600)) / (cell_capacity / 1000)) * 100;
	            cout << "Phase " << i+1 << ": Current = " << current[i] << " A, Time = " << times[i] << " s, SOC end = " << fixed << setprecision(2) << soc_end << " %" << endl;

			}
            else if (i < times.size())
			{
				soc_end = soc_end - ((current[i] * (times[i] - times[i-1]) / 3600) / (cell_capacity / 1000)) * 100;
	            cout << "Phase " << i+1 << ": Current = " << current[i] << " A, Time = " << times[i] << " s, SOC end = " << fixed << setprecision(2) << soc_end << " %" << endl;
		}
			else
            {
                soc_end = soc_end - ((current[i-1] * (TotalTime - times[i-1]) / 3600) / (cell_capacity / 1000)) * 100;
                cout << "Phase " << i+1 << ": Current = " << current[i-1] << " A, Time = " << TotalTime << " s, SOC end = " << fixed << setprecision(2) << soc_end << " %" << endl;
            }
        }       
    
        // Charge
        charge_current = charge_rate * (cell_capacity / 1000); 
        charge_time = ((SOCinit - soc_end) / 100) * (cell_capacity / 1000) * 3600 / charge_current;

        cout << "TAT: Current = " << fixed << setprecision(2) << charge_current << " A, Duration = " << (charge_time / 60) << " min, SOC end = " << fixed << setprecision(2) << SOCinit << " %" << endl;

        SimTotalTime = TotalTime + charge_time;
    }
    else
    {
        SimTotalTime = TotalTime;
    }
}

double get_max_batt_temp(int o, int* pp, double* T, int* R)
{
	int k;
	double max_batt_temp = 0.0;
	
	for (int i = 0; i < o; i++)
	{
		k = pp[i];
		if (T[k] > max_batt_temp && R[k] == 1)
		{
			max_batt_temp = T[k];
		}	
	}
	return max_batt_temp;
}