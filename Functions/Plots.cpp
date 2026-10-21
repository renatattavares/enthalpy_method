#include "..\Headers\Include.h"

void plot_mesh(int o, int Nx, int Ny, double time, double dx, double dy, int* pp, int* ww, int* ee, int* nn, int* ss, int* R, double* kx, double* ky, double* rho, double* cp, double* L, int* batt_pos, string dir) {

	double rx{}, ry{};
	struct stat sb;
	string iteration = to_string(int(time));
	string filename;
	string simulation_folder = dir + "\\" + "Mesh";

	// Create Results Folder if it does not exist
	if (stat(simulation_folder.c_str(), &sb) != 0)
	{
		if (_mkdir(simulation_folder.c_str()) != 0)
		{
			cout << "Failed to create folder " << simulation_folder << endl;
		}
	}
	
	// Create filename
	filename = simulation_folder + "\\" + "MeshPlots" + "_" + iteration + ".vtk";

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting data files..." << endl;
	
	// Build
	ofstream fout;
	fout.open(filename);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "Physical Properties plot" << endl;
	fout << "ASCII" << endl;
	fout << endl;
	fout << "DATASET RECTILINEAR_GRID" << endl;
	fout << "DIMENSIONS " << Ny << " " << Nx << " 1" << endl;
	fout << endl;
	fout << "Y_COORDINATES " << Nx << " float" << endl;


	for (int i = 0; i < Nx; i++)
	{
		rx = i * dx + 0.5 * dx;
		fout << rx << " ";
	}

	fout << endl;
	fout << endl;
	fout << "X_COORDINATES " << Ny << " float" << endl;

	for (int j = 0; j < Ny; j++)
	{
		ry = j * dy + 0.5 * dy;
		fout << ry << " ";
	}

	fout << endl;
	fout << endl;
	fout << "Z_COORDINATES 1 float" << endl;
	fout << "0" << endl;
	fout << endl;
	fout << "POINT_DATA " << Nx * Ny << endl;

	// Write PPs
	fout << "SCALARS pp int" << endl;
	fout << "LOOKUP_TABLE pp" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << pp[a] << endl;
	}

	// Write PPs
	fout << "SCALARS ww int" << endl;
	fout << "LOOKUP_TABLE ww" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << ww[pp[a]] << endl;
	}

	// Write PPs
	fout << "SCALARS ee int" << endl;
	fout << "LOOKUP_TABLE ee" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << ee[pp[a]] << endl;
	}

	// Write PPs
	fout << "SCALARS nn int" << endl;
	fout << "LOOKUP_TABLE nn" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << nn[pp[a]] << endl;
	}

	// Write PPs
	fout << "SCALARS ss int" << endl;
	fout << "LOOKUP_TABLE ss" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << ss[pp[a]] << endl;
	}

	// Write Region Information
	fout << "SCALARS Region int" << endl;
	fout << "LOOKUP_TABLE Region" << endl;

	for (int a = 0; a < o; a++) 
	{
		fout << R[pp[a]] << endl;
	}

	// Write Conductivity Information
	fout << "SCALARS Conductivity_X float" << endl;
	fout << "LOOKUP_TABLE Conductivity_X" << endl;

	for (int a = 0; a < o; a++) 
	{
		fout << kx[pp[a]] << endl;
	}

	// Write Conductivity Information
	fout << "SCALARS Conductivity_Y float" << endl;
	fout << "LOOKUP_TABLE Conductivity_Y" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << ky[pp[a]] << endl;
	}

	// Write Density Information
	fout << "SCALARS Density float" << endl;
	fout << "LOOKUP_TABLE Density" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << rho[pp[a]] << endl;
	}

	// Write Specific Heat Information
	fout << "SCALARS Specific_Heat float" << endl;
	fout << "LOOKUP_TABLE Specific_Heat" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << cp[pp[a]] << endl;
	}

	// Write Latent Heat Information
	fout << "SCALARS Latent_Heat float" << endl;
	fout << "LOOKUP_TABLE Latent_Heat" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << L[pp[a]] << endl;
	}

	// Write Cell Batt Position Information
	fout << "SCALARS Batt_Pos int" << endl;
	fout << "LOOKUP_TABLE Batt_Pos" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << batt_pos[pp[a]] << endl;
	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

}

void plot_coef(int o, int Nx, int Ny, double dx, double dy, int* pp, double* aw, double* ae, double* an, double* as, double* ap, double* b, double time, string results_folder) {

	double rx{}, ry{};
	struct stat sb;
	string iteration = to_string(int(time));
	string title = "Coefficients �C plot";
	string info = "SCALARS Coefficients float";
	string filename;
	string property = "Coef";
	string property_folder = results_folder + "\\" + property;

	// Create Results Folder if it does not exist
	if (stat(property_folder.c_str(), &sb) != 0)
	{
		if (_mkdir(property_folder.c_str()) != 0)
		{
			cout << "Failed to create folder " << property_folder << endl;
		}
	}

	// Create filename
	filename = property_folder + "\\" + property + "_" + iteration + ".vtk";

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting coefficients' files..." << endl;

	// Build
	ofstream fout;
	fout.open(filename);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "Physical Properties plot" << endl;
	fout << "ASCII" << endl;
	fout << endl;
	fout << "DATASET RECTILINEAR_GRID" << endl;
	fout << "DIMENSIONS " << Ny << " " << Nx << " 1" << endl;
	fout << endl;
	fout << "Y_COORDINATES " << Nx << " float" << endl;


	for (int i = 0; i < Nx; i++)
	{
		rx = i * dx + 0.5 * dx;
		fout << rx << " ";
	}

	fout << endl;
	fout << endl;
	fout << "X_COORDINATES " << Ny << " float" << endl;

	for (int j = 0; j < Ny; j++)
	{
		ry = j * dy + 0.5 * dy;
		fout << ry << " ";
	}

	fout << endl;
	fout << endl;
	fout << "Z_COORDINATES 1 float" << endl;
	fout << "0" << endl;
	fout << endl;
	fout << "POINT_DATA " << Nx * Ny << endl;

	fout << "SCALARS pp int" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << pp[a] << endl;
	}

	fout << "SCALARS aw float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << aw[pp[a]] << endl;
	}

	fout << "SCALARS ae float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << ae[pp[a]] << endl;
	}

	fout << "SCALARS an float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << an[pp[a]] << endl;
	}

	fout << "SCALARS as float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << as[pp[a]] << endl;
	}

	fout << "SCALARS ap float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << ap[pp[a]] << endl;
	}

	fout << "SCALARS b float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << b[pp[a]] << endl;
	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

}

void plot_sim(int o, int Nx, int Ny, double dx, double dy, int* pp, double* T, double* f, double time, string results_folder) {

	double rx{}, ry{};
	struct stat sb;
	string iteration = to_string(int(time));
	string filename;
	string simulation_folder = results_folder + "\\" + "Simulation";
	
	// Create Results Folder if it does not exist
	if (stat(simulation_folder.c_str(), &sb) != 0)
	{
		if (_mkdir(simulation_folder.c_str()) != 0)
		{
			cout << "Failed to create folder " << simulation_folder << endl;
		}
	}

	// Create filename
	filename = simulation_folder + "\\" + "SimPlots" + "_" + iteration + ".vtk";

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting simulation plots..." << endl;

	// Build VTK
	ofstream fout;
	fout.open(filename);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "Simulation plots" << endl;
	fout << "ASCII" << endl;
	fout << endl;
	fout << "DATASET RECTILINEAR_GRID" << endl;
	fout << "DIMENSIONS " << Ny << " " << Nx << " 1" << endl;
	fout << endl;
	fout << "Y_COORDINATES " << Nx << " float" << endl;

	for (int i = 0; i < Nx; i++)
	{
		rx = i * dx + 0.5 * dx;
		fout << rx << " ";
	}

	fout << endl;
	fout << endl;
	fout << "X_COORDINATES " << Ny << " float" << endl;

	for (int j = 0; j < Ny; j++)
	{
		ry = j * dy + 0.5 * dy;
		fout << ry << " ";
	}

	fout << endl;
	fout << endl;
	fout << "Z_COORDINATES 1 float" << endl;
	fout << "0" << endl;
	fout << endl;
	fout << "POINT_DATA " << Nx * Ny << endl;

	// Write Property Information
	fout << "SCALARS Temperature float" << endl;
	fout << "LOOKUP_TABLE Temperature" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << T[pp[a]] << endl;
	}

	// Write Property Information
	fout << "SCALARS Liquid_Fraction float" << endl;
	fout << "LOOKUP_TABLE Liquid_Fraction" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << f[pp[a]] << endl;
	}

	// Exit
	fout.close();
	cout << "\nDone." << endl << endl;
	
	return;
}

string create_results_folder() {

	// Get the current time
	time_t t = time(nullptr);
	tm* now = localtime(&t);

	// Create a string stream to format the date and time
	ostringstream oss;
	oss << (now->tm_year + 1900) << '_'
		<< setw(2) << setfill('0') << (now->tm_mon + 1) << '_'
		<< setw(2) << setfill('0') << now->tm_mday << '_'
		<< setw(2) << setfill('0') << now->tm_hour << '_'
		<< setw(2) << setfill('0') << now->tm_min << '_'
		<< setw(2) << setfill('0') << now->tm_sec;

	string dateTimeStr = oss.str();

	// Get current folder
	char cwd[1024];
	getcwd(cwd, 1024);
	string current_folder(cwd);
	
	// Construct results folder
	current_folder = current_folder + "\\Results\\Results_";
	string results_folder = current_folder + dateTimeStr;

	if (_mkdir(results_folder.c_str()) != 0) 
	{
		cout << "Failed creating results folder..." << endl;
		return results_folder;
	}
	return results_folder;
}

void log(double time, int o, int* pp, int* R, double* T, double* f, string results_folder) {

	double delta_temp = 0.0;
	double max_temp = 0.0;
	double min_temp = 1000.0;
	double total_f = 0.0;
	double tmax_pcm = 0.0;
	double tmed_pcm = 0.0;
	int can_melt = 0;
	int l = 0;

	// Log file
	ofstream fout;
	string filename;
	filename = results_folder + "/Log.dat";

	for (int i = 0; i < o; i++)
	{
		l = pp[i];
		if (T[l] > max_temp && R[l] == 1)
		{
			max_temp = T[l];
		}
		if (T[l] < min_temp && R[l] == 1)
		{
			min_temp = T[l];
		}
		if (R[l] == 0 || R[l] == 5)
		{
			total_f = total_f + f[l];
			tmed_pcm = tmed_pcm + T[l];
			can_melt++;
			if (T[l] > tmax_pcm)
			{
				tmax_pcm = T[l];
			}
		}
	}

	delta_temp = max_temp - min_temp;
	total_f = total_f / can_melt;
	tmed_pcm = tmed_pcm / can_melt;

	fout.open(filename, std::ios_base::app);
	fout << time << '\t' << max_temp << '\t' << min_temp << '\t' << delta_temp << '\t' << total_f << '\t' << tmax_pcm << '\t' << tmed_pcm << endl;
	fout.close();

	//cout << "Max temperature:" << max_temp << endl;
	//cout << "Min temperature:" << min_temp << endl;
	//cout << "Delta temperature:" << delta_temp << endl;
	//cout << "Liquid Fraction:" << total_f << endl;

}

void plot_res(int o, int Nx, int Ny, double dx, double dy, int* pp, double* P, int it, double time, string results_folder) {

	double rx{}, ry{};

	string sim_time = to_string(0);//to_string(time);
	string iteration = to_string(it);
	string title = "Residue �C plot";
	string info = "SCALARS Residue float";
	string filename;

	// Create Results Folder
	string property_folder = results_folder + "\\" + "Res";

	if (_mkdir(property_folder.c_str()) != 0)
	{
		cout << "Failed to create folder " << property_folder << endl;
	}

	// Create filename
	filename = property_folder + "\\Res_" + sim_time + "_" + iteration + ".vtk";

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting file..." << endl;

	// Build VTK
	ofstream fout;
	fout.open(filename);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << title << endl;
	fout << "ASCII" << endl;
	fout << endl;
	fout << "DATASET RECTILINEAR_GRID" << endl;
	fout << "DIMENSIONS " << Ny << " " << Nx << " 1" << endl;
	fout << endl;
	fout << "Y_COORDINATES " << Nx << " float" << endl;

	for (int i = 0; i < Nx; i++)
	{
		rx = i * dx + 0.5 * dx;
		fout << rx << " ";
	}

	fout << endl;
	fout << endl;
	fout << "X_COORDINATES " << Ny << " float" << endl;

	for (int j = 0; j < Ny; j++)
	{
		ry = j * dy + 0.5 * dy;
		fout << ry << " ";
	}

	fout << endl;
	fout << endl;
	fout << "Z_COORDINATES 1 float" << endl;
	fout << "0" << endl;
	fout << endl;
	fout << "POINT_DATA " << Nx * Ny << endl;

	// Write Property Information
	fout << info << endl;
	fout << "LOOKUP_TABLE Residue" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << P[pp[a]] << endl;
	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

	return;
}

void log_plot(double time, int o, int* pp, int* R, double* T, double* f) {

	double delta_temp = 0.0;
	double max_temp = 0.0;
	double min_temp = 1000.0;
	double total_f = 0.0;
	int can_melt = 0;
	int l = 0;

	for (int i = 0; i < o; i++)
	{
		l = pp[i];
		if (T[l] > max_temp && R[l] == 1)
		{
			max_temp = T[l];
		}
		if (R[l] == 0 || R[l] == 5)
		{
			total_f = total_f + f[l];
			can_melt++;
		}
	}

	total_f = total_f / can_melt;

	cout << "Max temperature:" << max_temp << endl;
	cout << "Liquid Fraction:" << total_f << endl;

}

void log_res(double time, int o, int* pp, int* R, double* T, double* f, string results_folder) {

	double delta_temp = 0.0;
	double max_temp = 0.0;
	double min_temp = 1000.0;
	double total_f = 0.0;
	int can_melt = 0;
	int l = 0;

	// Log file
	ofstream fout;
	string filename;
	filename = results_folder + "/Log_Residue.dat";

	for (int i = 0; i < o; i++)
	{
		l = pp[i];
		if (T[l] > max_temp && R[l] == 1)
		{
			max_temp = T[l];
		}
		if (T[l] < min_temp && R[l] == 1)
		{
			min_temp = T[l];
		}
		if (R[l] == 0 || R[l] == 5)
		{
			total_f = total_f + f[l];
			can_melt++;
		}
	}

	delta_temp = max_temp - min_temp;
	total_f = total_f / can_melt;

	fout.open(filename, std::ios_base::app);
	fout << time << '\t' << max_temp << '\t' << min_temp << '\t' << delta_temp << '\t' << total_f << endl;
	fout.close();

	//cout << "Max temperature:" << max_temp << endl;
	//cout << "Min temperature:" << min_temp << endl;
	//cout << "Delta temperature:" << delta_temp << endl;
	//cout << "Liquid Fraction:" << total_f << endl;

}