#include "..\Headers\Include.h"

void read_data(string filename,
    int& type, int& pouch,
    double& TotalTime, double& dt,
    double& dx, double& dy,
    double& D, double& h,
    double& t, double& l, double& w, double& t_fin, double& t_pcm,
    double& Tinitial, double& Q,
    double& Tn, double& Ts, double& Te, double& Tw,
    double& qn, double& qs, double& qe, double& qw,
    double& rho_bat, double& kx_bat, double& ky_bat, double& cp_bat,
    double& rho_pcm, double& k_pcm, double& ks_pcm, double& kl_pcm,double& cp_pcm, double& L_pcm, double& Tmelt,
    double& rho_alu, double& k_alu, double& cp_alu,
    double& rho_gra, double& kx_gra, double& ky_gra, double& cp_gra,
    double& rho_cop, double& k_cop, double& cp_cop,
    double& porosity,
    double& h_cp, double& T_cp,
    double& h_air, double& T_air,
    double& rho_cpcm, double& k_cpcm, double& cp_cpcm, double& L_cpcm,
    int& lines_batt_mod, int& cols_batt_mod,
    double& m_dot, double& cp_liq,
    vector<double>& times, vector<double>& q_dots,
    bool& tr_active, int& tr_cell, double& tr_q_dot, double& tr_time, double& tr_duration,
    double& SOCinit, vector<double>& current, double& cell_capacity, bool& variable_resistance,
    bool& charge, double& charge_rate, double& max_charge_temp, bool& cooling_mission)
{
    YAML::Node config = YAML::LoadFile(filename);

    type = config["problem"]["type"].as<int>();
    pouch = config["problem"]["pouch"].as<int>();

    TotalTime = config["time"]["TotalTime"].as<double>();
    dt = config["time"]["dt"].as<double>();

    dx = config["mesh"]["dx"].as<double>();
    dy = config["mesh"]["dy"].as<double>();

    auto bg = config["battery_geometry"];
    D = bg["D"].as<double>();
    h = bg["h"].as<double>();
    t = bg["t"].as<double>();
    l = bg["l"].as<double>();
    w = bg["w"].as<double>();
    t_fin = bg["t_fin"].as<double>();
    t_pcm = bg["t_pcm"].as<double>();

    Tinitial = config["initial_conditions"]["Tinitial"].as<double>();
    Q = config["initial_conditions"]["Q"].as<double>();

    auto bt = config["boundary_temperatures"];
    Tn = bt["Tn"].as<double>();
    Ts = bt["Ts"].as<double>();
    Te = bt["Te"].as<double>();
    Tw = bt["Tw"].as<double>();

    auto bf = config["boundary_fluxes"];
    qn = bf["qn"].as<double>();
    qs = bf["qs"].as<double>();
    qe = bf["qe"].as<double>();
    qw = bf["qw"].as<double>();

    auto m = config["materials"];

    auto mb = m["battery"];
    rho_bat = mb["rho"].as<double>();
    kx_bat = mb["kx"].as<double>();
    ky_bat = mb["ky"].as<double>();
    cp_bat = mb["cp"].as<double>();

    auto mp = m["pcm"];
    rho_pcm = mp["rho"].as<double>();
    ks_pcm = mp["ks"].as<double>();
    kl_pcm = mp["kl"].as<double>();
    k_pcm = ks_pcm;                     // Initial conductivity is solid phase
    cp_pcm = mp["cp"].as<double>();
    L_pcm = mp["L"].as<double>();
    Tmelt = mp["Tmelt"].as<double>();

    auto ma = m["aluminum"];
    rho_alu = ma["rho"].as<double>();
    k_alu = ma["k"].as<double>();
    cp_alu = ma["cp"].as<double>();

    auto mg = m["graphite"];
    rho_gra = mg["rho"].as<double>();
    kx_gra = mg["kx"].as<double>();
    ky_gra = mg["ky"].as<double>();
    cp_gra = mg["cp"].as<double>();

    auto mc = m["copper"];
    rho_cop = mc["rho"].as<double>();
    k_cop = mc["k"].as<double>();
    cp_cop = mc["cp"].as<double>();

    porosity = config["cpcm"]["porosity"].as<double>();

    auto conv = config["convection"];
    h_cp = conv["cold_plate"]["h_cp"].as<double>();
    T_cp = conv["cold_plate"]["T_cp"].as<double>();
    m_dot = conv["cold_plate"]["m_dot"].as<double>();
    cp_liq = conv["cold_plate"]["cp_liq"].as<double>();
    h_air = conv["natural"]["h_air"].as<double>();
    T_air = conv["natural"]["T_air"].as<double>();

    auto mod = config["module"];
    lines_batt_mod = mod["lines"].as<int>();
    cols_batt_mod = mod["columns"].as<int>();

    auto mission = config["mission"];
    variable_resistance = mission["variable_resistance"].as<bool>();
    max_charge_temp = mission["max_charge_temp"].as<double>();
    charge = mission["charge"].as<bool>();
    current = mission["current"].as<vector<double>>();
    SOCinit = mission["SOCinit"].as<double>();
    times = mission["times"].as<vector<double>>();
    cell_capacity = mission["cell_capacity"].as<double>();
    charge_rate = mission["charge_rate"].as<double>();
    q_dots = mission["q_dot"].as<vector<double>>();
    cooling_mission = mission["cooling"].as<bool>();
    int size_times = static_cast<int>(times.size());
    int size_qdots = static_cast<int>(q_dots.size());
    assert(size_times + 1 == size_qdots);

    auto tr = config["thermal_runaway"];
    tr_cell = tr["cell"].as<int>();
    tr_active = tr["active"].as<bool>();
    tr_q_dot = tr["q_dot"].as<double>();
    tr_time = tr["time"].as<double>();
    tr_duration = tr["duration"].as<double>();


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
}