# Dynamic_network_poisson
Asai and So (2024), "Dynamic Network Poisson Autoregression with Application to
COVID-19 Count Data", Journal of Data Science
Matlab codes

The data in the analysis is provided by:
https://www.sciencedirect.com/science/article/pii/S1201971220324541?via%3Dihub

confirmed1.mat: mY0 (T times 51 matrix, confirmed cases of covid-19 by state)

linkavar1.mat: mD0 (T times 51*51 matrix, dynamic network among states)

%%%%%%%% Section 3

Figure 1: main_check_lyap.m

Monte Carlo experiments

DN-PAR model: main_mc_tdnp_p1q1_dram.m

CN-PAR model: main_mc_cnp_p1q1.m


%%%%%%%% Section 4

Figure 2: main_conf_fig.m

Figure 3: main_conf_net_heatmap.m

Estimation via DRAM

CN-PAR model: main_conf_cnp_p1q1_dram.m

EN-PAR model: main_conf_ednp_p1q1_dram.m

DN-PAR model: main_conf_tdnp_p1q1_dram.m

CN-PAR-I model: main_conf_cnpi_p1q1_dram.m

EN-PAR-I model: main_conf_ednpi_p1q1_dram.m

DN-PAR-I model: main_conf_tdnpi_p1q1_dram.m

Forecasts

CN-PAR model: main_conf_cnp_p1q1_dram_for.m

EN-PAR model: main_conf_ednp_p1q1_dram_for.m

DN-PAR model: main_conf_tdnp_p1q1_dram_for.m

CN-PAR-I model: main_conf_cnpi_p1q1_dram_for.m

EN-PAR-I model: main_conf_ednp_p1q1_dram_for.m

TN-PAR-I model: main_conf_tdnp_p1q1_dram_for.m

Comparison of Forecasts: main_conf_for_com.m

Figure 4: main_conf_tdnp_p1q1_dram_for_fig.m


