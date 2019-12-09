This code is the companion to "Determining risk factors for triple whammy acute kidney injury: sex-specific modeling and analysis" by Jessica Leete, Francisco J. Lopez-Hernadez, and Anita T. Layton.

This code is the companion to "Sex-Specific Computational Models for Blood Pressure Regulation in the Rat" by Sameed Ahmed and Anita T. Layton.

This file contains information about how to run the code.

To create all files and figures from the first paper, run Leete2020.m
To run your own simulations, call solve_ss_numerical or run_sim_treatments with the appropraite inputs.

Leete2020.m
Running this will create all data files and figures for the above mentioned paper.
It reads from files in Human_IG_Data/. It writes to folders called Human_Data/ and Human_Figures/.
The files in Human_IG_Data/ are sufficient initial guesses for the solvers.

solve_ss_numerical.m
This solves the steady state system of equations. Saves output to Human_Data/
Required input parameters:
human - 1 to run a human simulation, 0 for the rat simulation
gg - 1 for male, 2 for female
IG - string containing the filename of the initial guess for the solver
Optional input parameters: (each should be a string followed by a value)
'ACEi', # - 0 for no treatment, 1 for normal dose, 2 for high dose. Default is 0.
'furosemide', # - 0 for no treatment, 1 for normal dose, 2 for high dose. Default is 0.
'NSAID', # - 0 for no treatment 1 for normal dose, 2 for high dose. Default is 0.
'Impaired Myogenic Response', # - 0 for normal, 1 for impaired. Default is 0.
'Low Water Intake', # - 0 for normal, 1 for low water intake. Default is 0.
'RSNA', # - value by which N_rsna is multiplied by to induce hypertension. In my simulations it is 1 for normotensive simulations, 2.5 for hypertensive. Default is 1.


run_sim_treatments.m
This solves the time course simulations. It saves output to Human_Data/.
Required input parameters:
human - 1 to run a human simulation, 0 for the rat simulation
gg - 1 for male, 2 for female
IC - string containing the filename of the initial condition for the solver
Optional input parameters: (each should be a string followed by a value)
'ACEi', # - value between 0 and 1 for the variable kappa_ACEi. Default is 0.
'furosemide', # - array of length 2 containing the values between 0 and 1 for the variables kappa_f and kappa_f_md.
'NSAID', # - 0 for no treatment 1 for normal dose, 2 for high dose. Default is [0 0].
'Impaired Myogenic Response', # - 0 for normal, 1 for impaired. Default is 0.
'Low Water Intake', # - 0 for normal, 1 for low water intake. Default is 0.
'RSNA', # - value by which N_rsna is multiplied by to induce hypertension. In my simulations it is 1 for normotensive simulations, 2.5 for hypertensive. Default is 1.

If you want to change the length of simulations, change the variable num_days on line 80.
If any simulation fails at a time value other than 0 or 100 (tchange), try changing the time vector tspan on line 84.


get_pars.m
function that returns parameter vector required by bp_reg_mod.m
Required input parameters:
species - string containing species to simulate. 'human' or 'rat'
sex - string containing sex of individual to simulate. 'male' or 'female'
Optional input parameters:
'Normal' - does not change the simulation
'm_Reab', boolean - if true, runs a simulation with male sodium and water reabsorption
'm_RAS', boolean - if true, runs a simulation with male renin angiotensin system
'm_RSNA', boolean - if true, runs a simulation with male renal sympathetic nervous activity
'm_RSNA_m_Reab', boolean - if true, runs a simulation with male sodium and water reabsorption and male renal sympathetic nervous activity
'm_RAS_m_Reab', boolean - if true, runs a simulation with male renin angiotensin system and male sodium and water reabsorption
'RSNA', # - value by which N_rsna is multiplied by to induce hypertension. In my simulations it is 1 for normotensive simulations, 2.5 for hypertensive. Default is 1.


bp_reg_mod.m
Model equations file
Required input parameters:
t - time value
x - variable array. Length 84 for human, 93 for rat simulations
x_p - variable derivative array. Length 84 for human, 93 for rat simulations
pars - parameter array from get_pars.m
tchange - time value at which to start applying drug treatments. A value of 100 was used.
Optional input parameters:
'AngII', # - value of Ang II infusion
'ACEi', # - value between 0 and 1 for ACE inhibition
'ARB', # - value between 0 and 1 for angiotensin receptor blocker treatment
'furosemide', [# #] - array of length 2 containing the values between 0 and 1 for the variables kappa_f and kappa_f_md.
'NSAID', # - value indicating NSAID treatement level. 0 for no treatment, 1 for normal dose, 2 for high dose.
'RPP', # - value to set renal perfusion pressure to.
'Denerve', boolean - if true, run simulation with unilateral renal denervation by fixing rsna = 1, which is baseline
'No TGF', boolean - if true, run simulation with no tubuloglomerular feedback.
'No Myo', boolean - if true, run simulation with no the myogenic response.
'Linear Myo', boolean - if true, run simulation with a linear myogenic response.
'Impaired Myogenic Response', boolean - if true, run simulation with an impaired myogenic response.
'Fixed Water Intake', boolean - if true, run simulation with fixed water intake
'Low Water Intake', boolean - if true, run simulation with low water intake
'm_RSNA', boolean - if true, runs a simulation with male renal sympathetic nervous activity
'm_RSNA_m_Reab', boolean - if true, runs a simulation with male sodium and water reabsorption and male renal sympathetic nervous activity
'm_AT2R', boolean - if true, run simulation with male angiotensin type 2 receptor


plot_scatter.m
Creates plots from Figure 4 in Leete 2020 paper.
Saves matlab file to Human_Figures/.
Requires hline from https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline


plot_riskfactors.m
Creates plots from Figure 5 in Leete 2020 paper.
Saves matlab file to Human_Figures/.
Saves .png file to Human_Figures/.
Requires hline from https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline


To create all files and figures from the second paper, run Ahmed2020.m
To run your own simulations, call run the individual scripts listed in Ahmed2020.m with your inputs.

Ahmed2020.m
Running this will create all data files and figures for the above mentioned paper.
It reads from and writes to folders called Rat_Data/ and Rat_Figures/.

run_sim.m
This runs the time course simulation and saves the figures to Rat_Figures/.

run_sim_AngII.m
This runs the time course simulation for Ang II infusion as in Ref. Sampson 2008.
It saves the figures to Rat_Figures/.

run_sim_RPP.m
This runs the time course simulation for manipulating renal perfusion pressure as in Ref. Hilliard 2011.
The rat is unilaterally denervated and anesthetized and body fluids are maintained.
It saves the figures to Rat_Figures/.

solve_ss_Phisodin.m
This solves the steady state solution for varing sodium intake.
It saves the figures to Rat_Figures/.

solve_ss_scenario.m
This solves the steady state solution for different scenarios.
It saves the data to Rat_Data/.

solve_ss_sen_anal.m
This solves the steady state solution to conduct a sensitivity analysis for varing certain parameters.
It saves the figures to Rat_Figures/.

