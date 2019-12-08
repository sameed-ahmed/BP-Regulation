% This script runs any of the simulations presented in Ahmed and Layton
% 2020 and saves the resulting figures in the folder "Rat_Figures". 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncomment the desired simulation.

% experiment = 'pressure-natriuresis';
% experiment = 'pressure-natriuresis whole range';
% experiment = 'sodium intake';
% experiment = 'Ang II infusion';
experiment = 'sensitivity analysis';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if     strcmp(experiment, 'pressure-natriuresis')
    run_sim_RPP
elseif strcmp(experiment, 'pressure-natriuresis whole range')
    run_sim_RPP_whole
elseif strcmp(experiment, 'sodium intake')
    solve_ss_Phisodin
elseif strcmp(experiment, 'Ang II infusion')
    run_sim_AngII
elseif strcmp(experiment, 'sensitivity analysis')
    solve_ss_sen_anal
end