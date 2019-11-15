

% experiment = 'baseline steady state';
experiment = 'baseline dynamic';
% experiment = 'pressure-natriuresis';
% experiment = 'pressure-natriuresis whole range';
% experiment = 'sodium intake';
% experiment = 'Ang II infusion';

if     strcmp(experiment, 'baseline steady state')
elseif strcmp(experiment, 'baseline dynamic')
    run_sim;
end