% Leete and Layton 2019 Triple Whammy AKI
% Run all simulations and create figures 4 and 5 for the paper

% It reads from files in Human_IG_Data/. It writes to folders called Human_Data/ and Human_Figures/.
% The files in Human_IG_Data/ are sufficient initial guesses for the steady state solver.
% Steady state solutions are then used as initial conditions in the time
% course simulations

% species indicator
human = 1;
% sexes to model
gender   = {'male', 'female'};
health_labels = {'normotensive','hypertensive'};
% factor by which to multiply rsna baseline (N_rsna)
rsna = [1,2.5]; %1 for normotensive, 2.5 for hypertensive
% drug treatments. Columns are ACEi, furosemide, and NSAID (normal and high doses)
models = [0,0,0;... %normal dose
    1,0,0;...
    0,1,0;...
    0,0,1;...
    1,1,0;...
    0,1,1;...
    1,0,1;...
    1,1,1;...
    2,0,0;...% high dose ACEi 9
    0,2,0;...% high dose furosemide 10
    0,0,2;...% high dose NSAID 11
    2,2,0;...% high dose combos 12
    0,2,2;...
    2,0,2;...
    2,2,2];

% Find steady state solutions. These will be the initial conditions for the
% time course model.
% This will loop through available initial guesses for the solver. Will
% display if a solution was found for each simulation or not.
disp('----------Computing steady state solutions----------')
disp({'sex','normo/hyp','impaired myogenic'})
for gg= 1:2
    for AA = 1:2
        for c =  0:1 % loop through healthy and impaired myogenic response
            disp({gender{gg}, health_labels{AA},num2str(c)})
            for i = [1,4,11] %loop through models to simulate
                for j=1%[1,4,11]%:length(models) % loop through data for initial guesses for solver
                    IG = sprintf('Human_IG_Data/human_%s_ss_%s_%s_%s_rsna%s.mat',gender{gg},num2str(models(j,1)),num2str(models(j,2)),num2str(models(j,3)),num2str(rsna(AA)));
                    %disp(IG)
                    [exitflag, imag] = solve_ss_numerical(human,gg,IG,'ACEi',models(i,1),'furosemide',models(i,2),'NSAID',models(i,3),'RSNA',rsna(AA),'Impaired Myogenic Response',c);
                    if (exitflag >0) && (imag == 0)
                        disp(['found solution for ', num2str(i)])
                        break
                    end
                    if (j==1)  && ((exitflag == 0) || (imag == 1))
                        disp(['no proper IG found for ', num2str(i)])
                    end   
                end
            end
        end
    end
end 
        
% % % Run time course simulations  
disp('----------Computing time course simulations----------')
disp({'sex','normo/hyp','impaired myogenic','low water intake'})
myo_ending = {'','_impairedmyo'};
for gg= 1:2 %loop through sex
    for AA = 1:2 %loop through normo/hypertension
        for c = 0:1 % loop through healthy and impaired myogenic response
            for w = 0:1 % loop through normal and low water intake
                disp({gender{gg}, health_labels{AA},num2str(c),num2str(w)})
                for i = 1:length(models) %loop through drug treatments to simulate
                    disp({char(9),'ACEi',num2str(models(i,1)),'furosemide',num2str(models(i,2)),'NSAID',num2str(models(i,3))})
                    IG = sprintf('Human_Data/human_%s_ss_%s_%s_%s_rsna%s%s.mat',gender{gg},num2str(0),num2str(0),num2str(models(i,3)),...
                                 num2str(rsna(AA)),myo_ending{c+1});
                    %disp(IG)
                    run_sim_treatments(human,gg,IG,'ACEi',models(i,1),'furosemide',models(i,2),'NSAID',models(i,3),...
                                                          'RSNA',rsna(AA),'Impaired Myogenic Response',c,'Low Water Intake', w);
                end
            end
        end
    end
end 
% % 
% % % Plot 2 day scatter plots (Figure 4)
plot_scatter();
% % 
% % % Plot risk factor bar charts (Figure 5)
plot_riskfactors();
        
        