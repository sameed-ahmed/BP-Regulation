% This is the nephron submodel for segmented water reabsorption. This 
% script is used to adjust different fractional water reabsorption in the 
% proximal tubule, distal tubule, and collecting duct, after inputting the
% glomerular filtration rate. Baseline fractional water reabsorption in 
% each segment is given by Rui's model of the transporters in each segment 
% of the nephron. 

% Rui's model is validated by experimental data on male and female rat
% water and sodium intake and excretion (as well as other quantities). Then
% there is a 'physiologically reasonable guess and check' for the 
% transporters along each nephron segment. As a result, the reabsorption 
% and flow of water, sodium, and other quantities are calculated in each 
% segment.

% Fixed quantities
% GFR, Frac-PT-Wreab, Frac-DT-Wreab, Frac-CD-Wreab

% Computed quantities
% PT-Wreab, DT-Wreab, CD-Wreab
% MD-UF, DT-UF; UF

function adjust_seg_wreab

% Gender.
gender = {'male','female'};
% Variables
vars = zeros(10,2);

for gg = 1:2 % gender

% Input GFR ml/min. Data from Munger - 1988.
if     strcmp(gender{gg},   'male')
    Phi_gfilt = 1.22;
elseif strcmp(gender{gg}, 'female')
    Phi_gfilt = 0.84;
end

if     strcmp(gender{gg},   'male')

% Fractional water reabsorption in each segment. Values from Layton - 2016.
eta_ptwreab = 0.86;
eta_dtwreab = 0.60;
eta_cdwreab = 0.78;

% Compute varying quantities.
Phi_ptwreab = Phi_gfilt * eta_ptwreab;
Phi_mdu = Phi_gfilt - Phi_ptwreab;
Phi_dtwreab = Phi_mdu * eta_dtwreab;
Phi_dtu = Phi_mdu - Phi_dtwreab;
Phi_cdwreab = Phi_dtu * eta_cdwreab;
Phi_u = Phi_dtu - Phi_cdwreab

vars(:,1) = [eta_ptwreab; eta_dtwreab; eta_cdwreab; ...
             Phi_ptwreab; Phi_dtwreab; Phi_cdwreab; ...
             Phi_gfilt  ; Phi_mdu    ; Phi_dtu    ; Phi_u];

% % Save quantities.
% save_data_name = sprintf('%s_seg_wreab_vars.mat', gender{gg});
% % save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'vars')

elseif strcmp(gender{gg}, 'female')

% % Compute varying quantities.
% Phi_u       = Phi_u; % same as male
% Phi_mdu     = Phi_mdu; % same as male
% Phi_ptwreab = Phi_gfilt - Phi_mdu;
% eta_ptwreab = Phi_ptwreab / Phi_gfilt
% eta_dtwreab = 0.60; % same as male
% Phi_dtwreab = Phi_mdu * eta_dtwreab;
% Phi_dtu     = Phi_mdu - Phi_dtwreab;
% Phi_cdwreab = Phi_dtu - Phi_u;
% eta_cdwreab = Phi_cdwreab / Phi_dtu

% Compute quantities with calibratiion.

% % male values
% eta_ptwreab = 0.86; % layton
% eta_dtwreab = 0.60;
% eta_cdwreab = 0.78;

% % Set pt, dt and compute cd
% Phi_u      = Phi_u; % same as male
% eta_ptwreab = 0.5 % calibrate
% Phi_ptwreab = Phi_gfilt * eta_ptwreab;
% Phi_mdu = Phi_gfilt - Phi_ptwreab;
% eta_dtwreab = 0.6 % calibrate
% Phi_dtwreab = Phi_mdu * eta_dtwreab;
% Phi_dtu     = Phi_mdu - Phi_dtwreab;
% Phi_cdwreab = Phi_dtu - Phi_u;
% eta_cdwreab = Phi_cdwreab / Phi_dtu

% Set dt, cd and compute pt
Phi_u      = Phi_u; % same as male
eta_dtwreab = 0.6 % calibrate
eta_cdwreab = 0.78 % calibrate 
eta_ptwreab = 1 - Phi_u / ( Phi_gfilt * (1 - eta_dtwreab) * (1 - eta_cdwreab) )

% vars(:,2) = [eta_ptwreab; eta_dtwreab; eta_cdwreab; ...
%              Phi_ptwreab; Phi_dtwreab; Phi_cdwreab; ...
%              Phi_gfilt  ; Phi_mdu    ; Phi_dtu    ; Phi_u];

% % Save quantities.
% save_data_name = sprintf('%s_seg_wreab_vars.mat', gender{gg});
% % save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'vars')

end % gender conditional

end % gender loop

end % function

























