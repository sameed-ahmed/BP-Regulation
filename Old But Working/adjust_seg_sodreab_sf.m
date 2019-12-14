% This is the nephron submodel for segmented sodium reabsorption. This 
% script is used to adjust different fractional sodium reabsorption in the 
% proximal tubule, distal tubule, and collecting duct, after inputting the
% glomerular filtration rate and sodium concentration.

% Fixed quantities
% GFR, [Sod], Frac-PT-Sodreab, Frac-DT-Sodreab

% Computed quantities
% Frac-CD-Sodreab
% PT-Sodreab, DT-Sodreab, CD-Sodreab
% MD-SodF, DT-SodF; U-SodF

function adjust_seg_sodreab_sf

% Gender.
gender = {'male','female'};
% Variables
vars = zeros(11,2);

for gg = 1:2 % gender

% Input GFR ml/min. Data from Munger 1988.
% Fractional sodium reabsorption in each segment.
if     strcmp(gender{gg},   'male')
    Phi_gfilt = 1.22;
    eta_ptsodreab = 0.8; % Karaaslan 2005
    eta_dtsodreab = 0.5; 
    eta_cdsodreab = 0.93;
elseif strcmp(gender{gg}, 'female')
    Phi_gfilt = 0.84;
    eta_ptsodreab = 0.5; % calibrated
    eta_dtsodreab = 0.6; 
    eta_cdsodreab = 0.95;
end
% Input [Sod] micro Eq/ml.
C_sod = 143;

if     strcmp(gender{gg},   'male')

% Compute quantities.
Phi_filsod    = Phi_gfilt * C_sod
Phi_ptsodreab = Phi_filsod * eta_ptsodreab;
Phi_mdsod     = Phi_filsod - Phi_ptsodreab
Phi_dtsodreab = Phi_mdsod * eta_dtsodreab;
Phi_dtsod     = Phi_mdsod - Phi_dtsodreab
Phi_cdsodreab = Phi_dtsod * eta_cdsodreab;
Phi_usod      = Phi_dtsod - Phi_cdsodreab;

% SF_filsod = Phi_filsod / 18 
% SF_mdsod  = Phi_mdsod  / 3.6
% SF_dtsod  = Phi_dtsod  / 1.8

elseif strcmp(gender{gg}, 'female')

% Compute quantities.
Phi_filsod    = Phi_gfilt * C_sod
Phi_ptsodreab = Phi_filsod * eta_ptsodreab;
Phi_mdsod     = Phi_filsod - Phi_ptsodreab
Phi_dtsodreab = Phi_mdsod * eta_dtsodreab;
Phi_dtsod     = Phi_mdsod - Phi_dtsodreab
Phi_cdsodreab = Phi_dtsod * eta_cdsodreab;
Phi_usod      = Phi_dtsod - Phi_cdsodreab;

% SF_filsod = Phi_filsod / 18 
% SF_mdsod  = Phi_mdsod  / 3.6
% SF_dtsod  = Phi_dtsod  / 1.8

end % gender conditional

end % gender loop

end % function

























