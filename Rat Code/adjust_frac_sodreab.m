% This is the nephron submodel sodium handling. This script is used to 
% adjust different fractional sodium reabsorption in the proximal tubule, 
% and distal tubule while keeping the urine sodium excretion the same. As a 
% result, the rest of the quantities related to the nephron change.

% This script also calculates the factor by which the quanities change.
% This is needed to reparametrize any multiplicative factors that involve
% these variables so that they are still 1 at baseline.

function adjust_frac_sodreab

% Fix fixed quantities

% Phi_filsod = 0.0606880034847096; % female
% eta_ptsodreab = 0.5; % female
% Phi_usod = 4.200e-4; % female

Phi_filsod = 0.0904082029353203; % male
eta_ptsodreab = 0.8; % male
Phi_usod = 6.300e-4; % male

eta_dtsodreab = 0.5;

% Compute varying quantities

Phi_ptsodreab = Phi_filsod * eta_ptsodreab;
Phi_mdsod = Phi_filsod - Phi_ptsodreab;
Phi_dtsodreab = Phi_mdsod * eta_dtsodreab;
Phi_dtsod = Phi_mdsod - Phi_dtsodreab;
Phi_cdsodreab = Phi_dtsod - Phi_usod;
eta_cdsodreab = Phi_cdsodreab / Phi_dtsod;

% Adjust scaling factors to be 1 at baseline
% Scaling factors that are adjusted are equation numbers 10, 24, 56.

% All female. Male not neccessary.
Phi_ptsodreab_scale = Phi_ptsodreab/0.0485530172515515;
Phi_mdsod_scale = Phi_mdsod/0.0121349862331580;
Phi_dtsodreab_scale = Phi_dtsodreab/0.00607635498966443;
Phi_dtsod_scale = Phi_dtsod/0.00605863124349362;
Phi_cdsodreab_scale = Phi_cdsodreab/0.00563863124349362;

end

































