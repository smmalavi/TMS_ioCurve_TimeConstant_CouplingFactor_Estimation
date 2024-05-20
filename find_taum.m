%==========================================================================
% find time constant 
%==========================================================================
% Main Matlab file for 
% “Identifiability analysis and noninvasive online estimation of 
% the first-order neural activation dynamics in the brain with 
% closed-loop transcranial magnetic stimulation,” 
% IEEE Trans on Biomedical Engineering, 70(9), 2564-2572, 2023.
%==========================================================================
%
% Seyed Mohammad Mahdi Alavi+, Stellantis (Chrysler), Canada 
% Fidel Vila-Rodriguez, Unitverisyt of British Columbia, Canada 
% Adam Mahdi, University of Oxford, UK
% Stefan M. Goetz, University of Cambridge (UK), Duke University (USA)
% +: code written by
% e-mail: mahdi.alavi.work@gmail.com
%
% April 2022
%==========================================================================


function taum_sol = find_taum(Tp1, Tp2, t31, t32, taum, mu, sigma, w)

taum_sol=abs((((mu*taum-sigma)*sin(w*Tp1)+w*cos(w*Tp1))*exp(-sigma*Tp1)-w*exp(-Tp1/taum))*t31-...
    (((mu*taum-sigma)*sin(w*Tp2)+w*cos(w*Tp2))*exp(-sigma*Tp2)-w*exp(-Tp2/taum))*t32);

end
