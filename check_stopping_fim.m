%==========================================================================
% Stopping Rule
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


if n>no_ini_pulses+4
%     rel_er_1= abs((t_est_f(n,:)-t_est_f(n-1,:))./t_est_f(n-1,:));
%     rel_er_2= abs((t_est_f(n-1,:)-t_est_f(n-2,:))./t_est_f(n-2,:));
%     rel_er_3= abs((t_est_f(n-2,:)-t_est_f(n-3,:))./t_est_f(n-3,:));
%     rel_er_4= abs((t_est_f(n-4,:)-t_est_f(n-5,:))./t_est_f(n-5,:));
%     rel_er_5= abs((t_est_f(n,:)-t_est_f(n-1,:))./t_est_f(n-1,:));
%     

    rel_er_t11= abs((t1_est_f(n,:)-t1_est_f(n-1,:))./t1_est_f(n-1,:));
    rel_er_t12= abs((t1_est_f(n-1,:)-t1_est_f(n-2,:))./t1_est_f(n-2,:));
    rel_er_t13= abs((t1_est_f(n-2,:)-t1_est_f(n-3,:))./t1_est_f(n-3,:));
    rel_er_t14= abs((t1_est_f(n-4,:)-t1_est_f(n-5,:))./t1_est_f(n-5,:));
    rel_er_t15= abs((t1_est_f(n,:)-t1_est_f(n-1,:))./t1_est_f(n-1,:));
    
    rel_er_t21= abs((t2_est_f(n,:)-t2_est_f(n-1,:))./t2_est_f(n-1,:));
    rel_er_t22= abs((t2_est_f(n-1,:)-t2_est_f(n-2,:))./t2_est_f(n-2,:));
    rel_er_t23= abs((t2_est_f(n-2,:)-t2_est_f(n-3,:))./t2_est_f(n-3,:));
    rel_er_t24= abs((t2_est_f(n-4,:)-t2_est_f(n-5,:))./t2_est_f(n-5,:));
    rel_er_t25= abs((t2_est_f(n,:)-t2_est_f(n-1,:))./t2_est_f(n-1,:));
    
    r_e_tm1=abs(taum_est_f(n)-taum_est_f(n-1))/taum_est_f(n-1); 
    r_e_tm2=abs(taum_est_f(n-1)-taum_est_f(n-2))/taum_est_f(n-2); 
    r_e_tm3=abs(taum_est_f(n-2)-taum_est_f(n-3))/taum_est_f(n-3); 
    r_e_tm4=abs(taum_est_f(n-3)-taum_est_f(n-4))/taum_est_f(n-4); 
    r_e_tm5=abs(taum_est_f(n-4)-taum_est_f(n-5))/taum_est_f(n-5); 
    r_e_tm=[r_e_tm1 r_e_tm2 r_e_tm3 r_e_tm4 r_e_tm5];
    
    if rel_er_t11<er_tol & rel_er_t12<er_tol & rel_er_t13<er_tol & rel_er_t14<er_tol & rel_er_t15<er_tol
        if rel_er_t21<er_tol & rel_er_t22<er_tol & rel_er_t23<er_tol & rel_er_t24<er_tol & rel_er_t25<er_tol
            if r_e_tm < er_tol_taum.*ones(1,5)
                n_conv_f=[n_conv_f n];
            end
        end
    end
    
end
    