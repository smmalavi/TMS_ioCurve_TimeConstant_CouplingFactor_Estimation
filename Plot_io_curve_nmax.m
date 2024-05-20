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
  

figure 
hold on

% ycurve_true=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
%     (1+(Vc_val/true_theta(3)).^true_theta(4)));


ycurve1_est_fim_nmax=real(t1_est_f(n,2)+(t1_est_f(n,1)-t1_est_f(n,2))./(1+(Vc_val/t1_est_f(n,3)).^t1_est_f(n,4)));
ycurve2_est_fim_nmax=real(t2_est_f(n,2)+(t2_est_f(n,1)-t2_est_f(n,2))./(1+(Vc_val/t2_est_f(n,3)).^t2_est_f(n,4)));




plot(Vc_base, 10.^y_base, 'dc','MarkerSize',8)

plot(Vc1_f(1:end), 10.^y1_f(1:end), 'ob','MarkerSize',8)

plot(Vc2_f(1:end), 10.^y2_f(1:end), 'or','MarkerSize',8)

plot(Vc_val,10.^ycurve_true1,'Color','[0 0 0.7]')
plot(Vc_val,10.^ycurve_true2,'Color','[0.7 0 0]')

plot(Vc_val, 10.^(ycurve1_est_fim_nmax), '--b','LineWidth',1)
plot(Vc_val, 10.^(ycurve2_est_fim_nmax), '--r','LineWidth',1)


% plot(Vc_u_matrix(end,1:end), 10.^y_u_matrix(end,1:end), '^b')
% plot(Vc_val, 10.^(ycurve_est_u_nmax), '--b','LineWidth',1)


xlabel('$V_C$ (normalized)', 'interpreter','latex')
ylabel('$y (V_{p-p})$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'log')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on
grid on


% strmax = ['c'];
% text(0.05,9e-3,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);




