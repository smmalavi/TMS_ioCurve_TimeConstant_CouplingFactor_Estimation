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


%close all
%figure
hold on

x_dg=linspace(0,1,10000);

myy1=real(true_theta1(2)+(true_theta1(1)-true_theta1(2))./...
    (1+((x_dg+sigma_x*randn(1,length(x_dg)))/true_theta1(3)).^true_theta1(4))+...
    sigma_y*randn(1,length(x_dg)));

plot(x_dg, 10.^myy1,'x','Color', 'b' ,'MarkerSize',8)


myy2=real(true_theta2(2)+(true_theta2(1)-true_theta2(2))./...
    (1+((x_dg+sigma_x*randn(1,length(x_dg)))/true_theta2(3)).^true_theta2(4))+...
    sigma_y*randn(1,length(x_dg)));
plot(x_dg, 10.^myy2,'x','Color', 'r','MarkerSize',8)


plot(Vc_val,10.^ycurve_true2,'Color','[0.7 0 0]','LineWidth',1)
plot(Vc_val,10.^ycurve_true1,'Color','[0 0 .7]','LineWidth',1)




%xlabel('$V_C$ (normalized)', 'interpreter','latex')
ylabel('$y (V_{p-p})$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'log')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on
grid on

 
%fig_font_size=16;
% legend
strmax = ['a'];
text(0.05,6e-3,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);






%% % legend
strmax = ['Ref. data'];
text(0.72,1.75e-6,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['Ref. IO curve'];
text(0.72,8.2e-7,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['Baseline'];
text(0.72,3.85e-7,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['Samples'];
text(0.72,1.8e-7,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['Estimation'];
text(0.72,8.1e-8,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);




line([0.61 0.71],[8.2e-7 8.2e-7],'Color','k','LineStyle','-');
%line([0.6 0.72],[4.2e-6 4.2e-6],'Color','b','LineStyle','--');
line([0.61 0.71],[8.1e-8 8.1e-8],'Color','k','LineStyle','--');

plot(0.66, 1.75e-6, 'Marker','x', 'MarkerSize',8,'MarkerEdgeColor','k');
plot(0.66, 3.85e-7, 'Marker','d', 'MarkerSize',8,'MarkerEdgeColor','k');
plot(0.66, 1.72e-7, 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor','k');

