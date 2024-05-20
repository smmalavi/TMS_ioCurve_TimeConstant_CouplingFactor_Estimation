%==========================================================================
% main file for plot
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




close all

fig_font_size=16;
leg_font_size=16;


%% io
ylim=[4*10^-8 10^-2];;

xlim=[0 1];

fig=figure 

Plot_io_curve; 
ax=gca;
ax.XLim=xlim;
ax.YLim=ylim;
strmax = ['a']; 
text(0.05,6e-3,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);



% 
strmax = ['Pulse Width:'];
text(0.08,1.8e-7,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'k','FontSize', leg_font_size);


strmax = ['$$120~\mu s$$'];

text(0.155,8e-8,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'r','FontSize', leg_font_size);
strmax = [','];
text(0.3,8e-8,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'k','FontSize', leg_font_size);


strmax = ['$$54~\mu s$$'];

text(0.33,8e-8,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'b','FontSize', leg_font_size);


xxmin=0.07;
xxmax=0.47;
yymin=5e-8;
yymax=3e-7;
plot([xxmin xxmax xxmax xxmin xxmin],[yymin yymin yymax yymax yymin],'k')



file_path_name=['fig-io-curve-data'];
saveas(gcf,file_path_name,'fig')
%saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')

%%


Plot_io_curve_nmax;
ax=gca;
ax.XLim=xlim;
ax.YLim=ylim;

%line([0 1],[50e-6 50e-6])
strmax = ['c']; 
text(0.05,6e-3,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);

file_path_name=['fig-io-nmax'];
saveas(gcf,file_path_name,'fig')
%saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')

%%
Plot_io_curve_nf;
ax=gca;
ax.XLim=xlim;
ax.YLim=ylim;

strmax = ['b']; 
text(0.05,6e-3,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);

file_path_name=['fig-io-nf'];
saveas(gcf,file_path_name,'fig')
%saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')


%% io parameters 
figure
hold on
plot([0 n],[true_yl true_yl],'k','LineWidth',1)
%plot(t_est_u(:,1),'--b','LineWidth',1)
plot(t1_est_f(:,1),'--b','LineWidth',1.5)
plot(t2_est_f(:,1),'--r','LineWidth',0.5)

%legend('True', 'Estimation','Location','northwest','Box','off')
xlabel('$n$', 'interpreter','latex')
ylabel('$\theta_1$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.YLim=[-7 0];
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

strmax = ['a'];
text(465,-.5,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);


strmax = ['True'];
text(75,-.5,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'k','FontSize', leg_font_size);

strmax = ['Estimation'];
text(75,-1,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'k','FontSize', leg_font_size);

line([15 70],[-0.5 -0.5],'Color','k','LineStyle','-');
line([15 70],[-1 -1],'Color','k','LineStyle','--');

strmax = ['Pulse Width:'];
text(20,-1.75,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'k','FontSize', leg_font_size);

strmax = ['$$120~\mu s$$'];
text(39,-2.25,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'r','FontSize', leg_font_size);
strmax = [','];
text(110,-2.25,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'k','FontSize', leg_font_size);

strmax = ['$$54~\mu s$$'];
text(125,-2.25,strmax,'interpreter','latex','HorizontalAlignment','left',...
    'FontName', 'Times New Roman', 'Color', 'b','FontSize', leg_font_size);

xxmin=15;
xxmax=200;
yymin=-2.5;
yymax=-1.4;
plot([xxmin xxmax xxmax xxmin xxmin],[yymin yymin yymax yymax yymin],'k')


file_path_name=['fig-yl'];
saveas(gcf,file_path_name,'fig')
%saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')

%%

figure
hold on
plot([0 n],[true_yh true_yh],'k')
plot(t1_est_f(:,2),'--b','LineWidth',1)
plot(t2_est_f(:,2),'--r','LineWidth',1)

%plot(t_est_u(:,2),'--b','LineWidth',1)

xlabel('$n$', 'interpreter','latex')
ylabel('$\theta_2$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

strmax = ['b'];
text(465,-.24,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);


file_path_name=['fig-yh'];
saveas(gcf,file_path_name,'fig')
%saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')


%%
figure
hold on
plot([0 n],[true_theta1(3) true_theta1(3)],'b')
plot([0 n],[true_theta2(3) true_theta2(3)],'r')

plot(t1_est_f(:,3),'--b','LineWidth',1)
plot(t2_est_f(:,3),'--r','LineWidth',1)

%plot(t_est_u(:,3),'--b','LineWidth',1)

xlabel('$n$', 'interpreter','latex')
ylabel('$\theta_3$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.YLim=[0 1];
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

strmax = ['c'];
text(465,0.92,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);


file_path_name=['fig-m'];
saveas(gcf,file_path_name,'fig')
%saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')


%%
figure
hold on
plot([0 n],[true_s1 true_s1],'b','LineWidth',1)
plot([0 n],[true_s2 true_s2],'r','LineWidth',1)

plot(t1_est_f(:,4),'--b','LineWidth',1)
plot(t2_est_f(:,4),'--r','LineWidth',1)

%plot(t_est_u(:,4),'--b','LineWidth',1)

%legend('True','Estimate')
xlabel('$n$', 'interpreter','latex')
ylabel('$\theta_4$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

strmax = ['d'];
text(465,92,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);


file_path_name=['fig-s'];
saveas(gcf,file_path_name,'fig')
%saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')

%%
figure
hold on
plot([0 n],[true_taum true_taum],'k','LineWidth',1)
%plot(taum_est_u,'--b','LineWidth',1)
plot(taum_est_f,'--k','LineWidth',1)
legend('True', 'Estimation','Location','northeast','Box','off')

%legend('True','Uniform Sampling', 'FIM SPE','Location','southeast','Box','off')
xlabel('$n$', 'interpreter','latex')

yt = get(gca, 'YTick');                                 % 'XTick' Values
set(gca, 'YTick', yt, 'YTickLabel', yt/1e-6) 
ylabel('$\tau~ (\mu s)$','interpreter','latex')
%ylabel('$\tau_m$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

file_path_name=['fig-taum'];
saveas(gcf,file_path_name,'fig')
%saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')

%%
figure
hold on
plot([0 n],[true_gr true_gr],'k','LineWidth',1)
%plot(taum_est_u,'--b','LineWidth',1)
plot(gr_est_f,'--k','LineWidth',1)
legend('True', 'Estimation','Location','northeast','Box','off')

%legend('True','Uniform Sampling', 'FIM SPE','Location','southeast','Box','off')
xlabel('$n$', 'interpreter','latex')

%yt = get(gca, 'YTick');                                 % 'XTick' Values
%set(gca, 'YTick', yt, 'YTickLabel', yt/1e-6) 
ylabel('$g$','interpreter','latex')
%ylabel('$\tau_m$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.YLim=[0 50];

ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

file_path_name=['fig-g'];
saveas(gcf,file_path_name,'fig')
%saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')


