%==========================================================================
% Main file
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

%
clear all
close all
clc

rng('shuffle'); 

% create a directry
% folder_index=randi([1 100],1,1);
% parent_folder=['gfx' num2str(folder_index)]
% mkdir(parent_folder)

%
n_max=500;


% cTMS pulse parameters 
R=0.1;%Ohm
r=20e-3;%mOhm
L=16e-6;%H
C=716e-6;%F
delta=3.2e-6;%(V/m)(A/s)
sigma=r/(2*L);
w=sqrt(1/(L*C)-sigma^2);
mu=(w^2+sigma^2);
k1=delta/(L*w);
min_Vc=0.01;
max_Vc=1;
Vc_val = linspace(0, 1, 100);
size_Vc_val=length(Vc_val);

min_taum=90e-6; %1e-6;
max_taum=220e-6;
x_dg=linspace(0,1,2000);


sigma_y=0.1;
sigma_x=0.05;
true_yl=-6;%-6.5+ (-5.5-(-6.5))*rand(1,1);
true_yh=-2.62;%-3+ (-2-(-3))*rand(1,1);
true_s1=8.84;%1+ (15-1)*rand(1,1);
true_s2=14.61;%(true_s1+1)+ (25-(true_s1+1))*rand(1,1);

true_taum=185e-6;%min_taum+(max_taum-min_taum)*rand(1,1);

min_Tp=10e-6;
%max_Tp=120e-6;

a =       97.54 ;% (95.69, 99.4)
b =    0.001206;%  (0.001108, 0.001304)
c =      -80.57;%  (-82.06, -79.08)
d =    -0.025 ;% (-0.02608, -0.02358)
max_Tp= (a*exp(b*true_taum*1e6) + c*exp(d*true_taum*1e6))
% using critical Tp

true_gr=37.12;%30+ (50-30)*rand(1,1);


% ObjFunc_Tp = @(z) find_next_Tp(z,[],true_taum,true_gr, k1, mu, sigma, w);% 
%     optsTp = optimoptions(@fmincon,'Algorithm','interior-point');
%     problem = createOptimProblem('fmincon','x0',max_Tp*rand,...
%                     'objective',ObjFunc_Tp,'lb',min_Tp,'ub',max_Tp,'options',optsTp);
% [true_Tp,fval_fimopt,flagm_fimopt,outptm_fimopt,manyminsm_fimopt] = run(GlobalSearch,problem);

% true_Tp1=round(10+ (70-10)*rand(1,1))*1e-6;
% true_Tp2=round((true_Tp1+10)+ (120-(true_Tp1+10))*rand(1,1))*1e-6;

% true_Tp1=round(10+ (max_Tp-10)*rand(1,1))*1e-6;
% true_Tp2=round(10+ (max_Tp-10)*rand(1,1))*1e-6;

true_Tp1=54e-6;%round(10+ (50-10)*rand(1,1))*1e-6;%
true_Tp2=120e-6;%round(true_Tp1*1e6+ (max_Tp-(true_Tp1*1e6))*rand(1,1))*1e-6;

true_theta1(1)=true_yl;
true_theta1(2)=true_yh;
true_theta1(4)=true_s1;

true_theta2(1)=true_yl;
true_theta2(2)=true_yh;
true_theta2(4)=true_s2;

%
true_tilde_rp1=k1/(mu*true_taum^2-2*sigma*true_taum+1)*...
    (((mu*true_taum-sigma)*sin(w*true_Tp1)+w*cos(w*true_Tp1))*exp(-sigma*true_Tp1)-...
    w*exp(-true_Tp1/true_taum));

true_tilde_rp2=k1/(mu*true_taum^2-2*sigma*true_taum+1)*...
    (((mu*true_taum-sigma)*sin(w*true_Tp2)+w*cos(w*true_Tp2))*exp(-sigma*true_Tp2)-...
    w*exp(-true_Tp2/true_taum));

%
true_theta1(3)=1/(true_gr*true_tilde_rp1);
true_theta2(3)=1/(true_gr*true_tilde_rp2);

% IO curve Tp1
myy1=real(true_theta1(2)+(true_theta1(1)-true_theta1(2))./...
    (1+((x_dg+sigma_x*randn(1,length(x_dg)))/true_theta1(3)).^true_theta1(4))+...
    sigma_y*randn(1,length(x_dg)));


ycurve_true1=real(true_theta1(2)+(true_theta1(1)-true_theta1(2))./...
    (1+(Vc_val/true_theta1(3)).^true_theta1(4)));

% IO curve Tp2
myy2=real(true_theta2(2)+(true_theta2(1)-true_theta2(2))./...
    (1+((x_dg+sigma_x*randn(1,length(x_dg)))/true_theta2(3)).^true_theta2(4))+...
    sigma_y*randn(1,length(x_dg)));


ycurve_true2=real(true_theta2(2)+(true_theta2(1)-true_theta2(2))./...
    (1+(Vc_val/true_theta2(3)).^true_theta2(4)));

%end
%
hold on
plot(x_dg, 10.^myy1,'xb')
plot(Vc_val,10.^ycurve_true1,'b')
plot(x_dg, 10.^myy2,'xr')
plot(Vc_val,10.^ycurve_true2,'r')

xlabel('$V_C$ (normalized)', 'interpreter','latex')
ylabel('$y (V_{p-p})$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'log')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = 18;
box on

grid on


%%
% Get basline 

n_base=50;% number of baseline data
Vc_base=zeros(1,n_base);
%y_base = virtualsubjectEIVStimulate_01(Vc_base, myp);% original data
%plot(Vc_base, y_base, 'dc')

y_base=real(true_theta1(2)+(true_theta1(1)-true_theta1(2))./...
    (1+((Vc_base+sigma_x*randn(1,length(Vc_base)))/true_theta1(3)).^true_theta1(4))+...
    sigma_y*randn(1,length(Vc_base)));



hold on
plot(Vc_base, 10.^y_base, 'dc')



%

% intialization
no_ini_pulses=3;
%Tp_f_ini=repmat(Tp_min+(Tp_max-Tp_min)*rand(1,1),1,3);
ini_Vc=min_Vc+(max_Vc-min_Vc)*rand(1,no_ini_pulses);
%ini_Tp=100e-6*ones(1,3);%min_Tp+(max_Tp-min_Tp)*rand(1,no_ini_pulses);
ini_Tp1=repmat(true_Tp1,1,3);%true_Tp1,true_Tp1];
ini_Tp2=repmat(true_Tp2,1,3);
ini_taum=min_taum+(max_taum-min_taum)*rand(1,1);


ini_tilde_rp1=k1/(mu*true_taum^2-2*sigma*true_taum+1)*...
    (((mu*true_taum-sigma)*sin(w*ini_Tp1)+w*cos(w*ini_Tp1)).*exp(-sigma*ini_Tp1)-...
    w*exp(-ini_Tp1/true_taum));


ini_y1=real(true_theta1(2)+(true_theta1(1)-true_theta1(2))./...
    (1+(true_gr*ini_tilde_rp1.*(ini_Vc+sigma_x*randn(1,length(ini_Vc)))).^true_theta1(4))+...
    sigma_y*randn(1,length(ini_Vc)));

ini_tilde_rp2=k1/(mu*true_taum^2-2*sigma*true_taum+1)*...
    (((mu*true_taum-sigma)*sin(w*ini_Tp2)+w*cos(w*ini_Tp2)).*exp(-sigma*ini_Tp2)-...
    w*exp(-ini_Tp2/true_taum));


ini_y2=real(true_theta2(2)+(true_theta2(1)-true_theta2(2))./...
    (1+(true_gr*ini_tilde_rp2.*(ini_Vc+sigma_x*randn(1,length(ini_Vc)))).^true_theta2(4))+...
    sigma_y*randn(1,length(ini_Vc)));


figure(1)
hold on
plot(ini_Vc, 10.^ini_y1, 'sg')
plot(ini_Vc, 10.^ini_y2, 'sg')


n=no_ini_pulses;
Vc1_f=ini_Vc;
Vc2_f=ini_Vc;
Tp1_f=repmat(true_Tp1,1,n_max);%Tp1 is known and fixed a priori
Tp2_f=repmat(true_Tp2,1,n_max);%Tp2 is known and fixed a priori
tilde_rp1_f=ini_tilde_rp1;
tilde_rp2_f=ini_tilde_rp2;
y1_f=ini_y1;
y2_f=ini_y2;

% IO curve 1 fitting 
[xData, yData] = prepareCurveData( [Vc_base Vc1_f], [y_base y1_f]);
cf_func1='b+(a-b)./(1+(x/c)^d)';
ft1 = fittype( cf_func1, 'independent', 'x', 'dependent', 'y' );
paramLB1=[-7 -3 0  1];%yl, yh, m, s
paramUB1=[-5 -2 1  100];%yl, yh, m, s
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Algorithm = 'Trust-Region';
opts.Lower = paramLB1;
opts.Upper = paramUB1;
opts.Robust ='LAR';%'LAR';% 'Bisquare';
ini_guess=paramLB1+ (paramUB1-paramLB1).*rand(1,4) ; 
opts.StartPoint =ini_guess;
         
[fitresult, gof] = fit( xData, yData, ft1, opts );

t1_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d];

ycurve1_est_fim_ini=real(t1_est_f(n,2)+(t1_est_f(n,1)-t1_est_f(n,2))./(1+(Vc_val/t1_est_f(n,3)).^t1_est_f(n,4)));

% IO curve 2 fitting 
[xData, yData] = prepareCurveData( [Vc_base Vc2_f], [y_base y2_f]);
cf_func1='b+(a-b)./(1+(x/c)^d)';
ft1 = fittype( cf_func1, 'independent', 'x', 'dependent', 'y' );
% paramLB1=[-7   -3 0  1];%yl, yh, m, s
% paramUB1=[-5 -2 1  50];%yl, yh, m, s
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Algorithm = 'Trust-Region';
opts.Lower = paramLB1;
opts.Upper = paramUB1;
opts.Robust ='LAR';%'LAR';% 'Bisquare';
ini_guess=paramLB1+ (paramUB1-paramLB1).*rand(1,4) ; 
opts.StartPoint =ini_guess;
         
[fitresult, gof] = fit( xData, yData, ft1, opts );

t2_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d];

ycurve2_est_fim_ini=real(t2_est_f(n,2)+(t2_est_f(n,1)-t2_est_f(n,2))./(1+(Vc_val/t2_est_f(n,3)).^t2_est_f(n,4)));

% initial estimate of tau

taum_func = @(taum) find_taum(Tp1_f(no_ini_pulses),Tp2_f(no_ini_pulses),...
    t1_est_f(n,3),t2_est_f(n,3), taum, mu, sigma, w);% 
    opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
%     problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
%         'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
% 
% [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
% 

   
    x0=min_Tp+(max_Tp-min_Tp)*rand(1,1);
    [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = fmincon(taum_func,x0,[],[],[],[],min_taum,max_taum,[],opts_taum);

% initial estimate of g_r

tilde_rp1_f(n)=k1/(mu*taum_est_f(n)^2-2*sigma*taum_est_f(n)+1)*...
    (((mu*taum_est_f(n)-sigma)*sin(w*Tp1_f(n))+w*cos(w*Tp1_f(n))).*exp(-sigma*Tp1_f(n))-...
    w*exp(-Tp1_f(n)/taum_est_f(n)));

tilde_rp2_f(n)=k1/(mu*taum_est_f(n)^2-2*sigma*taum_est_f(n)+1)*...
    (((mu*taum_est_f(n)-sigma)*sin(w*Tp2_f(n))+w*cos(w*Tp2_f(n))).*exp(-sigma*Tp2_f(n))-...
    w*exp(-Tp2_f(n)/taum_est_f(n)));

gr_est_f(n)=1/(t1_est_f(n,3)*tilde_rp1_f(n));

%
n_conv_f=[];
% n_conv_u=[];
er_tol=.01*ones(1,4);
er_tol_taum=.01*ones(1,5);

% Vc_u_matrix=NaN(n_max,n_max);
% y_u_matrix=NaN(n_max,n_max);

%bad_fit_flag=0;

for n=no_ini_pulses+1:n_max
    n
    
%     ObjFunc_Tp = @(z) find_next_Tp(z,Tp1_f,taum_est_f(n-1),true_gr, k1, mu, sigma, w);% 
%     optsTp = optimoptions(@fmincon,'Algorithm','interior-point');
%     problem = createOptimProblem('fmincon','x0',max_Tp*rand,...
%                     'objective',ObjFunc_Tp,'lb',min_Tp,'ub',max_Tp,'options',optsTp);
%     [Tp1_f(n),fval_fimopt,flagm_fimopt,outptm_fimopt,manyminsm_fimopt] = run(GlobalSearch,problem);

    % 1st IOC fitting
    ObjFunc_FIM = @(x) find_next_Vc(x,Vc1_f,t1_est_f(n-1,:));% 
    optsFIM = optimoptions(@fmincon,'Algorithm','interior-point');
    % problem = createOptimProblem('fmincon','x0',max_Vc*rand,...
    %                 'objective',ObjFunc_FIM,'lb',min_Vc,'ub',max_Vc,'options',optsFIM);
    % [Vc1_f(n),fval_fimopt,flagm_fimopt,outptm_fimopt,manyminsm_fimopt] = run(GlobalSearch,problem);
    % 
   
    x0=max_Vc*rand;
    [Vc1_f(n),fval_fimopt,flagm_fimopt,outptm_fimopt,manyminsm_fimopt] = fmincon(ObjFunc_FIM,x0,[],[],[],[],min_Vc,max_Vc,[],optsFIM);



    y1_f(n)=real(true_theta1(2)+(true_theta1(1)-true_theta1(2))./...
        (1+(true_gr*true_tilde_rp1*(Vc1_f(n)+sigma_x*randn(1,1))).^true_theta1(4))+...
        sigma_y*randn(1,1));

   
    [xData, yData] = prepareCurveData( [Vc_base Vc1_f], [y_base y1_f]);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Algorithm = 'Trust-Region';
    opts.Robust ='LAR';%'LAR';% 'Bisquare';
    opts.Lower = paramLB1;
    opts.Upper = paramUB1;
    opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
    [fitresult, gof] = fit( xData, yData, ft1, opts );
    t1_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
    
    % check bad-fit two times
    if n>no_ini_pulses+4
        % first time
        rel_er_1= abs((t1_est_f(n,:)-t1_est_f(n-1,:))./t1_est_f(n-1,:));
        if ~isempty(find(rel_er_1 > 20*er_tol))            
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t1_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
        end           
        % 2nd time
        rel_er_1= abs((t1_est_f(n,:)-t1_est_f(n-1,:))./t1_est_f(n-1,:));
        if ~isempty(find(rel_er_1 > 20*er_tol)) 
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t1_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
        end
        % 3rd time
        rel_er_1= abs((t1_est_f(n,:)-t1_est_f(n-1,:))./t1_est_f(n-1,:));
        if ~isempty(find(rel_er_1 > 20*er_tol)) 
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t1_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
        end
    end
        
    
    % 2nd IOC fitting
    ObjFunc_FIM = @(x) find_next_Vc(x,Vc2_f,t2_est_f(n-1,:));% 
    optsFIM = optimoptions(@fmincon,'Algorithm','interior-point');
    % problem = createOptimProblem('fmincon','x0',max_Vc*rand,...
    %                 'objective',ObjFunc_FIM,'lb',min_Vc,'ub',max_Vc,'options',optsFIM);
    % [Vc2_f(n),fval_fimopt,flagm_fimopt,outptm_fimopt,manyminsm_fimopt] = run(GlobalSearch,problem);
    % 


    x0=max_Vc*rand;
    [Vc2_f(n),fval_fimopt,flagm_fimopt,outptm_fimopt,manyminsm_fimopt] = fmincon(ObjFunc_FIM,x0,[],[],[],[],min_Vc,max_Vc,[],optsFIM);


   
    y2_f(n)=real(true_theta2(2)+(true_theta2(1)-true_theta2(2))./...
        (1+(true_gr*true_tilde_rp2*(Vc2_f(n)+sigma_x*randn(1,1))).^true_theta2(4))+...
        sigma_y*randn(1,1));

   
    [xData, yData] = prepareCurveData( [Vc_base Vc2_f], [y_base y2_f]);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Algorithm = 'Trust-Region';
    opts.Robust ='LAR';%'LAR';% 'Bisquare';
    opts.Lower = paramLB1;
    opts.Upper = paramUB1;
    opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
    [fitresult, gof] = fit( xData, yData, ft1, opts );
    t2_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
    
    % check bad-fit two times
    if n>no_ini_pulses+4
        % first time
        rel_er_1= abs((t2_est_f(n,:)-t2_est_f(n-1,:))./t2_est_f(n-1,:));
        if ~isempty(find(rel_er_1 > 20*er_tol))            
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t2_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
        end        
        % 2nd time
        rel_er_1= abs((t2_est_f(n,:)-t2_est_f(n-1,:))./t2_est_f(n-1,:));
        if ~isempty(find(rel_er_1 > 20*er_tol)) 
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t2_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
        end
        
        % 3rd time
        rel_er_1= abs((t2_est_f(n,:)-t2_est_f(n-1,:))./t2_est_f(n-1,:));
        if ~isempty(find(rel_er_1 > 20*er_tol)) 
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t2_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
        end
    end
    
    
    % find taum(n)
%     taum_func = @(taum) find_taum_fixed_Tp(Tp1_f(n),taum, t1_est_f(n,3), true_gr, k1, mu, sigma, w);% 
%     opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
%     problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
%         'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
%     [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
  
%     pts = min_taum+ (max_taum-min_taum).*rand(200,1);
%     tpoints = CustomStartPointSet(pts);
%     allpts = {tpoints};
%                 
%     [taum_est_f(n),fval_f(n),flagm_ms,outptm_ms,manyminsm_ms] = run(MultiStart,problem,allpts);

%  estimate of tau

    taum_func = @(taum) find_taum(Tp1_f(n),Tp2_f(n),t1_est_f(n,3),t2_est_f(n,3),taum,mu,sigma,w);% 
    opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
    % problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
    %     'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
    % [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
    % 
  
    x0=min_Tp+(max_Tp-min_Tp)*rand(1,1);
    [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = fmincon(taum_func,x0,[],[],[],[],min_taum,max_taum,[],opts_taum);



%     options = optimoptions('particleswarm','SwarmSize',100);
%     [taum_est_f(n),fval(n),exitflag] = particleswarm(taum_func,1,120e-6,max_taum,options);
    
    % Check bad-fit for tau
    if n>no_ini_pulses+4
        % first time
        r_e_tm1=abs(taum_est_f(n)-taum_est_f(n-1))/taum_est_f(n-1); 
        if  r_e_tm1 > 20*er_tol_taum(1)           
            taum_func = @(taum) find_taum(Tp1_f(no_ini_pulses),Tp2_f(no_ini_pulses),...
                t1_est_f(n,3),t2_est_f(n,3), taum, mu, sigma, w);% 
            opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
            % problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
            %     'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
            % [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);


            x0=min_Tp+(max_Tp-min_Tp)*rand(1,1);
            [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = fmincon(taum_func,x0,[],[],[],[],min_taum,max_taum,[],opts_taum);

        end
        
        % 2nd time 
        r_e_tm1=abs(taum_est_f(n)-taum_est_f(n-1))/taum_est_f(n-1); 
        if  r_e_tm1 > 20*er_tol_taum(1)           
            taum_func = @(taum) find_taum(Tp1_f(no_ini_pulses),Tp2_f(no_ini_pulses),...
                t1_est_f(n,3),t2_est_f(n,3), taum, mu, sigma, w);% 
            opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
            % problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
            %     'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
            % [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);

                        x0=min_Tp+(max_Tp-min_Tp)*rand(1,1);
            [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = fmincon(taum_func,x0,[],[],[],[],min_taum,max_taum,[],opts_taum);

        end
        
        % 3rd time 
        r_e_tm1=abs(taum_est_f(n)-taum_est_f(n-1))/taum_est_f(n-1); 
        if  r_e_tm1 > 20*er_tol_taum(1)           
            taum_func = @(taum) find_taum(Tp1_f(no_ini_pulses),Tp2_f(no_ini_pulses),...
                t1_est_f(n,3),t2_est_f(n,3), taum, mu, sigma, w);% 
            opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
            % problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
            %     'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
            % [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
                        x0=min_Tp+(max_Tp-min_Tp)*rand(1,1);
            [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = fmincon(taum_func,x0,[],[],[],[],min_taum,max_taum,[],opts_taum);

        end
    end
 
    % Estimate of g_r
    tilde_rp1_f(n)=k1/(mu*taum_est_f(n)^2-2*sigma*taum_est_f(n)+1)*...
        (((mu*taum_est_f(n)-sigma)*sin(w*Tp1_f(n))+w*cos(w*Tp1_f(n))).*exp(-sigma*Tp1_f(n))-...
        w*exp(-Tp1_f(n)/taum_est_f(n)));

    tilde_rp2_f(n)=k1/(mu*taum_est_f(n)^2-2*sigma*taum_est_f(n)+1)*...
        (((mu*taum_est_f(n)-sigma)*sin(w*Tp2_f(n))+w*cos(w*Tp2_f(n))).*exp(-sigma*Tp2_f(n))-...
        w*exp(-Tp2_f(n)/taum_est_f(n)));

    gr_est_f(n)=1/(t1_est_f(n,3)*tilde_rp1_f(n));
    
    % check bad-fit first time
    
%     if n>no_ini_pulses+4
%         rel_er_1= abs((t1_est_f(n,:)-t1_est_f(n-1,:))./t1_est_f(n-1,:));
%         rel_er_2= abs((t2_est_f(n,:)-t2_est_f(n-1,:))./t2_est_f(n-1,:));
%         r_e_tm1=abs(taum_est_f(n)-taum_est_f(n-1))/taum_est_f(n-1); 
%            
%         if ~isempty(find(rel_er_1 > 20*er_tol)) || ~isempty(find(rel_er_2 > 20*er_tol)) || r_e_tm1 > 20*er_tol_taum(1)           
%             opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
%             [fitresult, gof] = fit( xData, yData, ft1, opts );
%             t1_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
%     
%             taum_func = @(taum) find_taum_fixed_Tp(Tp1_f(n),taum, t1_est_f(n,3), true_gr, k1, mu, sigma, w);% 
%             opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
%             problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
%                 'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
%                 
%             %[taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
%   
%             pts = min_taum+ (max_taum-min_taum).*rand(200,1);
%             tpoints = CustomStartPointSet(pts);
%             allpts = {tpoints};
%                 
%             [taum_est_f(n),fval_f(n),flagm_ms,outptm_ms,manyminsm_ms] = run(MultiStart,problem,allpts);    
%         end
%         
%         
%         % check bad-fit for the 2nd time
%         
%         rel_er_1= abs((t1_est_f(n,:)-t1_est_f(n-1,:))./t1_est_f(n-1,:));
%         r_e_tm1=abs(taum_est_f(n)-taum_est_f(n-1))/taum_est_f(n-1); 
%            
%         if ~isempty(find(rel_er_1 > 20*er_tol)) || r_e_tm1 > 20*er_tol_taum(1)
%             
%             opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
%             [fitresult, gof] = fit( xData, yData, ft1, opts );
%             t1_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
%     
%     
%             taum_func = @(taum) find_taum_fixed_Tp(Tp1_f(n),taum, t1_est_f(n,3), true_gr, k1, mu, sigma, w);% 
%             opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
%             problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
%                 'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
%                 
%             %[taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
%   
%             pts = min_taum+ (max_taum-min_taum).*rand(200,1);
%             tpoints = CustomStartPointSet(pts);
%             allpts = {tpoints};
%                 
%             [taum_est_f(n),fval_f(n),flagm_ms,outptm_ms,manyminsm_ms] = run(MultiStart,problem,allpts);
%         end
%         
%   end
    
    
    %
    
    
    
    check_stopping_fim;
    
    
    %main_tau_est_fixed_Tp_uni;
    
end


%
plot_figures;

