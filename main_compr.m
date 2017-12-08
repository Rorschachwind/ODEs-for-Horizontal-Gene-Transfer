%{ 
------------12/06/2017-----------
-----------main function---------

------------Resistance-----------
         G is sensitive to Cm
         R is resistant to Cm
         Y is resistant to Cm

         G is resistant to Kan
         R is sensitive to Kan
         Y is resistant to Kan
-----Combination of eta_Hill-----
             1: G a R r
             2: G r R a
             3: G a R a
             4: G r R r
%}

%% global term
clear; clc;
global K1 K2 alpha1 alpha2 beta1 beta2 n m Nm mu_G_max mu_R_max mu_Y_max  initConc tspan eta0 A Atype

K1 = 0.01;
K2 = 0.01;
alpha1 = 3.4409e-07;
alpha2 = 0.0017;
beta1 = 25;
beta2 = 25;
n = 2;
m = 2;
Nm = 1E9; % carrying capacity
xLimit = [0 60];
 
% antibiotics type 
Atype = 'none';
A=0.0;
%Dilution rate
D=0.0;
%segregation error
kappa=0;

mu_Y_max = 0.31;
tspan= 0.1:0.1:60;
initConc = [1E-3*Nm; 1E-3*Nm; 0];
etaC_prime = 0.0375;
eta0 = (etaC_prime * mu_Y_max / Nm);
%% set up function
DD=0.0:0.005:0.15;
AA=0:0.05:1;
aa=length(AA);
dd=length(DD);
y0=[];
y1=[];
y2=[];
y3=[];
y4=[];

%this part should be commented if explore the effects of K1(K2)
% mu_G_eff=0:0.01:mu_G_max;
% mu_R_eff=0:0.01:mu_R_max;
% mu_Y_eff=0:0.01:mu_Y_max;
% Hill_GA = @(mu_G_eff) alpha1 + alpha2 * mu_G_eff.^n./(K1^n+mu_G_eff.^n);
% Hill_GR = @(mu_G_eff) alpha1 + alpha2 * K1^n./(K1^n+mu_G_eff.^n);
% Hill_RA = @(mu_R_eff) beta1 + beta2 * mu_R_eff.^m./(K2^m+mu_R_eff.^m);
% Hill_RR = @(mu_R_eff) beta1 + beta2 * K2^m./(K2^m+mu_R_eff.^m);
% Hill_YA=@(mu_Y_eff) alpha1 + alpha2 * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);
% Hill_YR=@(mu_Y_eff) alpha1 + alpha2 * K1.^n./(K1^n+mu_Y_eff.^n);
% 
% q_GA=integral(Hill_GA,0,mu_G_max)/mu_G_max;
% q_GR=integral(Hill_GR,0,mu_G_max)/mu_G_max;
% q_RA=integral(Hill_RA,0,mu_R_max)/mu_R_max;
% q_RR=integral(Hill_RR,0,mu_R_max)/mu_R_max;
% q_YA=integral(Hill_YA,0,mu_Y_max)/mu_Y_max;
% q_YR=integral(Hill_YR,0,mu_Y_max)/mu_Y_max;

%% heat map for a pair of parameters
% for i=1:aa
%     A=AA(i);
%     for j = 1:dd
%         
%         D=DD(j);
% 
%         %constant etaC with dilution D
% %         Version = 0;
% %         [tv,Fv0]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version),tspan,initConc);
%     %     Y0_1=Fv0(:,1);
%     %     Y0_2=Fv0(:,2);
%     %     Y0_3=Fv0(:,3);
%     %     y0(i)=(Y0_3(jj))/(Y0_1(jj)+Y0_3(jj)+Y0_2(jj));
% 
%         % G a R r
%         Version = 1;
%         [tv,Fv1]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version,q_GA,q_RR,q_YA),tspan, initConc);
%         [mu_eff_update1,etaGR1,etaYR1]= calcE (Fv1,Version,q_GA,q_RR,q_YA);
%         jj=length(Fv1(:,1));
%         Y1_1=Fv1(:,1);
%         Y1_2=Fv1(:,2);
%         Y1_3=Fv1(:,3);
%         y1(i,j)=Y1_3(jj)/(Y1_1(jj)+Y1_2(jj)+Y1_3(jj));
% 
%         %G r R a
%         Version = 2;
%         [tv,Fv2]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version,q_GA,q_RR,q_YA),tspan, initConc);
%         [mu_eff_update2,etaGR2,etaYR2]= calcE (Fv2,Version,q_GA,q_RR,q_YA);
%         Y2_1=Fv2(:,1);
%         Y2_2=Fv2(:,2);
%         Y2_3=Fv2(:,3);
%         y2(i,j)=Y2_3(jj)/(Y2_1(jj)+Y2_2(jj)+Y2_3(jj));

%         %G a R a
%         Version = 3;
%         [tv,Fv3]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version,q_GA,q_RR,q_YA),tspan, initConc);
%         [mu_eff_update3,etaGR3,etaYR3]= calcE (Fv3,Version,q_GA,q_RR,q_YA);
%         Y3_1=Fv3(:,1);
%         Y3_2=Fv3(:,2);
%         Y3_3=Fv3(:,3);
%         y3(i,j)=Y3_3(jj)/(Y3_1(jj)+Y3_2(jj)+Y3_3(jj));


%         %G r R r
%         Version = 4;
%         [tv,Fv4]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version,q_GA,q_RR,q_YA),tspan, initConc);
%         [mu_eff_update4,etaGR4,etaYR4]= calcE (Fv4,Version,q_GA,q_RR,q_YA);
%         Y4_1=Fv4(:,1);
%         Y4_2=Fv4(:,2);
%         Y4_3=Fv4(:,3);
%         y4(i,j)=Y4_3(jj)/(Y4_1(jj)+Y4_2(jj)+Y4_3(jj));

%     end
% end
% figure;
% heatmap(AA,DD,y1);
% xlabel('[Cm]');
% ylabel('Dilution rate/Hr');
% title('Fraction of transconjugants');
% figure;
% heatmap(AA,DD,y2);
% xlabel('[Cm]');
% ylabel('Dilution rate/Hr');
% title('Fraction of transconjugants');
% figure;
% heatmap(AA,DD,y3);
% xlabel('[Cm]');
% ylabel('Dilution rate/Hr');
% title('Fraction of transconjugants');
% figure;
% heatmap(AA,DD,y4);
% xlabel('[Cm]');
% ylabel('Dilution rate/Hr');
% title('Fraction of transconjugants');

%% this part can generate fraction of plasmid-carrying cell or transconjugants versus different parameters
%initial parameters set
DD=0.0:0.005:0.3;
AA=0:0.02:3;
aa=length(AA);
dd=length(DD);
KK=0:0.01:1;
kk=length(KK);
KAPA=0:0.001:0.1;
pp=length(KAPA);

for i= 1:kk
%this part can be changed for paramter scanning
%       D=DD(i);
%       A=AA(i);
%       kappa=KAPA(i);
    K1=KK(i);K2=KK(i);
    mu_G_eff=0:0.01:mu_G_max;
    mu_R_eff=0:0.01:mu_R_max;
    mu_Y_eff=0:0.01:mu_Y_max;
    Hill_GA = @(mu_G_eff) alpha1 + alpha2 * mu_G_eff.^n./(K1^n+mu_G_eff.^n);
    Hill_GR = @(mu_G_eff) alpha1 + alpha2 * K1^n./(K1^n+mu_G_eff.^n);
    Hill_RA = @(mu_R_eff) beta1 + beta2 * mu_R_eff.^m./(K2^m+mu_R_eff.^m);
    Hill_RR = @(mu_R_eff) beta1 + beta2 * K2^m./(K2^m+mu_R_eff.^m);
    Hill_YA=@(mu_Y_eff) alpha1 + alpha2 * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);
    Hill_YR=@(mu_Y_eff) alpha1 + alpha2 * K1.^n./(K1^n+mu_Y_eff.^n);

    q_GA=integral(Hill_GA,0,mu_G_max)/mu_G_max;
    q_GR=integral(Hill_GR,0,mu_G_max)/mu_G_max;
    q_RA=integral(Hill_RA,0,mu_R_max)/mu_R_max;
    q_RR=integral(Hill_RR,0,mu_R_max)/mu_R_max;
    q_YA=integral(Hill_YA,0,mu_Y_max)/mu_Y_max;
    q_YR=integral(Hill_YR,0,mu_Y_max)/mu_Y_max;
    
    %constant etaC with dilution D
    Version = 0;
    [tv,Fv0]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,[],[],[]),tspan,initConc);
     jj=length(tv);
    Y0_1=Fv0(:,1);
    Y0_2=Fv0(:,2);
    Y0_3=Fv0(:,3);
    y0(i)=(Y0_3(jj)+Y0_1(jj))/(Y0_1(jj)+Y0_3(jj)+Y0_2(jj));
    
    % G a R r
    Version = 1;
    [tv,Fv1]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
    [mu_eff_update1,etaGR1,etaYR1]= calcE (Fv1,Version,q_GA,q_RR,q_YA);
    Y1_1=Fv1(:,1);
    Y1_2=Fv1(:,2);
    Y1_3=Fv1(:,3);
    y1(i)=(Y1_3(jj)+Y1_1(jj))/(Y1_1(jj)+Y1_2(jj)+Y1_3(jj));

    %G r R a
    Version = 2;
    [tv,Fv2]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
    [mu_eff_update2,etaGR2,etaYR2]= calcE (Fv2,Version,q_GA,q_RR,q_YA);
    Y2_1=Fv2(:,1);
    Y2_2=Fv2(:,2);
    Y2_3=Fv2(:,3);
    y2(i)=(Y2_1(jj)+Y2_3(jj))/(Y2_1(jj)+Y2_2(jj)+Y2_3(jj));

    %G a R a
    Version = 3;
    [tv,Fv3]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
    [mu_eff_update3,etaGR3,etaYR3]= calcE (Fv3,Version,q_GA,q_RR,q_YA);
    Y3_1=Fv3(:,1);
    Y3_2=Fv3(:,2);
    Y3_3=Fv3(:,3);
    y3(i)=(Y3_1(jj)+Y3_3(jj))/(Y3_1(jj)+Y3_2(jj)+Y3_3(jj));
    
    %G r R r
    Version = 4;
    [tv,Fv4]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
    [mu_eff_update4,etaGR4,etaYR4]= calcE (Fv4,Version,q_GA,q_RR,q_YA);
    Y4_1=Fv4(:,1);
    Y4_2=Fv4(:,2);
    Y4_3=Fv4(:,3);
    y4(i)=(Y4_1(jj)+Y4_3(jj))/(Y4_1(jj)+Y4_2(jj)+Y4_3(jj));
    
    % eta_total
    xx(:,i) = etaYR1+etaGR1;
    yy(:,i) =  etaYR2+etaGR2;
    zz(:,i) = etaYR3+etaGR3;
    ww(:,i) =  etaYR4+etaGR4;
end
%% plot All figures
% Popfig(tv,FvS,Fv0,Fv1,Fv2,Fv3,Fv4,xx,yy,zz,ww)
%% dilution compare figure
% figure;
% YYY=[y0;y1;y2;y3;y4];
% Name=["G: activation | R: repression","G: repression | R: activation","G: activation | R: activation","G: repression | R: repression"];
% for i = 1 :4
%     subplot(1,4,i);
%     plot(DD,YYY(i+1,:),DD,YYY(1,:),'--','LineWidth',2);
%      ylim([0,0.25]);set(gca,'LineWidth',2,'Fontsize',25);
%      xlabel('Dilution rate/Hr','Fontsize',25);
%      ylabel('Fraction of Y cells','Fontsize',25);
%      title(Name(i),'Fontsize',25);
%      
% end
% subplot(1,4,4);
% legend('\eta(t)','\eta_0');
% H=legend('\eta(t)','\eta_0');
% set(H,'Fontsize',30);
% suptitle('Fraction of transconjugants under different dilution rate');
% h=suptitle('Fraction of transconjugants under different dilution rate');
% set(h,'Fontsize',25);
%% cm antibiotics compare figure
% figure;
% YYY=[y0;y1;y2;y3;y4];
% Name=["G: activation | R: repression","G: repression | R: activation","G: activation | R: activation","G: repression | R: repression"];
% for i = 1 :4
%     subplot(1,4,i);
%     plot(AA,YYY(i+1,:),AA,YYY(1,:),'--','LineWidth',2);
%       ylim([0,0.2]);
%      set(gca,'LineWidth',2,'Fontsize',25);
%      xlabel('[Kan]','Fontsize',25);
%      ylabel('Fraction of Y cells','Fontsize',25);
%      title(Name(i),'Fontsize',25);
% end
% subplot(1,4,4);
% legend('\eta(t)','\eta_0');
% H=legend('\eta(t)','\eta_0');
% set(H,'Fontsize',30);
% suptitle('Fraction of transconjugants under different [Kan]');
% h=suptitle('Fraction of transconjugants under different [Kan]');
% set(h,'Fontsize',25);
%% K1(K2) compare figure
figure;
YYY=[y0;y1;y2;y3;y4];
Name=["G: activation | R: repression","G: repression | R: activation","G: activation | R: activation","G: repression | R: repression"];
for i = 1 :4
    subplot(1,4,i);
    plot(KK,YYY(i+1,:),KK,YYY(1,:),'--','LineWidth',2);
        ylim([0.5,1]);
     set(gca,'LineWidth',2,'Fontsize',25);
     xlabel('Kappa','Fontsize',25);
     ylabel('Fraction of G+Y cells','Fontsize',25);
     title(Name(i),'Fontsize',25);
end
subplot(1,4,4);
legend('\eta(t)','\eta_0');
H=legend('\eta(t)','\eta_0');
set(H,'Fontsize',30);
suptitle('Fraction of plasmid-containing cells under different K1(K2)');
h=suptitle('Fraction of plasmid-containing cells under different K1(K2)');
set(h,'Fontsize',25);
