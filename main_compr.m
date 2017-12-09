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

% K1 = 0.01;
% K2 = 0.01;
alpha1 = 3.4409e-07;
alpha2 = 0.0017;
beta1 = 25;
beta2 = 25;
n = 2;
m = 2;
Nm = 1E9; % carrying capacity
xLimit = [0 60];
 
% antibiotics type 
Atype = 'cm';
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
DD=0.0:0.005:0.2;
AA=0.0:0.05:2.0;
KK=0.0:0.005:0.15;
kk=length(KK);
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
for i=1:aa
    A=AA(i);
%     D=DD(i);
    for j = 1:kk
%          A=AA(j);
        K1=KK(j);
        K2=KK(j);
         
        mu_G_eff=0:0.01:mu_G_max;
        mu_R_eff=0:0.01:mu_R_max;
        mu_Y_eff=0:0.01:mu_Y_max;
        Hill_GA = @(mu_G_eff) alpha1 + alpha2 * mu_G_eff.^n./(K1^n+mu_G_eff.^n);
        Hill_GR = @(mu_G_eff) alpha1 + alpha2 * K1^n./(K1^n+mu_G_eff.^n);
        Hill_RA = @(mu_R_eff) beta1 + beta2 * mu_R_eff.^m./(K2^m+mu_R_eff.^m);
        Hill_RR = @(mu_R_eff) beta1 + beta2 * K2^m./(K2^m+mu_R_eff.^m);
        Hill_YA=@(mu_Y_eff) alpha1 + alpha2 * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);
        Hill_YR=@(mu_Y_eff) alpha1 + alpha2 * K1^n./(K1^n+mu_Y_eff.^n);

        q_GA=integral(Hill_GA,0,mu_G_max)/mu_G_max;
        q_GR=integral(Hill_GR,0,mu_G_max)/mu_G_max;
        q_RA=integral(Hill_RA,0,mu_R_max)/mu_R_max;
        q_RR=integral(Hill_RR,0,mu_R_max)/mu_R_max;
        q_YA=integral(Hill_YA,0,mu_Y_max)/mu_Y_max;
        q_YR=integral(Hill_YR,0,mu_Y_max)/mu_Y_max;
        
        %constant etaC with dilution D
        Version = 0;
        [tv,Fv0]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,[],[],[]),tspan,initConc);
        jj=length(Fv0(:,1));
        Y0_1=Fv0(:,1);
        Y0_2=Fv0(:,2);
        Y0_3=Fv0(:,3);
        y0(i,j)=(Y0_3(jj)+Y0_1(jj))/(Y0_1(jj)+Y0_3(jj)+Y0_2(jj));
        y00(i,j)=(Y0_3(jj))/(Y0_1(jj)+Y0_3(jj)+Y0_2(jj));
% 
        % G a R r
        Version = 1;
        [tv,Fv1]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
        [mu_eff_update1,etaGR1,etaYR1]= calcE (Fv1,Version,q_GA,q_RR,q_YA);
        jj=length(Fv1(:,1));
        Y1_1=Fv1(:,1);
        Y1_2=Fv1(:,2);
        Y1_3=Fv1(:,3);
        y1(i,j)=(Y1_1(jj)+Y1_3(jj))/(Y1_1(jj)+Y1_2(jj)+Y1_3(jj));
        y11(i,j)=(Y1_3(jj))/(Y1_1(jj)+Y1_2(jj)+Y1_3(jj));
        %G r R a
        Version = 2;
        [tv,Fv2]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
        [mu_eff_update2,etaGR2,etaYR2]= calcE (Fv2,Version,q_GA,q_RR,q_YA);
        Y2_1=Fv2(:,1);
        Y2_2=Fv2(:,2);
        Y2_3=Fv2(:,3);
        y2(i,j)=(Y2_1(jj)+Y2_3(jj))/(Y2_1(jj)+Y2_2(jj)+Y2_3(jj));
        y22(i,j)=(Y2_3(jj))/(Y2_1(jj)+Y2_2(jj)+Y2_3(jj));

        %G a R a
        Version = 3;
        [tv,Fv3]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
        [mu_eff_update3,etaGR3,etaYR3]= calcE (Fv3,Version,q_GA,q_RR,q_YA);
        Y3_1=Fv3(:,1);
        Y3_2=Fv3(:,2);
        Y3_3=Fv3(:,3);
        y3(i,j)=(Y3_1(jj)+Y3_3(jj))/(Y3_1(jj)+Y3_2(jj)+Y3_3(jj));
        y33(i,j)=(Y3_3(jj))/(Y3_1(jj)+Y3_2(jj)+Y3_3(jj));


        %G r R r
        Version = 4;
        [tv,Fv4]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
        [mu_eff_update4,etaGR4,etaYR4]= calcE (Fv4,Version,q_GA,q_RR,q_YA);
        Y4_1=Fv4(:,1);
        Y4_2=Fv4(:,2);
        Y4_3=Fv4(:,3);
        y4(i,j)=(Y4_1(jj)+Y4_3(jj))/(Y4_1(jj)+Y4_2(jj)+Y4_3(jj));
        y44(i,j)=(Y4_3(jj))/(Y4_1(jj)+Y4_2(jj)+Y4_3(jj));

    end
end

%% heat map used to find maximum fraction of cells
figure;
heatmap(KK,AA,y1);
xlabel('KK');
% xlabel('[Cm]');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of GY');
figure;
heatmap(KK,AA,y2);
xlabel('KK');
% xlabel('[Cm]');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of GY');
figure;
heatmap(KK,AA,y3);
xlabel('KK');
% xlabel('[Cm]');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of GY');
figure;
heatmap(KK,AA,y4);
xlabel('KK');
% xlabel('[Cm]');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of GY');

figure;
heatmap(KK,AA,y11);
xlabel('KK');
% xlabel('[Cm]');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of Y');
figure;
heatmap(KK,AA,y22);
xlabel('KK');
% xlabel('[Cm]');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of Y');
figure;
heatmap(KK,AA,y33);
xlabel('KK');
% xlabel('[Cm]');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of Y');
figure;
heatmap(KK,AA,y44);
xlabel('KK');
% xlabel('[Cm]');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of Y');

figure;
heatmap(KK,AA,y0);
xlabel('KK');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of GY');
figure;
heatmap(KK,AA,y00);
xlabel('KK');
ylabel('[Cm]');
% ylabel('Dilution rate/Hr');
title('Fraction of Y');
%% imagesc figures comparison
figure;
subplot(2,2,1);
imagesc(AA,DD,y1);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel('[Kan]','Fontsize',20);
ylabel('Dilution rate/Hr','Fontsize',20);
% ylabel('[Cm]','Fontsize',20);
title('GA | RR','Fontsize',20);
% colorbar;
%title('Fraction of transconjugants');
subplot(2,2,2);
imagesc(AA,DD,y2);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel('[Kan]','Fontsize',20);
ylabel('Dilution rate/Hr','Fontsize',20);
% ylabel('[Cm]','Fontsize',20);
title('GR | RA','Fontsize',20);
% colorbar;
% title('Fraction of transconjugants');
subplot(2,2,3);
imagesc(AA,DD,y3);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel('[Kan]','Fontsize',20);
ylabel('Dilution rate/Hr','Fontsize',20);
% ylabel('[Cm]','Fontsize',20);
title('GA | RA','Fontsize',20);
% colorbar;
% title('Fraction of transconjugants');
subplot(2,2,4);
imagesc(AA,DD,y4);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel('[Kan]','Fontsize',20);
ylabel('Dilution rate/Hr','Fontsize',20);
% ylabel('[Cm]','Fontsize',20);
title('GR | RR','Fontsize',20);
colorbar;
h=suptitle('Franction of Plasmid-carrying Cells');
set(h,'Fontsize',25);

% title('Fraction of transconjugants');

figure;
subplot(2,2,1);
imagesc(AA,DD,y11);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel('[Kan]','Fontsize',20);
ylabel('Dilution rate/Hr','Fontsize',20);
% ylabel('[Cm]','Fontsize',20);
title('GA | RR','Fontsize',20);
% colorbar;
%title('Fraction of transconjugants');
subplot(2,2,2);
imagesc(AA,DD,y22);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel('[Kan]','Fontsize',20);
ylabel('Dilution rate/Hr','Fontsize',20);
% ylabel('[Cm]','Fontsize',20);
title('GR | RA','Fontsize',20);
% colorbar;
% title('Fraction of transconjugants');
subplot(2,2,3);
imagesc(AA,DD,y33);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel('[Kan]','Fontsize',20);
ylabel('Dilution rate/Hr','Fontsize',20);
% ylabel('[Cm]','Fontsize',20);
title('GA | RA','Fontsize',20);
% colorbar;
% title('Fraction of transconjugants');
subplot(2,2,4);
imagesc(AA,DD,y44);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel('[Kan]','Fontsize',20);
ylabel('Dilution rate/Hr','Fontsize',20);
% ylabel('[Cm]','Fontsize',20);
title('GR | RR','Fontsize',20);
colorbar;
h=suptitle('Franction of Transconjugants');
set(h,'Fontsize',25);


%% this part can generate fraction of plasmid-carrying cell or transconjugants versus different parameters
%initial parameters set
% DD=0.0:0.005:0.3;
% DD=[0.0 0.1 0.2 0.3];
% AA=0:0.5:2;
% aa=length(AA);
% dd=length(DD);
% KK=0.0:0.005:0.15;
% kk=length(KK);
% KAPA=0:0.001:0.1;
% pp=length(KAPA);

% for i= 1:kk
% %this part can be changed for paramter scanning
%     
% %        D=DD(i);
% %       A=AA(i);
% %       kappa=KAPA(i);
%      K1=KK(i);K2=KK(i);
%     mu_G_eff=0:0.01:mu_G_max;
%     mu_R_eff=0:0.01:mu_R_max;
%     mu_Y_eff=0:0.01:mu_Y_max;
%     Hill_GA = @(mu_G_eff) alpha1 + alpha2 * mu_G_eff.^n./(K1^n+mu_G_eff.^n);
%     Hill_GR = @(mu_G_eff) alpha1 + alpha2 * K1^n./(K1^n+mu_G_eff.^n);
%     Hill_RA = @(mu_R_eff) beta1 + beta2 * mu_R_eff.^m./(K2^m+mu_R_eff.^m);
%     Hill_RR = @(mu_R_eff) beta1 + beta2 * K2^m./(K2^m+mu_R_eff.^m);
%     Hill_YA=@(mu_Y_eff) alpha1 + alpha2 * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);
%     Hill_YR=@(mu_Y_eff) alpha1 + alpha2 * K1.^n./(K1^n+mu_Y_eff.^n);
% 
%     q_GA=integral(Hill_GA,0,mu_G_max)/mu_G_max;
%     q_GR=integral(Hill_GR,0,mu_G_max)/mu_G_max;
%     q_RA=integral(Hill_RA,0,mu_R_max)/mu_R_max;
%     q_RR=integral(Hill_RR,0,mu_R_max)/mu_R_max;
%     q_YA=integral(Hill_YA,0,mu_Y_max)/mu_Y_max;
%     q_YR=integral(Hill_YR,0,mu_Y_max)/mu_Y_max;
%     
% %     constant etaC with dilution D
%     Version = 0;
%     [tv,Fv0]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,[],[],[]),tspan,initConc);
%      jj=length(tv);
%     Y0_1=Fv0(:,1);
%     Y0_2=Fv0(:,2);
%     Y0_3=Fv0(:,3);
%     y0(i)=(Y0_3(jj)+Y0_1(jj))/(Y0_1(jj)+Y0_3(jj)+Y0_2(jj));
%     y00(i)=(Y0_3(jj))/(Y0_1(jj)+Y0_3(jj)+Y0_2(jj));
%     
%     % G a R r
%     Version = 1;
%     [tv,Fv1]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
%     [mu_eff_update1,etaGR1,etaYR1]= calcE (Fv1,Version,q_GA,q_RR,q_YA);
%     Y1_1=Fv1(:,1);
%     Y1_2=Fv1(:,2);
%     Y1_3=Fv1(:,3);
%     y1(i)=(Y1_3(jj)+Y1_1(jj))/(Y1_1(jj)+Y1_2(jj)+Y1_3(jj));
%     y11(i)=(Y1_3(jj))/(Y1_1(jj)+Y1_2(jj)+Y1_3(jj));
%     
%     %G r R a
%     Version = 2;
%     [tv,Fv2]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
%     [mu_eff_update2,etaGR2,etaYR2]= calcE (Fv2,Version,q_GA,q_RR,q_YA);
%     Y2_1=Fv2(:,1);
%     Y2_2=Fv2(:,2);
%     Y2_3=Fv2(:,3);
%     y2(i)=(Y2_1(jj)+Y2_3(jj))/(Y2_1(jj)+Y2_2(jj)+Y2_3(jj));
%     y22(i)=(Y2_3(jj))/(Y2_1(jj)+Y2_2(jj)+Y2_3(jj));
%     
%     %G a R a
%     Version = 3;
%     [tv,Fv3]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
%     [mu_eff_update3,etaGR3,etaYR3]= calcE (Fv3,Version,q_GA,q_RR,q_YA);
%     Y3_1=Fv3(:,1);
%     Y3_2=Fv3(:,2);
%     Y3_3=Fv3(:,3);
%     y3(i)=(Y3_1(jj)+Y3_3(jj))/(Y3_1(jj)+Y3_2(jj)+Y3_3(jj));
%     y33(i)=(Y3_3(jj))/(Y3_1(jj)+Y3_2(jj)+Y3_3(jj));
%     
%     %G r R r
%     Version = 4;
%     [tv,Fv4]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
%     [mu_eff_update4,etaGR4,etaYR4]= calcE (Fv4,Version,q_GA,q_RR,q_YA);
%     Y4_1=Fv4(:,1);
%     Y4_2=Fv4(:,2);
%     Y4_3=Fv4(:,3);
%     y4(i)=(Y4_1(jj)+Y4_3(jj))/(Y4_1(jj)+Y4_2(jj)+Y4_3(jj));
%     y44(i)=(Y4_3(jj))/(Y4_1(jj)+Y4_2(jj)+Y4_3(jj));
%     
%     % eta_total
%     xx(:,i) = etaYR1+etaGR1;
%     yy(:,i) =  etaYR2+etaGR2;
%     zz(:,i) = etaYR3+etaGR3;
%     ww(:,i) =  etaYR4+etaGR4;
% end
% %% plot All figures
% % Popfig(tv,FvS,Fv0,Fv1,Fv2,Fv3,Fv4,xx,yy,zz,ww)
%% dilution compare figure
% figure;
% YYY=[y0;y1;y2;y3;y4];
% Name=["G: activation | R: repression","G: repression | R: activation","G: activation | R: activation","G: repression | R: repression"];
% 
% for i = 1 :4
%     subplot(1,4,i);
%     plot(DD,YYY(i+1,:),DD,YYY(1,:),'--','LineWidth',2);
%      ylim([0.5,1.0]);
%      set(gca,'LineWidth',2,'Fontsize',25);
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


%plot a bar for plasmid-carrying cells 
% figure;
% % x=0:0.1:0.3;
% Y=[y0(1), y1(1),y2(1),y3(1),y4(1);y0(2), y1(2),y2(2),y3(2),y4(2);y0(3), y1(3),y2(3),y3(3),y4(3);y0(4), y1(4),y2(4),y3(4),y4(4);y0(5), y1(5),y2(5),y3(5),y4(5);y0(6), y1(6),y2(6),y3(6),y4(6)];
% b=bar(KK,Y,1,'LineWidth',2);
% % set(b,'Facecolor',
% ylim([0 1.0]);
% set(gca, 'XTick',[0 0.01 0.02 0.03 0.04 0.05]); 
% set(gca,'LineWidth',2,'Fontsize',18);
% xlabel('K1(K2)','Fontsize',20);
% ylabel('Fraction of Plasmid-carrying Cells','Fontsize',20);
% title('Ranking of 4 Scenarios under Different K1(K2)');
% legend({'\eta_0','GA|RR','GR|RA','GA|RA','GR|RR'},'Fontsize',18);
% 
% %plot bar for Y cell
% YYY=[y00;y11;y22;y33;y44];
% figure;
% Y=[y00(1), y11(1),y22(1),y33(1),y44(1);y00(2), y11(2),y22(2),y33(2),y44(2);y00(3), y11(3),y22(3),y33(3),y44(3);y00(4), y11(4),y22(4),y33(4),y44(4);y00(5), y11(5),y22(5),y33(5),y44(5);y00(6), y11(6),y22(6),y33(6),y44(6)];
% b=bar(KK,Y,1,'LineWidth',2);
% % set(b,'Facecolor',
% % ylim([0 0.2]);
% set(gca, 'XTick',[0 0.01 0.02 0.03 0.04 0.05]); 
% set(gca,'LineWidth',2,'Fontsize',18);
% xlabel('K1(K2)','Fontsize',20);
% ylabel('Fraction of Transconjugants','Fontsize',20);
% title('Ranking of 4 Scenarios under Different K1(K2)');
% legend({'\eta_0','GA|RR','GR|RA','GA|RA','GR|RR'},'Fontsize',18);

%% cm antibiotics compare figure
% % figure;
% % YYY=[y0;y1;y2;y3;y4];
% % Name=["G: activation | R: repression","G: repression | R: activation","G: activation | R: activation","G: repression | R: repression"];
% % for i = 1 :4
% %     subplot(1,4,i);
% %     plot(AA,YYY(i+1,:),AA,YYY(1,:),'--','LineWidth',2);
% %       ylim([0,0.2]);
% %      set(gca,'LineWidth',2,'Fontsize',25);
% %      xlabel('[Kan]','Fontsize',25);
% %      ylabel('Fraction of Y cells','Fontsize',25);
% %      title(Name(i),'Fontsize',25);
% % end
% % subplot(1,4,4);
% % legend('\eta(t)','\eta_0');
% % H=legend('\eta(t)','\eta_0');
% % set(H,'Fontsize',30);
% % suptitle('Fraction of transconjugants under different [Kan]');
% % h=suptitle('Fraction of transconjugants under different [Kan]');
% % set(h,'Fontsize',25);
%% K1(K2) compare figure
% figure;
% YYY=[y0;y1;y2;y3;y4];
% Name=["G: activation | R: repression","G: repression | R: activation","G: activation | R: activation","G: repression | R: repression"];
% for i = 1 :4
%     subplot(1,4,i);
%     plot(KK,YYY(i+1,:),KK,YYY(1,:),'--','LineWidth',2);
%         ylim([0.5,1]);
%      set(gca,'LineWidth',2,'Fontsize',25);
%      xlabel('Kappa','Fontsize',25);
%      ylabel('Fraction of G+Y cells','Fontsize',25);
%      title(Name(i),'Fontsize',25);
% end
% subplot(1,4,4);
% legend('\eta(t)','\eta_0');
% H=legend('\eta(t)','\eta_0');
% set(H,'Fontsize',30);
% suptitle('Fraction of plasmid-containing cells under different K1(K2)');
% h=suptitle('Fraction of plasmid-containing cells under different K1(K2)');
% set(h,'Fontsize',25);
