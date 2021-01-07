clc;
clear all; 
close all; 

%% Example 4.1 in:
% Versteeg, H.K., Malalasekera, W., 2007. An introduction to computational 
% fuid dynamics: the finite volume method. Pearson Education. pp. 118-121

%% Notes:
% The central differencing scheme have been used to discretized the equations
% while the Gauss-Siedel iteration method to solve the the set of algebraic
% equations. 

%% Inputs

N=5;             % Number of nodes
k=1000;          % Thermal conducitivity (diffusion coefficient)
ConvCrit=1e-6;   % Convergence criteria (for the Gauss-Seidel Scheme)
Area=10*10^(-3); % Cross sectional area [m^2]
j=1:N;           % Vector index
L=0.5;           % Length [m]
dx=L/N;          % Grid size [m]

%% Analytical solution %%

Xcal=0.5*dx:dx:0.5; 
X(1)=0; X(2:N+1)=Xcal; X(N+2)=0.5; % Distance vector
T_exact=800*X+100;                 % Analytical Solution

%%  Numerical solver using the FVM %%

%% Creating matrix A

Ta=100; Tb=500; % Boundary Conditions

% Inner node coefficients:

aw=k*Area/dx;
ae=k*Area/dx;
ap=aw+ae;

A=eye(N,N)*ap-diag(ones(1,N-1)*(aw),-1)-diag(ones(1,N-1)*(ae),1);

%% Boundary nodes

% First node:

Sp=-(2*k*Area)/dx;
Su_A=Ta*(2*k*Area)/dx;
ap=ae-Sp;
A(1,1)=ap;   % change in matrix A

% Last node:

Su_B=Tb*(2*k*Area)/dx;
ap=aw-Sp;
A(end,end)=ap;   % change in matrix A

%% Creating vector b:

b=zeros(N,1);
b(1,1)=Su_A;   % Assign source term (such that Eq. 4.18 is correct)
b(end,1)=Su_B; % Assign source term (such that Eq. 4.21 is correct)

%% Numerical Solution Using the FVM  %%

x0=zeros(N,1); % Initial guess of phi for the internal nodes

% Gauss-Siedel Method for Solving Ax=b: 

[x, residual, numItr] = gauss_seidel(A, b, x0, ConvCrit);

T=x;  % Temperature field

distance_num=[dx/2:dx:L-dx/2];

%% Plot data

figure(1);
plot (X, T_exact,'-k',distance_num,T,':sqk','LineWidth',1.5,'MarkerFaceColor','k');
set(gcf,'Units','centimeters');
afFigurePosition = [15 10 10 7.5];       % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); 
set(gca,'xlim',[0 0.5],'xtick',[0:0.1:0.5],'FontSize',8,'FontWeight','normal');
set(gca,'ylim',[0 600],'ytick',[0:100:600],'FontSize',8,'FontWeight','normal');
set(gcf,'color','w');
xlabel('Distance (m)','Fontsize',10); 
ylabel('Temperature (^oC)','Fontsize',10); 
legend('Exact solution','Numerical solution (UD)','Location','northwest');
title(['Example 4.1'],'FontWeight','normal','fontsize',10); 
