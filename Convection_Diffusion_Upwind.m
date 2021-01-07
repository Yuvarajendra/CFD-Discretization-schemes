clc; clear; close all;

%% Example 5.2, Case 1 in:
% Versteeg, H.K., Malalasekera, W., 2007. An introduction to computational 
% fuid dynamics: the finite volume method. Pearson Education. pp. 147-148

%% Notes:
% The upwind differencing scheme have been used to discretized the equations
% while the Gauss-Siedel iteration method to solve the the set of algebraic
% equations. 

%% Inputs

N=10;            % Number of nodes
ConvCrit=1e-2;  % Convergence criteria (for the Gauss-Seidel Scheme)
L=1.0;          % Length [m]
dx=L/N;         % Grid size [m]
rho=1.0;        % Density [kg m^-3]
u=0.1;          % Velocity [m s^-1]
u=2.5;          % Velocity [m s^-1]
F=rho*u;        % Convective flux term [kg m^-2 s^-1]
Gamma=0.1;      % Diffusion coefficient [kg m^-1 s^-1]
D=Gamma/dx;     % Diffusion conductance at cell faces [kg m^-2 s^-1]
Pe=F/D;         % Peclet number

disp (['Peclet number = ', num2str(Pe,2)]); % Display Pe number

%% Analytical solution, Case 1 %%

% Ua=1; Ub=0; % Boundary Conditions

N2=100; % Number of nodes analytical solution

distance_ana=zeros(N2+2,1);
distance_ana(1,1)=0; 
distance_ana(end,1)=L;

phi_exact=zeros(N2+2,1);
phi_exact(1,1)=1;
phi_exact(end,1)=0;

% Inner Points:

Xa=L/N2;
dxa=L/N2;

for r=2:N2+1 % loop over a inner points

phi_exact(r,1)=(2.7183-exp(Xa))/(1.7183);
distance_ana(r,1)=Xa;
Xa=Xa+dxa;

end

%%%%  Numerical solver using the FVM %%%%

%% Creating matrix A

% Inner nodes:

Sp=0;
Su=0;
ae=D+max(0,-F); % Note, Fw=Fe=F
aw=D+max(F,0);
ap=aw+ae-Sp;

A=eye(N,N)*ap+diag(ones(1,N-1)*(-aw),-1)+diag(ones(1,N-1)*(-ae),1);

%% Assign boundary conditions at first and last node

% First node:

Sigma_A=1; % at x=0 (boundary condition)

% source terms due to boundary conditions

Sp=-(2*D+F); 
Su_A=(2*D+F)*Sigma_A;

% coefficients: 

aw=0; 
ap=aw+ae-Sp;
A(1,1)=ap; % change in matrix A

% Last node:

Sigma_B=0;
Sp=-(2*D);
Su_B=(2*D)*Sigma_B;
ae=0;
aw=D+F;
ap=aw+ae-Sp;
A(N,N)=ap; % change in matrix A

%% Creating vector b:

b=zeros(N,1);
b(1,1)=Su_A; % Assign source term (such that Eq. 5.34 is correct)

% Note that b(N,N)=Su_B=0

%% Numerical Solution Using the FVM  %%

x0=zeros(N,1); % Initial guess of phi for the internal nodes

% Gauss-Siedel Method for Solving Ax=b: 

[x, residual, numItr] = gauss_seidel(A, b, x0, ConvCrit);

phi=x; % The transported scalar 

distance_num=[dx/2:dx:L-dx/2];

%% Plot data

figure(1);
plot (distance_ana, phi_exact,'-k',distance_num,phi,':sqk','LineWidth',1.5,'MarkerFaceColor','k');
set(gcf,'Units','centimeters');
afFigurePosition = [15 10 10 7.5];       % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); 
set(gca,'xlim',[0 1],'xtick',[0:0.2:1.0],'FontSize',8,'FontWeight','normal');
set(gca,'ylim',[0 1.05],'ytick',[0:0.2:1.0],'FontSize',8,'FontWeight','normal');
set(gcf,'color','w');
xlabel('Distance (m)','Fontsize',10); 
ylabel('$\phi$','interpreter','latex','FontSize',10);
legend('Exact solution','Numerical solution (UD)');
title(['Example 5.2 (Case 1)'],'FontWeight','normal','fontsize',10); 

%% Write data to text file (csv):

T=([distance_num', phi]); % setup output matrix

dlmwrite([pwd,'/file_name.csv'],T,'delimiter',',', 'precision', 6);

% For more details type "help dlmwrite" in the Command Window

