clc;
clear all;
close all;

%% Example 5.2, Case 1 in:
% Versteeg, H.K., Malalasekera, W., 2007. An introduction to computational 
% fuid dynamics: the finite volume method. Pearson Education. pp. 147-148

%% Notes:
% The upwind differencing scheme have been used to discretized the equations
% while the Gauss-Siedel iteration method to solve the the set of algebraic
% equations. 

%% Inputs
tic

N=5;            % Number of nodes
ConvCrit=10e-6;  % Convergence criteria (for the Gauss-Seidel Scheme)
L=1.0;          % Length [m]
dx=L/N;         % Grid size [m]
rho=1.0;        % Density [kg m^-3]
u=0.1;          % Velocity [m s^-1]
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

phi_exact(r,1)=(2.7183-exp(Xa))/(1.7183);           %%for u=0.1m/s
%phi_exact(r,1)=1+((1-exp(25*Xa))/(7.20*10^10));      %%for u=2.5m/s
distance_ana(r,1)=Xa;
Xa=Xa+dxa;

end

%%%%  Numerical solver using the FVM %%%%
Fe=F;
Fw=F;

if F>0
    alpha=1;
else
    alpha=0;
end

%% Creating matrix A

% Inner nodes:
Sp=0;
Su=0;
aww=-(1/8)*alpha*F;
aw= D + ((6/8)*alpha*F) + ((1/8)*alpha*F) + ((3/8)*(1-alpha)*F);
ae=D - ((3/8)*alpha*F) - ((6/8)*(1-alpha)*F) - ((1/8)*(1-alpha)*F);
aee=((1/8)*(1-alpha)*F);
ap=aww+aw+ae+aee+(Fe-Fw)-Sp;

A=eye(N,N)*ap+diag(ones(1,N-1)*(-aw),-1)+diag(ones(1,N-2)*(-aww),-2)+diag(ones(1,N-1)*(-ae),1);

%% Assign boundary conditions at first and last node

% First node:

Sigma_A=1; % at x=0 (boundary condition)

% source terms due to boundary conditions

Sp=-((8/3)*D + (2/8)*F + F); 
Su_A=((8/3)*D + (2/8)*F + F)*Sigma_A;

% coefficients: 

aw=0;
aww=0;
ae= D + ((1/3)*D) - ((3/8)*F);
aee=0;
ap=aww+aw+ae+aee+(Fe-Fw)-Sp;
A(1,1)=ap; % change in matrix A
A(1,2)=-ae;

% Second node:
Sp=(1/4)*F;
Su_B=-(1/4*F)*Sigma_A;
aww=0;
aw= D + ((7/8)*F) + ((1/8)*F);
ae=D-((3/8)*F);
aee=0;
ap=aww+aw+ae+aee+(Fe-Fw)-Sp;

A(2,2)=ap;
A(2,1)=-aw;
% Last node:

Sigma_B=0;
Sp=-((8/3)*D - F);
Su_C=((8/3)*D - F)*Sigma_B;
aww=-(1/8)*F;
aw= D + ((1/3)*D) + ((6/8)*F);
ae=0;
aee=0;
ap=aww+aw+ae+aee+(Fe-Fw)-Sp;
A(N,N)=ap; % change in matrix A
A(N,N-1)=-aw

%% Creating vector b:

b=zeros(N,1);
b(1,1)=Su_A;     % Assign source term 
b(2,1)=Su_B;     % Assign source term (such that Eq. 5.34 is correct)

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
title(['Example 5.2 (Case 1) Mesh 4N*4N'],'FontWeight','normal','fontsize',10); 

%% Write data to text file (csv):

T=([distance_num', phi]); % setup output matrix

dlmwrite([pwd,'/file_name.csv'],T,'delimiter',',', 'precision', 6);

% For more details type "help dlmwrite" in the Command Window
disp(A)
toc