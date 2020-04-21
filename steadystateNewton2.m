%SteadystateNewton2

% NLS_newton1D_complex.m:
% newton method to caluclate a complex steady stae method for the
% particular PDE iu_t = -1/2 u_xx + g*|u|^2 u
% Code copied from other resources mainly.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%START

%Definitions:

mu = 1; %time frequency of steady state
g = -1; % g=-1,1 are important focusing cases we consider for the PDE
L=-10; R=10; %left and right bounds of the interval
Nx=201; %number of discrete points in 1D mesh (x-direction)
x=linspace(L,R,Nx); %discrete space we solve the steady stae over
dx=x(2)-x(1); % mesh size
ONE=ones(Nx,1);
D2=spdiags([ONE,-2*ONE,ONE],-1:1,Nx,Nx); %discrete laplacian matrix (1D)
D2(1,Nx)=1; D2(Nx,1)=1; %adding in periodic BCs
D2=D2/(dx^2); %suitable normalising

% Defining Perturbation:
if(exist('V')==0) V=0;end; % potential is zero if not already defined
A=sqrt(2*mu); %amplitude of initial guess
if(exist('offset')==0) offset=0;end;
u0=A*sech(A*(x-offset-(R+L)/2)); % unpeturbed initial guess in centre of domain
u0=u0';
rng('default');
pert=1*0.5; %size of perturbation
up=u0+pert*(rand(Nx,1)-0.5); %peturbed IC. rand(a,b) returns a a-by-b array of random numbers in [0,1].
U=[real(up);imag(up)]; % total initial vector.
it=0; err=1; % initializing error

while(err>1e-6) %1e-6 is our error tolerance.
    it=it+1;
    Ur=U(1:Nx);
    Ui=U(Nx+1:end); % real and imagionary parts of U
    
    J11= -0.5*D2+diag(g*(3*Ur.^2 + Ui.^2) + mu + V);
    J12= diag(2*g*Ur.*Ui);
    J22= -0.5*D2+diag(g*(Ur.^2+3*Ui.^2)+mu+V);
    J = [J11,J12;J12,J22]; % Jacobian
    
    U2=Ur.^2+Ui.^2; % mod square of U
    Fr = -0.5*D2*Ur+(g*U2+V+mu).*Ur; %real RHS/F
    Fi = -0.5*D2*Ui+(g*U2+V+mu).*Ui; % imag RHS/F
    F=[Fr;Fi]; % total F / RHS
    DU=-J\F; % Newton correction
    U1=U+DU; % new newton step
    err=norm(U-U1); %convergence error
    
    U=U1; % update solution
end
u=U(1:Nx)+1i*U(Nx+1:end); % wrapping solution into complex vector
