% steadystateNewton2D

%Setting up the space:

mu = 0.5;  %time frequency of steady state
g = -1; % g=-1,1 are important focusing cases we consider for the PDE
L=-10; R=10; %left and right bounds for x variable
N=201; %number of discrete points in 1D mesh (square mesh)
x=linspace(L,R,N); %discrete x-space
y=linspace(L,R,N); % discrete y-space
h=x(2)-x(1); % = y(2)-y(1): mesh size

% Creating Laplacian:


    




