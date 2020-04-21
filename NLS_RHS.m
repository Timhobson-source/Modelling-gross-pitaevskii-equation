% NLS_RHS (function used in ODE_RK4)
function RHS = NLS_RHS(u,N,dx,g,V)
up=[u(N);u(1:N-1)]; %periodic BCs
um=[u(2:N);u(1)]; %periodic BCs
% up=[0;u(1:N-1)]; %Zero BCs
% um=[u(2:N);0]; %Zero BCs
% up=[u(1)-(u(2)-u(1));u(1:N-1)]; %linear interp
% um=[u(2:N);u(N)-(u(N-1)-u(N))]; %linear interp
% up=[u(1);u(1:N-1)]; %Laplace zero BCs
% um=[u(2:N);u(N)]; %Lapce zero BCs
RHS=1i*((0.5/dx^2)*(up-2*u+um)-(g*u.*conj(u)+V).*u);
return;
