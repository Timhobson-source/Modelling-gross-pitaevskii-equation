% generating solution plots
L=-10; R=10; %left and right bounds of the interval
Nx=201; %number of discrete points in 1D mesh (x-direction)
x=linspace(L,R,Nx); %discrete space we solve the steady stae over
Omega=0.075;
B=0.3; beta=0.5;
V=0.5*(Omega^2)*x.^2 + B*sech(beta*x).^2;
V=V';
offset=4;
steadystateNewton2;
figure(1)
solu=real(u); %since solution is fully real.
plot(x,u0,'.',x,solu,x,V,'--')
legend('u_0(x)','u(x)','V(x)')
xlabel('x')
ylabel('u(x), V(x)')

NLS_stability1D;