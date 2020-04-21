% 2D NLS_FD_EE
g=-1; %g<0 for bright solitons and g>0 for dark solitons
N=201;
L=-10; R=10;
x=linspace(L,R,N)';
y=linspace(L,R,N)';
dx=x(2)-x(1);
dy=y(2)-y(1);
mu=0.5;
A=sqrt(2*mu); % initial amplitude, A
c=[0,0]; % initial velocity
[X,Y]=meshgrid(L:0.5:R,L:0.5:R);
u0=A*sech(A*sqrt((X-(L+R)/2).^2+(Y-(L+R)/2).^2)).*exp(1i*(c(1)*X+c(2)*Y));
