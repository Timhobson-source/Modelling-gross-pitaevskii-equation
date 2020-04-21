% 2D Explicit Euler
g=-1;
L=-10; R=10;
N=100;
x=linspace(L,R,N); y=linspace(L,R,N);
dx=(R-L)/N; dy=(R-L)/N;
[x,y]=meshgrid(x,y);
w=0.5; A=sqrt(2*w);
ONE=ones(N,1);
D=spdiags([ONE,-2*ONE,ONE],-1:1,N,N);

if(exist('V')==0) V=0;end;

r=[0,0]; % initial position
c=[1,2]; % initial velocity
u0=A*sech(A*sqrt((x-r(1)).^2+(y-r(2)).^2)).*exp(1i*(c(1)*x+c(2)*y));

%Plotting u0:
%surf(x,y,abs(u0))

%colorbar

% EE method:
dt=0.1;
T=10; %max time
steps=floor(T/dt); %
figure(1);clf;
t0=0;
t=t0;
u=u0;
for k=1:steps
   u=u+dt*EE_2D_F(u,V,g,dx,dy,D);
   t=t+dt;
   surf(x,y,abs(u))
   xlabel("x")
   ylabel("y")
   zlabel("u")
   title("t="+num2str(t))
   drawnow;fprintf('Press any key to continue...\n'); pause
end
print("done!")
