% NLS_FD_RK4 : Integrating the pde

g=1; %g<0 for bright solitons and g>0 for dark solitons
N=201;
L=-10; R=10;
x=linspace(L,R,N)';
dx=x(2)-x(1);

mu=1;
A=sqrt(mu); % initial amplitude, A
c=0; % IC veclocity, c
offset=3.8;

Omega=0.1;B=1; beta=0.6;
V=0.5*(Omega^2)*x.^2; %+ B*sech(beta*x).^2;
%V=zeros(N,1);
if(exist('V')==0) V=0;end; %V=0 if not already defined
f0=sqrt(max(mu-V,0));

if(exist('offset')==0) offset=0;end;
u0=(A*tanh(A*(x-(L+R)/2))).*(A*tanh(A*(x-offset-(L+R)/2)).*exp(1i*c*x)).*(A*tanh(A*(x+offset-(L+R)/2)).*exp(1i*(-c)*x)).*f0; %IC for u. tanh for dark solitons, sech for bright
u=u0;allu=u;t=0;allt=t; isave=1; %store position and time

maxtime=40; %final time
dt=0.01; %timestep
nsave=200;ndisp=40; %snapshots to save
maxstep=fix(maxtime/dt);
stopdisp=fix(maxstep/ndisp);
stopsave=fix(maxstep/nsave);
if(dt>2*dx^2/sqrt(2)) error('Need smaller dt.'); end;

figure(1);clf; %maintimeloop
for k=1:maxstep
    t=t+dt;
    u=ODE_RK4(u,N,g,V,dx,dt); %update u using FD+RK4
    
    if(fix(k/stopsave)==k/stopsave) %saving progress
        isave=isave+1;
        allt=[allt,t];
        allu(:,isave)=u;
    end
    if(fix(k/stopdisp)==k/stopdisp) %plottingprogress
        plot(x,abs(u),x,V,"green",x,f0,"red") %real(u),x,imag(u),x,abs(u))
        title(['t=',num2str(t)]);
        axis([L R -0.2*A 1.1*A]);
        legend('abs(u)','Potential','f0','Location','Northeast'); drawnow; %'Re(u)','Im(u)',
    end
end

contourf(abs(allu)');
xlabel("x");
ylabel("t");