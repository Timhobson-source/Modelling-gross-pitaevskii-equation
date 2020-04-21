% NLS_stability1D
% we are computing the spectrum for a steady state of:
% i*u_t = -0.5*u_xx +g*(|u|^2)*u
% code based on that of Panos Kevrekidis et al.

steadystateNewton2; % getting a steady state solution using this program

%Jacobian:
M11 = -(-0.5*D2+diag(2*g*u.*conj(u)+V+mu));
M12 = -diag(g*u.*u);
M21 = -conj(M12);
M22 = -conj(M11);
M = 1i*[M11,M12;M21,M22];
neigs=0; % if neigs>0 compute neigs evals
if(neigs>0)
    optionseigs.disp=0;
    z0=0.5+0*1i;
    [vvv,eee]=eigs(sparse(M),neigs,z0,optionseigs);
else
    [vvv,eee]=eig(M);
end

ee=diag(eee); %eigenvalues
vv=vvv(1:Nx,:)+conj(vvv(Nx+1:end,:)); %eigenvectors
[temp,bb]=(sort(real(ee)));
bb=flipud(bb);
ee=ee(bb);
vv=vv(:,bb);
figure(2);clf
plot(ee,'o')
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')


