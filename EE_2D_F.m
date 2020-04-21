% F function
function F = EE_2D_F(u,V,g,dx,dy,D)
F = g*conj(u).*u.^2 + V.*u +0.5*((D*u)*(dx^-2) + (u*D)*(dy^-2));
return;