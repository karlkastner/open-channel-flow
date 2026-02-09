% test geometry for RK-flow
L=1024; n=L; x = (0:n-1)'*L/n; y = x'; s=L/6; zb = 0-(1-cos(2*pi*x/L)).*(1-cos(2*pi*y/L))/4; imagesc(zb); axis equal; colorbar; %plot(zb(end/2,:))
