function plot_spiral(uout,par,numPar)

th = linspace(0,2*pi,numPar.nx+1);
r = linspace(par.r1,par.r2,numPar.ny);

uu = reshape(uout(1:numPar.nx*numPar.ny),numPar.nx,numPar.ny)';
vv = reshape(uout(numPar.nx*numPar.ny+1:2*numPar.nx*numPar.ny),numPar.nx,numPar.ny)';

[X,Y] = meshgrid(th,r); u_new = [uu uu(:,1)];
[th, r, z] = pol2cart(X,Y,u_new);
figure; pcolor(th,r,z); shading interp;  % U equation
title('U-equation','FontSize',16);
colorbar; set(gca,'fontsize',16);
drawnow;

v_new = [vv vv(:,1)];
[th, r, z2] = pol2cart(X,Y,v_new);
figure; pcolor(th,r,z2); shading interp;  % U equation
title('V-equation','FontSize',16);
colorbar; set(gca,'fontsize',16);
drawnow;
