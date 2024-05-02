function J = Barkley_jacobian_weighted_operator(U,V,L1,L2,L1r,par,numPar)
% Linearization about spiral wave in a weighted radial frame. 
% Jacobian for spectral computations
% Assumes differentiation operators and spiral wave in short grid format. 
% No-flux boundary conditions
% Boundary conditions implicitly built into the finite difference operators
% i.e., still apply the PDE at the boundary
%% Set up

% Numerical parameters
nx = numPar.nx;
ny = numPar.ny;

% System parameters
a = par.a;
b = par.b;
ep = 1./par.ep;
delta = par.delta;  % Diffusion coefficient of v-equation
w = par.w;          % Weight
omega = par.omega;  % Angular frequency

% Define radial mesh for weight
hxn = par.r2/(ny-1);
R = (1:ny-1)*hxn;                    % radial mesh (not including origin)
R = repmat(R,nx,1);
R = R(:); R = [1;R];                 % Add one more point for the origin
R = w./R;   % Already includes weight


%% Jacobian              
nonlin_jacob1 = ep*(2*(1+b/a)*U-b/a-(1/a)*V-3*U.^2+(2/a)*U.*V); % d(nonlin)/dU
nonlin_jacob2 = ep*(-(1/a)*U + (1/a)*U.^2);                     % d(nonlin)/dV


dU = [ L2  + 2.*w.*L1r + omega*L1  + spdiags(nonlin_jacob1 + R + w.^2, 0,nx*(ny-1)+1,nx*(ny-1)+1);
    speye(nx*(ny-1)+1,nx*(ny-1)+1)];


dV = [ spdiags(nonlin_jacob2,0,nx*(ny-1)+1,nx*(ny-1)+1);
    delta.*L2 + 2.*delta.*w.*L1r  + omega*L1 + spdiags( delta.*R -1 + delta.*w.^2, 0,nx*(ny-1)+1,nx*(ny-1)+1)];


J = [dU, dV];
J = sparse(J);

    
    
    
    
    
    
    

