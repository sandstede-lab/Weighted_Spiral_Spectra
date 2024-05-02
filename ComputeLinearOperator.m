function [D1x,D2x, D1r] = ComputeLinearOperator(par,numPar)

%% rename parameters
	nx = numPar.nx;
	ny = numPar.ny;
    order = numPar.order;
	r1 = par.r1;
	r2 = par.r2;
    
% theta direction

switch numPar.thgrid
    case 'F' % Fourier
        [~,D2x] = fourdif(nx,2);    % 2nd derivative matrix        
        [~,D1x] = fourdif(nx,1);    % 1st derivative matrix
            
    case 'FD' % Finite Differences
        switch order 
            case '2'
            hx = 2*pi/nx;  % Periodic so 0 = 2*pi
            D1x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx); 
            D1x = D1x - D1x';
            D1x(1,end) = -1; D1x(end,1) = 1; % Periodic boundary conditions
            D1x = D1x/(2*hx);  % 2nd order FD centered difference first derivative matrix

            ex = ones(nx,1);
            D2x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx)- sparse(1:nx,1:nx,ex,nx,nx);
            D2x = D2x + D2x';
            D2x(end,1) = 1;
            D2x(1,end) = 1;
            D2x = D2x/(hx^2);  % 2nd order FD centered difference 2nd derivative matrix

            case '4' % 4th order 
            hx = 2*pi/nx;  % Periodic so 0 = 2*pi
            D1x = sparse(1:Nrh-1,[2:Nrh-1 Nrh],8*ones(Nrh-1,1),Nrh,Nrh) - sparse(1:Nrh-2,[3:Nrh-1 Nrh],ones(Nrh-2,1),Nrh,Nrh);
            D1x = (D1x - D1x');
            D1x(1,end-1:end) = [1, -8];D1x(2,end) = 1; D1x(end-1,1) = -1; D1x(end,1:2) = [8,-1]; % Periodic boundary conditions
            D1x = D1x/(12*hx); % First derivative matrix

            D2x = sparse(1:Nrh-1,[2:Nrh-1 Nrh],16*ones(Nrh-1,1),Nrh,Nrh) - sparse(1:Nrh-2,[3:Nrh-1 Nrh],ones(Nrh-2,1),Nrh,Nrh);
            D2x = (D2x + D2x' - 30*speye(Nrh)); 
            D2x(1,end-1:end) = [-1, 16]; D2x(2,end) = -1; D2x(end-1,1) = -1; D2x(end,1:2) = [16,-1]; % Periodic boundary conditions
            D2x = D2x/(12*hx^2); % 2nd derivative matrix
        end
        
end
Ix  = speye(nx,nx);

% radial direction
switch numPar.rgrid
    case 'FD' % No inner hole
        [R,D2y,D1y] = Compute_2D_radial_Laplacian_finite_difference(ny,r2,order);
        hy = r2/(ny-1);
        
    case 'FD_hole' % Inner hole present
        [r,D2y,D1y] = Compute_2D_radial_Laplacian_finite_difference_non_zero_radius(ny,r1,r2,order);
        R = sparse(1:ny,1:ny,1./r,ny,ny);
    
    case 'Cheby'  % Chebyshev grid 
        [r,D2y,D1y] = Compute_2D_radial_Laplacian_cheb_non_zero_radius(ny,r1,r2);
        R = sparse(1:ny,1:ny,1./r,ny,ny);
        
    case 'FD_symm' % No inner hole, assumes radial symmetry (don't use for spiral)
        [R,D2y,D1y] = Compute_2D_radial_Laplacian_finite_difference_symm(ny,r2);
        
	
end
Iy = speye(ny,ny);

% First derivative: d/d(theta)
D1x = kron(Iy,D1x);  % Works for all the cases. Just contains theta derivatives. 

% Assemble the Laplacian: 
D2x  =  kron(R.^2,D2x) + kron(D2y,Ix);

if strcmp(numPar.rgrid,'FD') % Need to add boundary conditions at r = 0
    D2x(1:nx,:) = 0;
    D2x(1:nx,nx+1:2*nx) = (4/nx)/(hy^2);  % Averaging method
    D2x(1:nx,1:nx) = -4/(hy^2).*speye(nx);
   
end

if nargout > 2
   D1r = kron(D1y,Ix); % First derivative matrix for radius
    
end



