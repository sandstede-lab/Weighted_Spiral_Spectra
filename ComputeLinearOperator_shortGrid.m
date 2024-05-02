function [D1x,D2x,D1r] = ComputeLinearOperator_shortGrid(par,numPar)
% Compute differentiation matrices on the short grid (one grid point at the
% origin). 

%% rename parameters
	nx = numPar.nx;
	ny = numPar.ny;
    order = numPar.order;
	r2 = par.r2;
    
% angular/theta direction

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
    case 'FD' % Standard FD on disk. Zero Neumann boundary conditions at r = 0 and r = r2
        [R,D2y,D1y] = Compute_2D_radial_Laplacian_finite_difference(ny,r2,order);
        
       
    case 'FD_weighted' % On a weighted grid - used for spectral computations with eigS
        
        hxn = r2/(ny-1);
        rn = (0:ny-1)'*hxn; rn(1) = 1;                      % radial mesh
        R = sparse(1:ny,1:ny,1./rn,ny,ny); 

        switch order
            case '2'
                % 2nd order 
                ex = ones(ny,1);
                D1y = sparse(1:ny-1,[2:ny-1 ny],ones(ny-1,1)/2,ny,ny); % 1st derivative matrix
                D1y = (D1y - D1y')/hxn; 
                D1y(1,:) = 0; D1y(ny,:) = 0; 
                D1y(ny,ny) = -par.w;  % Weighted boundary condition
        
                D2xn = sparse(1:ny-1,[2:ny-1 ny],ones(ny-1,1),ny,ny) - sparse(1:ny,1:ny,ex,ny,ny); % 2nd derivative matrix
                D2xn = (D2xn + D2xn'); 
                D2xn(1,:) = 0; 
                D2xn(1,1)=-2 ; D2xn(1,2) = 2;
                D2xn(ny,:) = 0 ;     % Neumann boundary conditions - alter r = 0 conditions after kron product
                D2xn(ny,ny-1) = 2;   % Weighted boundary conditions at outer edge
                D2xn(ny, ny) = -2*(hxn*par.w + 1); % Weighted boundary conditions at outer edge
                D2xn = D2xn./hxn^2;
        
        
            case '4'
                % 4th order 
                D1y = sparse(1:ny-1,[2:ny-1 ny],8*ones(ny-1,1),ny,ny) - sparse(1:ny-2,[3:ny-1 ny],ones(ny-2,1),ny,ny);
                D1y = (D1y - D1y')/12;  % r = 0 will be altered after kron
        
                % Boundary conditions in D1y
                D1y(1:2,:) = 0;
                D1y(2,1:3) = [-1/2,0,1/2];
                D1y(ny-1:ny,:) = 0; 
                D1y(ny-1,ny-2:ny) = [-1/2,0,1/2];  % second order at the boundary
                D1y = D1y/hxn;  % First derivative matrix
                
                D1y(ny,ny) = -par.w;  % weighted boundary condition at r= R
            
        
                D2xn = sparse(1:ny-1,[2:ny-1 ny],16*ones(ny-1,1),ny,ny) - sparse(1:ny-2,[3:ny-1 ny],ones(ny-2,1),ny,ny);
                D2xn = (D2xn + D2xn' - 30*speye(ny))/12; 
                D2xn(end-1:end,:) = 0;  D2xn(1:2,:) = 0; % Neumann boundary conditions

                % Use 2nd order at the boundary
                D2xn(1,1:2) = [-2,2];
                D2xn(2,1:3) = [1,-2,1];
                D2xn(end-1,end-2:end) = [1, -2, 1];  
                D2xn(end,end-1:end) = [2,-2*(hxn*par.w + 1)]; % weighted boundary condition
        
                % Use one-sided 4th order at outer boundary - alternate possibility
                %D2xn(end-1,end-6:end-1) = [-5/6, 61/12, -13, 107/6, -77/6, 15/4];  
                %D2xn(end,end-5:end) = [-5/6, 61/12, -13, 107/6, -77/6, 15/4];
                
                
                D2xn = D2xn./(hxn.^2); % 2nd derivative matrix   
        end

        %% 2D radial
        D2y = D2xn + R*D1y; % d_rr + d_r/r
  
end

%% Output matrices

hy = r2/(ny-1);
Iy = speye(ny,ny);

% First angular derivative: d/d(theta)
D1x = kron(Iy,D1x);  % Works for all the cases. Just contains theta derivatives. 
D1x = D1x(nx:end,nx:end); % Remove r=0 except for 1 row and comlumn
D1x(:,1) = 0; D1x(1,:) = 0; % Redundant - don't take the theta derivative at the origin (just one grid point). 


%% Assemble the Laplacian: 
% Angular portion
L1 = kron(R.^2,D2x);    % Components of the Laplacian
L1 = L1(nx:end,nx:end); % r1, r2 do not need to be updated
L1(:,1) = 0; L1(1,:) = 0; % r0 component - don't take the angular derivatives at the origin (just one point). 

% Radial components: 
switch order
    case '2'  % Confirmed this case works
        L2 = kron(D2y,Ix);      % r1 blocks need to be updated. 
        L2 = L2(nx:end,nx:end);
        L2(:,1) = 0; L2(1,:) = 0;   % r0 component

        L2(2:nx+1,1) = 1./(2.*hy.^2); % r1 update: 2nd order

        D2x = L1 + L2;

        % 5 point stencil at origin: 2nd order!!!
        D2x(1,:) = 0;
        stencil_idx = [2,floor(nx/4)+2, floor(nx/2)+2, floor(3*nx/4) + 2]; % gridpoints for the 5 point stencil (in r1) 
        D2x(1,stencil_idx) = 1./(hy.^2);
        D2x(1,1) = -4./(hy.^2);

    case '4' % Has 2nd order at boundaries

        L2 = kron(D2y,Ix);      % r1 blocks need to be updated. 
        L2 = L2(nx:end,nx:end);
        L2(:,1) = 0; L2(1,:) = 0;   % r0 component

        % r1 update: use 2nd order
        L2(2:nx+1,1) = 1./(2.*hy.^2); % r1 update: 2nd order

        % r2 update: use 2nd order
        L2(nx+2:2*nx+1,1) = -1./(24.*hy.^2); % r0 component of laplacian

        D2x = L1 + L2;

        % 5 point stencil at origin: 2nd order!!!
        D2x(1,:) = 0;
        stencil_idx = [2,floor(nx/4)+2, floor(nx/2)+2, floor(3*nx/4) + 2]; % gridpoints for the 5 point stencil (in r1) 
        D2x(1,stencil_idx) = 1./(hy.^2);
        D2x(1,1) = -4./(hy.^2);

end


if nargout > 2

     D1r = kron(D1y,Ix);              % First derivative matrix for radius
     D1r = D1r(nx:end,nx:end);        % On short grid.
      
     % 2nd order at the boundary near the origin (i.e., need to update
     % radius 1)
     D1r(:,1) = 0; % Set first column to 0
     D1r(2:nx+1,1) = -1./(2.*hy);  

     if order == '4'

        % Need to update the radius 2 portion too

        D1r(nx+2:2*nx+1,1) = 1./(12.*hy);


     end

end








