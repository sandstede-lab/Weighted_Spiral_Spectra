%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a numerical spiral wave solution for given model
% Assumes you have an initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear;

file_names.spiral = 'Barkley_spiral_r25_h0p05_delta0p2_b0p001_a0p7_ep0p02_positiveOmega.mat'; % Initial data
file_names.out_name = '';  % Out file 

file_names.problem = 'Barkley_2D_rotating';  % Name of system to solve

u0 = load(file_names.spiral); % Load the data
par = u0.par;
numPar = u0.numPar;
U0 = u0.U; 

par.Free = 'omega';


%% Additional set up and variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...
                   'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',70);
      
fh = str2func(file_names.problem); 

% Phase condition: picks out R/2
pc = zeros(numPar.nx,numPar.ny);
pc(:,ceil(numPar.ny/2)) = 1; 
phase_cond.pc = pc(:);

% Compute linear operator
R2 = par.r2;
par.r2 = 1;  % Compute the linear operators assuming normalized disk 
[L1, L2] = ComputeLinearOperator(par,numPar); 
par.r2 = R2; % Reset

%% Solve spiral

% Phase condition
us = (phase_cond.pc).*U0(1:numPar.nx*numPar.ny);
phase_cond.u_star = us(phase_cond.pc>0);

us_th = (L1*U0(1:numPar.nx*numPar.ny)); % d/d(theta) of phase condition
phase_cond.u_star_th = us_th(phase_cond.pc > 0);

initial_sol = [U0;par.(par.Free)];
uout = fsolve(@(y) fh(y,L1,L2,par,numPar,phase_cond),initial_sol,options); 
par.(par.Free) = uout(end);


plot_spiral(uout,par,numPar);

% save final data
U = uout(1:end-1);


%save(file_names.out_name ,'U','par','numPar');













