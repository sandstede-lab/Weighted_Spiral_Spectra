%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find pseudospectra of weighted spiral wave operator
% Options: pseudospectra (via svds) or condition number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clear;

comp_type = 'pseudoSpec'; % options: 'condNum','pseudoSpec'

file = 'Barkley_spiral_r25_h0p05_delta0p2_b0p001_a0p7_ep0p02_positiveOmega.mat';
out_file = '';

load(file);

spectra_problem = 'Barkley_jacobian_weighted_operator'; % Defines the spectral operator 
fh = str2func(spectra_problem); 


% Seed points for computations: See points should not be eigenvalues
loc_real = -3:0.2:0.2; % Real parts of grid coordinates
loc_imag = [linspace(0.01,2*par.omega+0.01,25), 2*par.omega+0.1];   % Imaginary parts of grid coordinates
[LR,LI] = meshgrid(loc_real,loc_imag); % grid 
LR = LR(:); LI = LI(:);
loc = LR + 1i.*LI;

weights = [0,0.5,1.0,1.5,2]; % Radial weight values

% Put spiral on short grid (one grid point at the origin)
% Method: Wheeler & Barkley SIADS (2006)
V = U(numPar.nx*numPar.ny+1:end);
U = U(1:numPar.nx*numPar.ny);

V = V(numPar.nx:end);
U = U(numPar.nx:end);

% Differentiation matrices options
numPar.rgrid = 'FD_weighted'; 

for m = 1:length(weights)
    
    par.w = weights(m);
    [L1,L2,L1r] = ComputeLinearOperator_shortGrid(par,numPar);

    J = fh(U,V,L1,L2,L1r,par,numPar); % Jacobian operator (weighted)

    II = speye(2*(numPar.nx*(numPar.ny-1)+1),2*(numPar.nx*(numPar.ny-1)+1));

    [LR,LI] = meshgrid(loc_real,loc_imag);
    
    switch comp_type

        case 'condNum'
            condition_number = zeros(length(loc),1); 

            parfor (j = 1:length(loc),6) % Loop over the seed points
    
                disp(j/length(loc)) % percentage done with computations
                 
                L_shift = J - loc(j).*II ;
                condition_number(j) = condest(L_shift);
            end

            condition_number = reshape(condition_number, size(LR));
            save([out_file num2str(m) '.mat'],'condition_number','LR','LI','par','numPar','file')

        case 'pseudoSpec'

            pseudospec = zeros(length(loc),1);

            parfor (j = 1:length(loc),6) % Loop over the seed points

                disp(j/length(loc)) % percentage done with computations
                pseudospec(j) =  svds( loc(j).*II - J, 1,'smallest') ;

                %%% Alternate version if sparse computations don't work well. Only works with small number of grid points. 
                %if pseudospec(j) == 0 
                    %pseudospec(j) = min(svd(full(loc(j).*II - J), 'econ'));
               % end
            end

            pseudospec = reshape(pseudospec, size(LR));
            save([out_file num2str(m) '.mat'],'pseudospec','LR','LI','par','numPar','file')

    end

end





