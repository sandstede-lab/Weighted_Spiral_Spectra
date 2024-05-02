%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find eigenvalues of spiral 
% Uses sparse methods (eigs)
% Uses radial weight w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear;

file = 'Barkley_spiral_r25_h0p05_delta0p2_b0p001_a0p7_ep0p02_positiveOmega.mat';
out_file = '';

spectra_problem = 'Barkley_jacobian_weighted_operator'; % Defines the spectral operator 
fh = str2func(spectra_problem); 

load(file);

eval_option = 'smallestabs';

switch eval_option
    case 'locs'

        % Seed points for eigs: define grid. Helpful to have seed points not be an eigenvalue
        loc_real = -2.9:0.5:0.1;
        loc_imag = -1:0.1:-0.1;
        [LR,LI] = meshgrid(loc_real,loc_imag); % grid
        LR = LR(:); LI = LI(:);
        loc = LR + 1i.*LI;
    
        numVals = 10;   % Number of eigenvalues to compute around each seed point


    case 'smallestabs'

        numVals = 400;


end

weights = [0, 0.5,1,1.5]; % Radial weights


% Put spiral on short grid (one grid point at the origin)
% Method: Wheeler & Barkley SIADS (2006)
V = U(numPar.nx*numPar.ny+1:end);
U = U(1:numPar.nx*numPar.ny);

V = V(numPar.nx:end);
U = U(numPar.nx:end);

numPar.rgrid = 'FD_weighted'; % construct weighted operators (outer boundary conditions modified)


for k = 1:length(weights)

    disp(k);
    par.w = weights(k); % Weight


    [L1,L2,L1r] = ComputeLinearOperator_shortGrid(par,numPar);
    J2 = fh(U,V,L1,L2,L1r,par,numPar);

    switch eval_option
        case 'locs'
            vals = cell(length(loc),1);

            for j = 1:length(loc)

               disp(j)
               vals{j} = eigs(J2,numVals,loc(j));
            
            end

            vals = cell2mat(vals);

             save([out_file num2str(k) '.mat'],'numPar','par','vals','file','eval_option','loc');

        case 'smallestabs'

            vals = eigs(J2,numVals,'smallestabs');

            save([out_file num2str(k) '.mat'],'numPar','par','vals','file','eval_option');

    end


   

end



