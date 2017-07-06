% Function to sample parameter space and calculate coefficients according
% to SKM algorithm, adapted from Murabito et al. 2014
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-05-11, last modified: 2017-07-06

function [DataOut] = MJanasch_SKM_Sampling_Code(iterations,InputDataStructure,MetConcDataIn)

% [CJ_rec,CS_rec,E_rec,MaxRealEigens,Parameters] = Sampling_code(M,Fluxes1,FullSto,RedSto,Link,seed)
%
% INPUT:
%
% * N       - Matlab structure containing the specifications of the model.
%             Such data structure is generated by the function
%             N = TranslateSBML(); N must be enriched with the cmetabolite
%             concentrations at steady-state.
%
% * Fluxes  - Vector containing the fluxes at theady state of the different
%             reactions in the model N. The entries in 'Fluxes' correspond
%             to the reactions in M, i.e. Fluxes(i) is the flux of the ith
%             reaction in N (N.reaction(i))
%
% * SFull   - Full stoichiometry matrix
%
% * SRed    - Reduced stoichiometry matrix
% 
% * L       - Link matrix
% 
% * seed    - (optional) seed to initialize the random number generator
%
%
% OUTPUTS:
%
% * CJ_rec  - 3-dimensional matrix of the scaled flux control coefficients
%             for all parameter sets resulting in the steady-state being
%             stable. CJ_rec(i,j,z) is the control exerted by enzyme j upon
%             flux i given the parameters sampled at iterazion z.
%
% * CS_rec  - 3-dimensional matrix of the scaled concentration control
%             coefficients for all parameter sets resulting in the steady-
%             state being stable.
%
% * E_rec   - 3-dimensional matrix of the scaled elasticity coefficients
%             for all parameter sets resulting in the steady-state being
%             stable.
%
% * MaxRealEigens -
%             Vector containing the maximal real part of the eigenvalues
%             for each parameter sampling iteration.
%
% * Parameters - 
%             2-dimensional matrix containing all the system's parameters
%             for each sampling iteration. Each row correspond to a
%             sampling iteration and each column to a parameter.
%


clc;

%--------------- VARIABLES DECLARATION ---------------%

% iterations = 1000; % No. of sampling iterations (Orig iterations=20000)
% Here as input to the function, defined in command line

load(InputDataStructure); % Load N, Fluxes, SFull, SRed and L

% Load metabolite concentration
for h = 1:length(N.species)
        N.species(h).initialConcentration = MetConcDataIn(1,h);
end



nOfSS = 0; % No. of sampling iterations resulting in stable steady-states.

F1 = 0.1; % Multiplicative factor defining the lower bound of the sampling intervals
F2 = 10;  % Multiplicative factor defining the upper bound of the sampling intervals

MaxRealEigens = zeros(iterations,1); % At iteration i, the maximal real
                                     % part of the eigenvalues is stored in
                                     % MaxRealEigens(i,1)
                                     
[m,n] = size(SFull); % Size of the full stoichiometric matrix
[m1,n] = size(SRed); % Size of the reduced stoichiometric matrix

CS = zeros(m,n); % Concentration Control Coefficients
CJ = zeros(n,n); % Flux Control Coefficients
                 
CJ_rec = zeros(n,n,iterations); % Flux control coefficients of stable systems
CS_rec = zeros(m,n,iterations); % Concentration control coefficients of stable systems
E_rec = zeros(n,m,iterations);  % Elasticity coefficients of stable systems
dfodc = {};
dfodc = zeros(n,m); % Matrix containing the numeric value of the derivative
                    % of the fluxes with respect to the metabolite
                    % concentrations.

RedJac1 = zeros(m1,m1); % Reduced Jacobian Matrix

R = {N.reaction.name}; % Names of the reactions
S = {N.species.id};    % Names of the metabolites



ID = eye(n);  % Unitary matrix with n. of columns/rows = n. of reactions

[ParID,ParVal] = GetAllParameters(N);
Parameters = zeros(iterations,length(ParID)); % At iteration i, all the
                                              % parameters (sampled and 
                                              % const) are stored in row i
                                              % of this matrix.


Conc = [];                                    % Initialize concentration 
                                              % vector
for i=1:length(N.species)
    if ~N.species(i).boundaryCondition        % if the boundary condition=0
                                              % (=internal metabolite, = no
                                              % constant set concentration)
        Conc(i) = N.species(i).initialConcentration; % fill in initial conc
    end
end
Conc = Conc(:);                                % align all columns after 
                                               % each other (example:
                                               % 6x5-matrix to 30x1 matrix)
Fluxes = Fluxes(:);                          % align all columns


%-- COMPUTE THE NORMALIZATION MATRICES FOR CJ and CS --%
nrmlz1 = compute_normalization_matrix(Fluxes);% create rxr matrix of fluxes
                                               % normalized by fluxes
                                               % every flux normalized/divided
                                               % by every flux
nrmlz3 = compute_normalization_matrix2(Fluxes,Conc);
                                               % create rxr matrix of conc
                                               % normalized by fluxes


%-- RETRIEVE THE FUNCTIONAL FORM OF THE RATE EQUATIONS --%
for j=1:length(R)
    REQ{j} = N.reaction(j).kineticLaw.formula;
end
[DER,Ssym] = compSymDeriv(N); % Contains the symbolic partial derivatives of
                       % the rate equations with respect to the 
                       % metabolites (concentrations that change, no 
                       % external metabolites)
                 
                                              
%-- COMPUTE THE SAMPLING INTERVAL OF EACH PARAMETER --%
%-- ACCORDING TO THE CONCENTRATIONS OF METABOLITES  --%
[ParValMin,ParValMax,vmax_indeces] = compParInt(N,F1,F2);

eval('default = 1;');       % Express that the default value for almost all
                            % parameters (except Ki, Ka, Keq)

tic                         % start timer

if ~exist('seed','var'); seed = 1; end; % if no seed defined, set it to 1
rng(seed);                              % define control of random number 
                                        % generater

for c = 1:iterations                    % for every interation
    
    %-- SAMPLE/EVALUATE PARAMETERS --%
    for p=1:length(ParID)               % for every parameter
        ParAux = log(ParValMin(p)) + rand(1)*(log(ParValMax(p)) - log(ParValMin(p)));
                                        % create a value, can be negative
            
        ParAux = exp(ParAux);           % e^alue, now always positive!
        eval([ParID{p},' = ParAux;']);  % SET PARAMETER VARIABLE TO
                                        % VALUE
        % eval is used here, since it literally creates the variables and
        % associates the related value (ParAux).
                                                                                
        Parameters(c,p) = ParAux;       % save the p-th parameter for the 
                                        % c-th iteration in matrix as well
    end

    %-------------------------------------------------------------
    %-- DETERMINING THE PARTIAL DERIVATIVES FOR STATE 1 --%
    %--  Assigning concentration values to metabolites  --%
    for i = 1:length(S)
        eval([S{i},' = N.species(i).initialConcentration;']); 
        % Like with parameters above, eval literally creates all the
        % substrate variables while looping through it and here associating
        % the initial concentrations 
    end

    %-- Computing Vmaxs --%
    for j = 1:length(vmax_indeces)
        eval([ParID{vmax_indeces(j)},' = 1;']);     % Set Vmax to 1
        % create all the Vmax variable and set them to 1
        
    end
    for j=1:length(R)                               % loop through reactions
        pippo = eval([REQ{j},';']);                            
        % Here the use of eval, with the creation of all the parameters
        % above, calculates the result(=pippo) of every rate equation with these
        % parameters, but Vmax=1!!! E.g. result of rate equation with all
        % parameters, but without Vmax!!! synthetic variable 'pippo'
        
        Fluxes(j);                                 % Not sure what this is
                                                    % about
                                                    
        Vmax = Fluxes(j)/pippo;                   
        
        % Connection between pippo and real fluxes is Fluxes = Vmax * pippo
        % therefore Vmax's for certain parameter set is Vmax=Fluxes/pippo
        
        
       
        eval([ParID{vmax_indeces(j)},' = Vmax;']);  % Create Vmax variable 
                                                    % and save new Vmax
        Parameters(c,vmax_indeces(j)) = Vmax;       % Save new Vmax in the 
                                                    % Parameters-matrix for
                                                    % the c-th iteration 
    end

    %-- Evaluating numerical value of partial derivatives --%
    for q=1:size(DER,1)
        for w=1:size(DER,2)
            dfodc(w,q) = eval([DER{q,w}]); % calculating results of derived
                                           % equations with sampled
                                           % parameters, saved in 'dfodc'
        end
    end
    %-------------------------------------------------------%

    % ---- COMPUTING THE REDUCED JACOBIANS ----%
    RedJac1(:,:) = SRed(:,:)*dfodc(:,:)*L(:,:);
    
    %---- COMPUTING THE EIGENVALUES ----%
    eigenvalues = eig(RedJac1);             % Calculate eigenvalues
    MaxRealEigen = max(real(eigenvalues));  % Find the maximum of the real
                                            % part of the eigenvalues
    MaxRealEigens(c,1) = MaxRealEigen;      % Save these for c-th iteration
                                            % in matrix
    
    %---- CHECKING FOR STABILITY ----%
    ok = 1;                                 % Initialize binary checking-
                                            % variable
    if(MaxRealEigen >= 0)                   % If the maximum of the real 
                                            % part of the eigenvalues is
                                            % larger or equal to 0
        ok = 0;                             % set checking-variable to 0
    end

    if(ok == 1)                             % If the checking-variable = 1
        
        nOfSS = nOfSS + 1;                  % Increase number of stable SS
        
        %---- COMPUTING AND RECORDING ELASTICITY COEFFICIENTS ----%
        E_rec(:,:,nOfSS) = dfodc .* (nrmlz3(:,:))';
        
        %---- COMPUTING THE CONCENTRATION CONTROL COEFFICIENTS ----%
        CS(:,:) = -(L*(inv(RedJac1(:,:)))*SRed);

        %---- COMPUTING THE FLUX CONTROL COEFFICIENTS ----%
        CJ(:,:) = (ID + dfodc(:,:)*CS(:,:))./nrmlz1(:,:);
        CJ_rec(:,:,nOfSS) = CJ; 
        CS_rec(:,:,nOfSS) = CS./nrmlz3;   
    end
% %% Output every 10 iterations about the % of stable steady states
%     ratio = (nOfSS/c)*100;
%     if(~mod(c,10) || c == iterations)
%     	fprintf('%d%s%d%s%6.2f%s\n',nOfSS,' -> ',c,' (',ratio,'%)');
%     end
end

toc

CJ_rec(:,:,(nOfSS+1):end) = [];
CS_rec(:,:,(nOfSS+1):end) = [];
E_rec(:,:,(nOfSS+1):end) = [];

%% Define output-data
DataOut.CJ_rec      = CJ_rec;
DataOut.CS_rec      = CS_rec;
DataOut.E_rec       = E_rec;
DataOut.Parameters  = Parameters;
DataOut.Conc        = Conc;
DataOut.dfodc       = dfodc;
DataOut.ParID       = ParID;


%==========================================================================
function nrmlz = compute_normalization_matrix(Fluxes)

nrmlz = Fluxes * ((1./Fluxes)');

%==========================================================================
function nrmlz = compute_normalization_matrix2(Fluxes,Conc)

nrmlz = Conc * ((1./Fluxes)');

%==========================================================================
function [DER,Ssym] = compSymDeriv(N)

R = {N.reaction.name};
S = {N.species.id};

% Transform S cell array into symbolic variables
Ssym=sym(S);


DER = {}; % This structure will contain the DERIVATIVES
          % of reaction j with respect to metabolite i.
    
% Computation of the derivatives
% fprintf('%s', 'Computing symbolic derivatives... ');
for j = 1:length(R)
    for i = 1:length(S)
        if ~(N.species(i).boundaryCondition) % find the ones with no constant
                                             % concentration,
                                             % diffenrentiate only if
                                             % changing metabolite conc,
                                             % 
                                             % When no changing conc., 
                                             % no sense in differentiating
                                             
            isitthere = IsSpeciesInReaction(N.species(i),N.reaction(j));
            if isitthere ~= 0
                %eval(['syms ', S{i}]);
                %DER{i,j} = diff(M.reaction(j).kineticLaw.formula,S{i});
                DER{i,j} = diff(N.reaction(j).kineticLaw.formula,Ssym(i));
            else
                DER{i,j}='0';
            end
        end
    end
end
% fprintf('%s\n', 'Done :-)');


%==========================================================================
function [ParValMin,ParValMax,vmax_indeces] = compParInt(N,F1,F2)

% Take in M and the boundary factors for the sampling intervals 

% Give out min and max values for parameters and indeces for all Vmax



% fprintf('%s', 'Evaluating parameter sampling intervals... ');

Conc = [];
S = {N.species.id};
for i=1:length(S)
    Conc(i) = N.species(i).initialConcentration;    % Get all initial conc.
end                                                 % not just internal

[ParID,ParVal] = GetAllParameters(N);   % Get all reaction parameters

ParValMin = ParVal;                     % Set parameter minimum bound to 
                                        % parameter value
ParValMax = ParVal;                     % Set parameter maximum bound to
                                        % parameter value

count=1;                                % initialize counting variable
for p = 1:length(ParID)                 % for every parameter
    
    par_aux = ParID{p};                 % for parameter p
    if(par_aux(end-1) == 'v')           % if it has a 'v' as second last 
                                        % character:
        par_aux2 = par_aux(3:end-2);    % save only NAME or 'max' or 'eq' 
    else
        par_aux2 = par_aux(3:end-3);    % if v as third last, meaning two
                                        % digit number reactions, save also
                                        % only NAME or 'max' or 'eq'
    end
    
    if ((~strcmp(par_aux2,'eq')) && (~strcmp(par_aux2,'max'))) 
        % if parameter p is NOT an equilibrium constant or a Vmax value
                                                                

           for i=1:length(S)               % for every metabolite
               if(strcmp(par_aux2,S{i}))   % if the parameter belongs to 
                                            % that metabolite (Km or Ki or
                                            % Ka)
                   ParValMin(p) = Conc(i)*F1; % set lower bound with conc
                   ParValMax(p) = Conc(i)*F2; % set upper bound with conc
                   break;
               end
           end
        
    elseif (strcmp(par_aux2,'max'))     % if parameter is a Vmax, set all 
                                        % to 1 (should be already?),
                                        % supposingly following of 
                                        % normalization, max V set to 1
        ParValMin(p) = 1;
        ParValMax(p) = 1;    
        vmax_indeces(count) = p;        % get Vmax indeces back (lost in 
                                        % par_aux-definition in beginning
                                        % of this subfunction (which
                                        % position there are in "Parameter"
        count = count+1;      
    end
end
% fprintf('%s\n', 'Done :-)');
