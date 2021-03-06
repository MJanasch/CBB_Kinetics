% Function to sample metabolite concentrations for CBB-model from NET-ranges
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-05-04, last modified: 2018-08-26

function [Y,MetConcDataSet,dGDataSet] = MJanasch_CBB_Metabolite_Sampling(NrSampling,InputDataStructure,InputNET)
%% function MJanasch_CBB_Metabolite_Sampling
% This function samples metabolite concentrations in the range of the NET-
% analysis and checks their thermodynamic consistency in the CBB network

%%%---Input---
% Nr of Samplings
% Concentration ranges from NET-analysis
% InputDataStructure, which contains:
%   - the proper Model-Structure (N), including Keq
%   - Stoichiometric matrices

%%%---Output---
% Matlab-files containing metabolite concentrations


tic % Start Timer

%% Parameter
T           = 303.15;       % Temperature in [K]
R           = 0.0083144598; % Universal gas constant in [kJ/(mol*K)]
%NrSampling  = 1E+7;        % Number of Sampling rounds (1E6 = 420 sec)

%% Initialize 

MetConcDataSet=[];

%% Read in concentration data from NET-analysis

% The function 'X=importdata('filename.txt')' separates numerical values 
% from text, stored in X.data and X-textdata,respectively
NET_Raw = importdata(InputNET); % Read in raw data from file
                                             

for i = 1:length(NET_Raw.data)  % For every metabolite
    
    Y.RangeNames(i,1)   = NET_Raw.textdata(i+1);  % Save the name of the met,
                                                % (i+1) because the first
                                                % entry is "Name"
    Y.MinBound(i,1) = NET_Raw.data(i,1);      % Save lower bound value, in mM!
    Y.MaxBound(i,1) = NET_Raw.data(i,2);      % Save upper bound value, in mM!
end
%% Read in Stoichiometry and Keq from N-structure File

load(InputDataStructure);

for j = 1:length(N.reaction)    % For every reaction
    Y.RxnNames{j,1}   = N.reaction(j).name;             % Save the name
    
    % Not Necessary and not part of standard SBML-structure
    %Y.RxnStoic{j,1}   = N.reaction(j).rxnStoichiometry; % Save stoich. 
                                                        % formula
                                                        
    % For each reaction, loop through parameters
    for k = 1:length(N.reaction(j).kineticLaw.localParameter)
        % If the first 4 characters of the parameter name are 'K_eq', it
        % means it is the equilibrium constant. If that is true, save the
        % value for that parameter
        if strncmp(N.reaction(j).kineticLaw.localParameter(k).name,'K_eq',4)
            Y.RxnKeq(j,1) = N.reaction(j).kineticLaw.localParameter(k).value;
        end
    end    
end

%% Remove Sink Reactions and Pi_Supply from S-matrix
% No thermodynamic driving force for these reactions

% First remove Biomass-metabolites from S-matrix (no metabolite
% concentrations)
S_No_Biomass = [];
for v = 1:length(N.species)
    if ~strncmp(N.species(v).name,'PPool',5) && ~strncmp(N.species(v).name,'BioM',4)
        S_No_Biomass(v,:) = SFullFull(v,:);
    end
end

% Remove Biomass reactions (no Keq value)
SFull_Mod = [];
for u = 1:length(N.reaction)
    if ~strncmp(N.reaction(u).name,'Sink',4) && ~strncmp(N.reaction(u).name,'Supply',6)
        SFull_Mod(:,u) = S_No_Biomass(:,u);
    end
end 

%% Sample metabolite concentrations

% Determine a random number in the range (max to min) of the NET-analysis
% concentrations, using 
% 'value = e^(ln(max) - ln(min)) * random(btwn 0 and 1) + ln(min)'
l=1;

for n = 1:NrSampling
MetConc = [];
    
    %%%---Create set of metabolite concentrations
    for m = 1:length(Y.RangeNames) % loop through metabolites
        
                % get random concentration in NET-range, logarithmic
                % distributed for concentrations go in dG formula with
                % logarithmic value!
        if ~any(regexp(Y.RangeNames{m},'/')) && ~any(regexp(Y.RangeNames{m},'Pool'))
            MetConc(m,1)=(exp((log(Y.MaxBound(m))-log(Y.MinBound(m)))*...
            rand(1)+log(Y.MinBound(m))))/1000; % Divide by 1000 to get from mM to M!
        end
        
    %%%---Check if sampled concentrations are within the ratios, applies to
    %%%---      NADPH/NADP, ATP/ADP, NADH/NAD
        % If there are speciefied rations, such as "NADPH/NADP"
        if any(regexp(Y.RangeNames{m},'/'))     
            % Split the ratio into the two species
            RatioBetween=regexp(Y.RangeNames{m},'/','split');
            
            Species1 = RatioBetween(1); % Name of species 1
            Species2 = RatioBetween(2); % Name of species 2
            
            % Get indexes to get corresponding concentration ranges
            SpeciesIndex1 = FindIndex(Y.RangeNames,RatioBetween(1));
            SpeciesIndex2 = FindIndex(Y.RangeNames,RatioBetween(2));
            
            
            n = 0; % Condition variable for while-loop
            
            while n < 1 % for as long as n = 0, see above
            
                % Calculate the ratio between the species from already
                % sampled concentration set
                Ratio = MetConc(SpeciesIndex1)/MetConc(SpeciesIndex2);
                
                % If the already sampled concentrations do not violate the
                % set ratio constraint
                if (Ratio >= Y.MinBound(m)) && (Ratio <= Y.MaxBound(m))
                  n = 1; % leave while loop
                
                % If the ratio between the 2 species is outside the allowed
                % range for the ratio, sample new concentrations for the
                % given species
                % This is repeated until the ratio is in the range and the
                % while loop can be left
                else
                    % New sampling for species 1
                    MetConc(SpeciesIndex1,1)=(exp((log(Y.MaxBound(...
                    SpeciesIndex1))-log(Y.MinBound(SpeciesIndex1)))*...
                    rand(1)+log(Y.MinBound(SpeciesIndex1))))/1000; % Divide by 1000 to get from mM to M!
                    % New sampling for species 2
                    MetConc(SpeciesIndex2,1)=(exp((log(Y.MaxBound(...
                    SpeciesIndex2))-log(Y.MinBound(SpeciesIndex2)))*...
                    rand(1)+log(Y.MinBound(SpeciesIndex2))))/1000; % Divide by 1000 to get from mM to M!
                end
            end
        end
    end

    %%%---Check for thermodynamic consistency of set from above
   
    % transpose S-matrix(watch out about full matrix or constant 
    % metabolites missing), loop through reactions, calculate 
    % dG = (S^T*lnMetConc+ln(Keq))*RT*-1
    % if any dG value >0, reject sample set!
    % Else save as part of 'MetConcDataSet' :)
    dG=[];
    SFullT = transpose(SFull_Mod);  % Transpose full S-matrix
    [r c] = size(SFullT);
    
    for p = 1:r                    % Loop through reactions, calculate dG 
                                   % for each reaction
        dG(p) = (SFullT(p,:)*-1*log(MetConc)+log(Y.RxnKeq(p)))*R*T*-1;     
    end

    % Find the indexes of PPool and Pi
    PPool_Index = FindIndex(Y.RangeNames,'PPool');
    P_i_Index = FindIndex(Y.RangeNames,'P_i');
    
    % Check feasibility by checking highest dG value of the set
    % Additionally check that the metabolite concentrations sum up to a
    % value below 100 mM
    if max(dG) < 0 && sum(MetConc*1000) < 100 % Multiply with 1000 to get from M to mM!
        MetConc(end+1,1)=MetConc(P_i_Index,1)*exp((log(Y.MaxBound(PPool_Index))-log(Y.MinBound(PPool_Index)))*rand(1)+log(Y.MinBound(PPool_Index)));
        
        dGDataSet(:,l) = dG;
        MetConcDataSet(:,l) = 1000*MetConc; % Multiply with 1000 to get from M to mM!
        l=l+1;    
    end
end 

toc % End Timer
%save('outfile.mat','MetConcDataSet')
end

function [n] = FindIndex(haystack, needle)
    n=find(ismember(haystack,needle));
end