%% Function to create dfodc-calculation file
% Markus Janasch, Ph.D. Student, KTH
% Created: 2018-02-22, last modified: -

% This function creates a function to calculate the dfodc for parameter
% sampling
function MJanasch_CreateFunction_CalDFODC(InputDataStructure)

    load(InputDataStructure)
    [ParID,ParVal] = ExtractParameters(N);
    
    
    fileID = fopen('MJanasch_CalDFODC.m','w');
    fprintf(fileID,'%10s\n\n','%% %% Function to calculate dfodc');
    fprintf(fileID,'%10s\n\n','function [dfodc,ParameterSet] = MJanasch_CalDFODC(N,ParameterSet,V_K_Indeces,Fluxes,SFull)');

    fprintf(fileID,'%10s\n\n','%% Define Metabolite Concentrations');


    for Species=1:length(N.species)
        fprintf(fileID,'\t%11s\t %1s%1.0f%1s\n',N.species(Species).name,'= N.species(',Species,').initialConcentration;');
    end    

    fprintf(fileID,'\n%10s\n','%% Define Parmeters');
    fprintf(fileID,'\n%10s\n','% Set Vmaxs and K to 1, if not already');

    fprintf(fileID,'\n\t%1s\n','for j = 1:length(V_K_Indeces)');
    fprintf(fileID,'\t\t%1s\n','ParameterSet(V_K_Indeces(j)) = 1;');
    fprintf(fileID,'\t%1s\n\n','end');

    RunVar = 1;
    for ParaNum = 1:length(ParID)
        fprintf(fileID,'\t%12s %1s%1.0f%1s\n',ParID{ParaNum},'= ParameterSet(',ParaNum,');');
        if (strncmp(ParID(ParaNum),'V_max',5)) || (strncmp(ParID(ParaNum),'K_PPool',7))
            V_K_Indeces(RunVar) = ParaNum;
            RunVar = RunVar + 1;
        end
    end    

    fprintf(fileID,'\n\n%10s\n\n','%% Rate Equations');
    for RxnNum = 1:length(N.reaction)
        fprintf(fileID,'\t%6s%1.0f%1s\t%1s%1s%1s\n','pippo(',RxnNum,',1)','= ',N.reaction(RxnNum).kineticLaw.formula,';');
    end    

    fprintf(fileID,'\n\n%10s\n','%% Calculate Vmax and K');
    fprintf(fileID,'\n\t%1s\n','for j = 1:length(N.reaction)');
    fprintf(fileID,'\t\t%1s\n','Fluxes(j);');
    fprintf(fileID,'\t\t%1s\n','Vmax_K = Fluxes(j)/pippo(j,1);');
    fprintf(fileID,'\t\t%1s\n','ParameterSet(V_K_Indeces(j)) = Vmax_K;');
    fprintf(fileID,'\t%1s\n\n','end');

    
    
    
    
    
    
    for VMaxNum = 1:length(V_K_Indeces)
        fprintf(fileID,'\t%10s%2s%1.0f%1s\n',ParID{V_K_Indeces(VMaxNum)},' = ParameterSet(V_K_Indeces(',VMaxNum,'));');   
    end

    fprintf(fileID,'\n\n%10s\n','%% Derivatives for DFODC matrix (Numerical values of the derivatives)');
    fprintf(fileID,'\t%1s\n','S = size(SFull,1);');
    fprintf(fileID,'\t%1s\n','R = size(SFull,2);');

    fprintf(fileID,'\n%10s\n','% Initialize dfodc matrix');
    fprintf(fileID,'\t%1s\n\n','clear dfodc;');


%% Calculate derivatives
    R = {N.reaction.name};
    S = {N.species.id};

    % Transform S cell array into symbolic variables
    Ssym=sym(S);
    DER = {}; 
    % This structure will contain the DERIVATIVES of reaction j with respect to metabolite i.
    % Computation of the derivatives

for j = 1:length(R)
    for i = 1:length(S)
        if ~(N.species(i).boundaryCondition)
            isitthere = IsSpeciesInReaction(N.species(i).name,N.reaction(j));
            if isitthere ~= 0
                DER{i,j} = diff(N.reaction(j).kineticLaw.formula,Ssym(i));
            else
                DER{i,j}='0';
            end
        end
    end
end

fprintf(fileID,'\n%10s\n','% Fill dfodc');


for RNum = 1:size(DER,2)
    for SNum = 1:size(DER,1)
        DER_Entry=char(DER{SNum,RNum});
       fprintf(fileID,'%12s%1.0f%1s%1.0f%1s\t%1s%1s%1s\n','dfodc(',RNum,',',SNum,')',' = ',DER_Entry,';'); 
    end
end


fprintf(fileID,'%1s','end');
fclose(fileID);
end

function [ParID_New,ParVal_New] = ExtractParameters(N)
    Running_Variable=1;
    for a = 1:length(N.reaction)
        for b=1:length(N.reaction(a).kineticLaw.parameter)
            ParID_New{Running_Variable}     = N.reaction(a).kineticLaw.parameter(b).name;
            ParVal_New(Running_Variable)    = N.reaction(a).kineticLaw.parameter(b).value;
            Running_Variable=Running_Variable+1;
        end
    end
end

% find species in this reaction
function [y] = IsSpeciesInReaction(Species,Reaction)
    y = 0;
    for c = 1:length(Reaction.product)
        if (strcmp(Species,Reaction.product(c).species))
            y = y + 1;
        end;
    end;

    for c = 1:length(Reaction.reactant)
        if (strcmp(Species, Reaction.reactant(c).species))
            y = y + 1;
        end;
    end;

    for c = 1:length(Reaction.modifier)
        if (strcmp(Species, Reaction.modifier(c).species))
            y = y + 1;
        end;
    end;
end