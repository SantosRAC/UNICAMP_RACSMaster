%
% Philip suggested to try optimization of demand reaction for each of the
% biomass precursors. Then, check for each of these precursors which
% reactions have significant flux. Delete each of these reactions to try to
% have a better idea about the impact of these reactions in a possible loop.
%

% Clear workspace before doing anything else
clear

% Initialize COBRA Toolbox
initCobraToolbox();

% Import K. brasiliensis GHG001 metabolic model
model_SBML = readCbModel('KbrasiliensisGHG1_MetabolicModel.xml',1000,'SBML');

% Fixing a weird number of metaolite formula vector (formula vector has one missing element)
newVectorMetFormula = cell(1559,1);
model_SBML.metFormulas = newVectorMetFormula;

% Add exchange reactions (talk with Philip, maybe add exchange and sink reactions only after
% solving the current problem with infeasible loops in model and simulation fluxes)


% Select a solver (GLPK)
changeCobraSolver('glpk');

% Find column with biomass pseudo-reaction
colS_biomass = find(model_SBML.c);

% List metabolites in biomass pseudo-reaction
biomassMetabs = model_SBML.mets(model_SBML.S(:,colS_biomass)<0);

% Size of array/ vector with biomass precursors
%numel(biomassMetabs); % Number of array elements
%length(biomassMetabs); % Lenght of the vector or largest array dimension

k=1;
m=1;

% Add demand reactions corresponding to biomass precursors
% and create a new model with the incorporated demand reactions
[model_SBML_newDemand,addedRxns] = addDemandReaction(model_SBML,biomassMetabs);

%
for i=1:length(biomassMetabs)
    model_SBML_newDemand.c = zeros(length(model_SBML_newDemand.c),1); % CHANGE OBJECTIVE FUNCTION TO NEW DEMAND RXN
    model_SBML_newDemand.c(strmatch(addedRxns{i},model_SBML_newDemand.rxns)) = 1;         
    FBAsolution = optimizeCbModel(model_SBML_newDemand); % OPTIMIZE
    if FBAsolution.f < 1E-9 % MAKE LIST OF WHICH BIOMASS PRECURSORS ARE SYNTHESIZED (THEY SHOULD NOT)
        missingMets(k) = biomassMetabs(i);
        k = k+1;
    else
        % Getting flux vector for all reaction in model
        o=1;
        MetsWithWeirdFlux=[];
        for a=1:length(FBAsolution.x)
            if FBAsolution.x(a) > 1E-9
                display(model_SBML_newDemand.rxns(a));
                display(FBAsolution.x(a));
                MetsWithWeirdFlux(o) = FBAsolution.x(a);
                o = o+1;
            end
        end
        for b=1:length(MetsWithWeirdFlux)
            mod_model_singleDeletion = model_SBML_newDemand;
            singleRxnDeletion(mod_model_singleDeletion,'FBA',MetsWithWeirdFlux(o));
            % Run FBA here again for this biomass precursor and check other
            % fluxes (affected by reaction deletion?)
            
        end
        presentMets(m) = biomassMetabs(i);
        m = m+1;
        return
    end
end