function [fvector,num_rxns_vector] = getNonZeroFluxReactions(inmodel)

model_SBML = inmodel;

% Choose solver (GLPK)
changeCobraSolver('glpk');

% Fix length of vector with metabolite formulas, if it different from the
% vector with metabolites
if length(model_SBML.mets) ~= length(model_SBML.metFormulas)
    newVectorMetFormula = cell(length(model_SBML.mets),1);
    model_SBML.metFormulas = newVectorMetFormula;
end

% Add demand reactions and close reactions (upper and lower bounds)
initVectorDM = length(model_SBML.rxns) + 1;
model_DM = addDemandReaction(model_SBML,model_SBML.mets);
endVectorDM = length(model_DM.rxns);
model_DM_closed = closeRxn(model_DM,model_DM.rxns(initVectorDM:end));

% Run FBA for each of the created demand reactions
fvector = zeros(endVectorDM,1);
num_rxns_vector = zeros(endVectorDM,1);

% Create file to output results of optimizations (demand reaction and other
% reactions with flux (that should not, given that the system has not
% exchange reactions)
fid=fopen('OutFileNonZeroFluxReactions_fixedAbsValue02262017.txt','w');

for i = initVectorDM:endVectorDM
    model_DM_test = model_DM_closed;
    model_DM_test.ub(i) = 1000;
    model_DM_test = changeObjective(model_DM_test,model_DM_test.rxns(i));
    soln = optimizeCbModel(model_DM_test,'max','one');
    num_rxns_vector(i) = length(find(abs(soln.x) > 0.001));
    vectorSolutionWithFluxes = find(abs(soln.x) > 0.001);
    if isempty(vectorSolutionWithFluxes)
       %sprintf('%s DMreactionWithoutReactionsWithFluxInOptm', char(model_DM_test.rxns(i)))
       fprintf(fid, '%s DMreactionWithoutReactionsWithFluxInOptm\n', char(model_DM_test.rxns(i)));
    else
        for nonzerofluxelm = 1:length(vectorSolutionWithFluxes)
            %sprintf('%s %s', char(model_DM_test.rxns(i)), char(model_DM_test.rxns(nonzerofluxelm)))
            %display(model_DM_test.rxns(vectorSolutionWithFluxes(nonzerofluxelm)))
            fprintf(fid, '%s %s\n', char(model_DM_test.rxns(i)), char(model_DM_test.rxns(vectorSolutionWithFluxes(nonzerofluxelm))));
        end
    end
    fvector(i) = soln.f;
end

fclose(fid);

end