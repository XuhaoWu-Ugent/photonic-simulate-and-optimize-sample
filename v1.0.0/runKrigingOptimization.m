%> @file "runKrigingOptimization.m"
%> @authors SUMO Lab Team
%> @version 2014a
%> @date Copyright 2006-2014
%>
%> This file is part of the Surrogate Modeling Toolbox ("SUMO Toolbox")
%> and you can redistribute it and/or modify it under the terms of the
%> GNU Affero General Public License version 3 as published by the
%> Free Software Foundation.  With the additional provision that a commercial
%> license must be purchased if the SUMO Toolbox is used, modified, or extended
%> in a commercial setting. For details see the included LICENSE.txt file.
%> When referring to the SUMO Toolbox please make reference to the corresponding
%> publication:
%>   - A Surrogate Modeling and Adaptive Sampling Toolbox for Computer Based Design
%>   D. Gorissen, K. Crombecq, I. Couckuyt, T. Dhaene, P. Demeester,
%>   Journal of Machine Learning Research,
%>   Vol. 11, pp. 2051-2055, July 2010.
%>
%> Contact : sumo@sumo.intec.ugent.be - http://sumo.intec.ugent.be
%> Signature
%>	 runKrigingOptimization()
%
% ======================================================================
%> @brief EGO optimization (kriging + EGO)
%>
%> Adapt the CHANGEME lines for your own problem.
% ======================================================================
function runKrigingOptimization()

close all;

%% configuration

% simulator
p = 'sumo-toolbox/src/matlab/playground/CoKrigingExamples/'; % CHANGEME: Path to simulator

if 0
    % Matlab
    %addpath(p); % add p to matlab path

    func = @cheapMath1D; % function handle to your matlab function
else
    % Python
    if count(py.sys.path,'') == 0
        insert(py.sys.path,int32(0),'');
    end
    
    mod = py.importlib.import_module('pythonSimulator');
    py.reload(mod);

    %func = @py.pythonSimulator.expensive; % function handle to your function
    func = @py.pythonSimulator.cheap; % function handle to your function

end

%x = [8.33268724207, 9.60086611679, 0.081202828722, 0.308432233019]*1e-6;
%the result of optimization by 3 wavelengths 
x = [8.3 17.6867186861 0.226970213952 0.531262189346]*1e-6;
y = func(x)
bounds =  [4.3,13.7,0.2,0.4; 12.3,21.7,0.3,0.6]*1e-6;
%bounds =  [8.3,9.6,0.1,0.3; 30.0,35.0,0.3,0.6]*1e-6;
%the BCs of optimization by 3 wavelengths

inDim = size(bounds,2); % number of input variables

transl = (bounds(2,:) + bounds(1,:))/2.0;
scale = (bounds(2,:) - bounds(1,:))/2.0;
[inFunc, outFunc] = calculateTransformationFunctions( [transl; scale] );

% general options
debug = true; % generate some extra plots
outputDir = 'results_egoKriging/';
nrMaxSamples = 250; % CHANGEME: maximum number of samples
nrInitialSamples = 50; % CHANGEME: number of initial points
distanceThreshold = 2.*eps;

% setup kriging model options
type = 'Kriging';
opts = feval([type '.getDefaultOptions'] );
opts.type = type;

theta0 = repmat(0.25,1,inDim);

lb = repmat(-2,1,inDim);
ub = repmat(2,1,inDim);

% CHANGEME: correlation function to use
%bf = BasisFunction( 'corrgauss', inDim, lb, ub, {'log'});
%bf = BasisFunction( 'correxp', inDim, lb, ub, {'log'});
%bf = BasisFunction( 'corrmatern32', inDim, lb, ub, {'log'});
bf = BasisFunction( 'corrmatern52', inDim, lb, ub, {'log'});

optimopts.derivativecheck = 'off';
optimopts.diagnostics = 'off';
optimopts.algorithm = 'active-set';
optimopts.maxFunEvals = 10000;
optimopts.maxIter = 500;
optimopts.gradobj = 'on';
opts.hpOptimizer = MatlabOptimizer( inDim, 1, optimopts );

% generate initial data sets
d = LatinHypercubeDesign( nrInitialSamples, inDim, 1 );
samples = d.generate();
simulatorSamples = outFunc(samples);
for i=1:size(simulatorSamples)
    values(i,:) = func( simulatorSamples(i,:) );
end

% select optimizer to use
candidatesPerSample = 25; % monte carlo

optimopts.gradobj = 'off';
optimopts.derivativecheck = 'off';
optimopts.diagnostics = 'off';
optimopts.algorithm = 'active-set';
optimizer = MatlabOptimizer(inDim, 1, optimopts);

%optimizer = MatlabPatternSearch(inDim, 1);
%optimizer = PSOtOptimizer(inDim, 1, 'maxiters', 200, 'popSize', 30, 'seedPSO', 0);

%optimopts.ep = 1e-6;
%optimopts.maxevals = 20000;
%optimopts.maxits = 1000;
%optimizer = DirectOptimizer(inDim, 1, optimopts);

optimizer = optimizer.setBounds(-ones(1,inDim), ones(1,inDim));

% candidateRankers to use
rankers = {expectedImprovement(inDim, 1), 'scaling', 'none' ...
    modelVariance(inDim, 1), 'scaling', 'none' };

%rankers = {knowledgeGradient(optimizer, inDim, 1) ...
%            maxvar(inDim, 1) };

%% main loop
nrSamples = size( samples, 1 );
nrIter = nrMaxSamples - sum(nrSamples, 1); % number of iterations

mkdir(outputDir);

modelPlot = [];
rankersPlots = [];

for i=1:nrIter
    
    fprintf('Iterion %i of %i (samples %i).\n', i, nrIter, nrSamples );
    
    % build and fit Kriging object
    k = KrigingModel( opts, theta0, 'regpoly0', bf, 'useLikelihood' );
    k = k.constructInModelSpace( samples, values );
    
    % optimize it
    %state.lastModels{1}{1} = k;
    state.lastModels{1}{1} = OutputFilterWrapper(k, 1);
    state.samples = samples;
    state.values = values;
        
    for j=1:length(rankers)
        
        rankers{j} = rankers{j}.initNewSamples(state);
        
        %% monte carlo
        nCandidates = nrSamples*candidatesPerSample;
        initialPopulation = rand(nCandidates, inDim) .* 2 - 1;
        
        %foundvalues = k.evaluateInModelSpace(initialPopulation);
        foundvalues = rankers{j}.score(initialPopulation, state);
        %[dummy, idx] = sort( foundvalues, 1, 'ascend' );
        [dummy, idx] = sort( foundvalues, 1, 'descend' );
        
        fprintf('\t- Generated %i candidates.\n', nCandidates );
        
        %% optimize best candidate
        
        % set initial population
        maxPopSize = optimizer.getPopulationSize();
        initialPopulation = initialPopulation(idx(1:maxPopSize,:),:);
        
        fprintf('\t- Optimizing best candidate...\n');
        %initialPopulation
        
        optimizer = optimizer.setInitialPopulation(initialPopulation);
        
        % give the state to the optimizer - might contain useful info such as # samples
        optimizer = optimizer.setState(state);
        
        optimFunc = @(x) rankers{j}.scoreMinimize(x, state);
        [dummy, xmin, fmin] = optimizer.optimize(optimFunc);
        
        dups = buildDistanceMatrix( xmin, samples, 1 );
        index = find(all(dups > distanceThreshold, 2));
        xmin = xmin(index,:);
        fmin = fmin(index,:);
        
        if ~isempty( xmin )
            break;
        end
        
    end
    
    if isempty( xmin )
        xmin = 2.*(rand(1,inDim) - 0.5);
        
        fprintf('\t- No unique samples found, randomizing.\n' );
    else
        fprintf('\t- Samples found using %s ranker.\n', class(rankers{j}));
    end
    
    %% evaluate new samples and add to set
    xmin = xmin(1,:);
    fmin = func( outFunc(xmin) );
    samples = [samples ; xmin(1,:)];
    values = [values ; fmin(1,:)];
    
    nrSamples = size( samples, 1 );
    
    %% nice plotjes
    if 0 %debug
        % sumo model
        plotOpts = struct('state', state, 'fixedSamples', xmin); % 'plotDerivatives', true);
        modelPlot = sliceContourPlot( modelPlot, k, plotOpts );
        title('Model plot');
        
        rankersPlots = sliceContourPlot( rankersPlots, rankers{j}, plotOpts );
        title('Ranker plot');
        pause(1); % wait 1 second
    end
    
    samplesSimulator = outFunc(samples);
    
    fprintf('\t- Best optimum so far:\n' );
    [fmin, idx] = min( values );
    xmin = samplesSimulator(idx,:)
    fmin
    
    %% save data
    save(fullfile(outputDir, 'samplesSimulator.mat'), 'samplesSimulator' );
    save(fullfile(outputDir, 'samples.mat'), 'samples' );
    save(fullfile(outputDir, 'values.mat'), 'values' );
    save(fullfile(outputDir, sprintf('sumo_model%i', nrSamples-1)), 'k');
end

disp('Finished');

%% final model + samples plot
%sliceContourPlot( modelPlot, k );

%% clean up
% remove p from matlab path
rmpath(p);

end