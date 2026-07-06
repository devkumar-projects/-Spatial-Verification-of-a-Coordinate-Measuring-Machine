function results = Metrology_Correction(config)
% METROLOGY_CORRECTION
% Robust correction and validation pipeline for CMM / Machine Checking Gauge data.
%
% The input files must contain at least six numeric columns:
%   X  Y  Z  NX  NY  NZ
%
% where [X Y Z] are point coordinates and [NX NY NZ] are probe directions.
%
% This implementation deliberately avoids correcting the data with a simple
% moving average. Instead, it:
%   1. validates and normalises the measurement data;
%   2. estimates physically interpretable geometric errors:
%        - centre offset [dX dY dZ];
%        - arm/radius error dR;
%   3. uses robust iteratively reweighted least squares (Huber weights);
%   4. separates systematic bias from directional hysteresis when a reverse
%      measurement file is provided;
%   5. reconstructs corrected points along the theoretical probe normals;
%   6. refits the corrected circle or sphere;
%   7. reports form error, RMS residual, bias, repeatability and hysteresis;
%   8. writes corrected point files, a metrics CSV file and diagnostic plots.
%
% Example:
%   cfg = struct;
%   cfg.theoreticalFile = 'theorique.txt';
%   cfg.forwardFile     = 'experimentale.txt';
%   cfg.reverseFile     = ''; % e.g. 'experimentale_retour.txt'
%   cfg.outputDirectory = 'metrology_output';
%   results = Metrology_Correction(cfg);
%
% Notes:
% - Units are preserved from the input files. The default label is millimetres.
% - The point order must be identical in theoretical and experimental files.
% - A reverse file is strongly recommended for a rigorous hysteresis analysis.
%
% Compatible with MATLAB R2019b or later.

    if nargin < 1 || isempty(config)
        config = struct;
    end

    config = applyDefaults(config);
    ensureOutputDirectory(config.outputDirectory);

    [Pth, Nth] = loadPointSet(config.theoreticalFile);
    [Pfw, Nfw] = loadPointSet(config.forwardFile);

    validateCorrespondence(Pth, Pfw, ...
        config.theoreticalFile, config.forwardFile);

    Nth = normaliseRows(Nth);
    Nfw = normaliseRows(Nfw);

    hasReverse = ~isempty(config.reverseFile) && isfile(config.reverseFile);

    if hasReverse
        [Prv, Nrv] = loadPointSet(config.reverseFile);
        validateCorrespondence(Pth, Prv, ...
            config.theoreticalFile, config.reverseFile);
        Nrv = normaliseRows(Nrv);
    else
        Prv = [];
        Nrv = [];
    end

    geometryType = classifyGeometry(Pth, config.planarityThreshold);

    % Normal-direction deviations relative to the theoretical trajectory.
    eForward = rowDot(Pfw - Pth, Nth);

    if hasReverse
        eReverse = rowDot(Prv - Pth, Nth);

        % Direction-independent part and direction-dependent hysteresis.
        eIdentification = 0.5 * (eForward + eReverse);
        hysteresis = 0.5 * (eForward - eReverse);
    else
        eReverse = [];
        eIdentification = eForward;
        hysteresis = [];
    end

    % Physical geometric model:
    % e_i ≈ n_i' * [dX dY dZ]' + dR
    A = [Nth, ones(size(Nth, 1), 1)];
    [beta, robustWeights, modelResidual] = robustLinearFit( ...
        A, eIdentification, config.huberK, ...
        config.maxIterations, config.tolerance);

    centreOffset = beta(1:3);
    radiusError = beta(4);
    systematicCorrection = A * beta;

    % Correct only the identified systematic geometric component.
    PfwCorrected = Pfw - systematicCorrection .* Nth;

    if hasReverse
        PrvCorrected = Prv - systematicCorrection .* Nth;
    else
        PrvCorrected = [];
    end

    % Rebuild normals from corrected geometry rather than filtering normal
    % components independently.
    if strcmp(geometryType, 'circle')
        fitTheory = fitCircle3DRobust(Pth, config);
        fitForwardBefore = fitCircle3DRobust(Pfw, config);
        fitForwardAfter = fitCircle3DRobust(PfwCorrected, config);
        NfwCorrected = circleNormals(PfwCorrected, fitForwardAfter);

        if hasReverse
            fitReverseBefore = fitCircle3DRobust(Prv, config);
            fitReverseAfter = fitCircle3DRobust(PrvCorrected, config);
            NrvCorrected = circleNormals(PrvCorrected, fitReverseAfter);
        else
            fitReverseBefore = [];
            fitReverseAfter = [];
            NrvCorrected = [];
        end
    else
        fitTheory = fitSphereRobust(Pth, config);
        fitForwardBefore = fitSphereRobust(Pfw, config);
        fitForwardAfter = fitSphereRobust(PfwCorrected, config);
        NfwCorrected = sphereNormals(PfwCorrected, fitForwardAfter.centre);

        if hasReverse
            fitReverseBefore = fitSphereRobust(Prv, config);
            fitReverseAfter = fitSphereRobust(PrvCorrected, config);
            NrvCorrected = sphereNormals(PrvCorrected, fitReverseAfter.centre);
        else
            fitReverseBefore = [];
            fitReverseAfter = [];
            NrvCorrected = [];
        end
    end

    eForwardCorrected = rowDot(PfwCorrected - Pth, Nth);

    if hasReverse
        eReverseCorrected = rowDot(PrvCorrected - Pth, Nth);
    else
        eReverseCorrected = [];
    end

    metrics = buildMetrics( ...
        eForward, eForwardCorrected, ...
        eReverse, eReverseCorrected, ...
        hysteresis, fitForwardBefore, fitForwardAfter, ...
        geometryType, config.units);

    % Write corrected files.
    forwardOutput = fullfile(config.outputDirectory, ...
        [config.outputPrefix '_forward_corrected.txt']);
    writematrix([PfwCorrected, NfwCorrected], forwardOutput, ...
        'Delimiter', 'tab');

    if hasReverse
        reverseOutput = fullfile(config.outputDirectory, ...
            [config.outputPrefix '_reverse_corrected.txt']);
        writematrix([PrvCorrected, NrvCorrected], reverseOutput, ...
            'Delimiter', 'tab');
    else
        reverseOutput = '';
    end

    metricsOutput = fullfile(config.outputDirectory, ...
        [config.outputPrefix '_metrics.csv']);
    writetable(metrics, metricsOutput);

    parameters = table( ...
        centreOffset(1), centreOffset(2), centreOffset(3), radiusError, ...
        rmsValue(modelResidual), max(abs(modelResidual)), ...
        mean(robustWeights), ...
        'VariableNames', { ...
        'CentreOffsetX', 'CentreOffsetY', 'CentreOffsetZ', ...
        'RadiusError', 'ModelResidualRMS', ...
        'ModelResidualMaximumAbsolute', 'MeanRobustWeight'});

    parametersOutput = fullfile(config.outputDirectory, ...
        [config.outputPrefix '_identified_parameters.csv']);
    writetable(parameters, parametersOutput);

    if config.generateFigures
        generateDiagnosticFigures( ...
            eForward, eForwardCorrected, ...
            eReverse, eReverseCorrected, hysteresis, ...
            fitForwardBefore, fitForwardAfter, ...
            geometryType, config);
    end

    results = struct;
    results.configuration = config;
    results.geometryType = geometryType;
    results.centreOffset = centreOffset;
    results.radiusError = radiusError;
    results.robustWeights = robustWeights;
    results.modelResidual = modelResidual;
    results.theoreticalFit = fitTheory;
    results.forwardFitBefore = fitForwardBefore;
    results.forwardFitAfter = fitForwardAfter;
    results.reverseFitBefore = fitReverseBefore;
    results.reverseFitAfter = fitReverseAfter;
    results.metrics = metrics;
    results.forwardCorrectedFile = forwardOutput;
    results.reverseCorrectedFile = reverseOutput;
    results.metricsFile = metricsOutput;
    results.parametersFile = parametersOutput;

    fprintf('\nMetrology correction completed successfully.\n');
    fprintf('Detected geometry : %s\n', geometryType);
    fprintf('Centre offset     : [%+.6g  %+.6g  %+.6g] %s\n', ...
        centreOffset(1), centreOffset(2), centreOffset(3), config.units);
    fprintf('Radius error      : %+.6g %s\n', radiusError, config.units);
    fprintf('Corrected file    : %s\n', forwardOutput);
    fprintf('Metrics file      : %s\n\n', metricsOutput);
end


function config = applyDefaults(config)
    defaults = struct;
    defaults.theoreticalFile = 'theorique.txt';
    defaults.forwardFile = 'experimentale.txt';
    defaults.reverseFile = '';
    defaults.outputDirectory = 'metrology_output';
    defaults.outputPrefix = 'mcg';
    defaults.units = 'mm';
    defaults.huberK = 1.345;
    defaults.maxIterations = 100;
    defaults.tolerance = 1e-12;
    defaults.planarityThreshold = 1e-5;
    defaults.generateFigures = true;

    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(config, names{k}) || isempty(config.(names{k}))
            config.(names{k}) = defaults.(names{k});
        end
    end
end


function ensureOutputDirectory(directory)
    if ~exist(directory, 'dir')
        mkdir(directory);
    end
end


function [P, N] = loadPointSet(filename)
    if ~isfile(filename)
        error('MetrologyCorrection:MissingFile', ...
            'Input file not found: %s', filename);
    end

    M = readmatrix(filename);
    M = M(any(isfinite(M), 2), :);

    if size(M, 2) < 6
        error('MetrologyCorrection:InvalidFormat', ...
            ['File "%s" must contain at least six numeric columns: ' ...
             'X Y Z NX NY NZ.'], filename);
    end

    M = M(:, 1:6);

    if any(~isfinite(M), 'all')
        error('MetrologyCorrection:InvalidData', ...
            'File "%s" contains NaN or infinite values.', filename);
    end

    P = M(:, 1:3);
    N = M(:, 4:6);

    if size(P, 1) < 6
        error('MetrologyCorrection:InsufficientPoints', ...
            'At least six measurement points are required.');
    end
end


function validateCorrespondence(P1, P2, file1, file2)
    if size(P1, 1) ~= size(P2, 1)
        error('MetrologyCorrection:PointCountMismatch', ...
            ['Files "%s" and "%s" do not contain the same number ' ...
             'of measurement points.'], file1, file2);
    end
end


function N = normaliseRows(N)
    norms = sqrt(sum(N.^2, 2));

    if any(norms < eps)
        error('MetrologyCorrection:ZeroNormal', ...
            'At least one probe normal has zero magnitude.');
    end

    N = N ./ norms;
end


function values = rowDot(A, B)
    values = sum(A .* B, 2);
end


function geometryType = classifyGeometry(P, threshold)
    centred = P - mean(P, 1);
    singularValues = svd(centred, 'econ');

    if numel(singularValues) < 3
        geometryType = 'circle';
        return;
    end

    relativePlanarity = singularValues(end) / max(singularValues(1), eps);

    if relativePlanarity < threshold
        geometryType = 'circle';
    else
        geometryType = 'sphere';
    end
end


function [beta, weights, residual] = robustLinearFit(A, b, huberK, maxIter, tol)
    beta = pinv(A) * b;
    weights = ones(size(b));

    for iteration = 1:maxIter
        residual = b - A * beta;
        scale = robustScale(residual);

        u = abs(residual) / max(huberK * scale, eps);
        weights = ones(size(u));
        mask = u > 1;
        weights(mask) = 1 ./ u(mask);

        Aw = A .* sqrt(weights);
        bw = b .* sqrt(weights);
        betaNew = pinv(Aw) * bw;

        if norm(betaNew - beta) <= tol * (1 + norm(beta))
            beta = betaNew;
            break;
        end

        beta = betaNew;
    end

    residual = b - A * beta;
end


function value = robustScale(x)
    centred = x - median(x);
    value = 1.4826 * median(abs(centred));

    if value < eps
        value = std(x);
    end

    if value < eps
        value = 1;
    end
end


function fit = fitSphereRobust(P, config)
    % Algebraic initialisation.
    A = [2 * P, ones(size(P, 1), 1)];
    b = sum(P.^2, 2);
    q = pinv(A) * b;

    centre = q(1:3);
    radius = sqrt(max(q(4) + dot(centre, centre), eps));

    for iteration = 1:config.maxIterations
        D = P - centre.';
        distances = sqrt(sum(D.^2, 2));
        residual = distances - radius;

        weights = huberWeights(residual, config.huberK);
        J = [-D ./ max(distances, eps), -ones(size(P, 1), 1)];

        Jw = J .* sqrt(weights);
        rw = residual .* sqrt(weights);
        delta = -pinv(Jw) * rw;

        centreNew = centre + delta(1:3);
        radiusNew = radius + delta(4);

        if norm([centreNew - centre; radiusNew - radius]) <= ...
                config.tolerance * (1 + norm([centre; radius]))
            centre = centreNew;
            radius = radiusNew;
            break;
        end

        centre = centreNew;
        radius = radiusNew;
    end

    D = P - centre.';
    radialResidual = sqrt(sum(D.^2, 2)) - radius;

    fit = struct;
    fit.type = 'sphere';
    fit.centre = centre;
    fit.radius = radius;
    fit.residual = radialResidual;
    fit.formError = max(radialResidual) - min(radialResidual);
    fit.rmsResidual = rmsValue(radialResidual);
end


function fit = fitCircle3DRobust(P, config)
    centroid = mean(P, 1);
    centred = P - centroid;

    [~, ~, V] = svd(centred, 'econ');
    basis = V(:, 1:2);
    planeNormal = V(:, 3);

    Q = centred * basis;
    x = Q(:, 1);
    y = Q(:, 2);

    % Algebraic initialisation in the best-fit plane.
    A = [2 * x, 2 * y, ones(size(x))];
    b = x.^2 + y.^2;
    q = pinv(A) * b;

    centre2D = q(1:2);
    radius = sqrt(max(q(3) + dot(centre2D, centre2D), eps));

    for iteration = 1:config.maxIterations
        D = Q - centre2D.';
        distances = sqrt(sum(D.^2, 2));
        residual = distances - radius;

        weights = huberWeights(residual, config.huberK);
        J = [-D ./ max(distances, eps), -ones(size(Q, 1), 1)];

        Jw = J .* sqrt(weights);
        rw = residual .* sqrt(weights);
        delta = -pinv(Jw) * rw;

        centreNew = centre2D + delta(1:2);
        radiusNew = radius + delta(3);

        if norm([centreNew - centre2D; radiusNew - radius]) <= ...
                config.tolerance * (1 + norm([centre2D; radius]))
            centre2D = centreNew;
            radius = radiusNew;
            break;
        end

        centre2D = centreNew;
        radius = radiusNew;
    end

    centre3D = centroid.' + basis * centre2D;
    D = Q - centre2D.';
    radialResidual = sqrt(sum(D.^2, 2)) - radius;

    fit = struct;
    fit.type = 'circle';
    fit.centre = centre3D;
    fit.radius = radius;
    fit.planeNormal = planeNormal;
    fit.basis = basis;
    fit.centroid = centroid.';
    fit.residual = radialResidual;
    fit.formError = max(radialResidual) - min(radialResidual);
    fit.rmsResidual = rmsValue(radialResidual);
end


function weights = huberWeights(residual, huberK)
    scale = robustScale(residual);
    u = abs(residual) / max(huberK * scale, eps);

    weights = ones(size(u));
    mask = u > 1;
    weights(mask) = 1 ./ u(mask);
end


function N = sphereNormals(P, centre)
    N = P - centre.';
    N = normaliseRows(N);
end


function N = circleNormals(P, fit)
    relative = P - fit.centre.';
    projected = relative - (relative * fit.planeNormal) * fit.planeNormal.';
    N = normaliseRows(projected);
end


function metrics = buildMetrics( ...
        eForward, eForwardCorrected, ...
        eReverse, eReverseCorrected, ...
        hysteresis, fitBefore, fitAfter, geometryType, units)

    metricName = { ...
        'Forward mean normal error'; ...
        'Forward RMS normal error'; ...
        'Forward peak-to-valley normal error'; ...
        'Forward standard deviation'; ...
        'Forward repeatability-based expanded uncertainty (k=2)'; ...
        [upperFirst(geometryType) ' form error']; ...
        [upperFirst(geometryType) ' fit RMS residual']};

    before = [ ...
        mean(eForward); ...
        rmsValue(eForward); ...
        peakToValley(eForward); ...
        std(eForward, 0); ...
        2 * std(eForward, 0); ...
        fitBefore.formError; ...
        fitBefore.rmsResidual];

    after = [ ...
        mean(eForwardCorrected); ...
        rmsValue(eForwardCorrected); ...
        peakToValley(eForwardCorrected); ...
        std(eForwardCorrected, 0); ...
        2 * std(eForwardCorrected, 0); ...
        fitAfter.formError; ...
        fitAfter.rmsResidual];

    if ~isempty(eReverse)
        metricName = [metricName; { ...
            'Reverse mean normal error'; ...
            'Reverse RMS normal error'; ...
            'Reverse peak-to-valley normal error'; ...
            'Bidirectional hysteresis RMS'; ...
            'Bidirectional hysteresis maximum absolute'}];

        before = [before; ...
            mean(eReverse); ...
            rmsValue(eReverse); ...
            peakToValley(eReverse); ...
            rmsValue(hysteresis); ...
            max(abs(hysteresis))];

        after = [after; ...
            mean(eReverseCorrected); ...
            rmsValue(eReverseCorrected); ...
            peakToValley(eReverseCorrected); ...
            rmsValue(hysteresis); ...
            max(abs(hysteresis))];
    end

    improvementPercent = 100 * (before - after) ./ max(abs(before), eps);
    unitColumn = repmat({units}, numel(metricName), 1);

    metrics = table(metricName, before, after, improvementPercent, unitColumn, ...
        'VariableNames', { ...
        'Metric', 'BeforeCorrection', 'AfterCorrection', ...
        'ImprovementPercent', 'Unit'});
end


function generateDiagnosticFigures( ...
        eForward, eForwardCorrected, ...
        eReverse, eReverseCorrected, hysteresis, ...
        fitBefore, fitAfter, geometryType, config)

    fig1 = figure('Visible', 'off', 'Name', 'Normal residuals');
    plot(eForward, 'LineWidth', 1.1);
    hold on;
    plot(eForwardCorrected, 'LineWidth', 1.1);
    grid on;
    xlabel('Measurement point index');
    ylabel(['Normal deviation (' config.units ')']);
    title('Normal-direction residuals before and after correction');
    legend('Before correction', 'After correction', 'Location', 'best');
    saveFigure(fig1, fullfile(config.outputDirectory, ...
        [config.outputPrefix '_normal_residuals.png']));
    close(fig1);

    fig2 = figure('Visible', 'off', 'Name', 'Form residuals');
    plot(fitBefore.residual, 'LineWidth', 1.1);
    hold on;
    plot(fitAfter.residual, 'LineWidth', 1.1);
    grid on;
    xlabel('Measurement point index');
    ylabel([upperFirst(geometryType) ' radial residual (' config.units ')']);
    title([upperFirst(geometryType) ' fit residuals']);
    legend('Before correction', 'After correction', 'Location', 'best');
    saveFigure(fig2, fullfile(config.outputDirectory, ...
        [config.outputPrefix '_form_residuals.png']));
    close(fig2);

    if ~isempty(eReverse)
        fig3 = figure('Visible', 'off', 'Name', 'Bidirectional comparison');
        plot(eForward, 'LineWidth', 1.0);
        hold on;
        plot(eReverse, 'LineWidth', 1.0);
        plot(eForwardCorrected, '--', 'LineWidth', 1.0);
        plot(eReverseCorrected, '--', 'LineWidth', 1.0);
        grid on;
        xlabel('Measurement point index');
        ylabel(['Normal deviation (' config.units ')']);
        title('Forward and reverse measurement comparison');
        legend('Forward before', 'Reverse before', ...
               'Forward after', 'Reverse after', ...
               'Location', 'best');
        saveFigure(fig3, fullfile(config.outputDirectory, ...
            [config.outputPrefix '_bidirectional_comparison.png']));
        close(fig3);

        fig4 = figure('Visible', 'off', 'Name', 'Hysteresis');
        plot(hysteresis, 'LineWidth', 1.1);
        yline(0, '--');
        grid on;
        xlabel('Measurement point index');
        ylabel(['Half-difference hysteresis (' config.units ')']);
        title('Direction-dependent hysteresis');
        saveFigure(fig4, fullfile(config.outputDirectory, ...
            [config.outputPrefix '_hysteresis.png']));
        close(fig4);
    end
end


function saveFigure(fig, filename)
    try
        exportgraphics(fig, filename, 'Resolution', 220);
    catch
        saveas(fig, filename);
    end
end


function value = peakToValley(x)
    value = max(x) - min(x);
end


function value = rmsValue(x)
    value = sqrt(mean(x.^2));
end


function text = upperFirst(text)
    text(1) = upper(text(1));
end
