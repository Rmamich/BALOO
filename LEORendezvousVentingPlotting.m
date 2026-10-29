%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot BALOO -- Earth Centered Venting and Maneuvering Scanario
%
%
% Created: 07/29/25
% Last Modified: 07/29/25
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

saveFig = true;
% blueCB = [68 , 119, 170]/255;
redCB = [238, 102, 119;]/255;
% purpleCB = [170, 51 , 119]/255;
greenCB = [34 , 136, 51 ]/255;
magentaCB = [170, 68, 153]/255; % BFF
darkGreenCB =  [17, 119, 51]/255; % iBFF
cyanCB = [102, 204, 238]/255; %
figLoc = [10, 10, 1850, 900];
grayCB = [0.7, 0.7, 0.7];

flagPlotTitle = true;                                       % true to include plot title
plotCheck = false;                                          % plots created from truth data as a check
baseFilepath = fullfile('.', 'Figures');                    % Sets base file path for Figure saving
baseName = 'simpleSim';                                     % Sets base name for Figure saving
timeUnit = '[s]';                                           % Sets time units for plotting
angUnit = '[deg]';                                          % Sets angular units for plotting
distUnit = '[km]';                                          % Sets linear distance units for plotting
plotStates = true;                                                    % Sets flag to plot MC for states or not
flagCovPlotOff = false;                                     % Truns cov plotting Cov on when false
plotMeas = true;                                            % true plots measurements
fontSize = 38;                                              % Sets font size for plotting script plots
line_width = 4;                                             % Sets line width for plotting script plots
textInterpreter = 'latex';                                  % Sets text interpreter for plotting script plots
axesLineWidth = 2;                                          % Sets axes width for plotting script plots
gridLineWidth = 1;                                          % Sets grid line width for plotting script plots
markerSize = 15;                                             % Sets marker size for plotting script plots
set(0, 'DefaultTextInterpreter', 'latex')                   % Sets text interpreter for plotting script plots
set(0, 'DefaultLegendInterpreter', 'latex')                 % Sets text interpreter for plotting script plots
set(0, 'DefaultTextFontSize', fontSize)                     % Sets font size for plotting script plots
set(0, 'DefaultLineLineWidth', line_width)                  % Sets line width for plotting script plots
set(0, 'DefaultAxesTickLabelInterpreter', textInterpreter)  % Sets text interpreter for plotting script plots
% flagZoomSave = 'zoomAll';                                   % Sets x limits for plotting script plots -- 'none' is default, 'zoomAll' enforces both limits, 'zoomEnd' enforces the upper limit
% xlimSet = [1000, 3000]; % NRHO post burn window
flagZoomSave = 'zoomEnd';
% xlimSet = [480, 700]; % NRHO burn window
xlimSet = [100, 300];                                      % Sets x limits for plotting script plots
% xlimSet = [780, 1000]; % NRHO bigger burn window

% xlimSet = [330, 530]; % NRHO burn window

% TODO: Add this flag to all scenario IC files
% New Flag for vBar compensation w/o accel state
flagCompensateMean = false;

%% Load in Data

load(fullfile('./MAT/EarthRel_venting_burnRendezvous_fixed_interacting_2modelsNoiseDifference/InitCovFactor/none/EKF', 'timeMat.mat'))

fileNameSaveCells_BALOO = './MAT/EarthRel_venting_burnRendezvous_BALOO/InitCovFactor/none/EKF/accelOnBiasOnLidarBiasOffj2j3j4Onj2j3j4On_noisyAccelerometer/fullCells';
fileNameSaveCells_6model = './MAT/EarthRel_venting_burnRendezvous_fixed_interacting/InitCovFactor/none/EKF/accelOnBiasOnLidarBiasOffj2j3j4Onj2j3j4On_noisyAccelerometer/fullCells';
fileNameSaveCells_2model = './MAT/EarthRel_venting_burnRendezvous_fixed_interacting_2modelsNoiseDifference/InitCovFactor/none/EKF/accelOnBiasOnLidarBiasOffj2j3j4Onj2j3j4On_noisyAccelerometer/fullCells';

fileNameSaveCells = {fileNameSaveCells_BALOO, fileNameSaveCells_6model, fileNameSaveCells_2model};

xErrAll_new = cell(3, 1);
covAllMCPlus_new = NaN(18, 3001, 3);
errAll_new = NaN(18, 3001, 3);
sigmaAll_new = NaN(18, 3001, 3);

for i = 1:3
    load(fullfile(fileNameSaveCells{i}, 'xErr.mat')); % , "xErrAll")
    load(fullfile(fileNameSaveCells{i}, 'covAllMCPlus.mat')); % , "covAllMCPlus")
    load(fullfile(fileNameSaveCells{i}, 'errAll.mat')); %, "errAll")
    load(fullfile(fileNameSaveCells{i}, 'sigmaAll.mat')); % , "sigmaAll")

    xErrAll_new{i} = xErrAll; % All err lines
    covAllMCPlus_new(:, :, i) = covAllMCPlus; % avg est 3 sigma
    errAll_new(:, :, i) = errAll; % avg err
    sigmaAll_new(:, :, i) = sigmaAll;


end


titles_states = {'BALOO', '6 Model IMM', '2 Model IMM'};
axisUnits_states = {'[$m$]', '[$m/s$]', '[$m/s^2$]', '[$m$]', '[$m/s$]', '[$m/s^2$]'};
figNames_states = {'RelPos', 'RelVel', 'RelAccel', 'ChaserPos', 'ChaserVel', 'ChaserBias'};
accelName = 'accelOnBiasOn';

nStateFigs = 6;
% truthState = truthState(1:18, :);
flagPlotFrame = 'N';


%% Measurement Flag Handling

% [m, ~] = size(meas);



titlesRes = {'Range Residuals v Time', 'Relative Azimuth Residuals v Time', 'Relative Elevation Residuals v Time'};
titlesMeas = {'Range v Time', 'Relative Azimuth v Time', 'Relative Elevation v Time'};
figNameMeas = {'Rng', 'RelAz', 'RelEl'};
resUnits = {'Range [$m$]', 'Az ', 'El '};

if contains(flagPlotFrame, 'rtn')
    baseAxis = {'R', 'T', 'N'};
    baseAxis_accel = {'X', 'Y', 'Z'};
else
    baseAxis = {'X', 'Y', 'Z'};
end

figNameSave = './Figures/EarthRel_venting_burnRendezvous/';


if ~exist(figNameSave, 'dir')
    mkdir(figNameSave)
end



if strcmp(flagZoomSave, 'zoomEnd')
    figNameEnd = '_zoomedEnd.png';
elseif strcmp(flagZoomSave, 'zoomAll')
    figNameEnd  ='_zoomedAll.png';
else
    figNameEnd = '.png';
end


%% Data Proceeing for Plotting
% n = sum(n_elm);
N = length(time);
% fileNameSaveCells = fullfile(fileBase, strcat(accelName, flagGravModelTruth, flagGravModel, burnName), 'fullCells');
nMC = length(xErrAll);
plotList = 1:nMC;



targetStates = 1:9;
deputyStates = 10:18;



if time(1) ~= 0
    time = time - time(1);
end

% if contains(flagScenario, 'nonDim')
%     time = time*tCharacteristic;
% end

if strcmp(timeUnit, '[h]')
    time = time / 3600; % Convert to Hours
elseif strcmp(timeUnit, '[min]')
    time = time / 60; % Convert to min
end



%% Plot Results
countAx = 1;
figCount = 1;
count = 1;
count1 = 1;
fig = [];



%%

%% Relativ Position -- not linked y 
iiStatePlot = 1; % unmodeled target acceleration
fig(figCount) = figure(figCount);
iiState = 2; % Y axis
% ax = [];

for j = 1:3
    ax(countAx) = subplot(1, 3, j);
    ax(countAx).FontSize = fontSize;
    ax(countAx).GridLineWidth = gridLineWidth;
    ax(countAx).LineWidth = axesLineWidth;
    set(gcf, 'Position', figLoc)
    grid minor
    hold on
    box on
    countAx = countAx + 1;
    xErrAll = xErrAll_new{j};

    for iiMC = plotList
        xErr = xErrAll{iiMC};
        p1 = plot(time, xErr(iiStatePlot*3-(3-iiState), :), 'color', grayCB, 'LineWidth',1);
    end
    plotLabelLines = p1;
    if flagPlotTitle
        title(titles_states{j}, 'FontSize', fontSize)
    end
    xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)

    % grid minor

    ylabel(strcat('$\mathbf{r}_{t/c, y}$ Error', ' ', axisUnits_states{iiStatePlot}),  'FontSize',fontSize)
    % if ~flagCovPlotOff
        p4 = plot(time, errAll_new(iiStatePlot*3-(3-iiState), :, j), 'k', 'LineWidth', line_width + 0.5);
        p5 = plot(time, 3*sigmaAll_new(iiStatePlot*3-(3-iiState), :, j), ':', 'Color', darkGreenCB, 'LineWidth', line_width);
        plot(time, -3*sigmaAll_new(iiStatePlot*3-(3-iiState), :, j), ':', 'Color', darkGreenCB, 'LineWidth', line_width)
        p3 = plot(time, 3*covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j), '--', 'Color', magentaCB, 'LineWidth', line_width);
        plot(time, -3*covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j), '--', 'Color', magentaCB, 'LineWidth', line_width)


    % end
    if strcmp(flagZoomSave, 'zoomEnd')
        xlim([xlimSet(1), xlimSet(2)])
    elseif strcmp(flagZoomSave, 'zoomAll')
        xlim(xlimSet)
    end

    yl = ylim;
    p6 = plot([121, 121], yl, '--', 'color', greenCB,'LineWidth', line_width);
    p7 = plot([180, 180], yl, '--', 'color', redCB,'LineWidth', line_width);
    plotLabelLines = [p1,... % p2,...
        p3, p5, p4, p6, p7];
    if j == 3
        legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
            'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
    end

    ylim(yl);

    figCount = figCount + 1;
    linkaxes(ax,'x');
end
if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end

%%
fig(figCount) = figure(figCount);
% ax = [];
for j = 1:3
    ax1(countAx) = subplot(1, 3, j);
    ax1(countAx).FontSize = fontSize;
    ax1(countAx).GridLineWidth = gridLineWidth;
    ax1(countAx).LineWidth = axesLineWidth;
    set(gcf, 'Position', figLoc)
    grid minor
    hold on
    box on
    countAx = countAx + 1;
    xErrAll = xErrAll_new{j};

    for iiMC = plotList
        xErr = xErrAll{iiMC};
        p1 = plot(time, xErr(iiStatePlot*3-(3-iiState), :), 'color', grayCB, 'LineWidth',1);
    end
    plotLabelLines = p1;
    if flagPlotTitle
        title(titles_states{j}, 'FontSize', fontSize)
    end
    xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)

    % grid minor

    ylabel(strcat('$\mathbf{r}_{t/c, y}$ Error', ' ', axisUnits_states{iiStatePlot}),  'FontSize',fontSize)
    % if ~flagCovPlotOff
        p4 = plot(time, errAll_new(iiStatePlot*3-(3-iiState), :, j), 'k', 'LineWidth', line_width + 0.5);
        p5 = plot(time, 3*sigmaAll_new(iiStatePlot*3-(3-iiState), :, j), ':', 'Color', darkGreenCB, 'LineWidth', line_width);
        plot(time, -3*sigmaAll_new(iiStatePlot*3-(3-iiState), :, j), ':', 'Color', darkGreenCB, 'LineWidth', line_width)
        p3 = plot(time, 3*covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j), '--', 'Color', magentaCB, 'LineWidth', line_width);
        plot(time, -3*covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j), '--', 'Color', magentaCB, 'LineWidth', line_width)

        % yl = ylim;
        % p6 = plot([121, 121], yl, '--', 'color', greenCB,'LineWidth', line_width);
        % p7 = plot([180, 180], yl, '--', 'color', redCB,'LineWidth', line_width);
        % plotLabelLines = [p1,... % p2,...
        %     p3, p5, p4, p6, p7];
        % if j == 1
        %     legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
        %         'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
        % end
    % end
    if strcmp(flagZoomSave, 'zoomEnd')
        xlim([xlimSet(1), xlimSet(2)])
    elseif strcmp(flagZoomSave, 'zoomAll')
        xlim(xlimSet)
    end

    yl = ylim;
    p6 = plot([121, 121], [-0.3, 0.3], '--', 'color', greenCB,'LineWidth', line_width);
    p7 = plot([180, 180], [-0.3, 0.3], '--', 'color', redCB,'LineWidth', line_width);
    plotLabelLines = [p1,... % p2,...
        p3, p5, p4, p6, p7];
    if j == 3
        legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
            'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
    end

    ylim(yl);

    figCount = figCount + 1;
    linkaxes(ax1,'x');
end
linkaxes(ax1, 'y');
if j == 3
    legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
        'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
end

if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat('yLinked_',figNames_states{iiStatePlot}, figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end



%% target acceleration 
iiStatePlot = 3; % unmodeled target acceleration
fig(figCount) = figure(figCount);
iiState = 2; % Y axis
% ax = [];

for j = 1:3
    ax2(countAx) = subplot(1, 3, j);
    ax2(countAx).FontSize = fontSize;
    ax2(countAx).GridLineWidth = gridLineWidth;
    ax2(countAx).LineWidth = axesLineWidth;
    set(gcf, 'Position', figLoc)
    grid minor
    hold on
    box on
    countAx = countAx + 1;
    xErrAll = xErrAll_new{j};

    for iiMC = plotList
        xErr = xErrAll{iiMC};
        p1 = plot(time, xErr(iiStatePlot*3-(3-iiState), :), 'color', grayCB, 'LineWidth',1);
    end
    plotLabelLines = p1;
    if flagPlotTitle
        title(titles_states{j}, 'FontSize', fontSize)
    end
    xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)

    % grid minor

    ylabel(strcat('$\mathbf{a}_{t, y}$ Error', ' ', axisUnits_states{iiStatePlot}),  'FontSize',fontSize)
    if ~flagCovPlotOff
        p4 = plot(time, errAll_new(iiStatePlot*3-(3-iiState), :, j), 'k', 'LineWidth', line_width + 0.5);
        p5 = plot(time, 3*sigmaAll_new(iiStatePlot*3-(3-iiState), :, j), ':', 'Color', darkGreenCB, 'LineWidth', line_width);
        plot(time, -3*sigmaAll_new(iiStatePlot*3-(3-iiState), :, j), ':', 'Color', darkGreenCB, 'LineWidth', line_width)
        p3 = plot(time, 3*covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j), '--', 'Color', magentaCB, 'LineWidth', line_width);
        plot(time, -3*covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j), '--', 'Color', magentaCB, 'LineWidth', line_width)


    end
    if strcmp(flagZoomSave, 'zoomEnd')
        xlim([xlimSet(1), xlimSet(2)])
    elseif strcmp(flagZoomSave, 'zoomAll')
        xlim(xlimSet)
    end



    figCount = figCount + 1;
    linkaxes(ax2,'x');

    yl = ylim;
    p6 = plot([121, 121], yl, '--', 'color', greenCB,'LineWidth', line_width);
    p7 = plot([180, 180], yl, '--', 'color', redCB,'LineWidth', line_width);
    plotLabelLines = [p1,... % p2,...
        p3, p5, p4, p6, p7];
    ylim(yl)
    if j == 3
        legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
            'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
    end
end
if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end

%% target acceleration equal zoom ylim 
iiStatePlot = 3; % unmodeled target acceleration
fig(figCount) = figure(figCount);
iiState = 2; % Y axis
% ax = [];

for j = 1:3
    ax3(countAx) = subplot(1, 3, j);
    ax3(countAx).FontSize = fontSize;
    ax3(countAx).GridLineWidth = gridLineWidth;
    ax3(countAx).LineWidth = axesLineWidth;
    set(gcf, 'Position', figLoc)
    grid minor
    hold on
    box on
    countAx = countAx + 1;
    xErrAll = xErrAll_new{j};

    for iiMC = plotList
        xErr = xErrAll{iiMC};
        p1 = plot(time, xErr(iiStatePlot*3-(3-iiState), :), 'color', grayCB, 'LineWidth',1);
    end
    plotLabelLines = p1;
    if flagPlotTitle
        title(titles_states{j}, 'FontSize', fontSize)
    end
    xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)

    % grid minor

    ylabel(strcat('$\mathbf{a}_{t, y}$ Error', ' ', axisUnits_states{iiStatePlot}),  'FontSize',fontSize)
    if ~flagCovPlotOff
        p4 = plot(time, errAll_new(iiStatePlot*3-(3-iiState), :, j), 'k', 'LineWidth', line_width + 0.5);
        p5 = plot(time, 3*sigmaAll_new(iiStatePlot*3-(3-iiState), :, j), ':', 'Color', darkGreenCB, 'LineWidth', line_width);
        plot(time, -3*sigmaAll_new(iiStatePlot*3-(3-iiState), :, j), ':', 'Color', darkGreenCB, 'LineWidth', line_width)
        p3 = plot(time, 3*covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j), '--', 'Color', magentaCB, 'LineWidth', line_width);
        plot(time, -3*covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j), '--', 'Color', magentaCB, 'LineWidth', line_width)

        % plotLabelLines = [p1,... % p2,...
        %     p3, p5, p4];
        % if j == 1
        %     legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
        %         'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Interpreter','latex', 'location', 'SouthEast')
        % end
    end
    if strcmp(flagZoomSave, 'zoomEnd')
        xlim([xlimSet(1), xlimSet(2)])
    elseif strcmp(flagZoomSave, 'zoomAll')
        xlim(xlimSet)
    end



    figCount = figCount + 1;
    linkaxes(ax3,'x');

    yl = ylim;
    p6 = plot([121, 121], [-0.6, 0.6], '--', 'color', greenCB,'LineWidth', line_width);
    p7 = plot([180, 180], [-0.6, 0.6], '--', 'color', redCB,'LineWidth', line_width);
    plotLabelLines = [p1,... % p2,...
        p3, p5, p4, p6, p7];
    ylim(yl)
    if j == 3
        legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
            'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
    end
end
if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end
linkaxes(ax3, 'y');
if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat('yLinked_',figNames_states{iiStatePlot}, figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end




%% Covariance ratio plot
iiStatePlot = 3;
iiState = 2;

fig(figCount) = figure(figCount);
% figCount = figCount+1;
set(gcf, 'Position', figLoc)


set(gca, 'FontSize', fontSize);
% gcf.FontSize = fontSize;
set(gca, 'GridLineWidth', gridLineWidth);
set(gca, 'LineWidth', axesLineWidth);

colors = [grayCB; darkGreenCB; magentaCB];
lines = {'-', '--', ':'};

hold on
for j = 1:3
sampleSigma = sigmaAll_new(iiStatePlot*3-(3-iiState), :, j);
filterSigma = covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j);

ratioSigma = sampleSigma./filterSigma;

plot(time, ratioSigma, lines{j}, 'Color', colors(j, :), 'LineWidth', line_width)

end


plot([time(1), time(end)], [1, 1], 'k', 'LineWidth', line_width + 0.5)
plot([121, 121], [0, 1], '--', 'color', greenCB,'LineWidth', line_width);
plot([180, 180], [0, 1], '--', 'color', redCB,'LineWidth', line_width);

grid minor 
box on
ylabel('Filter $\sigma$/ Sample $\sigma$ [$m/s^2$]',  'Interpreter','latex', 'fontSize', fontSize)
xlabel('Time [s]',  'Interpreter','latex', 'fontSize', fontSize)
if strcmp(flagZoomSave, 'zoomEnd')
    xlim([xlimSet(1), xlimSet(2)])
elseif strcmp(flagZoomSave, 'zoomAll')
    xlim(xlimSet)
end

ylim([0, 1])

legend('BALOO', '6Model IMM', '2Model IMM', 'Perfect Match', 'Start Burn', 'End Burn', 'Interpreter','latex', 'fontSize', fontSize);





if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat('accelSigmaRatio', figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end



%%
iiStatePlot = 1;
iiState = 2;

figCount = figCount + 1;

fig(figCount) = figure(figCount);
% figCount = figCount+1;
set(gcf, 'Position', figLoc)


set(gca, 'FontSize', fontSize);
% gcf.FontSize = fontSize;
set(gca, 'GridLineWidth', gridLineWidth);
set(gca, 'LineWidth', axesLineWidth);

colors = [grayCB; darkGreenCB; magentaCB];
lines = {'-', '--', ':'};

hold on
for j = 1:3
sampleSigma = sigmaAll_new(iiStatePlot*3-(3-iiState), :, j);
filterSigma = covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j);

errRatioSigma = (sampleSigma-filterSigma)./filterSigma;
ratioSigma = (sampleSigma)./filterSigma;
plot(time, ratioSigma, lines{j}, 'Color', colors(j, :), 'LineWidth', line_width)
% plot(time, errRatioSigma, lines{j}, 'Color', colors(j, :), 'LineWidth', line_width)

end


plot([time(1), time(end)], [1, 1], 'k', 'LineWidth', line_width + 0.5)
plot([121, 121], [0, 2], '--', 'color', greenCB,'LineWidth', line_width);
plot([180, 180], [0, 2], '--', 'color', redCB,'LineWidth', line_width);

grid minor 
box on
ylabel('Filter $\sigma$/ Sample $\sigma$ [$m/s^2$]',  'Interpreter','latex', 'fontSize', fontSize)
xlabel('Time [s]',  'Interpreter','latex', 'fontSize', fontSize)
if strcmp(flagZoomSave, 'zoomEnd')
    xlim([xlimSet(1), xlimSet(2)])
elseif strcmp(flagZoomSave, 'zoomAll')
    xlim(xlimSet)
end

ylim([0, 2])

legend('BALOO', '6Model IMM', '2Model IMM', 'Perfect Match', 'Start Burn', 'End Burn', 'Interpreter','latex', 'fontSize', fontSize);





if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat('positionSigmaRatio', figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end

