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
figLoc = [10, 10, 1850, 1400];
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
% xlimSet = [100, 300];                                      % Sets x limits for plotting script plots
% xlimSet = [780, 1000]; % NRHO bigger burn window

xlimSet = [330, 540]; % NRHO burn window

% TODO: Add this flag to all scenario IC files
% New Flag for vBar compensation w/o accel state
flagCompensateMean = false;

flagPlotFrame = 'N';
mu = 0.012150584270; % mass ratio from paper
L = 384400; % km -- characteristic length 
L_m = L*1000;
tCharacteristic = 3.751902619517e5; % s -- characteristic time from paper  

%% Load in Data

load(fullfile('./MAT/NRHO_RelNav_nonDim_burnRendezvous_venting_BALOO_more/InitCovFactor/none/EKF', 'timeMat.mat'))
fileNameSaveCells_BALOO = './MAT/NRHO_RelNav_nonDim_burnRendezvous_venting_BALOO_more/InitCovFactor/none/EKF/accelOnBiasOnLidarBiasOffCRTBPCRTBP_noisyAccelerometer/fullCells';
% fileNameSaveCells_BALOO = './MAT/NRHO_RelNav_nonDim_burnRendezvous_venting_bigManeuver_BALOO/InitCovFactor/none/EKF/accelOnBiasOnLidarBiasOffCRTBPCRTBP_noisyAccelerometer/fullCells';
% load(fullfile('./MAT/NRHO_RelNav_nonDim_burnRendezvous_venting_bigManeuver_BALOO/InitCovFactor/none/EKF', 'timeMat.mat'))

% 'NRHO_RelNav_nonDim_burnRendezvous_venting_bigManeuver_BALOO';
fileNameSaveCells_6model = './MAT/NRHO_RelNav_nonDim_burnRendezvous_venting_fixed_interacting/InitCovFactor/none/EKF/accelOnBiasOnLidarBiasOffCRTBPCRTBP_noisyAccelerometer/fullCells';
fileNameSaveCells_2model = './MAT/NRHO_RelNav_nonDim_burnRendezvous_venting_fixed_interacting_2modelsNoiseDifference/InitCovFactor/none/EKF/accelOnBiasOnLidarBiasOffCRTBPCRTBP_noisyAccelerometer/fullCells';

fileNameSaveCells = {fileNameSaveCells_BALOO, fileNameSaveCells_6model, fileNameSaveCells_2model};

xErrAll_new = cell(3, 1);
covAllMCPlus_new = NaN(18, 3001, 3);
errAll_new = NaN(18, 3001, 3);
sigmaAll_new = NaN(18, 3001, 3);
avg_nodes = NaN(3, 3001);

for i = 1:3
    load(fullfile(fileNameSaveCells{i}, 'xErr.mat')); % , "xErrAll")
    load(fullfile(fileNameSaveCells{i}, 'covAllMCPlus.mat')); % , "covAllMCPlus")
    load(fullfile(fileNameSaveCells{i}, 'errAll.mat')); %, "errAll")
    load(fullfile(fileNameSaveCells{i}, 'sigmaAll.mat')); % , "sigmaAll")
    % load(fullfile(fileNameSaveCells{i}, 'numNodes.mat'));
    load(fullfile(fileNameSaveCells{i}, 'num_nodesPlot_all.mat'))


    xErrAll_new{i} = xErrAll; % All err lines
    covAllMCPlus_new(:, :, i) = covAllMCPlus; % avg est 3 sigma
    errAll_new(:, :, i) = errAll; % avg err
    sigmaAll_new(:, :, i) = sigmaAll;
    avg_nodes(i, :) = round(mean(num_nodesPlot_all));


end


titles_states = {'BALOO', '7 Model IMM', '2 Model IMM'};
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

figNameSave = './Figures/NRHORel_venting_burnRendezvous/';


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
time = time*tCharacteristic;

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

ylabel_position = {'$\mathbf{r}_{t/c, x}$', '$\mathbf{r}_{t/c, y}$', '$\mathbf{r}_{t/c, z}$'};
%% Relativ Position -- not linked y 
iiStatePlot = 1; % unmodeled target acceleration
fig(figCount) = figure(figCount);
iiState = 2; % Y axis
% ax = [];



for j = 1:3
    xErrAll = xErrAll_new{j};
    for iiState = 1:3
        ax(countAx) = subplot(4, 3, ((iiState-1)*3+j));
        ax(countAx).FontSize = fontSize;
        ax(countAx).GridLineWidth = gridLineWidth;
        ax(countAx).LineWidth = axesLineWidth;
        set(gcf, 'Position', figLoc)
        grid minor
        hold on
        box on
        countAx = countAx + 1;


        for iiMC = plotList
            xErr = xErrAll{iiMC};
            p1 = plot(time, xErr(iiStatePlot*3-(3-iiState), :), 'color', grayCB, 'LineWidth',1);
        end
        plotLabelLines = p1;
        if flagPlotTitle
            if iiState == 1
                title(titles_states{j}, 'FontSize', fontSize)
            end
        end
        % if iiState == 3
        %     xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)
        % end
        % grid minor
        if j == 1
            ylabel(strcat(ylabel_position{iiState}, ' ', axisUnits_states{iiStatePlot}),  'FontSize',fontSize)
        end
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
        p6 = plot([351, 351], yl, '--', 'color', greenCB,'LineWidth', line_width);
        p7 = plot([380, 380], yl, '--', 'color', redCB,'LineWidth', line_width);
        % plotLabelLines = [p1,... % p2,...
        %     p3, p5, p4, p6, p7];
        % if j == 3
        %     legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
        %         'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
        % end

        % ylim(yl);


        linkaxes(ax,'x');
    end
    ax(countAx) = subplot(4, 3, ((3)*3+j));
    ax(countAx).FontSize = fontSize;
    ax(countAx).GridLineWidth = gridLineWidth;
    ax(countAx).LineWidth = axesLineWidth;
    set(gcf, 'Position', figLoc)
    grid minor
    hold on
    box on
    countAx = countAx + 1;
    xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)
    ylabel('Nodes [-]',  'FontSize',fontSize)
    plot(time, avg_nodes(j, :), 'ok')
    if strcmp(flagZoomSave, 'zoomEnd')
        xlim([xlimSet(1), xlimSet(2)])
    elseif strcmp(flagZoomSave, 'zoomAll')
        xlim(xlimSet)
    end
    yl_node = ylim;
    plot([351, 351], yl_node, '--', 'color', greenCB,'LineWidth', line_width);
    plot([380, 380], yl_node, '--', 'color', redCB,'LineWidth', line_width);

end

figCount = figCount + 1;
if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, 'allAxes',figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end

%% position y linked axes 

fig(figCount) = figure(figCount);


for j = 1:3
    xErrAll = xErrAll_new{j};
    for iiState = 1:3
        ax1(j, iiState) = subplot(4, 3, ((iiState-1)*3+j));
        ax1(j, iiState).FontSize = fontSize;
        ax1(j, iiState).GridLineWidth = gridLineWidth;
        ax1(j, iiState).LineWidth = axesLineWidth;
        set(gcf, 'Position', figLoc)
        grid minor
        hold on
        box on


        for iiMC = plotList
            xErr = xErrAll{iiMC};
            p1 = plot(time, xErr(iiStatePlot*3-(3-iiState), :), 'color', grayCB, 'LineWidth',1);
        end
        plotLabelLines = p1;
        if flagPlotTitle
            if iiState == 1
                title(titles_states{j}, 'FontSize', fontSize)
            end
        end
        % if iiState == 3
        %     xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)
        % end
        % grid minor
        if j == 1
            ylabel(strcat(ylabel_position{iiState}, ' ', axisUnits_states{iiStatePlot}),  'FontSize',fontSize)
        end
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

        % yl = ylim;
        p6 = plot([351, 351], [-0.3, 0.3], '--', 'color', greenCB,'LineWidth', line_width);
        p7 = plot([380, 380], [-0.3, 0.3], '--', 'color', redCB,'LineWidth', line_width);
        plotLabelLines = [p1,... % p2,...
            p3, p5, p4, p6, p7];
        % if j == 3
        %     legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
        %         'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
        % end

        % ylim(yl);


        linkaxes(ax1,'x');
    end
    ax(countAx) = subplot(4, 3, ((3)*3+j));
    ax(countAx).FontSize = fontSize;
    ax(countAx).GridLineWidth = gridLineWidth;
    ax(countAx).LineWidth = axesLineWidth;
    set(gcf, 'Position', figLoc)
    grid minor
    hold on
    box on
    countAx = countAx + 1;
    xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)
    ylabel('Nodes [-]',  'FontSize',fontSize)
    plot(time, avg_nodes(j, :), 'ok')
    if strcmp(flagZoomSave, 'zoomEnd')
        xlim([xlimSet(1), xlimSet(2)])
    elseif strcmp(flagZoomSave, 'zoomAll')
        xlim(xlimSet)
    end
    yl_node = ylim;
    plot([351, 351], yl_node, '--', 'color', greenCB,'LineWidth', line_width);
    plot([380, 380], yl_node, '--', 'color', redCB,'LineWidth', line_width);
end
linkaxes(ax1(:, 1), 'y')
linkaxes(ax1(:, 2), 'y')
linkaxes(ax1(:, 3), 'y')

if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat('yLinked_',figNames_states{iiStatePlot}, 'allAxes',figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end



%% target acceleration  -- not linked
ylabel_accel = {'$\mathbf{a}_{t, x}$', '$\mathbf{a}_{t, y}$', '$\mathbf{a}_{t, z}$'};
figCount = figCount + 1;

iiStatePlot = 3; % unmodeled target acceleration
fig(figCount) = figure(figCount);
iiState = 2; % Y axis
% ax = [];

for j = 1:3
    xErrAll = xErrAll_new{j};
    for iiState = 1:3
        ax2(j, iiState) = subplot(4, 3, ((iiState-1)*3+j));
        ax2(j, iiState).FontSize = fontSize;
        ax2(j, iiState).GridLineWidth = gridLineWidth;
        ax2(j, iiState).LineWidth = axesLineWidth;
        set(gcf, 'Position', figLoc)
        grid minor
        hold on
        box on
        countAx = countAx + 1;


        for iiMC = plotList
            xErr = xErrAll{iiMC};
            p1 = plot(time, xErr(iiStatePlot*3-(3-iiState), :), 'color', grayCB, 'LineWidth',1);
        end
        plotLabelLines = p1;
        if flagPlotTitle
            if iiState == 1
                title(titles_states{j}, 'FontSize', fontSize)
            end
        end
        % if iiState == 3
        %     xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)
        % end
        % grid minor
        if j == 1
            ylabel(strcat(ylabel_accel{iiState}, ' ', axisUnits_states{iiStatePlot}),  'FontSize',fontSize)
        end
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
        p6 = plot([351, 351], yl, '--', 'color', greenCB,'LineWidth', line_width);
        p7 = plot([380, 380], yl, '--', 'color', redCB,'LineWidth', line_width);
        plotLabelLines = [p1,... % p2,...
            p3, p5, p4, p6, p7];
        % if j == 3
        %     legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
        %         'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
        % end

        % ylim(yl);


        linkaxes(ax2,'x');
    end
    ax(countAx) = subplot(4, 3, ((3)*3+j));
    ax(countAx).FontSize = fontSize;
    ax(countAx).GridLineWidth = gridLineWidth;
    ax(countAx).LineWidth = axesLineWidth;
    set(gcf, 'Position', figLoc)
    grid minor
    hold on
    box on
    countAx = countAx + 1;
    xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)
    ylabel('Nodes [-]',  'FontSize',fontSize)
    plot(time, avg_nodes(j, :), 'ok')
    if strcmp(flagZoomSave, 'zoomEnd')
        xlim([xlimSet(1), xlimSet(2)])
    elseif strcmp(flagZoomSave, 'zoomAll')
        xlim(xlimSet)
    end
    yl_node = ylim;
    plot([351, 351], yl_node, '--', 'color', greenCB,'LineWidth', line_width);
    plot([380, 380], yl_node, '--', 'color', redCB,'LineWidth', line_width);
end


if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, 'allAxes',figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end

%% target acceleration equal zoom ylim 
figCount = figCount + 1;

iiStatePlot = 3; % unmodeled target acceleration
fig(figCount) = figure(figCount);
iiState = 2; % Y axis
% ax = [];

for j = 1:3
    xErrAll = xErrAll_new{j};
    for iiState = 1:3
        ax3(j, iiState) = subplot(4, 3, ((iiState-1)*3+j));
        ax3(j, iiState).FontSize = fontSize;
        ax3(j, iiState).GridLineWidth = gridLineWidth;
        ax3(j, iiState).LineWidth = axesLineWidth;
        set(gcf, 'Position', figLoc)
        grid minor
        hold on
        box on


        for iiMC = plotList
            xErr = xErrAll{iiMC};
            p1 = plot(time, xErr(iiStatePlot*3-(3-iiState), :), 'color', grayCB, 'LineWidth',1);
        end
        plotLabelLines = p1;
        if flagPlotTitle
            if iiState == 1
                title(titles_states{j}, 'FontSize', fontSize)
            end
        end
        % if iiState == 3
        %     xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)
        % end
        % grid minor
        if j == 1
            ylabel(strcat(ylabel_accel{iiState}, ' ', axisUnits_states{iiStatePlot}),  'FontSize',fontSize)
        end
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
        p6 = plot([351, 351], [-0.26, 0.26], '--', 'color', greenCB,'LineWidth', line_width);
        p7 = plot([380, 380], [-0.26, 0.26], '--', 'color', redCB,'LineWidth', line_width);
        plotLabelLines = [p1,... % p2,...
            p3, p5, p4, p6, p7];
        % if j == 3
        %     legend(plotLabelLines, 'ith Run Error$^+$',... 'Avg $3\sigma^-$',
        %         'Avg Filter $3\sigma^+$', 'Sample 3$\sigma^+$', 'Sample Avg Err$^+$', 'Start Burn', 'End Burn', 'Interpreter','latex', 'location', 'SouthEast')
        % end

        % ylim(yl);


        linkaxes(ax3,'x');
    end
    ax(countAx) = subplot(4, 3, ((3)*3+j));
    ax(countAx).FontSize = fontSize;
    ax(countAx).GridLineWidth = gridLineWidth;
    ax(countAx).LineWidth = axesLineWidth;
    set(gcf, 'Position', figLoc)
    grid minor
    hold on
    box on
    countAx = countAx + 1;
    xlabel(strcat('Time ', timeUnit), 'FontSize',fontSize)
    ylabel('Nodes [-]',  'FontSize',fontSize)
    plot(time, avg_nodes(j, :), 'ok')
    if strcmp(flagZoomSave, 'zoomEnd')
        xlim([xlimSet(1), xlimSet(2)])
    elseif strcmp(flagZoomSave, 'zoomAll')
        xlim(xlimSet)
    end
    yl_node = ylim;
    plot([351, 351], yl_node, '--', 'color', greenCB,'LineWidth', line_width);
    plot([380, 380], yl_node, '--', 'color', redCB,'LineWidth', line_width);
end

linkaxes(ax3(:, 1), 'y')
linkaxes(ax3(:, 2), 'y')
linkaxes(ax3(:, 3), 'y')


if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat('yLinked_',figNames_states{iiStatePlot}, 'allAxes', figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end




%% Covariance ratio plot
iiStatePlot = 3;
iiState = 2;
figCount = figCount + 1;

ylabel_accel_ratio = {'Avg $\sigma_{a_{t_x}}$/ Sample $\sigma_{a_{t_x}}$ [$m/s^2$]', 'Avg $\sigma_{a_{t_y}}$/ Sample $\sigma_{a_{t_y}}$ [$m/s^2$]', 'Avg $\sigma_{a_{t_z}}$/ Sample $\sigma_{a_{t_z}}$ [$m/s^2$]'};

fig(figCount) = figure(figCount);
% figCount = figCount+1;
set(gcf, 'Position', figLoc)


set(gca, 'FontSize', fontSize);
% gcf.FontSize = fontSize;
set(gca, 'GridLineWidth', gridLineWidth);
set(gca, 'LineWidth', axesLineWidth);

colors = [grayCB; darkGreenCB; magentaCB];
lines = {'-', '--', ':'};

for j = 1:3
    subplot(1, 3, j)
    title(titles_states{j}, 'FontSize', fontSize)
    hold on
    set(gca, 'FontSize', fontSize);
    for iiState = 1:3   
        sampleSigma = sigmaAll_new(iiStatePlot*3-(3-iiState), :, j);
        filterSigma = covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j);

        ratioSigma = sampleSigma./filterSigma;

        plot(time, ratioSigma, lines{iiState}, 'Color', colors(iiState, :), 'LineWidth', line_width)

    end



plot([time(1), time(end)], [1, 1], 'k', 'LineWidth', line_width + 0.5)
plot([351, 351], [0, 1], '--', 'color', greenCB,'LineWidth', line_width);
plot([380, 380], [0, 1], '--', 'color', redCB,'LineWidth', line_width);



grid minor 
box on

    ylabel(ylabel_accel_ratio{iiState},  'Interpreter','latex', 'fontSize', fontSize)

    xlabel('Time [s]',  'Interpreter','latex', 'fontSize', fontSize)
    

if strcmp(flagZoomSave, 'zoomEnd')
    xlim([xlimSet(1), xlimSet(2)])
elseif strcmp(flagZoomSave, 'zoomAll')
    xlim(xlimSet)
end

ylim([0, 1])

if j == 3
    legend('X', 'Y', 'Z', 'Perfect Match', 'Start Burn', 'End Burn', 'Interpreter','latex', 'fontSize', fontSize, 'location', 'best');
end

end


if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat('accelSigmaRatio', 'allAxes', figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end



%%
iiStatePlot = 1;
figCount = figCount + 1;

ylabel_pos_ratio = {'Avg $\sigma_{r_{t/c_x}}$/ Sample $\sigma_{r_{t/c_x}}$ [$m/s^2$]', 'Avg $\sigma_{r_{t/c_y}}$/ Sample $\sigma_{r_{t/c_y}}$ [$m/s^2$]', 'Avg $\sigma_{r_{t/c_z}}$/ Sample $\sigma_{r_{t/c_z}}$ [$m/s^2$]'};

fig(figCount) = figure(figCount);
% figCount = figCount+1;
set(gcf, 'Position', figLoc)


set(gca, 'FontSize', fontSize);
% gcf.FontSize = fontSize;
set(gca, 'GridLineWidth', gridLineWidth);
set(gca, 'LineWidth', axesLineWidth);

colors = [grayCB; darkGreenCB; magentaCB];
lines = {'-', '--', ':'};

for j = 1:3
    subplot(1, 3, j)
    title(titles_states{j}, 'FontSize', fontSize)
    hold on
    set(gca, 'FontSize', fontSize);
    for iiState = 1:3   
        sampleSigma = sigmaAll_new(iiStatePlot*3-(3-iiState), :, j);
        filterSigma = covAllMCPlus_new(iiStatePlot*3-(3-iiState), :, j);

        ratioSigma = sampleSigma./filterSigma;

        plot(time, ratioSigma, lines{iiState}, 'Color', colors(iiState, :), 'LineWidth', line_width)

    end



plot([time(1), time(end)], [1, 1], 'k', 'LineWidth', line_width + 0.5)
plot([351, 351], [0, 1], '--', 'color', greenCB,'LineWidth', line_width);
plot([380, 380], [0, 1], '--', 'color', redCB,'LineWidth', line_width);



grid minor 
box on

    ylabel(ylabel_pos_ratio{iiState},  'Interpreter','latex', 'fontSize', fontSize)

    xlabel('Time [s]',  'Interpreter','latex', 'fontSize', fontSize)
    

if strcmp(flagZoomSave, 'zoomEnd')
    xlim([xlimSet(1), xlimSet(2)])
elseif strcmp(flagZoomSave, 'zoomAll')
    xlim(xlimSet)
end

ylim([0, 1])

if j == 3
    legend('X', 'Y', 'Z', 'Perfect Match', 'Start Burn', 'End Burn', 'Interpreter','latex', 'fontSize', fontSize, 'location', 'best');
end

end



if saveFig
    figNameCells{figCount} = fullfile(figNameSave, strcat('positionSigmaRatio', 'allAxes',figNameEnd));
    exportgraphics(gcf, figNameCells{figCount}, 'resolution',300)
    % saveas(fig(figCount), fullfile(figNameSave, strcat(figNames_states{iiStatePlot}, figNameEnd)))
end

