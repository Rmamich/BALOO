%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Truth Gen, Filter, and Plotting
%
%
% Created: 11/07/22
% Last Modified: 12/05/22
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%house keeping
clear
close all


%% Main Simulation Flags
% flagScenario ='EarthRel_gravDiff_burnRendezvous_meanAdjust_smallerBuff';%_meanAdjust_sub';%'EarthRel_procNoise_burnRendezvous';% 'EarthRel_gravDiff';%  _fixedQ_biggerQ % 'EarthRel_Venting'; % 'EarthRel_procNoise_burnRendezvous'; %  %'EarthRel_procNoise';%x_burnRendezvous_fixedQ_biggerQ   _fixedQ_biggerQ'; % '; % 'EarthRel_procNoise'; % 'EarthRel_Venting'; % ; % 'EarthRel_procNoise_burnRendezvous'; % 'EarthRel_Venting_lowerRW'; % 'EarthRel_procNoise'; %; %_noVWR 'EarthRel_Venting'; % Determiines which IC file to run 'EarthRel_procNoise_burnRendezvous' % _burnRendezvous_noProcNoise
% 'EarthRel_gravDiff_burnRendezvous_meanAdjust_sub_tryingSomething'
% flagScenario ='EarthRel_gravDiff_burnRendezvous_vBarAdjust';
flagScenario = 'NRHO_Venting_RelNav_nonDim'; % _direct
% flagScenario = 'EarthRel_Venting';
% flagScenario = 'EarthRel_gravDiff'; %
% flagScenario ='EarthRel_gravDiff_burnRendezvous_fixedQ_biggerQ';
% flagScenario = 'EarthRel_gravDiff_fixed';
% flagScenario = 'NRHO_Venting_RelNav_nonDim_fixed';
flagSetFlags = false;                                       % Flag to clear all or not for subscripts -- it's existence means all subsequent scripts DO NOT clear before running
runTruthGen = true;                                         % true runf truth generation scripts, false supresses this call and uses .mat files if they exist for scripts downstream 
runMC = true;                                               % true runs the filter, false supresses the call and uses .mat files if they exist for plotting
plotMC = true;                                              % true calls plotting script, if false, no plots from the filter will be generated
getStats = false;
saveFig = false;                                            % Saves plots when true
flagPlotFrame = 'rtn';                                      % plotting frame 'rtn' for radial, transverse, normal, or 'N' for inertial


flagGenCellsForPlotting = true;




%% Plotting Flags
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
flagZoomSave = 'zoomEnd';                                   % Sets x limits for plotting script plots -- 'none' is default, 'zoomAll' enforces both limits, 'zoomEnd' enforces the upper limit
% xlimSet = [250, 8528];                                      % Sets x limits for plotting script plots
xlimSet = [300, 5500];

% TODO: Add this flag to all scenario IC files
% New Flag for vBar compensation w/o accel state
flagCompensateMean = false;
%% Call Truth

if contains(flagScenario, 'Earth')
    if contains(flagScenario, 'Rendezvous')
        if contains(flagScenario, 'lagerBuff')
            IC_dualBurnRendezvous_largeBuffer;
        elseif contains(flagScenario, 'gravDiff')
            if contains(flagScenario, 'vBarAdjust')
                IC_dualBurnRendezvous_gravDiff_noAccelState_vBar;
            elseif contains(flagScenario, 'fixed')
                IC_dualBurnRendezvous_gravDiff_plainEKF;
            else
                IC_dualBurnRendezvous_gravDiff;
            end
        else
            if contains(flagScenario, 'fixed')
                IC_dualBurnRendezvous_plainEKF;
            else
                IC_dualBurnRendezvous;
            end
            
        end
    elseif contains(flagScenario, 'Venting')
        IC_venting;
    elseif contains(flagScenario, 'PeriodSet')
        IC_procNoiseOnly;
    elseif contains(flagScenario, 'longerBuffer')
        IC_procNoiseOnly_longerBuffer;
    elseif contains(flagScenario, 'gravDiff')
        if contains(flagScenario, 'fixed')
            IC_gravDiff_plainEKF_BALOO
        else
            IC_gravDiff;
        end
    else
        IC_procNoiseOnly;
    end
else
    if contains(flagScenario, 'nonDim')
        if contains(flagScenario, 'direct')
            IC_NRHOVenting_nonDim_directMeas;
        elseif contains(flagScenario, 'fixed')
            IC_NRHOVenting_nonDim_plainEKF;
        else
            IC_NRHOVenting_nonDim;
        end
    else
        IC_NRHOVenting;
    end
end

if nMCTruth > 1
    nMC = 1;
end


for iiTruth = firstPlot:nMCTruth
    warning('off','all');
    if runTruthGen
        if contains(flagScenario, 'Kayhan')
            % Should run for either Kayhan Rel Nav Scenario or Kayhan Ineertial
            KayhanSpaceTruthTraj;
        elseif contains(flagScenario, 'Earth')
            if contains(flagScenario, 'procNoise')
                if contains(flagScenario, 'burn')
                    if contains(flagScenario, 'Rendezvous')
                        Earth_ConstQ_DualBurnRendezvous
                    else
                        Earth_ConstQ_burns
                    end
                elseif contains(flagScenario, 'PeriodSet')
                    Earth_ConstQ_PeriodSet;
                else
                    Earth_ConstQ;
                end
            elseif contains(flagScenario, 'Venting')
                EarthVenting_scratch;
                EarthOrbitVenting;
            elseif contains(flagScenario, 'gravDiff')
                if contains(flagScenario, 'Rendezvous')
                    Earth_NoQ_DualBurnRendezvous;
                else
                    Earth_gravDiff;
                end
            else
                simpleCM_Earth;
            end
        elseif contains(flagScenario, 'Erros')
            simpleCMSim;
        else
            if contains(flagScenario, 'nonDim')
                if contains(flagScenario, 'fixed')
                    NRHO_Scratch_keepitNonDim_zeroVOutVents;
                    NRHOVenting_nonDim;
                else
                    NRHOVenting_Scratch_keepitNonDim; % randomize vents
                    NRHOVenting_nonDim;
                end
            else
                NRHOVenting;
            end
        end
    end


    if runMC
        runFilterMC;
    end

end

if nMCTruth > 1
    nMC = nMCTruth;
end

if plotMC
    fprintf('Running Plotting Tool\n')
    MCplotNew;
end

if getStats
    fprintf('Running Gen Stats\n')
    genMaxAvgStats;
end


clear flagSetFlags

% load gong.mat
% sound(y)

