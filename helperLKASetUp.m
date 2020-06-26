% Set up Script for the Lane Keeping Assist (LKA) Example
%
% This script initializes the LKA example model. It loads necessary control
% constants and sets up the buses required for the referenced model
%
%   This is a helper script for example purposes and may be removed or
%   modified in the future.

%   Copyright 2017-2018 The MathWorks, Inc.

%% General Model Parameters
Ts = 0.1;               % Simulation sample time                (s)

%% Ego Car Parameters
% Dynamics modeling parameters
m       = 1575;     % Total mass of vehicle                          (kg)
Iz      = 2875;     % Yaw moment of inertia of vehicle               (m*N*s^2)
lf      = 1.2;      % Longitudinal distance from c.g. to front tires (m)
lr      = 1.6;      % Longitudinal distance from c.g. to rear tires  (m)
Cf      = 19000;    % Cornering stiffness of front tires             (N/rad)
Cr      = 33000;    % Cornering stiffness of rear tires              (N/rad)

%% Controller parameter
PredictionHorizon = 30; % Number of steps for preview    (N/A)

%% Bus Creation
% Create buses for lane sensor and lane sensor boundaries
createLaneSensorBuses
% Create the bus of actors from the scenario reader
modelName = 'LKATestBenchExample';
wasModelLoaded = bdIsLoaded(modelName);
if ~wasModelLoaded
    load_system(modelName)
end
% load the bus for scenario reader
blk=find_system(modelName,'System','driving.scenario.internal.ScenarioReader');
s = get_param(blk{1},'PortHandles');
get(s.Outport(1),'SignalHierarchy');

%% Create scenario and road specifications
[scenario,roadCenters,laneSpecification] = createDoubleCurveScenario;

% You can use Driving Scenario Designer to explore the scenario
% drivingScenarioDesigner(scenario)
% drivingScenarioDesigner('LKATestBenchScenario')

%% Generate data for Simulink simulation  
[driverPath,x0_ego,y0_ego,v0_ego,yaw0_ego,simStopTime] = ...
    createDriverPath(scenario,6);