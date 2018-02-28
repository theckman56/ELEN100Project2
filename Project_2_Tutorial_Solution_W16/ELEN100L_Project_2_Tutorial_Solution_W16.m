%% ELEN 100L (Electric Circuits II): Project 2 Tutorial
%
% <<ELEN100L_Project_2_Tutorial_figure_01.png>>
%
% *Hard Copy Deliverables:*
%
% # Hard copy for hand calculations.
% # A MATLAB script and publish the solution using MATLAB's *publish*
% feature.
% # Turn in MATLAB scripts and a document of the run-time results.
% # Turn in oscilloscope images of the measured results.
%
% *Soft Copy Deliverables:*
%
% # Turn in all MATLAB files.
% # Turn in all LTSpice files.
% # Turn in all oscilloscope images of the measured results.
%

%% Initialize MATLAB Environment
%

clear; clc; clf; cla; close all;
format long; format compact;

%% Setup global variables
%

% These Ideal Design element values are fixed in the circuit.
VG  = 1;                   % Generator voltage

R1_ideal = 2.0;            % Ohms
R2_ideal = 2.0;            % Ohms
R3_ideal = 1.0;            % Ohms
C1_ideal = 0.5;            % Farads

% Build an array for the R elements.
R_ideal = [R1_ideal, R2_ideal, R3_ideal];

% Build an array for the C elements.
C_ideal = [(0), (0),  (0)        ; ...
           (0), (0),  (0)        ; ...
           (0), (0), -(C1_ideal)];

% Build an array for the source elements.
B = [VG; 0; 0];

% Build an array for the the time vector.
time = [0, 5];

% Build an array for the the initial conditions.
x0 = [0; 0; 0];     % Assume everything is zero to start

% These values are used for plotting purposes.
fignum = 1;
plot_left   = 0;       plot_right = time(2);    % x-axis range (seconds)
plot_bottom = 0;       plot_top   = VG+0.1;     % y-axis range (volts)

%% Problem 2
%
% <<ELEN100L_Project_2_Tutorial_figure_02.png>>
%

fignum = fignum+1;

%%
% Display the component values for the Ideal design.
%

display(' ');
display('The Ideal Design component values are:');
fprintf('    R1 = %+11.4f Ohms.\n',   R1_ideal);
fprintf('    R2 = %+11.4f Ohms.\n',   R2_ideal);
fprintf('    R3 = %+11.4f Ohms.\n',   R3_ideal);
fprintf('    C1 = %+11.4e Farads.\n', C1_ideal);

%%
% Calculate the transient response for the Ideal design.
%

% Update the resistor variables used in the proj2E100_transient function
% before calling the ode23t solver.
R1_circuit = R1_ideal;
R2_circuit = R2_ideal;
R3_circuit = R3_ideal;

options = odeset('mass', C_ideal, 'RelTol', 0.1e-9);
[t, x] = ode23t(@proj2E100_transient, time, x0, options);

% Capture a voltage point with index.
v3_point_of_interest = 400e-3; % This particular value is chosen only for
                               % plot annotation illustration purposes.

v3_shifted_mag = abs( x(:,3) - v3_point_of_interest );
v3_point_of_interest_index = ...
    find( min( v3_shifted_mag ) == v3_shifted_mag );
v3_point_of_interest = x(v3_point_of_interest_index, 3);

% Capture the time stamp at voltage point of interest.
t_point_of_interest  = t(v3_point_of_interest_index);

%%
% Generate the plot for the transient response.

fignum = 5;  % Bump the figure number to match tutorial documentation.

fignum = fignum+1; figObj = figure(fignum);  % Establish a figure number
set(fignum, 'Name', ...
    ['Transient Response Ideal Design']);    % Name the figure

Tr_ideal_Plot = plot(t, x);                  % Generate plot
grid on;                                     % Turn grid on
xlabel('Time (seconds)');                    % Label the x-axis
ylabel('Amplitude (volts)');                 % Label the y-axis
axis([plot_left, plot_right, ...
      plot_bottom, plot_top]);               % Bound plot
title(['Figure ',num2str(fignum,'%-2.u'),...
       ': Transient Response']);
legend('v_1(t)', 'v_2(t)', 'v_3(t)', 'Location', 'NorthEast');

% Add annotation to the plot.
str_curs = ['\leftarrow ', num2str(t_point_of_interest), 's, ', ...
            num2str(v3_point_of_interest),'V'];
text(t_point_of_interest, ...
     v3_point_of_interest, ...
     str_curs, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

%%
% Display the MATLAB point of interest values for the Ideal design.
%

display(' ');
display('The MATLAB point of interest values are:');
fprintf('    V3 p.o.i. = %+11.4f Volts.\n',   v3_point_of_interest);
fprintf('     t p.o.i. = %+11.4f seconds.\n', t_point_of_interest);

%% Problem 3
%
% <<ELEN100L_Project_2_Tutorial_figure_07.png>>
%

fignum = fignum+1;

%%
% The LTSpice model for the circuit is shown below.
%
% <<ELEN100L_Project_2_Tutorial_figure_08.png>>
%

fignum = fignum+1;

%%
% The LTSpice model voltage source setup, to correlate with the MATLAB unit
% step function definition in Appendix #3, is shown below.
%
% <<ELEN100L_Project_2_Tutorial_figure_09.png>>
%
% <<ELEN100L_Project_2_Tutorial_figure_10.png>>
%

fignum = fignum+2;

%%
% The LTSpice model transient analysis simulation setup is shown below.
%
% <<ELEN100L_Project_2_Tutorial_figure_11.png>>
%

fignum = fignum+1;

%%
% The LTSpice model for the simulation result is shown below.
%
% <<ELEN100L_Project_2_Tutorial_figure_12.png>>
%

fignum = fignum+1;

% Capture point of interest voltage from the plot.
ltspice_v3_point_of_interest = 400.00051e-3;

% Capture point of interest time stamp from the plot.
ltspice_t_point_of_interest  = 1.6094435;

%%
% Display the LTSpice voltage point of interest values for the Ideal
% design.
%

display(' ');
display('The LTSpice point of interest values are:');
fprintf('    V3 p.o.i. = %+11.4f Volts.\n', ...
    ltspice_v3_point_of_interest);
fprintf('     t p.o.i. = %+11.4f seconds.\n', ...
    ltspice_t_point_of_interest);

%%
% Calculate the percent difference at the point of interest values between
% MATLAB and LTSpice Ideal Designs.
%

diff_ideal_v3_point_of_interest  = ...
    (ltspice_v3_point_of_interest - v3_point_of_interest) ...
    /abs(v3_point_of_interest)*100;

diff_ideal_t_point_of_interest  = ...
    (ltspice_t_point_of_interest - t_point_of_interest) ...
    /abs(t_point_of_interest)*100;

display(' ');
display('The % difference between MATLAB and LTSpice at');
display('the point of interest:');
fprintf('    MATLAB  V3 p.o.i. = %+11.4f Volts.\n', ...
         v3_point_of_interest);
fprintf('    LTSpice V3 p.o.i. = %+11.4f Volts.\n', ...
         ltspice_v3_point_of_interest);
fprintf('        %% diff = %+8.4f (%%).\n', ...
         diff_ideal_v3_point_of_interest);

fprintf('    MATLAB   t p.o.i. = %+11.4f seconds.\n', ...
         t_point_of_interest);
fprintf('    LTSpice  t p.o.i. = %+11.4f seconds.\n', ...
         ltspice_t_point_of_interest);
fprintf('        %% diff = %+8.4f (%%).\n', ...
         diff_ideal_t_point_of_interest);

%% Program execution complete
%

display(' ');
disp('Program execution complete....');

%% MATLAB code listing
%