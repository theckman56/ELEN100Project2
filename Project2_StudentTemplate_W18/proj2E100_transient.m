function F = proj2E100_transient(t, x)
%
% This function provides the right-hand side of the differential equation
% for the MATLAB solver.

F = [ ? ];
% Vector F is initialized.

R1 = evalin('base','R1_circuit');
R2 = evalin('base','R2_circuit');
R3 = evalin('base','R3_circuit');
R4 = evalin('base','R4_circuit');
R5 = evalin('base','R5_circuit');
% These are the resistor values for the circuit. If you want to perform a
% simulation with different R values, set them in the Rx_circuit variables
% that are defined in the main MATLAB script (or workspace) before calling
% this function.

if (t >= ? ) && (t <= ? )
    h = ? ;
else
    h = ? ;
end;
% h(t) represents an approximation of the step function, with a rise time
% of 1 microsecond.

F(1) = ? ;
F(2) = ? ;
F(3) = ? ;
F(4) = ? ;
F(5) = ? ;
% This is the right-hand side of the DAE. It is written in terms of R, so
% that you don't need to rewrite the code every time your element values
% change.
