%{
Case Study #2 EMTH171  2020
E-Bike

Task 1 Description:   
For this task, assume that there is wind resistance. Assuming 
it starts from standstill, what final velocity and distance 
travelled will the e-bike achieve on a level road with no pedal 
assistance and a constant 4A motor current over a total time 
of 30 secs?  

Orginal code by: P.J. Bones   UCECE
Code edited by: John Elliott and Samuel Vallance
Last modified:	08/10/2020
%}


clear
clc
close all


%{
===========================================================================
                            | Variable Setup | 
===========================================================================
%}


% Constants
g = 9.81;		% Gravitational acceleration in m.s^-2

% E-bike constants
M = 95;		    % Mass of e-bike + rider in kg
Jw = 0.5;		% Motor/wheels moment of intertia in kg.m^2
fd0 = 2;		% Static (& low speed) drag force in N
beta = 0.25;	% Wind resistance factor in N.m^-2.s^2
km = 3.2;		% Motor factor in N.m per amp
rw = 0.35;		% Wheel radius in m

% Simulation parameters
Im = 4;		    % Fixed motor current in A for Task 1
tf = 30;		% Final time in seconds for Task 1
tdelta = 5;	% Time step in seconds

% Calculate derived quantities
Tm = Im * km;   % Fixed motor torque in N.m
Jef = Jw + rw^2 * M;
                % Effective moment of intertia in kg.m^2

% Arrays to store variables
N = tf / tdelta + 1;
uarray = zeros (1, N);
sarray = zeros (1, N);
domdtarray = zeros (1, N);

% Initial conditions (standing start)
om = 0;         % Initial omega_m ('om')
u = 0;          % Speed along road in m/sec
s = 0;          % Distance along road from start in m


%{
===========================================================================
                 | Euler's Method and Speed Calculations | 
===========================================================================
%}


% Iteratively compute u and s as the e-bike proceeds, advancing 
% time by tdelta sec each step
index = 1;
for t = 0:tdelta:(tf-tdelta)
   % Compute the motor acceleration at the start of the step
   domdt = (Tm - rw * (fd0 + beta*u^2)) / Jef;

   % Estimate the state at the end of the time step (n = 'next')
   omn = om + tdelta * domdt;   	% Estimate om by Euler's method
   un = omn * rw;     			    % Estimate u at end of step
   sn = s + tdelta * (un + u) / 2;	% Estimate s at end of step (Trap.)

   domdtarray(index) = domdt;       % Store current values
   uarray(index) = u;
   sarray(index) = s;
   index = index + 1;
   
   % Advance to next time step
   om = omn;
   u = un; 
   s = sn;  
end
% Store final values (As the caculated values only get asigned in the 
% arrays the next loop).
domdtarray(index) = domdt;   
uarray(index) = un;
sarray(index) = s;


%{
===========================================================================
                      | Display Results and Graphs | 
===========================================================================
%}


fprintf("Final Time: %d s\n", tf)
fprintf("Final Angular Speed (Motor): %.4f rads^-1\n", omn)
fprintf("Final Linear Speed (Motor): %.4f ms^-1\n", un)
fprintf("Final Distance: %.4f m\n", sn)

figure(1)
% Task 0 using Euler's method:  d\omega_m/dt
subplot (3,1,1), plot (0:tdelta:tf, domdtarray)
ylabel ('d\omega_m/dt for e-bike (rad.s^{-2})')
xlabel ('Time (sec)')
% Task 0 using Euler's method:  e-bike speed u
subplot (3,1,2), plot (0:tdelta:tf, uarray)
ylabel ('speed (m.s^{-1})')
xlabel ('Time (sec)')
% Task 0 using Trapezium method:  e-bike distance s
subplot (3,1,3), plot (0:tdelta:tf, sarray)
ylabel ('distance (m)')
xlabel ('Time (sec)')
