%% TEST SETUP
% this should include all of the solver-independent setup parameters to
% standardise across the various runs.

name = 'peakwave';
desc = 'Forced wave with single sharp peak';
b = pi/2;
R = 0.0;
C = 20;
L = 32;
Tmax = 64;
g = @(x) 1 - 0.0 * sin(2*pi*x/L);

v = 2.55;
f = @(x,u,t) 0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x-v*t)/L)-1)/0.001) - 1) + 1);

% name = 'wackywaves';
% desc = 'Forced waves multiple interacting peaks';
% b = 5*pi/8;
% R = 0.0;
% C = 3;
% L = 32;
% Tmax = 64;
% g = @(x) 1 - 0.0 * sin(2*pi*x/L);
% 
% v = 2.47;
% f = @(x,u,t) 0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x-v*t)/L)-1)/0.001) - 1) + 1);

% name = 'wackywaves2';
% desc = 'More forced waves multiple interacting peaks';
% b = pi/2;
% R = 0.1;
% C = 20;
% L = 32;
% Tmax = 64;
% g = @(x) 1 - 0.0 * sin(2*pi*x/L);
% 
% v = 2.55;
% f = @(x,u,t) 0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x-v*t)/L)-1)/0.001) - 1) + 1);

% name = 'peakwave';
% desc = 'Forced wave with single sharp peak';
% b = pi/2;
% R = 0.0;
% C = 20;
% L = 32;
% Tmax = 64;
% g = @(x) 1 - 0.0 * sin(2*pi*x/L);
% 
% v = 2.55;
% f = @(x,u,t) 0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x-v*t)/L)-1)/0.001) - 1) + 1);

% name = 'peakwaves';
% desc = 'Pair of forced wave with single sharp peaks';
% b = pi/2;
% R = 0.0;
% C = 20;
% L = 32;
% Tmax = 32;
% g = @(x) 1 - 0.0 * sin(2*pi*x/L);
% 
% v = 2.45;
% f_1 = @(x,u,t) 0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x-v*t)/L)-1)/0.001) - 1) + 1);
% f_2 = @(x,u,t) 0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x+L/2-v*t)/L)-1)/0.001) - 1) + 1);
% f = @(x,u,t) f_1(x,u,t) + f_2(x,u,t);

% name = 'wavefront';
% desc = 'Slow forcing drives a wavefront that eventually catches the origional (at about 32)';
% b = pi/2;
% R = 0.0;
% C = 2;
% L = 32;
% Tmax = 16;
% g = @(x) 1 - 0.0 * sin(2*pi*x/L);
% 
% v = 1;
% f = @(x,u,t) 0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x-v*t)/L)-1)/0.001) - 1) + 1);


% name = 'blockwave';
% desc = 'As for wavefront but with a larger wave';
% b = pi/2;
% R = 0.0;
% C = 2;
% L = 32;
% Tmax = 32;
% g = @(x) 1 - 0.0 * sin(2*pi*x/L);
% 
% v = 1.5;
% f = @(x,u,t) 0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x-v*t)/L)-1)/0.001) - 1) + 1);

% name = 'inversepeak';
% desc = 'Suction used to force a single narrow trough';
% b = pi/2;
% R = 0.0;
% C = 20;
% L = 32;
% Tmax = 32;
% g = @(x) 1 - 0.0 * sin(2*pi*x/L);
% 
% v = 1.5;
% f = @(x,u,t) -0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x+L/2-v*t)/L)-1)/0.001) - 1) + 1);

% name = 'trench';
% desc = 'Suction used to force a wide trench trough';
% b = pi/2;
% R = 0.0;
% C = 20;
% L = 32;
% Tmax = 32;
% g = @(x) 1 - 0.0 * sin(2*pi*x/L);
% 
% v = 2;
% f = @(x,u,t) -0.2 * (1.0127784694779528 * (exp((cos(2*pi*(x+L/2-v*t)/L)-1)/0.001) - 1) + 1);







