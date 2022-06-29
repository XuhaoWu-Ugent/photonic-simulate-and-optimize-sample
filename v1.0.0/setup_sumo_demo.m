% This Script configures Matlab for a demo of optimizing a component
% - using the SUMO toolbox
% - and 2 python simulator functions (defined in 'python_simulator')
% This module should have a function 'cheap' and a function 'expensive'
% and it should configure itself upon import

disp('Configuring Matlab for Optimization Demo');
S = dbstack('-completenames');
[current_path n e] = fileparts(S(1).file);

addpath(current_path);
addpath(strcat(current_path, '\sumo-toolbox'));

% startup sumo toolbox
startup;

% python config
luceda_path = 'C:\luceda\ipkiss_313\python\envs\ipkiss3\';
luceda_python = strcat(luceda_path, 'python.exe');
p = '.';

setenv('path',strcat(getenv('path'), ';', luceda_path, 'Scripts'));
setenv('pythonpath',p);
pyversion(luceda_python);

python_simulator = 'pythonSimulator';
