% This script simulates the model and returns solutions
addpath('matlabinclude');
addpath('Housetrade');
addpath('..');

clear
clc
%sum(sum(sol.ccp_tau{2},3),1)'
%% 2: Persistently heterogeneous consumer economy
close all;

% common
mp = setparams.default(); % parameters used for illustration
mp.ntypes = 8;
mp.nhousetypes=4; % switch to 1 property type
%mp.lbl_cartypes = {' '}; % no label for the only house
s = trmodel.index(mp);

% 1. no transaction cost
% 1.a set parameters
mp.transcost = 0;
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
lbl = table(mp.lbl_types);
writetable(lbl,['lbl4x8.xlsx']);
Ta = table(sol.h_tau);
writetable(Ta,'htau4x8.xlsx');
Ta = table(sol.p);
writetable(Ta, 'p4x8.xlsx')


