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
mp = setparams.simple(); % parameters used for illustration
mp.ntypes = 2;
mp.nhousetypes=1; % switch to 1 property type
mp.lbl_cartypes = {' '}; % no label for the only house
s = trmodel.index(mp);

% 1. no transaction cost
% 1.a set parameters
mp.transcost = 0;
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
lbl = table(mp.lbl_types);
writetable(lbl,['lbl.xlsx']);
Ta = table(sol.h_tau);
writetable(Ta,'htau.xlsx');
a = sol.p;

% 1.a set parameters
mp.ntypes = 1;
mp.tw = [1]';
mp.mum = {5}; % only rich
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
b = sol.p;

% 1.a set parameters
mp.mum = {15}; % only poor
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
c = sol.p;
Ta = table([a b c]);
writetable(Ta, 'p.xlsx')




%  .1 transaction cost
mp.ntypes = 2;
mp.tw = [.5 .5]';
mp.mum = {5; 15};
% 1.a set parameters
mp.transcost = .25;
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
a = sol.p;

Ta = table(sol.h_tau);
writetable(Ta,'htau_tc.xlsx');


mp.ntypes = 1;
mp.tw = [1]';
mp.mum = {5}; % only rich
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
b = sol.p;

mp.mum = {15}; % only poor
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
c = sol.p;
Ta = table([a b c]);
writetable(Ta,'p_tc.xlsx');




%  .1 transaction cost
mp.ntypes = 2;
mp.tw = [.5 .5]';
mp.mum = {5; 15};
% 1.a set parameters
mp.transcost = 0;
mp.tc_sale = .25;
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
a = sol.p;

Ta = table(sol.h_tau);
writetable(Ta,'htau_stc.xlsx');


mp.ntypes = 1;
mp.tw = [1]';
mp.mum = {5}; % only rich
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
b = sol.p;

mp.mum = {15}; % only poor
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
c = sol.p;
Ta = table([a b c]);
writetable(Ta,'p_stc.xlsx');


%  .1 transaction cost
mp.ntypes = 2;
mp.tw = [.5 .5]';
mp.mum = {5; 15};
% 1.a set parameters
mp.transcost = .25;
mp.tc_sale = .25;
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
a = sol.p;

Ta = table(sol.h_tau);
writetable(Ta,'htau_bstc.xlsx');


mp.ntypes = 1;
mp.tw = [1]';
mp.mum = {5}; % only rich
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
b = sol.p;

mp.mum = {15}; % only poor
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
c = sol.p;
Ta = table([a b c]);
writetable(Ta,'p_bstc.xlsx');


%  mp.pupgrade 2
mp.ntypes = 2;
mp.tw = [.5 .5]';
mp.mum = {5; 15};
% 1.a set parameters
mp.transcost = 0;
mp.pupgrade =.5;
mp.pupfix = 0;
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
a = sol.p;
Ta = table(sol.h_tau);
writetable(Ta,'htau_pup.xlsx');


mp.ntypes = 1;
mp.tw = [1]';
mp.mum = {5}; % only rich
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
b = sol.p;

mp.mum = {15}; % only poor
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline 
c = sol.p;
Ta = table([a b c]);
writetable(Ta,'p_pup.xlsx');

%  .1 energy cost
mp = setparams.higheng();
mp.ntypes = 2;
mp.nhousetypes=1; % switch to 1 property type
mp.lbl_cartypes = {' '}; % no label for the only house
s = trmodel.index(mp);
mp.ntypes = 2;
mp.tw = [.5 .5]';
mp.mum = {5; 15};
% 1.a set parameters

mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
a = sol.p;
Ta = table(sol.h_tau);
writetable(Ta,'htau_eng.xlsx');
mp.ntypes = 1;
mp.tw = [1]';
mp.mum = {5}; % only rich
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
b = sol.p;
mp.mum = {15}; % only poor
mp = trmodel.update_mp(mp);
% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
c = sol.p;
Ta = table([a b c]);
writetable(Ta,'p_eng.xlsx');



