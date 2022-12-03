% This script simulates the model and returns solutions
addpath('matlabinclude');
addpath('Housetrade');
addpath('..');

clear
clc

%% 2: Persistently heterogeneous consumer economy
close all;

% common
mp = setparams.default(); % parameters used for illustration
mp.ntypes = 2;
mp.nhousetypes=2; % switch to 1 property type
mp.lbl_cartypes = {' '}; % no label for the only house
s = trmodel.index(mp);

% 1. Normal transaction cost

% 1.a set parameters
mp.transcost = 0;
mp = trmodel.update_mp(mp);

% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline
Ta = table(sol.h_tau);
writetable(Ta,'htau.xlsx');
lbl = table(mp.lbl_types);
writetable(lbl,['lbl.xlsx']);

% 1.c graphs
%f = graphs.agg_holdings(mp,s,sol);
%ylim = get(gca, 'ylim'); 
%legend('off');
%saveas(f, 'results/illustration/example_holdings_tc_0.eps', 'epsc');
%f1 = graphs.prices_with_hom(mp,s,sol);
%legend('off');
%saveas(f1, 'results/illustration/example_prices_tc_0.eps', 'epsc');

%% 2. With high transaction costs
%
%% 2.a set parameters
%mp.transcost=10;
%mp.psych_transcost = {0}; 
%mp=trmodel.update_mp(mp); % standard model parameters
%
%% 2.b solve model
%[sol]=equilibrium.solve(mp, s); % solve model in baseline
%
%% 2.c graphs
%f = graphs.agg_holdings(mp,s,sol);
%set(gca, 'ylim', ylim);
%legend('location', 'north'); 
%saveas(f, 'results/illustration/example_holdings_tc_10.eps', 'epsc');
%f2 = graphs.prices_with_hom(mp,s,sol);
%s  aveas(f2, 'results/illustration/example_prices_tc_10.eps', 'epsc');



