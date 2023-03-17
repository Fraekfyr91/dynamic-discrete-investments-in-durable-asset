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

pup_list = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75};
for p=1:length(pup_list)
    mp.pupgrade = pup_list{p};
    mp = trmodel.update_mp(mp);
    [sol]=equilibrium.solve(mp, s); % solve model in baseline
    for tau=1:mp.ntypes 
      ta = table(sol.ccp_tau{tau});
      txt = sprintf('ccp%.0f.xlsx', tau);
      txt = strcat(string(p), txt);
      writetable(ta,txt);
    end
end
