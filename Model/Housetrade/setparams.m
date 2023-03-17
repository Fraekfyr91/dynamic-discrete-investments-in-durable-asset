% The main static class for setting paremeters 
% November 2022

classdef setparams
methods (Static)
  function mp = default()
    % Default parameters used for examples
    % SYNTAX :
    %   [mp] =setparams.default
    %
    % See also trmodel.setparams to get list of acctive parameters used in trmodel (all set to zero)
    
    mp.modeltype = 'reducedform'; % uses u_0 and u_a, not the "deep" structural parameters
    mp.pupfix = 0;
    mp.pupgrade = 0.25; % price of upgrading house 
    mp.lbl_housetypes={'LDH', 'SDH', 'LHP', 'SHP'}; % house type labels 
    mp.abar_j0 = {20, 20, 20, 20};                 % house specific value of abar 
    mp.nhousetypes=numel(mp.abar_j0);       % number of house types 

    mp.lbl_types={'R1A', 'P1A', 'R2A', 'P2A', 'R2A1C', 'P2A1C', 'R2A2C', 'P2A2C'};  % consumer type labels 
    mp.tw=[0.2702, 0.2702, 0.1422, 0.1422, 0.0392, 0.0392, 0.0484, 1 - sum(0.2702 + 0.2702 + 0.1422 + 0.1422+  0.0392 + 0.0392 + 0.0484)]';                 % distribution of the types in the population (tw vector must sum to 1)
    mp.ntypes=numel(mp.tw);         % number of types of consumers

    mp.sigma=1.;             % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.sigma_s=0.5; 
    %mp.acc_0={-5.}; % = log(0.01), i.e. 1% accidents per year 
    %mp.acc_a={0.0};  % 

    mp.mum = {3; 9; 2.5; 8; 2; 7; 1.5; 6};     % marginal utility of money
    mp.psych_transcost = {0}; 
    mp.transcost =  0;      % fixed transaction cost for buying car house
    mp.ptranscost=0.0;       % proportional component of transaction costs, as a fraction of the price of the house the consumer buys
    mp.u_0={4 4 4 4};         % intercept in utility function (ntypes x nhousetypes cell) 
    mp.u_a={-.3 -.3 -.3 -.3};      % slope coefficient on house-age in utility function (ntypes x nhousetypes cell) 
    mp.u_a_sq = {-.0, -.0 -.0, -.0 };

    % tax policy parameters (must be before specification of prices of fuel and cars)
    mp.vat       =   25/100;     % value added tax

    mp.housetax_lo =  100/100;   % registration tax (below kink, K_housetax_hi)
    mp.housetax_hi =  100/100;   % registration tax (above kink, K_housetax_hi)
    mp.tax_fuel  =  100/100;     % proportional fuel tax 
    mp.K_housetax_hi = 1;       % mp.K_housetax_hi before mp.housetax_hi tax kicks in
    
    %mp.housetax_lo =  105/100;   % registration tax (below kink, K_housetax_hi)
    %mp.housetax_hi =  180/100;   % registration tax (above kink, K_housetax_hi)
    %mp.tax_fuel  =  100/100;     % proportional fuel tax 
    %mp.K_housetax_hi = 81;       % mp.K_housetax_hi before mp.housetax_hi tax kicks in
        
    % Prices before taxes (prices after taxes are computed when solving the model)
    mp.pnew={3.2833, 2.6639, 3.4477, 2.8283};
    %mp.pnew={4.5, 6.801};% new house prices (housetype specific)
    mp.pscrap={1.3188,1.3188, 1.3188,1.3188};% scraphouse prices (housetype specific)
    %mp.p_fuel={.000000662, .000001391};       % price of fuel (DKK/liter) from average fuel price in 2008 from our dataset
    mp.p_fuel={.000000662, .000000662, .000001391, .000001391};
    [mp.p_fuel_notax, mp.pnew_notax, mp.pscrap_notax] = trmodel.price_notax(mp);
    
    % coefficent with price per kilometer by consumer type
    mp.fe       = {16000, 16000, 17000, 17000, 17500, 17500, 18000,18000};              % type specific fuel efficency (km/liter) 
    mp.size     = {165, 118, 165, 118};
    mp.pkm_par = 0.015;
    %mp.pkm_par = 0;
    mp.size_par = {0.04, 0.04, 0.05, 0.05, 0.06, 0.06, 0.07, 0.07}; 
    mp.db.specification = 'linear';
    %mp.db.pkm   = {-26.2549; -26.3131}; % coefficient on pkm; 
    %mp.db.house = {35.74};              % house type specific coefficient on car; 
    %mp.db.tau   = {-0.0968; -0.0723};   % coefficent with consumer type fixed effect by consumer type            
    %mp.db.a1    = {-0.3444};            % coefficient on a*1; 
    %mp.db.a2    = {0.00246};            % coefficient on a^2; 

    mp=trmodel.setparams(mp);

   
  end % end of setparams

  function mp = higheng()
    % Default parameters used for examples
    % SYNTAX :
    %   [mp] =setparams.default
    %
    % See also trmodel.setparams to get list of acctive parameters used in trmodel (all set to zero)
    
    mp.modeltype = 'reducedform'; % uses u_0 and u_a, not the "deep" structural parameters
    mp.pupfix = 0;
    mp.pupgrade = 0.25; % price of upgrading house 
    mp.lbl_housetypes={'Luxury', 'Normal'}; % house type labels 
    mp.abar_j0 = {20,10};                 % house specific value of abar 
    mp.nhousetypes=numel(mp.abar_j0);       % number of house types 

    mp.lbl_types={'Rich'; 'Poor'};  % consumer type labels 
    mp.tw=[.5 .5]';                 % distribution of the types in the population (tw vector must sum to 1)
    mp.ntypes=numel(mp.tw);         % number of types of consumers

    mp.sigma=1.;             % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.sigma_s=0.5; 
    %mp.acc_0={-5.}; % = log(0.01), i.e. 1% accidents per year 
    %mp.acc_a={0.0};  % 

    mp.mum = {5; 15};     % marginal utility of money
    mp.psych_transcost = {0}; 
    mp.transcost =  0;      % fixed transaction cost for buying car house
    mp.ptranscost=0.0;       % proportional component of transaction costs, as a fraction of the price of the house the consumer buys
    mp.u_0={4 4};         % intercept in utility function (ntypes x nhousetypes cell) 
    mp.u_a={-.25 -.25};      % slope coefficient on house-age in utility function (ntypes x nhousetypes cell) 
    mp.u_a_sq = {-.00, -.00};

    % tax policy parameters (must be before specification of prices of fuel and cars)
    mp.vat       =   25/100;     % value added tax

    mp.housetax_lo =  100/100;   % registration tax (below kink, K_housetax_hi)
    mp.housetax_hi =  100/100;   % registration tax (above kink, K_housetax_hi)
    mp.tax_fuel  =  100/100;     % proportional fuel tax 
    mp.K_housetax_hi = 1;       % mp.K_housetax_hi before mp.housetax_hi tax kicks in
    
    %mp.housetax_lo =  105/100;   % registration tax (below kink, K_housetax_hi)
    %mp.housetax_hi =  180/100;   % registration tax (above kink, K_housetax_hi)
    %mp.tax_fuel  =  100/100;     % proportional fuel tax 
    %mp.K_housetax_hi = 81;       % mp.K_housetax_hi before mp.housetax_hi tax kicks in
        
    % Prices before taxes (prices after taxes are computed when solving the model)
    mp.pnew={3.2833, 3.2833};
    %mp.pnew={4.5, 6.801};% new house prices (housetype specific)
    mp.pscrap={1.3188,1.3188};% scraphouse prices (housetype specific)
    %mp.p_fuel={.000000662, .000001391};       % price of fuel (DKK/liter) from average fuel price in 2008 from our dataset
    mp.p_fuel={.000000662*2, .000001391*2};
    [mp.p_fuel_notax, mp.pnew_notax, mp.pscrap_notax] = trmodel.price_notax(mp);
    
    % coefficent with price per kilometer by consumer type
    mp.fe       = {18000,18000};              % house specific fuel efficency (km/liter) 
    mp.size     = {165, 165};
    mp.pkm_par = 0.02;
    %mp.pkm_par = 0;
    mp.size_par = {0.01, 0.01};
    mp.db.specification = 'linear';
    %mp.db.pkm   = {-26.2549; -26.3131}; % coefficient on pkm; 
    %mp.db.house = {35.74};              % house type specific coefficient on car; 
    %mp.db.tau   = {-0.0968; -0.0723};   % coefficent with consumer type fixed effect by consumer type            
    %mp.db.a1    = {-0.3444};            % coefficient on a*1; 
    %mp.db.a2    = {0.00246};            % coefficient on a^2; 

    mp=trmodel.setparams(mp);
  end
  function mp = simple()
    % Default parameters used for examples
    % SYNTAX :
    %   [mp] =setparams.default
    %
    % See also trmodel.setparams to get list of acctive parameters used in trmodel (all set to zero)
    
    mp.modeltype = 'reducedform'; % uses u_0 and u_a, not the "deep" structural parameters
    mp.pupfix = 0;
    mp.pupgrade = 0.25; % price of upgrading house 
    mp.lbl_housetypes={'Luxury', 'Normal'}; % house type labels 
    mp.abar_j0 = {20,10};                 % house specific value of abar 
    mp.nhousetypes=numel(mp.abar_j0);       % number of house types 

    mp.lbl_types={'Rich'; 'Poor'};  % consumer type labels 
    mp.tw=[.5 .5]';                 % distribution of the types in the population (tw vector must sum to 1)
    mp.ntypes=numel(mp.tw);         % number of types of consumers

    mp.sigma=1.;             % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.sigma_s=0.5; 
    %mp.acc_0={-5.}; % = log(0.01), i.e. 1% accidents per year 
    %mp.acc_a={0.0};  % 

    mp.mum = {5; 15};     % marginal utility of money
    mp.psych_transcost = {0}; 
    mp.transcost =  0;      % fixed transaction cost for buying car house
    mp.ptranscost=0.0;       % proportional component of transaction costs, as a fraction of the price of the house the consumer buys
    mp.u_0={4 4};         % intercept in utility function (ntypes x nhousetypes cell) 
    mp.u_a={-.25 -.25};      % slope coefficient on house-age in utility function (ntypes x nhousetypes cell) 
    mp.u_a_sq = {-.00, -.00};

    % tax policy parameters (must be before specification of prices of fuel and cars)
    mp.vat       =   25/100;     % value added tax

    mp.housetax_lo =  100/100;   % registration tax (below kink, K_housetax_hi)
    mp.housetax_hi =  100/100;   % registration tax (above kink, K_housetax_hi)
    mp.tax_fuel  =  100/100;     % proportional fuel tax 
    mp.K_housetax_hi = 1;       % mp.K_housetax_hi before mp.housetax_hi tax kicks in
    
    %mp.housetax_lo =  105/100;   % registration tax (below kink, K_housetax_hi)
    %mp.housetax_hi =  180/100;   % registration tax (above kink, K_housetax_hi)
    %mp.tax_fuel  =  100/100;     % proportional fuel tax 
    %mp.K_housetax_hi = 81;       % mp.K_housetax_hi before mp.housetax_hi tax kicks in
        
    % Prices before taxes (prices after taxes are computed when solving the model)
    mp.pnew={3.2833, 3.2833};
    %mp.pnew={4.5, 6.801};% new house prices (housetype specific)
    mp.pscrap={1.3188,1.3188};% scraphouse prices (housetype specific)
    %mp.p_fuel={.000000662, .000001391};       % price of fuel (DKK/liter) from average fuel price in 2008 from our dataset
    mp.p_fuel={.000000662, .000001391};
    [mp.p_fuel_notax, mp.pnew_notax, mp.pscrap_notax] = trmodel.price_notax(mp);
    
    % coefficent with price per kilometer by consumer type
    mp.fe       = {18000,18000};              % house specific fuel efficency (km/liter) 
    mp.size     = {165, 165};
    mp.pkm_par = 0.02;
    %mp.pkm_par = 0;
    mp.size_par = {0.01, 0.01};
    mp.db.specification = 'linear';
    %mp.db.pkm   = {-26.2549; -26.3131}; % coefficient on pkm; 
    %mp.db.house = {35.74};              % house type specific coefficient on car; 
    %mp.db.tau   = {-0.0968; -0.0723};   % coefficent with consumer type fixed effect by consumer type            
    %mp.db.a1    = {-0.3444};            % coefficient on a*1; 
    %mp.db.a2    = {0.00246};            % coefficient on a^2; 

    mp=trmodel.setparams(mp);
  end

end %methods

end %class
