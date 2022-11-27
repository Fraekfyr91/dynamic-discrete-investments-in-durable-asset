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
    mp.pupgrade = 10; % price of upgrading house 
    mp.lbl_housetypes={'Luxury', 'Normal'}; % house type labels 
    mp.abar_j0 = {25,25};                 % house specific value of abar 
    mp.nhousetypes=numel(mp.abar_j0);       % number of house types 

    mp.lbl_types={'Rich'; 'Poor'};  % consumer type labels 
    mp.tw=[.5 .5]';                 % distribution of the types in the population (tw vector must sum to 1)
    mp.ntypes=numel(mp.tw);         % number of types of consumers

    mp.sigma=1.;             % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.sigma_s=0.5; 
    %mp.acc_0={-5.}; % = log(0.01), i.e. 1% accidents per year 
    %mp.acc_a={0.0};  % 

    mp.mum = {.1; .3};     % marginal utility of money
    mp.psych_transcost = {0}; 
    mp.transcost = 1.5;      % fixed transaction cost for buying car house
    mp.ptranscost=0.0;       % proportional component of transaction costs, as a fraction of the price of the house the consumer buys
    mp.u_0={6  6.5};         % intercept in utility function (ntypes x nhousetypes cell) 
    mp.u_a={-.5 -.475};      % slope coefficient on house-age in utility function (ntypes x nhousetypes cell) 

    % tax policy parameters (must be before specification of prices of fuel and cars)
    mp.vat       =   25/100;     % value added tax
    mp.housetax_lo =  105/100;   % registration tax (below kink, K_housetax_hi)
    mp.housetax_hi =  180/100;   % registration tax (above kink, K_housetax_hi)
    mp.tax_fuel  =  100/100;     % proportional fuel tax 
    mp.K_housetax_hi = 81;       % mp.K_housetax_hi before mp.housetax_hi tax kicks in
        
    % Prices before taxes (prices after taxes are computed when solving the model)
    mp.pnew={200,260};      % new house prices (housetype specific)
    mp.pscrap={1,1};        % scraphouse prices (housetype specific)
    mp.p_fuel=10.504;       % price of fuel (DKK/liter) from average fuel price in 2008 from our dataset
    [mp.p_fuel_notax, mp.pnew_notax, mp.pscrap_notax] = trmodel.price_notax(mp);
    
    % coefficent with price per kilometer by consumer type
    mp.fe       = {25,20};              % house specific fuel efficency (km/liter)
    mp.db.specification = 'linear';
    mp.db.pkm   = {-26.2549; -26.3131}; % coefficient on pkm; 
    mp.db.house = {35.74};              % house type specific coefficient on car; 
    mp.db.tau   = {-0.0968; -0.0723};   % coefficent with consumer type fixed effect by consumer type            
    mp.db.a1    = {-0.3444};            % coefficient on a*1; 
    mp.db.a2    = {0.00246};            % coefficient on a^2; 

    mp=trmodel.setparams(mp);
  end % end of setparams

end %methods

end %class
