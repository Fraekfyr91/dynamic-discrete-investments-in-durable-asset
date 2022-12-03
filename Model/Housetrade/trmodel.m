% The main static class for the HOUSE TRADE model.
% Contains all needed model parts.
% 
% November 2022
%

classdef trmodel

methods (Static)
  function [mp] =setparams(mp0)
    % Standard parameters of the model
    % SYNTAX :
    %   [mp] = trmodel.setparams       set all parameters to zeros
    %   [mp] = trmodel.setparams(mp0)  update patameters with parameters in mp0
    %
    % NOTE: When adding new parameters to trmodel, they must be added here 
    % DO NOT ENTER VALUES OF PARAMETERS - JUST ZEROS (or values that turn off the parameter)
    % See also setparams.default to get a default parameters

    % ******************************************
    % switches
    % ******************************************
    mp.es=1;          % use model with endogenous scrapage if 1 
    mp.fixprices=0;   % set to zero, to compute model solution at a given set of prices without solving for equilibrium 
    mp.ll_scrap=true; % add scrap decisions to likelihood if true
    % ******************************************
    % misc. parameters
    % ******************************************
    %mp.dk_population =2.5;     % number of danish households in millions

    % ******************************************
    % Consumer types and House types
    % ******************************************
    % Some parameters below are consumer type specific (in rows) and House type specific (in columns) 
    % If number of entries is smaller than ntypes, last entry is repeated)

    mp.lbl_housetypes={'Luxury', 'Normal'}; % house type labels 
    mp.abar_j0 = {25,25};                   % House specific value of abar 
    mp.nhousetypes=numel(mp.abar_j0);       % number of House types 

    mp.lbl_types={'Rich', 'Poor'};  % consumer type labels
    mp.tw=[.5 .5]';                 % distribution of the types in the population (tw vector must sum to 1)
    mp.ntypes=numel(mp.tw);         % number of types of consumers

    % ******************************************
    % discount factor (beta)
    % ******************************************
    mp.bet=.95;  
   
    % ******************************************
    % parameters of GEV distribution (theta)
    % ******************************************
    mp.sigma=1;         % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.sigma_s=1;       % the extreme value scale parameter for the idiosyncratic extreme value shocks affecting the sale vs scrappage decision
 
    % ******************************************
    % utility patameters (theta)
    % ******************************************
    mp.mum = {1};   % marginal utility of money
    mp.phi = {0.01}; % scale of heating utilty (set low to exclude heating)

    % transactions costs parameters 
    mp.transcost =0;      % fixed transaction cost for byuing a house
    mp.ptranscost=0;      % proportional component of transaction costs, as a fraction of the price of the house the consumer buys
    mp.psych_transcost      = {0};      % utility cost of bying a house
    mp.psych_transcost_nohouse= {0};    % additional utility cost of bying a house, when not owning a house
    mp.nohousesc   =0;  % additional search cost incurred by a consumer who has no house
    mp.tc_sale     =0;    % set the fixed component of the transaction cost to a seller (cleaning/painting/realestate cost to make a house qualified to be resold)
    mp.tc_sale_age =0;    % coefficient of age on the total sales costs
    mp.ptc_sale    =0;    % set the proportional component of the transaction cost to a seller (cleaning/painting/realestate cost to make a house qualified to be resold)
    %mp.tc_sale_even=0;    % incremental sales transactions costs during inspection years (even years) necessary to pass inspection before car can be sold

    % house utility parameters (see also trmodel.u_house) 
    % Reduced form house utility mp.u_0{tp, house}+mp.u_a{tp, house}*house_age + mp.u_a_sq{tp, house}*house_age.^2
    mp.u_0={0};         % intercept in utility function (ntypes x nhousetypes cell) 
    mp.u_a={0};         % slope coefficient on house-age in utility function (ntypes x nhousetypes cell) 
    mp.u_a_sq = {0};    % house age squared 

    mp.convexutility=0; % this is a switch, if 1 then mp.u.a_sq is forced to be positive via the transformation exp(mp.u_a_sq) in u_house

    mp.u_og={0};        % utilty of outside good  (ntypes x 1 cell) 
    %mp.u_even={0};      % utility during even inspection years (expect to be a disutility or negative value due to the disutility of having car inspected)        
    % tax policy parameters (must be before specification of prices of fuel and houses)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%FIX HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mp.vat           =  0;   % value added tax
    mp.housetax_lo   =  0;   % registration tax (below kink, K_housetax_hi) %
    mp.housetax_hi   =  0;   % registration tax (above kink, K_housetax_hi)
    mp.tax_fuel      =  0;   % proportional fuel tax 
    mp.K_housetax_hi =  0;   % mp.K_housetax_hi before mp.housetax_hi tax kicks in

    % ******************************************
    % Prices before taxes (prices after taxes are computed when solving the model)
    % if house or consumer type parameters have size smaller than nhousetypes and ntypes, the last element is repeated.
    % ******************************************
    mp.pnew_notax   = {100};  % new house prices (housetype specific)
    mp.pscrap_notax = {0};    % scraphouse prices (housetype specific)
    mp.p_fuel_notax = 5;      % pricer of fuel (DKK/liter) from average fuel price before tax
    
    % ******************************************
    % Parameters of reduced form heating equation 
    % x = mp.db.pkm{tp}*pkm{car} + mp.db.car{car}+ mp.db.tau{tp} +mp.db.a1{tp}*a+mp.db.a2{tp}*a.^2;
    % ******************************************

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%FIX HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Reduced form driving equation coefficients are stored in db structure 
    mp.fe={20};       % house specific fuel efficency (km/liter)
    mp.db.specification = 'linear'; % 'linear' or 'log-log' 
    mp.db.pkm = {0};  % coefficient on pkm; 
    mp.db.pkm = {0};  % coefficient on pkm; 
    mp.db.house = {0};  % house type specific coefficient on house; 
    mp.db.tau = {0};  % coefficent with house type fixed effect by consumer type            
    mp.db.a1  = {0};  % coefficient on a*1; 
    mp.db.a2  = {0};  % coefficient on a^2; 
    
    % ******************************************
    % Reduced form or structural form
    % ******************************************
    % Structural parameters can be identified from reduced from parameters
    % for house demand and heating 
    % by using trmodel.update_structural_par to obtain "structural" utility parameters mp.sp. 
    % Because of implicit dependence on type specific marginal utility of money, all these coefficients have to be consumer type specific
    % coefficent with price per kilometer by consumer type
    
    mp.modeltype = 'reducedform'; 
    % if mp.modeltype == 'reducedform': reudced form parameters mp.u_0, mp.u_a, and mp.u_a_sq are used in trmodel.u_house. 
    % if mp.modeltype == 'structuralform': structural parameters mp.sp,
    % mp.mum are used in trmodel.u_house
    % 
    % If the model allows for driving, the reduced form parameters mp.u_0, mp.u_a, and mp.u_a_sq 
    % are not policy invariant since they on fuel-prices (and marginal utility of money and other structural parameters)
    % So to run counterfactuals where fuelprices are changed you need to set mp.modeltype = 'structuralform';
    %
    % Estimation strategy: 
    %   First set mp.modeltype = 'reducedform' to estimate reduced form parameters during estimation. 
    %   Then set mp.modeltype = 'structuralform' for running counter-factuals that changes fuel prices. 

    % ********************************************************************************
    % update parameters
    % ********************************************************************************
    % update mp with mp0 (default values are over written by input values mp0)
    if nargin>0
      mp.db=combinestructs(mp.db, mp0.db);
      mp=combinestructs(mp, rmfield(mp0,'db'));
    end

    % update endogenous parameters
    mp=trmodel.update_mp(mp);
    [mp.p_fuel_notax, mp.pnew_notax, mp.pscrap_notax] = trmodel.price_notax(mp);
  end % end of setparams

  function [p_fuel_notax, pnew_notax, pscrap_notax] = price_notax(mp)

    % This function computes prices before taxes given after tax prices and tax rates in mp
    % When implementing new taxes/subsidies remember to update the corresponding after tax function below. 

    % fuel prices before taxes
    p_fuel_notax=mp.p_fuel/(1+mp.tax_fuel);  
          
    pnew_notax = cell(1, mp.nhousetypes); 
    for house=1:mp.nhousetypes 

      pnew_notax{house} = trmodel.phouse_notax(mp.pnew{house}, mp.K_housetax_hi, mp.housetax_lo, mp.housetax_hi, mp.vat);
    
      % scrap price before tax: not currently implemented 
      if nargout >=3 
        pscrap_notax{house} = mp.pscrap{house};
      end
    end 
    
  end

  function pnew_notax = phouse_notax(phouse, kink_threshold, housetax_low, housetax_high, vat)
    % phouse_notax(): new house prices before registration taxes and VAT. 
    %
    % INPUTS: (all scalar floats)
    %   phouse: consumer price of the house 
    %   kink_threshold: the threshold where the registration tax rate shifts 
    %   housetax_high: marginal tax rate above kink
    %   housetax_low: marginal tax rate below kink
    %   vat: Value Added Tax (should be in [0; 1).)
    % 
    % OUTPUT: 
    %   pnew_notax: house price without any taxes

    assert(housetax_low <= housetax_high, 'Low-bracket rate is above high-bracket rate: this sounds like an error (unless you are doing a fancy counterfactual, then comment this assertion out)');
    assert(vat < 1.0, 'Expecting VAT = 0.25 (or at least < 1.0). Delete this assertion if you are trying a crazy counterfactual and know what you are doing.')
    
    % cutoff = tax amount paid if a house price is precisely equal to the kink threshold value
    
    % how much is paid if the pre-tax house price (with VAT) puts it
    % precisely at the cutoff 
    price_paid_at_cutoff = (1+housetax_low)*kink_threshold;

    % compute the final inversion separately depending on whether we are
    % above or below the cutoff 
    if phouse <= price_paid_at_cutoff 
        pnew_notax = phouse / ((1+housetax_low)*(1+vat)); 
    else 
        numer = phouse + (housetax_high - housetax_low)*kink_threshold;
        denom = (1+vat)*(1+housetax_high); 
        pnew_notax = numer/denom; 
    end
        
  end

  function [p_fuel, pnew, pscrap] = price_aftertax(mp)
    % price_aftertax(): This function computes prices after taxes given before tax prices and tax rates in mp
    % When implemneting new taxes/subisdies remember to update the correspoding after tax function below. 

    p_fuel  = mp.p_fuel_notax*(1+mp.tax_fuel); % scalar 
          
    pnew    = cell(1, mp.nhousetypes); 
    pscrap  = cell(1, mp.nhousetypes); 
    
    for ihouse=1:mp.nhousetypes

      pnew_notax = trmodel.phouse_after_passthrough(mp.pnew_notax{ihouse}, mp, ihouse); 

      pnew{ihouse} = trmodel.phouse_aftertax(pnew_notax, mp.K_housetax_hi, mp.housetax_lo, mp.housetax_hi, mp.vat);

      if nargout >= 3
          pscrap{ihouse} = mp.pscrap_notax{ihouse}; % currently, no special treatment 
      end
    end
  end

  function passthrough = set_up_passthrough(mp_baseline, rate)

    passthrough = struct(); 
    passthrough.pnew_notax_baseline = mp_baseline.pnew_notax;
    passthrough.pnew_baseline = mp_baseline.pnew; 
    passthrough.rate = rate; 
    
    passthrough.housetaxes_baseline = struct(); 
    passthrough.housetaxes_baseline.K_housetax_hi   = mp_baseline.K_housetax_hi; 
    passthrough.housetaxes_baseline.housetax_lo     = mp_baseline.housetax_lo; 
    passthrough.housetaxes_baseline.housetax_hi     = mp_baseline.housetax_hi; 
    passthrough.housetaxes_baseline.vat             = mp_baseline.vat; 
    
  end

  function phouse = phouse_after_passthrough(phouse_raw, mp, ihouse, DOPRINT)
      % phouse_after_passthrough(): adds a mechanical firm-response to the
      % raw pre-tax price. The amount to add is set so that a target
      % passthrough rate is achieved 
      %
      % INPUTS: 
      %     phouse_raw (double): house price, pre-tax and pre-passthrough 
      %     mp (struct): model parameters 
      %     ihouse (int): index for the house 
      %     DOPRINT (bool, optional): print details 
      % 
      % OUTPUT: 
      %     phouse (double): house price after firm "markup" has adjusted 
      % 
      
      if nargin < 4
          DOPRINT = false;
      end
      
      if isfield(mp, 'passthrough')
          assert(ihouse>=1 && ihouse<=mp.nhousetypes);
          assert(isfield(mp.passthrough, 'pnew_notax_baseline'));
          assert(isfield(mp.passthrough, 'rate'));
          
          % 1. compute price that *would* have prevailed absent any changes
          % in firm behavior 
          phouse_at_full_passthrough = trmodel.phouse_aftertax(phouse_raw, mp.K_housetax_hi, mp.housetax_lo, mp.housetax_hi, mp.vat);
          mp0 = mp.passthrough.housetaxes_baseline; % tax rates in the baseline
          phouse_baseline_full_pasthrough = trmodel.phouse_aftertax(phouse_raw, mp0.K_housetax_hi, mp0.housetax_lo, mp0.housetax_hi, mp0.vat);
          delta_tax = phouse_at_full_passthrough - phouse_baseline_full_pasthrough;
          
          % 2. change in manufacturer price
          % NOTE: "-1" because manufacturers move *opposite* of the policy makers intended direction
          delta_firm_price = (-1) * (1-mp.passthrough.rate) * delta_tax;
          
          % 3. final price before taxes get applied
          phouse = phouse_raw + delta_firm_price;
          
          if DOPRINT
              fprintf('At full passthrough: p0 = %5.2f to p1 = %5.2f: implied delta tax = %5.2f\n', phouse_baseline_full_pasthrough, phouse_at_full_passthrough, delta_tax); 
              fprintf('Requested passthrough = %5.2f%% => delta p raw = %5.2f\n', 100.0*mp.passthrough.rate, delta_firm_price); 
              fprintf('Final result: p0 = %5.2f -> p with passthrough = %5.2f\n', phouse_raw, phouse); 
          end
         
      else
          % nothing to do
          phouse = phouse_raw;
      end
  end

  function phouse_incl_taxes = phouse_aftertax(phouse_notax, K_housetax_hi, housetax_lo, housetax_hi, vat)
      % phouse_notax(): new house prices before taxes
      %
      % INPUTS: 
      %   phouse_notax: consumer price of the house
      %   K_housetax_hi: the threshold where the registration tax rate shifts 
      %   housetax_hi: marginal tax rate above kink
      %   housetax_lo: marginal tax rate below kink
      %   vat: Value Added Tax 
      % 
      % OUTPUT: 
      %   phouse_incl_taxes: house price including all taxes

      assert(housetax_lo <= housetax_hi, 'hi/low bracket rates reversed: this could be an error!'); 
      assert(vat <= 1.0, 'VAT should be in [0;1] (unless you are doing a crazy large counterfactual)');
      
      if phouse_notax*1.25 <= K_housetax_hi
          % no top-tax will be paid 
          phouse_incl_taxes = (1+housetax_lo)*1.25*phouse_notax; 
      else
          % price is in the top-bracket 
          phouse_incl_taxes = (1+housetax_hi)*1.25*phouse_notax - (housetax_hi - housetax_lo)*K_housetax_hi; 
      end
      
  end

  function [s] = index(mp, abar_j)
    % This function set up the state space for the trmodel. 
    % 
    % SYNTAX: [s] = trmodel.index(mp, abar_j)
    % 
    % INPUT: 
    %   mp:     structure with model parameters (see trmodel.setparams)
    %   abar_j: mp.nhousetypes cell array with max age (forced scrap age) for each house j=1,..,mp,nhousetypes
    %
    % OUTPUT:  
    %   s:  a struct with the following elements:
    %
    %     id: struct that holds decision indexes groups of decisions. 
    %     is: struct that holds state indexes groups of states
    %     ip: mp.nhousetypes x 1 cell of indexes used house prices. 
    %     ns: number of states
    %     nd: number of decisions
    %     np: number of price parameters
    %     tr: struct with sub-indexes for transitions
    %     ipt: struct that holds indexes for the post trade distribution
    %
    %  Example: 
    %  s=trmodel.index(mp, {16,22}) gives
    %   s = 
    %     struct with fields:
    %       abar_j: {[16]  [22]}
    %           id: [1x1 struct]
    %           is: [1x1 struct]
    %           ip: {[1 2 ... 15]  [16 17 ... 36]}
    %           ns: 39
    %           nd: 40
    %           np: 36
    %           tr: [1x1 struct]
    %          ipt: [1x1 struct]
    %  where
    % 
    %   s.id has the following fields
    %                keep: 1
    %               trade: {[2 3 ... 17]  [18 19 ... 39]}
    %          trade_used: {[3 4 ...  17]  [19 20 ...  39]}
    %           trade_new: {[2]  [18]}
    %                 age: [NaN 0 1 ...  15 0 1 2 ...  21 NaN]
    %               purge: 40
    %    
    %   s.is has the following fields
    %                house: {[1 2  ... 16]  [17 18 ... 38]}
    %       house_ex_scrap: {[1 2  ... 15]  [17 18 ... 37]}
    %                scrap: {[16]  [38]}
    %                  age: [1 2 ... 16 1 2 ... 22 NaN]
    %              nohouse: 39
    %       
    %               choice: [3 4 ... 17 2 19 20 ... 39 18 40]
    %            choice_up: [3 4 ... 17 2 19 20 ... 39 18 40]
    %                state: [1 2 ...  39]
    %  
    %  s.ip mp.nhousetypes x 1 cell with the following elements
    %       s.ip{1} = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
    %       s.ip{2} = [16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36]
    %  
    %  s.ipt has the following fields
    %       house   : {[1 2 ... 16]  [17 18 ... 38]}
    %       nohouse : 39
    %       new     : {[1]  [17]}
    %       used    : {[2 3 ... 16]  [18 19 ... 38]}
    %       age     : [0 1 ... 15 0 1 2 ... 21 NaN]   
    % 
    %  usage: 
    %   s.id.keep, s.id.trade{j},s.trade_used{j}, s.trade_new{j}, and s.id.purge are choice indexes that 
    %   for example used for indexing columns in matrices such as ccp_tau, and decision specific utility
    %
    %   s.id.upgrade{j}, s.upgrade_used{j}, s.upgrade_new{j}, and s.id.purge are choice indexes that 
    %   for example used for indexing columns in matrices such as ccp_tau, and decision specific utility
    %
    %   s.is.house{j}, s.is.house_ex_scrap{j},s.is.scrap{j}, and s.is.nohouse are state indexes that 
    %   for example used for indexing rows in matrices such as ccp_tau, ev_tau, q, q_tau{tau}.
    %
    %   For example, ccp_tau{tau}(s.is.house_ex_scrap{j}, s.id.keep, s.id.keep) gives the vector of probabilities 
    %   of keeping and not upgrading a house of type j with ages 1,2,..,abar_j{j}-1. The age of such a house can also be accessed by 
    %   s.is.age(s.is.house_ex_scrap{j})  
    %
    % Some example descriptions of subindices: 
    % 
    % s.ip  for each house type, the indices for each traded house age (1,....abar_j{i}-1), i=1...nhousetypes sequentially ordered
    %       for that the union of these indices over all house types is the sequence 1,2,....,abar_j{1}+...abar_j{nhousetype}-nhousetypes, 
    %
    % s.is.houes (mp.nhousestypes cell array): 
    %   s.is.house{j} is a vector holding the indexes in the overall ordering of house states for house type j - except the no house state. 
    % s.is.house_ex_clunker (mp.nhousetypes cell array): 
    %   s.is.house_ex_clunker{j} is a vector holding the indices in the overall ordering of house states for house type j
    %    - except the oldest house (Scrap) and the no house state. 
    % s.is.clunker (mp.nhousetypes cell array): 
    %   s.is.scrap{j} is a scalar holding the index for the oldest house of type j
    % s.is.nohouse (scalar): 
    %   index for the state of having no house
    % s.is.age (s.ns vector): 
    %   s.is.age vector of house ages corresponding to state indexes. Age of no house is coded as NaN


    if nargin==1
        s.abar_j=mp.abar_j0; 
    else
        s.abar_j=abar_j;
    end
 
    % state and decision indexes
    s.id.keep=1;    % keep decision
    ei=0;
    eip=0;
    for j=1:mp.nhousetypes
      % indexes states and trade decisions (including new house)
      si=ei+1;
      ei=ei+s.abar_j{j};

      % house owner state
      s.is.house{j}=[si:ei];
      s.is.house_ex_scrap{j}=s.is.house{j}(1:end-1);
      s.is.scrap{j}=s.is.house{j}(end);

      % trade decision 
      s.id.trade{j}=1+(si:ei); % add 1 to account for keep beeing the first index
      s.id.trade_used{j}=s.id.trade{j}(2:end);
      s.id.trade_new{j}=s.id.trade{j}(1);

      % upgrade decision

      s.id.upgrade{j} = 1+ (si:ei);
      s.id.upgrade_used{j} = s.id.upgrade{j}(2:end);
      s.id.upgrade_new{j} = s.id.upgrade{j}(1);
      % corresponding house age values
      s.is.age(s.is.house{j})=(1:s.abar_j{j})';   % not possible to enter period with new house
      s.id.age(s.id.trade{j})=(0:s.abar_j{j}-1)'; % not possible to buy scrap
      
      %indexes for used house prices and excess demand
      sip=eip+1;
      eip=eip+s.abar_j{j}-1; 
      s.ip{j}=(sip:eip); 
    end

    s.ns=ei+1;             % number of states      (add one for no house state)
    s.nd=ei+2;             % number of decisions   (add one for keep and purge)
    s.np=eip;              % number prices          

    s.is.nohouse=s.ns;    % nohouse state
    s.id.purge=s.ns+1;    % purge decision

    s.is.age(s.is.nohouse)=nan;
    s.id.age(s.id.keep)= nan;
    s.id.age(s.id.purge)=nan;

    % Indices for trade transition matrix
    % build index for chocies correnspoding to columns in deltaT matrix (see definition of deltaK in paper)
    s.tr.choice=[];
    s.tr.choice_up = [];
    for j=1:mp.nhousetypes; % loop over house
        % first used, then new
        s.tr.choice=[s.tr.choice s.id.trade{j}(2:end) s.id.trade{j}(1)]; 
        s.tr.choice_up =[s.tr.choice_up s.id.upgrade{j}(2:end) s.id.upgrade{j}(1)];
    end
    s.tr.choice= [s.tr.choice s.id.purge];  % purge in last column 
    s.tr.choice_up= [s.tr.choice_up s.id.purge];
  
    % index for all rows and columns in delta
    s.tr.state=[s.is.house{:} s.is.nohouse]; 

    % index for next period house state 
    s.tr.next_keep=s.tr.state+1;   % state after keep (current house ages one period)
    s.tr.next_trade=s.tr.state;    % state after trading
    s.tr.next_upgrade = s.tr.state;

    for j=1:mp.nhousetypes;
      s.tr.next_keep([s.is.scrap{j}])=[s.is.house{j}(1)];
    end

    s.tr.next_keep(s.is.nohouse)=s.is.nohouse;

    % Age transition matrices
    s.Q  = sparse(s.tr.state,s.tr.next_keep,1,s.ns,s.ns);
    s.dQ = s.Q; 
    % Trade and upgrade transition matrices
    s.F  = sparse(s.tr.state,s.tr.state,1,s.ns, s.ns);
    s.U  = sparse( s.tr.state, s.tr.state,1, s.ns, s.ns);
    s.dF = s.F; 
    s.dU = s.U;

    % indexes for post trade distribution (before aging)
    s.ipt.house=s.is.house; 
    s.ipt.nohouse=s.is.nohouse; 
    for j=1:mp.nhousetypes;
      s.ipt.new{j}=s.ipt.house{j}(1);
      s.ipt.used{j}=s.ipt.house{j}(2:end);
    end 
    s.ipt.age=s.is.age -1; 

    % indexes for post upgrade distributions
    %s.ipu.house = s.is.house;
    %s.ipu.nohouse = s.is.nohouse;
    %for j=1:mp.nhousetypes;
    %    s.ipu.new{j}= s.ipu.house{j}(1);
    %    s.ipu.used{j} = s.ipu.house{j}(2:end);
    %end
  end % end of index

  function [u, ev_scrap, ccp_scrap]=utility(mp, s, tau, price_vec)
    %   utility:  Function to compute utility for each consumer in each state and for each decision 
    %             and its derivative wrt to used house prices. 
    %
    %   SYNTAX:     [u, ev_scrap, ccp_scrap]=trmodel.utility(mp, s, price_j)
    %
    %   INPUTS:
    %     mp          structure with model parameters (see trmodel.setparams)
    %     s:          structure with indexes for dimensions. see trmodel.index
    %     tau:        household type 
    %     price_j0    (optional) starting values for price_j. price_j0 is a mp.nhousetypes dimensional cell array 
    %
    %
    %   OUTPUTS          
    %   u:      a cell array of dimension mp.ntypes whose elements are matrices that store the utility
    %           received from trading for a house of type j and age a for a house of type j' and age d, or "purge"
    %           trade the current house for the outside good, or to keep the house, for each consumer type.
    %
    %           u_tau{t} has s.ns rows and s.nd columns and rows is indexed by s.is and columns by s.id. 
    %
    %           Rows correspond to the different states 
    %           i.e the type/age of the currently held house or the state of no house. 
    %           The last row corresponds to the state of having no house, ie. outside good.
    %
    %           The columns correspond to the different possible decisions a consumer has.
    %           The different possible decisions are: keep, trade and purge. 
    %
    %           The first column is the utility of the decision to keep (s.id.keep=1) 
    %           for each possible age of house (or outside good in last column)
    %
    %           NOTE: The utility of keep contain a missing value in in rows where the age equals the scrappage age (is.scrap{j})
    %           for each house type, since a consumer is not allowed to keep a house once it reaches the scrappage age, and for a consumer
    %           who has no house, the law row of the u_tau provides the utility of continuing to stay with the outside good (keeping outside good)
    % 
    %           The last column, s.id.purge=s.nd, is the utility of the decision to purge 
    %           (i.e. choose the outside good, i.e. move to having no house)
    %
    %           The columns in between, from 2 to s.nd-1 or [s.id.trade{:}], correspond to the decision to trade 
    %           for a house of age a, a=0,...,abar{j}-1 for each of the
    %           nhousetypes of houses
    % 
    %  ev_scrap  the extra option value from being able to scrap instead of sell a used house (only non-zero if mp.es=1)
    %            This is the s.ns x 1 vector of integrated value of scrap choice relative to selling choice
    %            i.e. logsum(mp.mum{tau}*(pscrap-psell), 0)  
    %  ccp_scrap: 
    %           s.ns x 1 vector of scrapping probabilities (conditional on not keeping)
    %
    %           Gillingham, Iskhakov, Munk-Nielsen, Rust and Schjerning, September 2020

    % update the prices in the vector price_vec to the cell array price_j
    for i=1:mp.nhousetypes;
      price_j{i}=price_vec(s.ip{i});
    end;

    % net sales price after transactions cost 
    [psell]=trmodel.psell(mp, s, price_j); 
    
    % net purchase price after transactions cost 
    [pbuy]=trmodel.pbuy(mp, s, price_j); 
    
    u=zeros(s.ns,s.nd,s.nd); % utility: number of states by number of decisions (see state_decision_index for information about of indexing)

    % utility of scrap option (when selling)
    [ccp_scrap, ev_scrap]=trmodel.p_scrap(mp, s, price_j, tau);
    
    
    
    for j=1:mp.nhousetypes
      % utility of keeping (when feasible)
      u(s.is.house_ex_scrap{j}, s.id.keep, s.id.keep) =  trmodel.u_house(mp, s.is.age(s.is.house_ex_scrap{j}), tau, j);
      
      % utility of trading and not upgrading
      u(: , s.id.trade{j}, s.id.keep)                 =  trmodel.u_house(mp, s.id.age(s.id.trade{j}), tau, j) ...
                                      + ev_scrap-mp.mum{tau}*(pbuy(:,s.id.trade{j})-psell)-mp.psych_transcost{tau};

      % utility of trading and upgrading
      B = trmodel.u_house(mp, s.id.age(s.id.upgrade{j}), tau, j) ...
                                      + ev_scrap-mp.mum{tau}*mp.pupgrade * flipud(s.id.age(s.id.upgrade{j})')' -mp.mum{tau}*(pbuy(:,s.id.trade{j})-psell+mp.pupgrade) ...
                                      -mp.psych_transcost{tau};
      % Add the utility to avoid a extra loop
      u(: , s.id.trade{j}, s.id.upgrade{j})       = u(:, s.id.trade{j}, s.id.upgrade{j})+B;

      % utilty of keeping and upgrading 
      u(: , s.id.keep, s.id.upgrade{j})           =  trmodel.u_house(mp, s.id.age(s.id.upgrade{j}), tau, j) ...
                                      + ev_scrap-mp.mum{tau}*mp.pupgrade* flipud(s.id.age(s.id.upgrade{j})')'-mp.psych_transcost{tau};
    end % end loop over house-types
 
    % utility of purging (i.e. selling/scrapping the current house and not replacing it, i.e. switching to the outside good)
    u(:, s.id.purge, s.id.purge)                  =  mp.u_og{tau} + mp.mum{tau}*psell + ev_scrap;   

    % additional psych transactions cost and monetary search cost incurred
    % by a consumer who has no house
    u(s.is.nohouse, [s.id.trade{:}], s.is.nohouse)   =  u(s.is.nohouse , [s.id.trade{:}], s.is.nohouse) ... 
                                             -  mp.psych_transcost_nohouse{tau} - mp.mum{tau}*mp.nohousesc; 
    
    u([s.is.scrap{:}] , s.id.keep, s.id.keep)  =  nan; % not possible to keep scrap

    u(s.is.nohouse , s.id.keep, s.id.keep)         =  nan; % not possible to keep no house

  end % end of utility

  function [ccp_scrap, ev_scrap]=p_scrap(mp, s, price_j, tau)
    %   p_scrap:  Function to compute expected utility of the option to scrap rather than sell  
    %             for each consumer in each house state, the probability of scrapping, and the
    %             derivatives with respect to mum{t} and the sell-side transaction cost parameters
    %
    %   SYNTAX:   [ccp_scrap, ev_scrap]=trmodel.p_scrap(mp, s, psell, dpsell, tau)

    psell=trmodel.psell(mp, s, price_j); 

    pscrap=nan(s.ns,1);
    for j=1:mp.nhousetypes
      pscrap(s.is.house_ex_scrap{j})=mp.pscrap{j};
    end

    ev_scrap = (mp.es==1)*logsum([mp.mum{tau}*(pscrap-psell), zeros(s.ns,1)], mp.sigma_s); 
    ccp_scrap=1-exp(-(ev_scrap)/mp.sigma_s);

    for j=1:mp.nhousetypes
      ccp_scrap(s.is.scrap{j})=1;
      ev_scrap(s.is.scrap{j})=0;
    end
    
  end % end of p_scrap
  
  function [uv]=u_house(mp, house_age, tau, house)
    % INPUTS: 
    %     house_age: vector 
    %     tau: household type, scalar index
    %     house: house type, scalar index 

    pkm = mp.p_fuel./mp.fe{house}*1000; % p_fuel is measured in 1000 DKK/l, mp.fe{j} is measured as km/l, but pkm was in DKK/km in regression)
    
    switch mp.modeltype
      case 'structuralform' % structural form - requires that mp.sp are set by trmodel.update_structural_par
        % evaluate utility 
        uv = mp.sp.alpha_0(tau,house) + mp.sp.alpha_a(tau,house) .* house_age + mp.sp.alpha_a_sq(tau,house) * house_age.^2 ... 
        - 1./(2*mp.sp.phi(tau)) .* (max(0,mp.sp.gamma_0(tau,house) + mp.sp.gamma_a(tau) .* house_age - pkm .* mp.mum{tau})).^2;
      case 'reducedform'
        if (mp.convexutility)
        uv=mp.u_0{tau, house}+mp.u_a{tau, house}*house_age + exp(mp.u_a_sq{tau, house})*house_age.^2;
        else
        uv=mp.u_0{tau, house}+mp.u_a{tau, house}*house_age + mp.u_a_sq{tau, house}*house_age.^2;
        end

      otherwise 
      error('Unexpected reduced form type, "%s".', mp.modeltype); 
    end

    % add the (dis)utility from car inspections in even years 
    % (after the first inspection at age 4)
    %inspection=(1-mod(car_age,2)).*car_age>=4; % dummy for inspection year
    %uv= uv+mp.u_even{tau,car}.*inspection; 
  end % end of u

  function [pbuy] = pbuy(mp, s, price_j)     
    pbuy=nan(1, s.nd); 
    for j=1:mp.nhousetypes
      pbuy(s.id.trade_new{j})=  mp.transcost + [mp.pnew{j}*(1+mp.ptranscost)];  
      pbuy(s.id.trade_used{j})= mp.transcost + [price_j{j}*(1+mp.ptranscost)];  
      pbuy(s.id.purge)=0;
    end
  end

  function [psell] = psell(mp, s, price_j)     
    % psell is the selling price net of seller-side transactions costs
    psell=nan(s.ns,1);
    for j=1:mp.nhousetypes
      % car_age=s.is.age(s.is.car{j})'; 
      % inspection=(1-mod(car_age,2)).*car_age>=4;
      % tc=[price_j{j}; mp.pscrap{j}] *mp.ptc_sale+mp.tc_sale+mp.tc_sale_age*car_age+mp.tc_sale_even*inspection;
      % psell(s.is.car{j})=[price_j{j}; mp.pscrap{j}]-tc;

      house_age=s.is.age(s.is.house_ex_scrap{j})';
      %inspection=(1-mod(car_age,2)).*car_age>=4; % dummy for inspection year
      %psell(s.is.car_ex_clunker{j})=[price_j{j}*(1-mp.ptc_sale)-mp.tc_sale-mp.tc_sale_age*car_age-mp.tc_sale_even*inspection];
      psell(s.is.house_ex_scrap{j})=[price_j{j}*(1-mp.ptc_sale)-mp.tc_sale-mp.tc_sale_age*house_age];
      psell(s.is.scrap{j})=[mp.pscrap{j}];
      psell(s.is.nohouse)=0;

    end
  end % end of trade_cost

  function [x]=heating(mp, a, tau, j)
    % Reduced form or structural form heating equation 

    % x = mp.db.pkm{tp}*pkm{car} + mp.db.car{car}+ mp.db.tau{tp} +mp.db.a1{tp}*a+mp.db.a2{tp}*a.^2;
    % trmodel.heating: heating in km/day 
    % Syntax: [x]=trmodel.heating(mp, a, tau, j)
    % OUTPUTS: 
    %   x: optimal driving 1000 km / year 
    %
    % NOTE: if mp.sp parameters are updated appropriately, the two
    % modeltypes should give identical driving (by construction). 

    pkm = mp.p_fuel./mp.fe{j}*1000; % p_fuel is measured in 1000 DKK/l, mp.fe{j} is measured as km/l, but pkm was in DKK/km in regression)
    switch mp.modeltype          
      case 'structuralform' % structural driving equation   
        x = -1./mp.sp.phi(tau).*max(0,mp.sp.gamma_0(tau,j) + mp.sp.gamma_a(tau) .* a - pkm .* mp.mum{tau}); 
      case 'reducedform' % use the estimated OLS equation directly 
        x =  mp.db.house{j} + mp.db.tau{tau} + mp.db.a1{tau}*a + mp.db.pkm{tau}*pkm; 
      otherwise 
        error('Unexpected value in mp.modeltype, "%s".', mp.modeltype); 
    end
  end

  function [F] = age_transition(mp, s)
    % age_transition.m: age transition probability matrices 
    % 
    %  SYNTAX:   [F]= trmodel.age_transition(mp, s)
    %
    % INPUTS: 
    %   mp:     structure with model parameters
    %   s:      structure with state and decision indexes (generated by trmodel.states) 
    %
    % OUTPUS: 
    %   F:      structure with age transition matrices 
    %
    %   The fields of F are block-diagonal matrices with a block for each house and 1 for no house  
    %     F.notrans:   s.ns x s.ns extended physical transition matrix 
    %     F.trans:     s.ns x s.ns state transition matrix conditional on trading a house or upgrading 

    F.notrans =s.Q;
    F.trans   =s.F;

  end

  function [delta, deltaK, deltaT, delta_scrap, deltaTU, deltaU] = trade_transition(mp, s, ccp, ccp_scrap)
    % Keeping transition probabilities
    deltaK=sparse(1:s.ns,1:s.ns,ccp(:,s.id.keep,s.id.keep),s.ns,s.ns);

    % Transition probabilities
    % trade, no upgrade
    deltaT     = ccp(s.tr.state,s.tr.choice, s.id.keep); 
    % trade and upgrade
    deltaTU1    = ccp(s.tr.state,s.tr.choice, s.tr.choice_up); 
    %no trade, upgrade
    deltaU     = ccp(s.tr.state,s.id.keep,s.tr.choice_up); 


    %deltaTU=permute(sum(deltaTU1,3),[1,2,3]);
    deltaTU=permute(sum(deltaTU1,2),[1,3,2]);
    
    deltaU = permute(deltaU,[1,3,2]);
    %disp(size(deltaTU))
    %disp(size(deltaT))
    %disp(size(deltaU))
    %disp(size(deltaK))
    %deltaTU = sum(deltaTU,1);
    %deltaU = permute(deltaU,[1,3,2]);
    
    %disp(size(permute(deltaU,[2,1,3])))
    %  transition probability matrix
    delta=deltaT + deltaU +deltaTU+ deltaK; 

    % sum of trades conditional on upgrading, used to add to calculate
    % total demand
    %deltaTU=permute(sum(deltaTU1,3),[1,2,3]);
    % raaranging upgrade probalilities, used for calculating keeping
    %deltaU = permute(deltaU,[1,3,2]);
  
    %  transition probability matrix
    %delta=deltaT + deltaU +deltaTU+ deltaK; 

    if nargout>3
      delta_scrap =(1-ccp(:, s.id.keep, s.id.keep)).*ccp_scrap;
    end

  end % end of trade_transition
 

  function [ev1, ccp, dev] = bellman(mp, s, util, F, ev)
    
    % bellman.m:    Implements the Bellman operator ev1=Gamma(ev) for the
    %               consumer's problem of trading houses 
    %               to maximize discounted utility in the presence of a
    %               secondary market for houses
    %               This function returns the *expected value function* ev and is a fast vectorized
    %               program that also calculates choice probabilities implied by current expected value ev
    %               and the derivative of the Bellman operator
    %
    % INPUTS: 
    %   mp:     structure with model parameters
    %   s:      structure with state and decision indexes (generated by trmodel.states) 
    %   util:   utility matrix (n.s x n.s+1) (states in rows, decision in columns)
    %   Q:      s.ns x s.ns extended physical transition matrix (block-diagonal matrix with a block for each house and 1 for no house)
    %   F:      s.ns x s.ns state transition matrix conditional on trading or upgrading a house (block-diagonal matrix with a block for each house and 1 for no house)
    %   ev:     s.ns x 1 vector of expected values 
    %
    % OUTPUS: 
    %   ev1:    s.ns  x 1 vector of updated expected values after evaluating the Bellman operator 
    %   ccp:    s.ns  x s.ns + 1 x s.ns +1  matrix of conditional choice probabilities  
    %   dev:    s.ns  x s.ns matrix of derivatives of Bellman operator
    % 
    %  Fixed point of Bellman can be found by calling dpsolver
    % 
    %  SYNTAX: 
    %     [ev,ccp,dev]= dpsolver.poly(@(ev) trmodel.bellman(mp, s, util, F, ev), ev0, mp, mp.bet);


    v=nan(s.ns,s.nd,s.nd);        % choice specific value functions (states by decision)
   
    % calculate the values for trading and each upgrade choice 
    for j=1:s.ns
        v(:, [s.id.trade{:}], j)         =util(:,[s.id.trade{:}], j) + (mp.bet*F.trans([s.is.house{:}],:)*ev)';
    end
    %  calculate the values for keeping
    v([s.is.house{:}],s.id.keep, s.id.keep)=    util([s.is.house{:}],s.id.keep, s.id.keep) + mp.bet*F.notrans([s.is.house{:}],:) *ev;

    % calculate the values for all options involving buying a house but not
    % upgrading 
    v(:, [s.id.trade{:}], s.id.keep)         =util(:,[s.id.trade{:}], s.id.keep) + (mp.bet* F.trans([s.is.house{:}],:)* ev)';

    % calculate the values for all options involving upgrading a house but not
    % trading  
    v(:, s.id.keep,[s.id.upgrade{:}])        =util(:, s.id.keep, [s.id.upgrade{:}]) + reshape((mp.bet* F.trans([s.is.house{:}],:)* ev), [1,1,size((mp.bet* F.trans([s.is.house{:}],:)* ev),1)]);
    
    % calculate the values for the purge decision (ev of no house)
    v(:,s.id.purge, s.id.purge)               =util(:,s.id.purge)   + mp.bet * ev(s.is.nohouse);                     
    
    v(isnan(util))=nan;
    % calculate expected values
    ev1= trmodel.logsum(v, mp.sigma); 
    
    if nargout>1
      % calculate choice probabilities across all choices
      ccp=exp((v-ev1)/mp.sigma);
      ccp(isnan(ccp))=0;  % restore the nans in the elements of ccp_cell corresponding to infeasible choices
    end
    
    if nargout>2
      %  compute trade transition probability matrix
      delta=trmodel.trade_transition(mp, s, ccp);

      % derivative of Bellman operator wrt ev      
      dev=mp.bet*delta*F.trans;
    
    end
  end % end of bellman


  function mp = update_mp(mp)

    %if (~iscell(mp.u_even))
    %     mp.u_even={mp.u_even};
    %end


    % update model dependent parameters when changing mp.ntypes or mp.nhousetypes
    % if cells are smaller than mp.ntypes or mp.nhousetypes the last element is repeated

    mp.db.house               =repcell(mp.db.house,1,mp.nhousetypes); 
    mp.db.pkm                 =repcell(mp.db.pkm,mp.ntypes,1); 
    mp.db.tau                 =repcell(mp.db.tau,mp.ntypes,1); 
    mp.db.a1                  =repcell(mp.db.a1,mp.ntypes,1); 
    mp.db.a2                  =repcell(mp.db.a2,mp.ntypes,1); 

    % house-specific parameters
    param_j       = {'lbl_housetypes', 'abar_j0','fe', 'pnew_notax', 'pscrap_notax'};
    for k=1:numel(param_j)
      mp.(param_j{k}) =repcell(mp.(param_j{k}),1, mp.nhousetypes);
    end  

    % household type-specific parameters
    param_tau     = {'lbl_types','mum','phi','u_og','psych_transcost','psych_transcost_nohouse'};
    for k=1:numel(param_tau)
      mp.(param_tau{k}) =repcell(mp.(param_tau{k}),mp.ntypes,1);
    end  

    % house-and-household specific parameters
    %param_tau_j   = {'u_0','u_a','u_a_sq','u_even'};
    param_tau_j   = {'u_0','u_a','u_a_sq'};

    for k=1:numel(param_tau_j)
      mp.(param_tau_j{k}) =repcell(mp.(param_tau_j{k}),mp.ntypes, mp.nhousetypes);
    end  

    % the fundamental prices are p_fuel_notax, pnew_notax, plus the tax rates 
    [mp.p_fuel, mp.pnew, mp.pscrap] = trmodel.price_aftertax(mp);
  
    if (mp.ntypes == 1)
       mp.tw=1;
    end

    if (abs(sum(mp.tw)-1.0) > 1e-12)
       fprintf('Error in trmodel.setup: mp.tw vector does not sum to 1\n');
    end
    
  end % end of setup

  %function sp = update_structural_par(mp)
  %    % SYNTAX: sp = trmodel.update_structural_par(mp)
  %    % NOTE: this should be called after estimation, but not whenever
  %    % changing e.g. fuel taxes (then the structural parameters are
  %    % fixed). 
  %    
  %    sp = struct(); 
  %    
  %    % 1. phi coefficients (on squared driving)
  %    mum = cell2mat(mp.mum);
  %    sp.phi = mum ./ cell2mat(mp.db.pkm); % ntypes*1 vector
  %    
  %    % 2. compute gamma coefficients
  %    d0 = cell2mat(mp.db.car) + cell2mat(mp.db.tau); % adding a row and column vector -> a ntypes*ncartypes matrix
  %    d1 = cell2mat(mp.db.a1);
  %    sp.gamma_0 = - d0 .* sp.phi; % ntypes*ncartypes
  %    sp.gamma_a = - d1 .* sp.phi; % ntypes*ncar
  %    
  %    % 3. compute alpha coefficients
  %    pkm=mp.p_fuel./cell2mat(mp.fe)*1000; % p_fuel is measured in 1000 DKK/l, mp.fe{j} is measured as km/l, but pkm was in DKK/km in regression)
  %    if (mp.convexutility)
  %    sp.alpha_a_sq = exp(cell2mat(mp.u_a_sq)) + 1./(2*sp.phi)      .* (sp.gamma_a.^2);
  %    else
  %    sp.alpha_a_sq = cell2mat(mp.u_a_sq) + 1./(2*sp.phi)      .* (sp.gamma_a.^2);
  %    end
  %    sp.alpha_a    = cell2mat(mp.u_a)    + sp.gamma_a./sp.phi .* (sp.gamma_0 - mum.*pkm);
  %    sp.alpha_0    = cell2mat(mp.u_0)    + 1./(2*sp.phi)      .* (sp.gamma_0 - mum.*pkm).^2;
      
      % 4. evaluate utility
      %uv = sp.alpha_0 + sp.alpha_a * car_age + sp.alpha_a_sq * car_age.^2 ...
      %    - 1/(2*sp.phi) * (sp.gamma_0 + sp.gamma_a * car_age - pkm * mum);

      % 5. check to make sure no negative driving:  issue warnings if this is found

  %    for t=1:mp.ntypes
  %     for c=1:mp.ncartypes
  %        if (sp.gamma_0+sp.gamma_a*(mp.abar_j0{c}-1) < mum(t)*pkm)
  %          fprintf('Warning: update_structural_par  negative driving predicted for household type %i and cartype %i\n',t,c);
  %        end
  %     end
  %    end
  %end

  function logsum=logsum(v, sigma)

	% logsum.m: Returns a mx1 vector of log-sum values of the m x n x n matrix of values v
	%
	% INPUTS: 
	% 	v      	m x n x n matrix of values v
	% 	sigma   scalar (scale fctor on iid ev error) 
	% 
	% OUTPUS: 
	% 	logsum: m x 1 lector of logsum values 
	%
	% SYNTAX:  
	% 	P=logsum(v, sigma);     
	% 	P=logsum(v):
	%
	% See also:
	% 	logit            

	if nargin ==1
		sigma=1;
	end
    % Max across all choises
	maxv=nanmax(v,[], [2 3]);
	if sigma==0
		logsum=maxv;
		return
	end
	v=v- maxv;
    % sum across all choises
	logsum=maxv + sigma*log(nansum(exp(v/sigma),'all'));

end % end logsum
  

end %methods

end %class
