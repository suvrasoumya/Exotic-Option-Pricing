% Clean the workspace

clear 

%%
% Setting the parameters


sig = 0.3;  % volatility
r = 0.05; % Risk free return (continuous compound)
div = 0.00; % dividend rate
s_0 = 100; % stock price at time 0
E = 95; % option strike price
X = 90; % Barrier price
steps = 30000; %Step in Binomial tree
T = 1; % Years till expiry
deltaT = T/steps; % one step time


% To be noted X < E (Barrier price < Strike )
if X > E
   error('The barrier X must be lesser than the strike of down and out European Call')
end


% taking care to include dividend yield in only the construction phase of tree
r_mod = r-div;


% Use the method where u*d=1 due to its recombinant nature
% Calculation of A
A = 0.5*(exp(-r_mod*deltaT) + exp((r_mod+sig^2)*deltaT));

% Calculation of u,d and p
u =  A + sqrt(A^2 -1);
d = A - sqrt(A^2 -1);
p = ( exp((r_mod)*deltaT)-d)/(u-d); % p for up 1-p for down

barrier_tree = Tree_Calculation(s_0,E,u,d,p,r,steps,deltaT,X)

% Using explicit method to compute price of the barrier
eu_call = blsprice(s_0,E, r, deltaT*steps, sig, div);
bls_bar = blsprice(X*X/s_0,E, r, deltaT*steps, sig, div);

barrier_explicit = eu_call - bls_bar*((s_0/X)^((2*r/(sig^2) - 1)*(-1)))

%% Function to simulate a binomial tree

function option_price = Tree_Calculation(S_0,E,u,d,p,r,steps,deltat,X)

discount = exp(-r*deltat);

% Initialize price tree
priceTree = nan(steps+1,steps+1);
priceTree(1,1) = S_0;
% Creating price tree
for idx = 2:steps+1
    priceTree(1:idx-1,idx) = priceTree(1:idx-1,idx-1)*u;
    priceTree(idx,idx) = priceTree(idx-1,idx-1)*d;
    
end
% Initialize option tree
optionTree = nan(size(priceTree));

% Calculating value for an arbitrary payoff function
optionTree(:,end) = payoff_func(priceTree(:,end),E);

steps = size(priceTree,2)-1;
% Calculating option price at every step
for idx = steps:-1:1
    optionTree(1:idx,idx) = discount*(p*optionTree(1:idx,idx+1) ...
        + (1-p)*optionTree(2:idx+1,idx+1));
    % Barrier option crtitcal condition
     optionTree((priceTree(:,idx) < X),idx) = 0;
    
end


option_price = [optionTree(1)];
end

%Barrier payoff
function payoff = payoff_func(share_price, strike)

payoff = max(share_price-strike,0);

end


