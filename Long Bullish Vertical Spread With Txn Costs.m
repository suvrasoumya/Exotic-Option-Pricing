%% Clean the workspace

clear 

%% Initialize parameters
s_0 = 0:1:100;
E_1 = 45;
E_2 = 55;
S_max=150; % maximum share for mesh
S_min=0; % minimum share for mesh
T = 0.5; % time to expiry
sig = 0.4; % Volatility
r=0.1; % risk free rate

% Compute portfolio values
basket_payoff = max(s_0-E_1,0)' - max(s_0-E_2,0)';
bkt_with_txn = vertical_spread(s_0,E_1,E_2,r,sig,T,S_max,S_min,1);
bkt_without_txn = vertical_spread(s_0,E_1,E_2,r,sig,T,S_max,S_min,0);


%% Plot for payoff of basket of options at expiry and with / without transaction costs

figure(1);
plot(s_0,basket_payoff,"r-");
hold on;
plot(s_0,bkt_with_txn',"b-");
plot(s_0,bkt_without_txn',"o--");
legend('payoff','with transaction costs','without transaction costs');
hold off;

%% Computing the delta of the portfolio by taking finite difference while delta_S = 1

delta_bkt = basket_payoff(2:end)'-basket_payoff(1:end-1)';
delta_txn = bkt_with_txn(2:end) - bkt_with_txn(1:end-1);
delta_without_txn = bkt_without_txn(2:end) - bkt_without_txn(1:end-1);

%% Plot for delta of basket of options at expiry and with / without transaction costs

figure(2);
plot(s_0(2:end),delta_bkt,"r-");
hold on;
plot(s_0(2:end),delta_txn,"b-");
plot(s_0(2:end),delta_without_txn,"o--");
legend('payoff','with transaction costs','without transaction costs');
hold off;

%% explicit scheme code referenced from Prof Aitor Bergara at EHU

% suitable changes have been made to include transaction cost functionality

function portfolio_price = vertical_spread(s_0,E_1,E_2,r,sigma,T,S_max,S_min,is_txn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% As we have explicit scheme, we need delta_s^2 > tdelta_t to maintain
% stability
M=10000; % Number of time points
N=150; % Number of share price points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=(T/M); % Time step = 10^-4
ds=(S_max-S_min)/N; % Price step = 1; Hence, the scheme is stable

% Initializing the matrix of the option value
v(1:N,1:M) = 0.0;

% Initial conditions at t=T V(S,T) = Payoff(S)
v(1:N,1) = max((S_min+(0:N-1)*ds -E_1),0) - max((S_min+(0:N-1)*ds -E_2),0);
% Boundary conditions at S=0; V(0,t)=0
v(1,2:M)=zeros(M-1,1)'; 
% V(S,t)=Payoff * exp[-r(T-t)] as S ->infininty.
v(N,2:M)=(max((N-1)*ds+S_min -E_1,0) - max((N-1)*ds+S_min -E_2,0))*exp(-r*(1:M-1)*dt);

% Determining the matrix coeficients of the explicit algorithm
aa=0.5*dt*(sigma*sigma*(1:N-2).*(1:N-2)-r*(1:N-2))';
bb=1-dt*(sigma*sigma*(1:N-2).*(1:N-2)+r)';
cc=0.5*dt*(sigma*sigma*(1:N-2).*(1:N-2)+r*(1:N-2))';
% Calculate the cost matrix
kj = dt*(0.25 * sigma*sigma*(1:N-2).*(1:N-2) * sqrt(2/pi))';
% Implementing the explicit algorithm and subtracts transaction cost from
% the value to the portfolio (if transaction cost exists)
for i=2:M
    if is_txn
        v(2:N-1,i)=bb.*v(2:N-1,i-1)+cc.*v(3:N,i-1)+aa.*v(1:N-2,i-1) - kj.*abs(v(1:N-2,i-1)+v(3:N,i-1)-2*v(2:N-1,i-1));
    else
        v(2:N-1,i)=bb.*v(2:N-1,i-1)+cc.*v(3:N,i-1)+aa.*v(1:N-2,i-1);
    end
end
% Reversal of the time components in the matrix as the solution of the BlackScholes
% equation was performed backwards
v=fliplr(v);
% Get call price by interpolating against option price
portfolio_price = interp1(ds*(1:N),v(1:N,1)',s_0);
end
