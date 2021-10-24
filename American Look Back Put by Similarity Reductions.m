%% Clean the workspace

clear 

%% Setting the parameters and obtain results

r = 0.1;
sig = 0.2;
T = 1;
dividend = 0;

sampling_t_A = 0.5:1:11.5;
sampling_t_B = [3.5 7.5 11.5];
ksi = [0.85 0.95 1.05];
am_lookback = LB_American_Iter(ksi,r, sig, 0.01, 0.01, T, dividend, 0, 0) % Calculation for 0 sampling
am_lookback_A = LB_American_Iter(ksi,r, sig, 0.01, 0.01, T, dividend, 1, sampling_t_A) % Calculation for case A
am_lookback_B = LB_American_Iter(ksi,r, sig, 0.01, 0.01, T, dividend, 1, sampling_t_B) % Calculation for case B


%% Define the functions for American Lookback Put


function put_Price = LB_American_Iter(S0,r,sig, del_ksi,del_t, time, div, is_disc, disc_samp)

% Function to calculate the price of a vanilla European

% Create mesh 

ksi_vector = 0:del_ksi:3;
t_vector = 0:del_t:time;
omega = 1;

% Taking effective interest rate
r = r - div;

% Specify the boundary conditions for call

% Get the number of grid points
M = length(ksi_vector)-1;
N = length(t_vector)-1;

% Pre-allocate the output
price_mesh(1:M+1,1:N+1) = nan;

% Constructing price mesh for call price and setting boundary conditions
price_mesh(:,end) = max(1-ksi_vector,0);
price_mesh(1,:) = exp(-r*t_vector(end:-1:1));
price_mesh(end,:) = 0;

% Calculate the coefficients of tridiagonal matrix
j = 0:M;
sig2 = sig*sig;
pj = (del_t/4)*(sig2*(j.^2) - r*j);
qj = -(del_t/2)*(sig2*(j.^2) + r);
rj = (del_t/4)*(sig2*(j.^2) + r*j);

% The tridiagonal matrix for LHS and RHS are C and D respectively
C = -diag(pj(3:M),-1) + diag(1-qj(2:M)) - diag(rj(2:M-1),1);
D = diag(pj(3:M),-1) + diag(1+qj(2:M)) + diag(rj(2:M-1),1);

% Solve at each node
% bound_S function is used to fetch the boundary condition for the stock
% at +/- infinity for option price

bound_S = zeros(size(D,2),1);

round_decimals = -log10(del_t);
disc_time = round(round(disc_samp/12,round_decimals)*(10^round_decimals),0);

for idx = N:-1:1
    if length(bound_S)==1
        bound_S = pj(2)*(price_mesh(1,idx)+price_mesh(1,idx+1)) + ...
            rj(end)*(price_mesh(end,idx)+price_mesh(end,idx+1));
    else
        bound_S(1) = pj(2)*(price_mesh(1,idx)+price_mesh(1,idx+1));
        bound_S(end) = rj(end)*(price_mesh(end,idx)+price_mesh(end,idx+1));
    end
    
    % Condition to check whether to do discrete sampling or not
    if is_disc
        %Applying jump condition across sampling dates
        if any(disc_time(:) == idx)
           x = round(min(ksi_vector,1)/del_ksi,0)+1;
           for jump = 2:size(x,2)-1
                price_mesh(jump,idx) = max(ksi_vector(jump),1)*price_mesh(x(jump),idx+1);
           end
           continue;
        end
    end
    
    price_mesh(2:M,idx) = proj_sor(C,(D*price_mesh(2:M,idx+1) + bound_S),1000,1.e-6,price_mesh(2:M,:),idx,omega);
end

% Calculate the option price
put_Price = interp1(ksi_vector,price_mesh(:,1),S0);
end

% Reference Taking Gauss-Siedel pseudo-code from University of Waterloo
% Publication - modified for SOR


function x = proj_sor( M, b, N, e ,price,idx,omega)
% Solve Mx = b
% The diagonal entries of M and their inverses
n = length( b );
d = diag( M );

if ~all( d )
    error 'at least one diagonal entry is zero';
end

invd = d.^-1;
% Matrix of off-diagonal entires of N
Moff = M - diag( d );

% Use d.^-1*b as the first approximation to x
invdb = invd.*b;
x = d.*b;

%              -1
% Iterate x = D  (b - M   *x)
%                      off
for k = 1:N
    xprev = x;
    for i = 1:n
        %Projected SOR
        % Here is the critical difference between European and American -
        % max used to satisfy linear complimentarity 
        x(i) = max((x(i) - omega*(x(i) - (invdb(i) - invd(i).*(Moff(i,:)*x)))), price(i,idx+1));
    end
    % Check for tolerance
    if norm( x - xprev, inf ) < e
        return;
    end
end
error 'the method did not converge';
end