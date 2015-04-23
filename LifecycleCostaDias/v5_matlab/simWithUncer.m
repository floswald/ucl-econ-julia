function [ y, c, a, v,e,logy1 ] = simWithUncer(policyA1,EV,startingA)

% This function takes the policy functions and value functions, along with
% starting assets and returns simulated paths of income, consumption,
% assets and value

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to
global mu sigma rho T r Tretire
global Agrid  Ygrid numSims interpMethod uncertaintyMethod hcIncPDF;
global normBnd

%% ------------------------------------------------------------------------
% Initialise arrays that will hold the paths of income consumption, value
% and assets

% Arguments for output
y = NaN(T, numSims);            % income
c = NaN(T, numSims);            % consumption
v = NaN(T, numSims);            % value
a = NaN(T + 1,numSims);         % this is the path at the start of each period, so we include the 'start' of death

% Other arrays that will be used below
e = NaN(T, numSims);            % the innovations to log income
logy1 = NaN(1, numSims);        % draws for initial income
ly = NaN(T, numSims);           % log income
ypathIndex = NaN(T, numSims);   % holds the index (location) in the vector 

%% ------------------------------------------------------------------------
% Obtain paths using the initial condition and the policy and value
% functions
%-------------------------------------------------------------------------%#
% Scenario where we have hardcoded the discrete
% number of income draws and their pdf

if (uncertaintyMethod == 0)     % we have hardcoded the income process        

    %First get random discrete draws using subroutine getDiscreteDraws.
    %More details given in that routine. Given a seed (which
    %we set equal to the time period), it draws randomly from vector Ygrid(t, :)
    %with pdf hcIncPDF
    for t = 1:1:T 
        seed = t;
        [ y(t, :), ypathIndex(t, :)] = getDiscreteDraws( Ygrid(t, :), hcIncPDF , 1, numSims, seed );
    end

    for s = 1:1: numSims              % loop through individuals
        a(1, s) = startingA;   
        for t = 1:1:T                 % loop through time periods for a particular individual
            
            if (t >= Tretire)
                ixY = 1;
            else 
                ixY = ypathIndex(t, s);
            end

            tA1 = policyA1(t, :, ixY);   %the relevant part of the policy function
            tV = EV(t, :, ixY);          %the relevant part of the value function
            a(t+1, s)   = interp1(Agrid(t, :), tA1 ,a(t, s),interpMethod, 'extrap');
            v(t, s)   = interp1(Agrid(t, :),tV,a(t, s),interpMethod, 'extrap');                                         
            c(t, s) = a(t, s)  + y(t, s) - (a(t+1, s)/(1+r));
        end   %t      
    end % s

%----------------------------------------%
% Scenario where income draws are normally distributed 

elseif (uncertaintyMethod == 1) 
    seed1 = 1223424;
    seed2 = 234636;
    sig_inc = sigma/ ((1-rho^2)^0.5);
    [ e ] = getNormalDraws( 0, sigma,  T, numSims, seed1);  % normally distributed random draws for the innovation
    [ logy1 ] =  getNormalDraws( mu, sig_inc,  1, numSims, seed2); % a random draw for the initial income       
    for s = 1:1: numSims                           % loop through individuals
        a(1, s) = startingA;                   
        ly(1, s) = truncate(logy1(1, s), -normBnd*sig_inc,normBnd*sig_inc );
        y(1, s) = exp(ly(1, s));
        for t = 1:1:T                              % loop through time periods for a particular individual               
            if (t >= Tretire)                      % set income to zero if the individual has retired
                y(t, s) = 0;
            end
            if (t < Tretire)                       % first for the before retirement periods
                clear tA1 tV;                      %necessary as the dimensions of these change as we wor through this file
                tA1(:,:) = policyA1(t, :, :);  % the relevant part of the policy function
                tV(:,:) = EV(t, :, :);         % the relevant part of the value function                
                a(t+1, s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tA1, a(t, s), y(t, s));                    
                v(t  , s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tV , a(t, s), y(t, s));
                
                if (t ~= T)  % Get next year's income
                    ly(t+1, s) = (1 -rho) * mu + rho * ly(t, s) + e(t + 1, s);
                    ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc );
                    y(t+1, s) = exp( ly(t+1, s) );                
                end % if (t ~= T)
            else                          % next for the post retirement periods
                clear tA1 tV;
                tA1 = policyA1(t, :, 1);  % the relevant part of the policy function
                tV = EV(t, :, 1);         % the relevant part of the value function                
                a(t+1, s) = interp1(Agrid(t,:)', tA1, a(t, s), 'linear', 'extrap');
                v(t,   s) = interp1(Agrid(t,:)', tV , a(t, s), 'linear', 'extrap');
                if (t ~= T)
                    y(t+1, s) = 0;                  
                end % if (t ~= T)
            end % if (t < Tretire)
            
            % Check whether next period's asset is below the lowest
            % permissable
            if ( a(t+1, s) < Agrid(t+1, 1))
               [ a(t+1, s) ] = checkSimExtrap( Agrid(t+1, 1),y(t, s), t ); 
            end
                
            c(t, s) = a(t, s)  + y(t, s) - (a(t+1, s)/(1+r));
        end   %t      
    end % s
end % elseif (uncertaintyMethod == 1)
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

  
 
 
end