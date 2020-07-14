% These functions are based on open source code available at: 
%       github.com/saucermanlab/netflux
% for simulating logic-based signaling networks as described originally in: 
%       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). 'Modeling 
%       cardiac B-adrenergic signaling with normalized-Hill differential
%       equations: comparison with a biochemical model.' BMC Systems 
%       Biology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE_master.m
% Solves the ODE system described in our accompanying publication: 
% Irons & Humphrey (2020): Cell signaling model for arterial mechanobiology, 
% PLOS Computational Biology.
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,y]= ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules,ODEfilename)

tspan=linspace(0,25,100);
 
[t,y]=ode45(@(t,y) ode_default(t,y,speciesNames, p0,tau,ymax,reactionRules,ODEfilename),tspan,y0);
y=real(y); 

end

function dydt=ode_default(t,y,speciesNames,p0,tau,ymax, reactionRules,ODEfilename, file_opt)
%Uniform reaction parameters
rpar=zeros(3,length(reactionRules));
rpar(1,:)=p0(1);
rpar(2,:)=p0(2);
rpar(3,:)=p0(3);

mapVar = cell(length(speciesNames),1);
for i = 1:length(speciesNames)
    mapVar{i} = sprintf('%s = %d;',speciesNames{i},i);
end
mapVar=horzcat(mapVar{1:end});
eval(mapVar); %Assign Stress=1 etc

load(ODEfilename); %Previously stored list of ODEs for our specific system
ODElist=horzcat(ODElist{1:end});
eval(ODElist);

end

%% Utility functions
% Implementation of logic gates, directly from open source code available at: 
%       github.com/saucermanlab/netflux
%   and as described originally in: 
%       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). 'Modeling 
%       cardiac B-adrenergic signaling with normalized-Hill differential
%       equations: comparison with a biochemical model.' BMC Systems 
%       Biology.
%
function fact = act(x,rpar)
% hill activation function with parameters w (weight), n (Hill coeff), K';
    w = rpar(1);
    n = rpar(2);
    EC50 = rpar(3); 
    beta = (EC50.^n - 1)./(2*EC50.^n - 1);
    K = (beta - 1).^(1./n);
   fact = min(w,w.*(beta.*x.^n)./(K.^n + x.^n));
    if fact>w                 % cap fact(x)<= 1';
        fact = w;
    end
end

function finhib = inhib(x,rpar)
% inverse hill function with parameters w (weight), n (Hill coeff), K (K0.5)';
    finhib = rpar(1) - act(x,rpar);
end

function z = OR(x,y)
% OR logic gate
    z = x + y - x*y;
end
    
function z = AND(rpar,varargin)
% AND logic gate, multiplying all of the reactants together
    w = rpar(1);
    if w == 0
        z = 0;
    else
        v = cell2mat(varargin);
        z = prod(v)/w^(nargin-2); % need to divide by w^(#reactants-1) to eliminate the extra w's
    end
end