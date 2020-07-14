% This script is partly based on open source code available at: 
%       github.com/saucermanlab/netflux
% for simulating logic-based signaling networks as described originally in: 
%       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). 'Modeling 
%       cardiac B-adrenergic signaling with normalized-Hill differential
%       equations: comparison with a biochemical model.' BMC Systems 
%       Biology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knockdown_Script.m
% Perturbs each node in turn by reducing Ymax from 1 to 0.1, whilst measuring
% changes in steady state outputs of each of species, as described in our 
% accompanying publication: Irons & Humphrey (2020): Cell signaling model 
% for arterial mechanobiology, PLOS Computational Biology. (Fig 3)
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clear global  

%% Set simulation options
numInputs=5; %#Model inputs
ival=0.2; 
Wss_IC=0.5;

load('reactions_final.mat');
ODEfilename='ODElist_final';

k_val=0.1; %Ymax for each species

%Uniform parameters
w=1; 
n=1.25;
EC50=0.55;
p0=[w;n;EC50];

if EC50^n>1/2
    warning('EC50^n>1/2: Negative B in normalised Hill function');
end         
       
%Initial conditions
y0=zeros(1,length(speciesNames));
y0(1:numInputs)=ival;  
y0(2)=Wss_IC; 

%% Generate baseline data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);

[Ordering,OrderingStr]=DefineOrdering(); %Reorder species for output matrix 

%% Perturb nodes and store absolute changes in steady states
StoreAbs=zeros(length(speciesNames));  
count=0;
for k_idx=Ordering
    count=count+1;
    y0(1:numInputs)=ival; 
    y0(2)=Wss_IC;

    ymax=ones(1,length(speciesNames));
    ymax(k_idx)=k_val;
    
    if k_idx<=numInputs
        y0(k_idx)=k_val;
    end

    [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);

    AbsDiff=transpose(y(end,:)-y_ref(end,:));
    AbsDiff=AbsDiff(Ordering,:);

    StoreAbs(count,:)=AbsDiff;
end 

StoreAbs(1:5,:)=[]; StoreAbs(:,[1:5])=[]; OrderingStr([1:5])=[]; 
StoreAbs=transpose(StoreAbs);
StoreAbs=StoreAbs(end:-1:1,:); %flip for imagesc
            
figure();
h=imagesc(StoreAbs);
yticks([1:1:1+size(StoreAbs,1)]);
yticklabels(OrderingStr(end:-1:1));
xticks(1:1:1+size(StoreAbs,2));
xticklabels(OrderingStr);xtickangle(270);
% colormap custom_cmap
colorbar
set(gca,'TickLength',[0 0])
set(gca,'Fontsize',8)
set(gcf,'Pos',[587   168   657   573])
suptitle('Absolute difference');
                
               