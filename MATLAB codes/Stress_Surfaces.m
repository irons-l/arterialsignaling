% This script is partly based on open source code available at: 
%       github.com/saucermanlab/netflux
% for simulating logic-based signaling networks as described originally in: 
%       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). 'Modeling 
%       cardiac B-adrenergic signaling with normalized-Hill differential
%       equations: comparison with a biochemical model.' BMC Systems 
%       Biology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stress_Surfaces.m
% Simulate steady states of the model equations under differing basal inputs
% and perturbed stress inputs, as described in our accompanying publication: 
% Irons & Humphrey (2020): Cell signaling model for arterial mechanobiology,
% PLOS Computational Biology. (Reproduces figures in S3 Appendix)
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clear global  

%% Set simulation options
ActiveInputs=[1,4,5]; %Stress, AngII, Integrins, SACs
ival=0.2; 
Wss_IC=0.5;

load('reactions_final.mat');
ODEfilename='ODElist_final';

%Uniform parameters
w=1; 
n=1.25;
EC50=0.55;
p0=[w;n;EC50];

if EC50^n>1/2
    warning('EC50^n>1/2: Negative B in normalised Hill function');
end
               
%% Set basal and perturbed levels (in the publication, a finer spacing of 0.05 was used)
bsam=0:0.1:0.5;
psam=0:0.1:0.5;
pertidx=1; %Stress perturbation

outputNames={'TGFb1';'TSP1';'TIMP';'MMP1';'MMP2';'MMP9';'NO';'ET1';'Col1mRNA';'Col3mRNA';'ActomyosinActivity';'SMCproliferation'};

mapVar = cell(length(speciesNames),1);
for i = 1:length(speciesNames)
    mapVar{i} = sprintf('%s = %d;',speciesNames{i},i);
end
mapVar=horzcat(mapVar{1:end});
eval(mapVar); %Assign Stress=1 etc

store=zeros(length(bsam)*length(outputNames),length(psam));

%% Calculate change in steady states for each (b,p) pair
for bidx=1:length(bsam)
    for pidx=1:length(psam)
    %Initial conditions
    y0=zeros(1,length(speciesNames));
    y0(2)=Wss_IC;
    y0(ActiveInputs)=bsam(bidx);
    y0(pertidx)=y0(pertidx)+psam(pidx);

    %Generate data
    [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);
    y=y(end,:); %Steady states

    for i=1:length(outputNames)
        store(bidx+length(bsam)*(i-1),pidx)=y(eval(outputNames{i}));
    end
        
    end
end
  
%% Plot surfaces for each species of interest
for i=1:length(outputNames)
    figure(1); subplot(4,3,i); hold on;
    temp_store=store(length(bsam)*(i-1)+1:length(bsam)*i,1:length(psam));    
    surf(bsam,psam,temp_store'); xlabel('b');ylabel('p'); zlabel(outputNames{i});
    view([60,30]); grid on; box on; zlim([0,1]);

    %normalise temp_store to be fold changes for each ival
    for j=1:size(temp_store,1) %for each i row
        temp_store(j,:)=temp_store(j,:)/temp_store(j,1); %normalise that row by the first value
    end

    figure(2); subplot(4,3,i); hold on; 
    surf(bsam,psam,temp_store','FaceAlpha',0.8); xlabel('b');ylabel('p'); zlabel(outputNames{i});
    view([60,30]); grid on; box on; %zlim([0,1]);
    hold on;
    surf(bsam,psam,ones(size(temp_store')),'FaceAlpha',0,'LineWidth',0.5,'Edgecolor','r');
    xlim([0.15,bsam(end)]); ylim([psam(1),psam(end)]);
end

% save(['b_p_combinations_pertidx',num2str(pertidx),'.mat']); 


  
