% This script is partly based on open source code available at: 
%       github.com/saucermanlab/netflux
% for simulating logic-based signaling networks as described originally in: 
%       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). 'Modeling 
%       cardiac B-adrenergic signaling with normalized-Hill differential
%       equations: comparison with a biochemical model.' BMC Systems 
%       Biology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DoseResponseSurfaces_StressAngII.m
% Simulate dose response surfaces to simultaneous Stress and exogenous AngII 
% perturbations, as described in our accompanying publication: 
% Irons & Humphrey (2020): Cell signaling model for arterial mechanobiology,
% PLOS Computational Biology. (Reproduces Fig 5A)
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clear global

%Define sampling points. Note, finer scale used in publication.
Stress_vec=[0:0.2:1];%0:0.1:1;
Ang_vec=[0:0.1:1];%0:0.05:1

Store2=[];
%% Set simulation options
ActiveInputs=[4,5]; %Integrins, SACs
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

mapVar = cell(length(speciesNames),1);
for i = 1:length(speciesNames)
    mapVar{i} = sprintf('%s = %d;',speciesNames{i},i);
end
mapVar=horzcat(mapVar{1:end});
eval(mapVar); %Assign Stress=1 etc

samples={'TGFb1','TSP1','TIMP','MMP1','MMP2','MMP9'};

%% Store values for each Stress and AngII combination
for i=1:length(Stress_vec)
    Stress_val=Stress_vec(i);
    for j=1:length(Ang_vec)
        Ang_val=Ang_vec(j);
        %Initial conditions
        y0=zeros(1,length(speciesNames));
        y0(ActiveInputs)=ival;  
        y0(2)=Wss_IC;
        y0(1)=Stress_val;
        y0(3)=Ang_val;
        %Generate data
        [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);
        Store=y(end,[TGFb1,TSP1,TIMP,MMP1,MMP2,MMP9]);     
        Store2=[Store2;Store];
    end
end

%% Plot dose response surfaces
for idx=1:length(samples)    
    surf_mat=[];
    count=0;
    for i=1:length(Stress_vec)
        Stress_val=Stress_vec(i);
        for j=1:length(Ang_vec)
            Ang_val=Ang_vec(j);
            count=count+1;        
            surf_mat(i,j)=Store2(count,idx);
        end
    end
    figure(23); subplot(2,3, idx);
    surf(Stress_vec,Ang_vec,surf_mat');
%     view([0,90])
    xlabel('Stress'); ylabel('AngII');
    zlim([0,1]); colorbar;caxis([0,1])
    title(samples{idx});
%     shading interp
end
set(gcf,'Pos',[308  326  974  258]);