% This script is partly based on open source code available at: 
%       github.com/saucermanlab/netflux
% for simulating logic-based signaling networks as described originally in: 
%       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). 'Modeling 
%       cardiac B-adrenergic signaling with normalized-Hill differential
%       equations: comparison with a biochemical model.' BMC Systems 
%       Biology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StressAngII_doses.m
% Simulate response to differing Stress and exogenous AngII combinations,
% as described in our accompanying publication: 
% Irons & Humphrey (2020): Cell signaling model for arterial mechanobiology,
% PLOS Computational Biology. (Reproduces Fig 4)
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clear global

%% Set simulation options
ActiveInputs=[4,5]; %Integrins, SACs - fixed throughout
Wss_IC=0.5;
   
load('reactions_final.mat');
ODEfilename='ODElist_final';

mapVar = cell(length(speciesNames),1);
for i = 1:length(speciesNames)
    mapVar{i} = sprintf('%s = %d;',speciesNames{i},i);
end
mapVar=horzcat(mapVar{1:end});
eval(mapVar); %Assign Stress=1 etc
  
%Uniform Parameters
w=1; 
n=1.25;
EC50=0.55;
p0=[w;n;EC50];

if EC50^n>1/2
    warning('EC50^n>1/2: Negative B in normalised Hill function');
end

samples=[0.1,0.1;
         0.2,0.2;
         0.4,0.2];
     
for sidx=1:size(samples,1)
    
    bval=samples(sidx,1);
    p=samples(sidx,2);
    Store2=[];

    %%Baseline
    pertval=p;
    pertidx=[];

    %Initial conditions
    y0=zeros(1,length(speciesNames));
    y0(ActiveInputs)=0.2;  
    y0(2)=Wss_IC;
    y0(1)=bval;
    y0(3)=bval; 
    %Generate baseline data
    [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);
    Store=y(end,[TGFb1,TSP1,TIMP,MMP1,MMP2,MMP9]);
    Store2=[Store2;Store];

    %%Low stress, high angII
    y0(1)=bval;
    y0(3)=bval+2*pertval; 
    %Generate data
    [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);
    Store=y(end,[TGFb1,TSP1,TIMP,MMP1,MMP2,MMP9]);
    Store2=[Store2;Store];

    %%Mid stress, low angII
    y0(1)=bval+pertval;
    y0(3)=bval; 
    %Generate data
    [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);
    Store=y(end,[TGFb1,TSP1,TIMP,MMP1,MMP2,MMP9]);
    Store2=[Store2;Store];

    %%Mid stress, high angII
    y0(1)=bval+pertval;
    y0(3)=bval+2*pertval; 
    %Generate data
    [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);
    Store=y(end,[TGFb1,TSP1,TIMP,MMP1,MMP2,MMP9]);
    Store2=[Store2;Store];

    %%High stress, low angII
    y0(1)=bval+2*pertval;
    y0(3)=bval; 
    %Generate data
    [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);
    Store=y(end,[TGFb1,TSP1,TIMP,MMP1,MMP2,MMP9]);
    Store2=[Store2;Store];

    %%High stress, high angII
    y0(1)=bval+2*pertval;
    y0(3)=bval+2*pertval; 
    %Generate data
    [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);
    Store=y(end,[TGFb1,TSP1,TIMP,MMP1,MMP2,MMP9]);
    Store2=[Store2;Store];

    sn={'TGFb1','TSP1','TIMP','MMP1','MMP2','MMP9'};
    figure();
    for i=1:6 
        subplot(6,1,i);
        bdata=Store2(:,i)/Store2(1,i);  
        c=categorical({'low \sigma','mid \sigma', 'high \sigma'});
        c = reordercats(c,{'low \sigma','mid \sigma', 'high \sigma'});
        b=bar(c,reshape(bdata,[2,3])');
        b(1).FaceColor=[0.9,0.9,0.9];
        b(2).FaceColor=[0.6,0.6,0.6];
        title(sn(i));
    end
    % l=legend( {'\fontsize{12} No AngII, ','\fontsize{12} AngII'});
    % l.Position=[  0.4729    0.9251    0.2582    0.0576];
    set(gcf,'Pos', [761    49   331   948]);
    suptitle(['\fontsize{13} b=',num2str(bval),', p=',num2str(p),', '])

end
