% This script is partly based on open source code available at: 
%       github.com/saucermanlab/netflux
% for simulating logic-based signaling networks as described originally in: 
%       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). 'Modeling 
%       cardiac B-adrenergic signaling with normalized-Hill differential
%       equations: comparison with a biochemical model.' BMC Systems 
%       Biology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IO_Script.m
% Simulates input-output relations for our optimal basal parameters and 
% quantifies qualitative matches to experimental observations. The values 
% of b, p, n, and EC50 were looped over to find this chosen/optimal combination,
% as described in our accompanying publication: 
% Irons & Humphrey (2020): Cell signaling model for arterial mechanobiology,
% PLOS Computational Biology. (Reproduces Fig 2B)
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clear global  

%% Set simulation options
ActiveInputs=[1,4,5]; %Stress, (AngII), Integrins, SACs

load('reactions_final.mat');
ODEfilename='ODElist_final';
Wss_IC=0.5;
bval=0.2;
pertval=0.3;
EC50=0.55;
n=1.25;

Diff_store=[];

%Uniform Parameters
w=1; 
n=1.25;
EC50=0.55;
p0=[w;n;EC50];               

%Initial conditions, reference case
y0=zeros(1,length(speciesNames));
y0(ActiveInputs)=bval;  
y0(2)=Wss_IC;

%Generate baseline data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);

for i_idx=1:3 %perturb Stress, Wss and AngII
    
    %Initial conditions, perturbed cases
    y0=zeros(1,length(speciesNames));
    y0(ActiveInputs)=bval;  
    y0(2)=Wss_IC;
    y0(i_idx)=y0(i_idx)+pertval;

    [~,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);

    Diff_store=[Diff_store;y(end,:)-y_ref(end,:)]; 
end

Diff_store2=Diff_store;
Qual_store=sign(Diff_store2);

%Reorder
[Ordering,OrderingStr]=DefineOrdering(); %Reorder species names

Diff_store=Diff_store(:,Ordering);
Diff_store2=Diff_store2(:,Ordering);
Qual_store=Qual_store(:,Ordering);

Exp_validation=CreateValidationMatrix(speciesNames, Ordering);

%%% Compare to experimental validation matrix
count=0; nullcount=0;
for i=1:size(Exp_validation,1)
    for j=1:size(Exp_validation,2)
        if Qual_store(i,j)==Exp_validation(i,j)
           count=count+1;
        else
            if abs(Exp_validation(i,j))~=0.1 %0.1 is default 'unknown' for experimental matrix
                  disp (['Mismatch for input: ', OrderingStr{i} ', output: ', OrderingStr{j}])
                  else %no experimental data for this pair
                      nullcount=nullcount+1;
             end
        end 
    end
end
        
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%']);
disp(['Basal input: ', num2str(bval)]);
disp(['Input perturbation: ', num2str(pertval)]);
disp(['n: ', num2str(n)]);
disp(['EC50: ', num2str(EC50)]);
disp(['Agreement on ', num2str(count), ' of ', num2str(numel(Exp_validation)-nullcount), ' experimental observations: ', num2str(100*count/((numel(Exp_validation)-nullcount))), '%'])
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%']);

%%Plot subset of interest
SpeciesIdx=[12,22,23,24,21,25,30,31,35,39,26,27,49,50];
Diff_store2=Diff_store2(:,SpeciesIdx);
Qual_store=Qual_store(:,SpeciesIdx);
SpeciesSubset={'TGFb1', 'MMP1', 'MMP2', 'MMP9',...
            'TSP1', 'TIMP', 'NO', 'ET1', 'Akt', 'p70S6K', 'Col1mRNA',...
            'Col3mRNA', 'ActomyosinActivity', 'SMCproliferation'};
           
figure();
h=imagesc(Diff_store2(3:-1:1,:));
xticks([1:1:length(SpeciesSubset)+1]);
xticklabels(SpeciesSubset);%(end:-1:1));
xtickangle(40);

yticks(1:1:3);
yticklabels(OrderingStr(3:-1:1)); 
caxis([-0.5,1]);
% colormap custom_cmap
colorbar
set(gca,'TickLength',[0 0])
set(gca,'Fontsize',8)
set(gcf,'Pos',[302   414   868   206])
