% This script is partly based on open source code available at: 
%       github.com/saucermanlab/netflux
% for simulating logic-based signaling networks as described originally in: 
%       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). 'Modeling 
%       cardiac B-adrenergic signaling with normalized-Hill differential
%       equations: comparison with a biochemical model.' BMC Systems 
%       Biology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BaseScript_master.m
% Simulate timecourses of the model equations described in our
% accompanying publication: Irons & Humphrey (2020): Cell signaling model
% for arterial mechanobiology, PLOS Computational Biology.
% Here, we use uniform reaction parameters and constant basal inputs as 
% described in the publication.
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clear global  

%% Set simulation options
load('reactions_final.mat'); %contains speciesNames, reactionRules, tau, ymax 

%Initial conditions
y0=zeros(1,length(speciesNames));
ActiveInputs=[1,4,5]; %Stress, Integrins, SACs
ival=0.2; 
Wss_IC=0.5;
y0(ActiveInputs)=ival;  
y0(2)=Wss_IC;

ODEfilename='ODElist_final'; %contains list of ODEs for our system

%Uniform parameters, reference case
w=1; 
n=1.25;
EC50=0.55;
p0=[w;n;EC50];

if EC50^n>1/2
    warning('EC50^n>1/2: Negative B in normalised Hill function');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Modify for specific knockdowns
% mapVar = cell(length(speciesNames),1);
% for i = 1:length(speciesNames)
%     mapVar{i} = sprintf('%s = %d;',speciesNames{i},i);
% end
% mapVar=horzcat(mapVar{1:end});
% eval(mapVar); %Assign Stress=1 etc
% ymax(ETBR)=0;
% ymax(NO)=0.5;% ymax(ET1)=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate timecourse data
[t,y]=ODE_master(p0,y0,speciesNames,tau,ymax, reactionRules, ODEfilename);

%Plot example timecourses
AngII=7;
AT1R=23;
figure();
plot(t,y(:,AngII),'LineWidth',1); hold on;
plot(t,y(:,AT1R),'LineWidth',1);
xlabel('Time'); ylabel('Species Activity');
legend('AngII','AT1R');
set(gcf,'Pos',[495  408  560  207])

%All species
figure(); imagesc(transpose(y(:,end:-1:1)));
yticks([1:1:1+size(y,1)]);
yticklabels(speciesNames(end:-1:1));
xlabel('Time');
colorbar
set(gcf,'Pos',[475  102 376  614])
