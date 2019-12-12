clc;tic;
load interaction.mat;
[nd,nc] = size(interaction); % nd:diseases number, nc:microbes number

disSim01  = GSD( interaction' ); % disease gaussian
circSim01  = GSM( interaction'); % microbe gaussian

disSim02  = cosSim( interaction); % disease cosine
circSim02  = cosSim( interaction');   % microbe cosine

dis_sim  = combineSim(disSim02,disSim01); % combine disease cosine and gaussian
circRNA_sim  = combineSim(circSim02,circSim01); % combine microbe cosine and gaussian

% dis_sim  = combineSim(disSim01,disSim02);
% circRNA_sim  = combineSim(circSim01,circSim02);

F = NCPLDA(circRNA_sim,dis_sim,interaction');

%%
microbe_list = importdata('microbe_list.txt','%s');
disease_list = importdata('disease_list.txt');

%%
[F,LP_rank_known] =Rank( F, interaction',microbe_list,disease_list);
Write_file( F );
toc;
