%% Fonction pour le Calcul de la cohérence Cortico-musculaire par FFT
%% Ahmed ADHAM - 2022

%%
% cette fonction va chercher les fichiers avec "previb" ou "postvib" dans leur nom et calculer la cohérence
% Les électrodes d'intérêt pour l'EEG sont les électrodes entre 7 et 70 dans le fichier
% 
%%
clear

load('EEG_Struc_64.mat') % charger noms électrodes + positions
%
%% Paramètres du signal et du calcul 

srate = 1024 %Hz - Fréquence d'échantillonage en Hertz 
Interval_1 = [4000:8000] % Interval de temps où calculer la cohérence dans le signal
repertoire_1 = '/Users/Ahmed/Downloads/wetransfer_s7j1contract2b_interpbad_notch_band_2022-07-21_1633/pre/'; % Répertoire où se trouvent les .mat Pré-vibration
repertoire_2 = '/Users/Ahmed/Downloads/wetransfer_s7j1contract2b_interpbad_notch_band_2022-07-21_1633/post/'; % Répertoire où se trouvent les .mat Post-vibration
EMG_chanel = 6; % Quel canal EMG utiliser : mettre 6 pour FCR, 5 pour ECR 
Freq_space = [1:0.5:50]; % Espace de fréquences à étudier (1 à 50Hz tous les 0.5Hz). Attention au temps de calcul ! 
Utiliser_laplacien = 0; %0 = pas de laplacien ; 1 = laplacien // je conseille de ne pas le mettre mais à toi de voire
condition{1} = 'pre-vibration' % nom condition_1
condition{2} = 'post-vibration'% nom condition_2
n_sujet = 1; % nombre de sujets dans le protocole
n_analyse = 10; % nombre d'analyses
Hanning_window = 1000; % La cohérence par FFT demande une fenêtre d'ordre du filtre : 
% c'est un paramètre un peu arbitraire qui détermine la résultion de
% fréquence ... (trop (ex 2000) sera illisible, pas assez (ex 100) sera tout
% mélangé; du coup j'ai mit 1000 mais libre à chacun de changer
% j'ai mit un petit truc en plus pour simuler une rythme à 10Hz sur C3 et sur le canal EMG ci dessous
% : 
% par défaut laisse sim à 0, sauf si tu veux simuler; 
% Tu peux changer la puissance du signal simulé en jouant sur le 0.000001
% essaie egalement de voire ce que ca donne avec et sans filtre laplacien
% et essaie ensuite de changer le hanning_window ;) 

sim = 0;
power_sim = 0.0001;

%% Paramètres d'affichage des résultats 

Electrode_diff = 32; 
%électrode sur laquelle calculer la différence de cohérence (32 = C3, cf tableau en bas) 

ylim_d = [-0.2 0.2]; % échelle cohérence figure de différence (figure 100)
ylim_c = [0 0.3]; % échelle cohérence figure individuelle (figure 1/2)
c_lim = [0 0.3]; % échelle des couleurs entre 0 et 1 pour les cartes globales (figure 1/2)

f1 = 8;  % Fréquence affichée pour la figure 1
f2 = 12; % Fréquence affichée pour la figure 2
f3 = 20; % Fréquence affichée pour la figure 3
f4 = 30; % Fréquence affichée pour la figure 4




%% c'est parti :) ! 

for sujet = 1:n_sujet

for pre_post = 1:2 
        
repertoire_1s = [repertoire_1 'sujet_' num2str(sujet) '/'] 
repertoire_2s = [repertoire_2 'sujet_' num2str(sujet) '/'] 

Name_file{sujet,1} = dir(fullfile(repertoire_1s, 'data*'))
Name_file{sujet,2} = dir(fullfile(repertoire_2s, 'data*'))


for file = 1:n_analyse 

disp(strcat({'Il reste : ' num2str(file),' fichier / ', num2str(n_analyse) 'sujet'  num2str(pre_post)}))

if pre_post == 1
file_to_load = strcat(repertoire_1s,Name_file{sujet,1}(file).name);
load(file_to_load);
end

if pre_post == 2
file_to_load = strcat(repertoire_2s,Name_file{sujet,2}(file).name);
load(file_to_load);
end

times = 1/1000:1/1000:size(F,2)/1000; % temps

sig = 0; 
if sim == 1
sig = cos(2*10*pi*times)*power_sim; % simulation 
end


EEG(:,:) = F(7:70,:); %% Si ca bug ici c'est que Brainstorm a changé le nom d'enregistement des données EEG, préviens moi ! 

if pre_post == 2
EEG(32,:) = EEG(32,:) + sig; % ajout (ou pas) de la simulation sur C3
end

EMG = F(EMG_chanel,:) + sig; 

if Utiliser_laplacien == 0 
EEG_Ana = EEG;
end 
if Utiliser_laplacien == 1
EEG_Ana = laplacian_perrinX(EEG,[EEG_64.chanlocs.X],[EEG_64.chanlocs.Y],[EEG_64.chanlocs.Z]);
end


%% Calcul de la cohérence

for chan = 1:64
[coh(file,chan,:),frex] = mscohere(EEG_Ana(chan,Interval_1),EMG(Interval_1),hann(Hanning_window),[],Freq_space,srate);
end

end

COH_save{pre_post}(:,:,:,sujet) = coh; % epoch % electrode % frex % subject

end 
end

% calculer la cohérence à partir de la transformée de Fourrier.
% J'ai fait une fonction qui la calcule à partir des ondelettes mais je ne
% suis pas sur du résultat et c'est beaaaucuoup plus long ... on pourra
% essayer mais à voir déjà si ça marche (la cohérence par ondelette n'est
% pas de base dans MatLab)

%% Affichage Différence

% Différence de Cohérence sur Electrode_diff entre post-pré
figure(100)
plot(frex,squeeze(mean(mean(COH_save{2}(:,Electrode_diff,:,:),1),4)) - squeeze(mean(mean(COH_save{1}(:,Electrode_diff,:,:),1),4)))
ylim(ylim_d)
title([ 'Différence de cohérence : Post-Pré , ' num2str(EEG_64.chanlocs(Electrode_diff).labels)]) % dans l'exemple entre 1 et 4 secondes
grid on
h_a = figure(100);
name_a = ['difference electrode' num2str(EEG_64.chanlocs(Electrode_diff).labels) ];
saveas(h_a,name_a,'jpg')


%% Affichage cohérence pour chaque groupe
% affichage de la cohérence 
% en haut : cohérence sur Electrode_diff
% en dessous = carte corticale pour chaque groupe aux fréquences spécifiées

for pre_post = 1:2
    
figure(pre_post)
subplot(3,2, [1 2])
plot(frex,squeeze(mean(mean(COH_save{pre_post}(:,Electrode_diff,:,:),1),4)))
title([ 'Cohérence ' condition{pre_post} ' ' num2str(EEG_64.chanlocs(Electrode_diff).labels)]) % dans l'exemple entre 1 et 4 secondes
ylim(ylim_c)


subplot(323)
topoplotIndie(squeeze(mean(mean(COH_save{pre_post}(:,:,dsearchn(frex',f1),:),1),4)),EEG_64.chanlocs,'numcontour',0)
title({num2str(f1) 'Hz'})
set(gca,'clim', c_lim)
colorbar

subplot(324)
topoplotIndie(squeeze(mean(mean(COH_save{pre_post}(:,:,dsearchn(frex',f2),:),1),4)),EEG_64.chanlocs,'numcontour',0)
title({num2str(f2) 'Hz'})
set(gca,'clim', c_lim)
colorbar

subplot(325)
topoplotIndie(squeeze(mean(mean(COH_save{pre_post}(:,:,dsearchn(frex',f3),:),1),4)),EEG_64.chanlocs,'numcontour',0)
title({num2str(f3) 'Hz'})
set(gca,'clim', c_lim)
colorbar

subplot(326)
topoplotIndie(squeeze(mean(mean(COH_save{pre_post}(:,:,dsearchn(frex',f4),:),1),4)),EEG_64.chanlocs,'numcontour',0)
title({num2str(f4) 'Hz'})
set(gca,'clim',  c_lim)
colorbar

end

h_b = figure(1);
name_b = ['Coherence Pre'];
saveas(h_b,name_b,'jpg')

h_c = figure(2);
name_c = ['Coherence Post'];
saveas(h_c,name_c,'jpg')
