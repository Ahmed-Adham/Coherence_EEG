%% Fonction pour le Calcul de la cohérence Cortico-musculaire par FFT
%% Ahmed ADHAM - 2022

%%
% Cette fonction va chercher les dossiers les fichiers contenant "data" dans leur nom calculer la cohérence.
% Il faut spécifier le chemin des fichier pré-vibration et des fichiers post-vibration ! 
% Les électrodes d'intérêt pour l'EEG sont les électrodes entre 7 et 70 dans le fichier
% Une fois la cohérence calculée elle est enregistrée dans la variable COH_save 
% que tu peux sauvegarder / exporter pour des stats par la suite 
%%

clear 

load('EEG_Struc_64.mat') % charger noms électrodes + positions
%
%% définition des paramètres


srate = 1024 %Hz - Fréquence d'échantillonage en Hertz 
Interval_1 = [4000:8000] % Interval de temps où calculer la cohérence dans le signal
repertoire_1 = '/Users/Ahmed/Downloads/wetransfer_s7j1contract2b_interpbad_notch_band_2022-07-21_1633/pre/'; % Répertoire où se trouvent les .mat Pré-vibration
repertoire_2 = '/Users/Ahmed/Downloads/wetransfer_s7j1contract2b_interpbad_notch_band_2022-07-21_1633/post/'; % Répertoire où se trouvent les .mat Post-vibration
EMG_chanel = 6; % Quel canal EMG utiliser : mettre 6 pour FCR, 5 pour ECR 
Freq_space = [1:0.5:50]; % Espace de fréquences à étudier (1 à 50Hz tous les 0.5Hz). Attention au temps de calcul ! 
Utiliser_laplacien = 1; %0 = pas de laplacien ; 1 = laplacien
condition{1} = 'pre-vibration' % nom condition_1
condition{2} = 'post-vibration'% nom condition_2
Hanning_window = 1000; % La cohérence par FFT demande une fenêtre d'ordre du filtre : 
% c'est un paramètre un peu arbitraire qui détermine la résultion de
% fréquence ... (trop (ex 2000) sera illisible, pas assez (ex 100) sera tout
% mélangé; du coup j'ai mit 1000 mais libre à chacun de changer


%% c'est parti :) ! 

for pre_post = 1:2 
    
Name_file{1} = dir(fullfile(repertoire_1, 'data*'))
Name_file{2} = dir(fullfile(repertoire_2, 'data*'))

for file = 1:length(Name_file{pre_post})

disp(strcat({'Il reste : ' num2str(file),' / ', num2str(length(Name_file{pre_post})) 'dans le groupe' num2str(pre_post)}))

if pre_post == 1
file_to_load = strcat(repertoire_1,Name_file{pre_post}(file).name);
load(file_to_load);
end

if pre_post == 2
file_to_load = strcat(repertoire_2,Name_file{pre_post}(file).name);
load(file_to_load);
end

EEG(:,:) = F(7:70,:); %% Si ca bug ici c'est que Brainstorm a changé le nom d'enregistement des données EEG, préviens moi ! 
EMG = F(EMG_chanel,:);

if Utiliser_laplacien == 0 
EEG_Ana = EEG;
end 
if Utiliser_laplacien == 1
EEG_Ana = laplacian_perrinX(EEG,[EEG_64.chanlocs.X],[EEG_64.chanlocs.Y],[EEG_64.chanlocs.Z]);
end

times = 0:1/1000:size(EEG_Ana,2)/1000; % temps

%% Calcul de la cohérence

for chan = 1:64
[coh(file,chan,:),frex] = mscohere(EEG_Ana(chan,Interval_1),EMG(Interval_1),hann(Hanning_window),[],Freq_space,srate);
end

end

COH_save{pre_post}(:,:,:) = coh;

end
% calculer la cohérence à partir de la transformée de Fourrier.
% J'ai fait une fonction qui la calcule à partir des ondelettes mais je ne
% suis pas sur du résultat et c'est beaaaucuoup plus long ... on pourra
% essayer mais à voir déjà si ça marche (la cohérence par ondelette n'est
% pas de base dans MatLab) 

%% Paramètres d'affichage des résultats 


ylim_d = [-0.2 0.2]; % échelle cohérence figure de différence (figure 100)
ylim_c = [0 0.3]; % échelle cohérence figure individuelle (figure 1/2)
c_lim = [0 0.3]; % échelle des couleurs entre 0 et 1 pour les cartes globales (figure 1/2)

f1 = 8;  % Fréquence affichée pour la figure 1
f2 = 12; % Fréquence affichée pour la figure 2
f3 = 20; % Fréquence affichée pour la figure 3
f4 = 30; % Fréquence affichée pour la figure 4

%% Affichage Différence


%électrode sur laquelle calculer la différence de cohérence (32 = C3, cf tableau en bas) 
Electrode_diff = 32; 

% Différence de Cohérence sur Electrode_diff entre post-pré
figure(100)
plot(frex,squeeze(mean(COH_save{2}(:,Electrode_diff,:),1))-squeeze(mean(COH_save{1}(:,Electrode_diff,:),1)))
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
plot(frex,squeeze(mean(COH_save{pre_post}(:,Electrode_diff,:),1)))
title([ 'Cohérence ' condition{pre_post} ' ' num2str(EEG_64.chanlocs(Electrode_diff).labels)]) % dans l'exemple entre 1 et 4 secondes
ylim(ylim_c)


subplot(323)
topoplotIndie(squeeze(mean(COH_save{pre_post}(:,:,dsearchn(frex',f1)),1)),EEG_64.chanlocs,'numcontour',0)
title({num2str(f1) 'Hz'})
set(gca,'clim', c_lim)
colorbar

subplot(324)
topoplotIndie(squeeze(mean(COH_save{pre_post}(:,:,dsearchn(frex',f2)),1)),EEG_64.chanlocs,'numcontour',0)
title({num2str(f2) 'Hz'})
set(gca,'clim', c_lim)
colorbar

subplot(325)
topoplotIndie(squeeze(mean(COH_save{pre_post}(:,:,dsearchn(frex',f3)),1)),EEG_64.chanlocs,'numcontour',0)
title({num2str(f3) 'Hz'})
set(gca,'clim', c_lim)
colorbar

subplot(326)
topoplotIndie(squeeze(mean(COH_save{pre_post}(:,:,dsearchn(frex',f4)),1)),EEG_64.chanlocs,'numcontour',0)
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

%% Channels 


%{ 
C3 = 32, C4 = 36

    1     2     3       4    5            6     7    8     9     10
0  'FP1';'FPZ'; FP2';'AF7';'AF3';    0  'Afz';'AF4';'AF8';'F9';'F7';
11 'F3'; 'F1';'Fz';'F2';'F4';       16  'F8';'F10';'FT9';'FT7';'FC5';
21 'FC3';'FC1';'FCz';'FC2';'FC4';   26  'FC6';'FT8';'FT10';'T9';'T7';
31 'C5'; 'C3';'C1';'Cz';'C2';       36  'C4';'C6';'T8';'T10';'TP9';
41 'TP7';'CP5';'CP3';'CP1';'CPz';   46  'CP2';'CP4';'CP6';'TP8';'TP10';
51 'P9'; 'P7';'P3';'P1';'Pz';       56  'P2';'P4';'P8';'P10';'PO7';
61 'POZ';'PO8';'O1';'O2'
}% 

%% 
