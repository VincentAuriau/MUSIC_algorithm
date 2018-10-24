close all; clear all; clc;


% SIGNAL RECU

% Caractérisiques du signal
frequences = [450;800]; % frequences des signaux

M = length(frequences); % Nombre de sources

NbrMesures = 1000; % Nombre de mesures réalisées par les capteurs
T=0:0.0001:0.0001*NbrMesures;%le temps
C=[2,2]; %Amplitude des signaux
B=zeros(2,NbrMesures);
for i=1:1:NbrMesures
    for k=1:1:2
        B(k,i)=C(k)*exp(1j*(2*pi*frequences(k)*T(i)));%signal reçu par le premier capteur à l'instant i issu du kième signal incident
    end
end

% Pulsation
w =2*pi*[frequences].';

NbrCapteurs = 10; % Nombre de capteurs
r = [(0:9)].'; % Réseau linéaire uniforme
% Matrice des vecteurs directeurs de l'antenne
o=r*w;
A=zeros(10,2);
for i=1:1:10
    for j=1:1:2    
        A (i,j)= cosd(o(i,j))+1j*sind(o(i,j));
    end
end

% Bruit
SNR =10; % SignaltoNoiseRatio
n = max(C)*(randn(NbrCapteurs,NbrMesures)*10^(-SNR/20) + 1j*randn(NbrCapteurs,NbrMesures)*10^(-SNR/20))/sqrt(2);

% Signal reçu
x = A*B + n;

% Algorithme MUSIC

% Matrice de covariance des signaux reçus
Rxx = x*x'/NbrMesures;

% Decompositions en vecteurs propres
[E,D] = eig(Rxx);
lambda=diag(D);
% Les vecteurs sont triés apr valeurs propres croissantes
En = E(:,1:end-M); % Vecteurs propres appartenant à l'espace bruit

% Recherche des frequences
f = (0:1:1000).'; % valeurs des frequences pour la recherche

 
% Point de recherche correspondant
VectF = 2*pi*[f].';
p=r*VectF;
VectDirecteurs=zeros(10,length(f));
for i=1:1:10
    for j=1:1:length(f)      
        VectDirecteurs (i,j)= cosd(p(i,j))+1j*sind(p(i,j));%Correspond au a(theta) du rapport
    end
end
% Spectre de MUSIC

Z = abs(VectDirecteurs'*En*En'*VectDirecteurs); %Z=1/P(thêta)fonction=zeros(1,length(Z));
for i=1:1:length(Z)
    Z(i,i)=1/Z(i,i);
    fonction(i)=Z(i,i);
end
% Affichage
figure();
plot(f,diag(Z));
title('Algorithme MUSIC');
xlabel('frequences (Hz)');
ylabel('Spectre Music ');
grid on; axis tight;
max1=0;
freqmax1=0;
max2=0;
freqmax2=0;
[pks,fr] = findpeaks(fonction);%Détection des pics avec MatLab

%Récupération des deux plus grands pics (max1 & max2) et leur fréquences
%correspondantes (freqmax1 & freqmax2)
for i=1:1:length(pks)
    if pks(i)>max1
        if pks(i)>max2 %max2 est le plus grand des pics
            max1=max2;%Si on doit remplacer max2, alors max1 prend la valeur de max2
            freqmax1=freqmax2;%On réalise les changements correspondants pour les fréquences
            max2=pks(i);
            freqmax2=fr(i);
        else
            max1=pks(i);
            freqmax1=fr(i);
        end
    end
end

