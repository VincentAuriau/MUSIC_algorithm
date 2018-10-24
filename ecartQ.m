function [ Q ] = ecartQ( STNR )
%Fonction qui calcule l'écart quadratique moyen pour 100 Simulations pour
%un SNR donné
%
close all; clc;
Q=0;%notre écart quadratique

% SIGNAL RECU

% directions des sources du signal
frequences = [800;450]; % frequences des signaux

M = length(frequences); % Nombre de sources

w=0.01;%pulsation
NbrMesures = 1000; % Nombre de mesures réalisées par les capteurs
T=0:0.0001:0.0001*NbrMesures;
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


G=zeros(100,1000);
%Maintenant qu'on a formé le signal émis, on réalise les cent simulations
%et on fait varier le bruit à chaque simulation
for h=1:1:100
    G(h)=[0];
    
    %Bruit
    n = max(C)*(randn(NbrCapteurs,NbrMesures)*10^(-(STNR)/20) + 1j*randn(NbrCapteurs,NbrMesures)*10^(-(STNR)/20))/sqrt(2);
    
    % Signal reçu
    x = A*B + n;
    
    % Algorithme MUSIC
    
    % Matrice de covariance des signaux reçus
    Rxx = x*x'/NbrMesures;
    
    % Decompositions en vecteurs propres
    [E,D] = eig(Rxx);
    lambda=diag(D);% les valeurs prorpres
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
            
            VectDirecteurs (i,j)= cosd(p(i,j))+1j*sind(p(i,j));
        end
    end
    % Spectre de MUSIC
    
    Z = abs(VectDirecteurs'*En*En'*VectDirecteurs);
    fonction=zeros(length(Z),1);
    for i=1:1:length(Z)
        Z(i,i)=1/Z(i,i);
        fonction(i)=Z(i,i);
    end
    max1=0;
    freqmax1=0;
    max2=0;
    freqmax2=0;
    [pks,fr] = findpeaks(fonction);
    for i=1:1:length(pks)
        if pks(i)>max1
            if pks(i)>max2
                max1=max2;
                freqmax1=freqmax2;
                max2=pks(i);
                freqmax2=fr(i);
            else
                max1=pks(i);
                freqmax1=fr(i);
            end
        end
    end
   
%Calcul de l'écart quadadratique moyen en fonction des fréquences attendues
%(ici on a juste considéré 450 Hz et 800Hz) :le carré de la différence
%entre les fréquences attendues et obtenues
Q=Q+(max(freqmax1,freqmax2)-800).^2+(min(freqmax1,freqmax2)-450).^2;
end
Q=Q/100;
end


