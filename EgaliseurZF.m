% TS217- TP Egalisation
% Pascal Vallet (IPB)
% 2014

clear all;
close all;
clc;

%% Paramètres

% Longueur de la séquence binaire transmise
N=5000; 

% SNR en dB
SNR=[0:20];
sigma=10.^(-SNR/20); % ecart-type du bruit correspondant

% Quelques canaux à tester ...

%h=[0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0.21 0.03 0.07]; % Canal Proakis A
%h=[0.407 0.815 0.407]; % Canal Proakis B
%h=[0.227 0.460 0.688 0.460 0.227]; % Canal Proakis C

%h=[1;0.5]; % CANAL TEST 1
%h=[1;0.1;0.9].'; % CANAL TEST 2
h=[0.5;0.8;0.5]; % CANAL TEST 3

K=length(h); % longueur du canal
P= 20; % ordre des filtres ZF/MMSE/DFE (filtre direct), A COMPLETER
Q= 20; % ordre du filtre de retour DFE, A COMPLETER

%% Simulation des signaux

% Bits
bits=rand(N,1) > 0.5;
s=2*(bits-0.5); % N symboles i.i.d. BPSK

% Quelques variables ...
zZF=zeros(N,length(SNR)); % sortie de l'égaliseur ZF
bitsZF=zeros(N,length(SNR)); % bits estimés après ZF
berZF=zeros(1,length(SNR)); % BER après ZF

% Generation du bruit

% Filtre ZF
% Matrice H
H = conv2(h, eye(P))';

% Retard optimal
vec=[];
for d=1:P+K-1
    e=zeros(P+K-1,1);
    e(d)=1;
    vec=[vec norm((H')*(pinv(H)')*e)];
end
[argvalue, argmax] = max(vec);

% Energie de l'IES
dZF= argmax; % retard optimal du filtre ZF
eZF=zeros(P+K-1,1);
eZF(dZF)=1;
Eies=norm(eZF)-norm((H')*(pinv(H)')*eZF);

fZF= (pinv(H)')*eZF;
%%
for i=1:length(SNR) 
    
    y =  filter(h, 1, s); % Observations en sortie du canal, A COMPLETER
    %y = y(K:end);
    
    % Egalisation ZF/MMSE
    for n=P+K-1:N 
        zZF(n,i)=real(fZF'*y(n:-1:n-P+1));
    end
    bitsZF(:,i)=real(zZF(:,i)) > 0; 
    berZF(:,i)=sum(abs(bits(P+K-1-(dZF-1):N-(dZF-1))-bitsZF(P+K-1:N,i)),1)/(N-(P+K-2)); 
end

%% Graphes
Nfft = 512;

% Fonction de transfert canal/ZF
plot((1:Nfft)/Nfft - 0.5, fftshift(abs(fft(h,Nfft))));
title('Reponse en frequence du filtre canal');
xlabel('Frequence normalisee'); ylabel('Amplitude');

%% Constellation des échantillons reçus y(n)
scatterplot(y);

%% Courbes de BER
% Probabilité d'erreur du canal AWGN = Q(sqrt(2 Eb/N0)), où Eb=1 et N0 = sigma^2
figure;
semilogy(SNR,berZF,'-b^',SNR,berMMSE,'-rd',SNR,berDFE,'-gv',SNR,1-normcdf(sqrt(2./sigma.^2),0,1),'k-','LineWidth',3);
grid on;
xlabel('SNR');
ylabel('BER')
legend('ZF','MMSE','DFE','AWGN','Location','SouthWest');




%%
