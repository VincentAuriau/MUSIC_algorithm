%Algorithlme qui permet de tracer l'écart quadratique moyen de MUSIC en
%fonction du SNR
U=0:1:40;
J=0:1:40;
for i=1:1:length(U)
    J(i)=ecartQ(U(i)); %Appel de la fonction ecartQ qui calcule l'écart quadaratique pour 100 simulations pour un SNR donné
end
plot(U,J);
title('Ecart Quadratique moyen en fonction du SNR');
xlabel('SNR');
ylabel('EQ moyen');
