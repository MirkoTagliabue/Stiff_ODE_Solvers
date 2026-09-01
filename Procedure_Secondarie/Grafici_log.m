
% I valori di N_vett ed Err_vett qui già in input sono quelli ottenuti con
% il metodo di Eulero Implicito con t_0=0, T=1000.

N_vett = [1000, 10000, 100000, 1000000];
Err_vett = [3.53628e+06, 4.68461e+05, 4.84397e+04, 4.86061e+03];


m = ( log(Err_vett(end))-log(Err_vett(1)) ) / (log(N_vett(end)) - log(N_vett(1)));


fprintf("Il coefficiente angolare vale: \n\t %.4f \n\n", m);


str_m = sprintf('coeff. angol. m = %.4f', m);
figure(1), loglog(N_vett,Err_vett), 
title('Andamento errore in scala logaritmica'),
subtitle(str_m);
xlabel('Numero nodi'), ylabel('Errore'),

