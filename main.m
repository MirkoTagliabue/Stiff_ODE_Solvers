

clearvars, close all

addpath('Metodi_Numerici');
addpath('Procedure_Secondarie');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ELENCO DEI METODI IMPLEMENTATI:

%   1) METODI CON MESH OMOGENEA:
% 'Radau_IIA_5'
% 'BDF_2'
% 'BDF_3'
% 'Pareschi_Russo'
% 'Eulero_Implicito'
% 'Runge_Kutta_4'


%   2) METODI CON MESH NON OMOGENEA:
% 'Gauss_Legendre_6_mesh_non_omogenea'
% 'Crank_Nicolson_mesh_non_omogenea'


%   3) METODI ADATTIVI:
% 'Gauss_Legendre_2_6'
% 'Runge_Kutta_Fehlberg_4_5'



            % SELEZIONARE IL METODO DA USARE:
metodo = 'Gauss_Legendre_2_6';

% NB: è possibile inoltre cambiare alcuni parametri (come h, T, passo_N,
% toll, ecc) nel menù sotto

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETRI INIZIALI:

t_0 = 0;
T = 1000;


% parametri per una mesh omogenea:
h = 0.01;


% parametri per una mesh NON omogenea (VINCOLI: necessario T > 100+passo_3 
% e necessario t_0 < 20-passo_1)
% Oss: se passo_1 = passo_2 = passo_3, si ottiene una mesh omogenea
passo_1 = 0.001;     % passo tra t_0 e 20
passo_2 = 0.01;      % passo tra 20 e 100
passo_3 = 10;        % passo tra 100 e T


% parametri per una mesh adattiva:
h_iniz = 0.001;  % h iniziale "di prova"
toll = 1e-3;     % stima dell'errore assoluto locale
r = 1/2;         % coefficiente cautelativo per il calcolo del nuovo passo h


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definiamo la matrice del sistema ed i parametri iniziali:
[M, y_0] = genera_dati_iniziali();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RISOLVO IL PROBLEMA CON UN METODO NUMERICO:
% si può modificare la variabile metodo nell'intestazione del main

switch metodo
    
            % Metodi su mesh omogenea:  
    
    case 'Radau_IIA_5' 
        t_vett = genera_mesh_omogenea(t_0,T,h);
        U = Radau_IIA_5(M,y_0,t_vett);

    case 'BDF_2'
        t_vett = genera_mesh_omogenea(t_0,T,h);
        U = BDF_2(M,y_0,t_vett);

    case 'BDF_3'
        t_vett = genera_mesh_omogenea(t_0,T,h);
        U = BDF_3(M,y_0,t_vett);

    case 'Pareschi_Russo'
        t_vett = genera_mesh_omogenea(t_0,T,h);
        U = Pareschi_Russo(M,y_0,t_vett);

    case 'Eulero_Implicito'
        t_vett = genera_mesh_omogenea(t_0,T,h);
        U = Eulero_Implicito(M,y_0,t_vett);

    case 'Runge_Kutta_4'
        t_vett = genera_mesh_omogenea(t_0,T,h);
        U = Runge_Kutta_4(M,y_0,t_vett);


            % Metodi su mesh non omogenea:

    case 'Gauss_Legendre_6_mesh_non_omogenea'
        t_vett = genera_mesh_non_omogenea(t_0,T,passo_1,passo_2,passo_3);
        U = Gauss_Legendre_6_mesh_non_omogenea(M,y_0,t_vett);

    case 'Crank_Nicolson_mesh_non_omogenea'
        t_vett = genera_mesh_non_omogenea(t_0,T,passo_1,passo_2,passo_3);
        U = Crank_Nicolson_mesh_non_omogenea(M,y_0,t_vett);    


            % Metodi su mesh adattiva:

    case 'Gauss_Legendre_2_6'
        [U, t_vett] = Gauss_Legendre_2_6(M, y_0, t_0, T, h_iniz, toll, r);

    case 'Runge_Kutta_Fehlberg_4_5'
        [U, t_vett] = Runge_Kutta_Fehlberg_4_5(M, y_0, t_0, T, h_iniz, toll, r);


            % Default case

    otherwise
        error(['Metodo non riconosciuto, controlla di aver scritto ' ...
            'correttamente il nome della variabile Metodo'])

        
end % end switch-case



% *************************************************************************
% DETERMINO LA SOLUZIONE ESATTA E L'ERRORE:

% Estrapolo tre funzioni soluzione dalla matrice soluzione: u_1, u_2, u_15
u_1 = U(1,:);
u_2 = U(2,:);
u_15 = U(15,:);

% Inizializzo le corrispettive soluzioni esatte:
y_1 = zeros(1,numel(t_vett));
y_2 = zeros(1,numel(t_vett));
y_15 = zeros(1,numel(t_vett));

Err_ass_glob = 0;
max_norm_y_n = 0;


for n=1:numel(t_vett)

    t = t_vett(n);
    sol_n = expm( M*(t-t_0) ) * y_0;

    y_1(n) = sol_n(1);
    y_2(n) = sol_n(2);
    y_15(n) = sol_n(15);

    e_n = norm( U(:,n) - sol_n );

    norm_y_n = norm( sol_n );

    if e_n > Err_ass_glob
        Err_ass_glob = e_n;
    end


    if norm_y_n > max_norm_y_n
        max_norm_y_n = norm_y_n;
    end


end  %end for


% Calcolo l'errore relativo globale:
Err_rel_glob = Err_ass_glob / max_norm_y_n;


% Stampo gli errori:
fprintf("L'errore assoluto globale vale: \n\t %.5e \n", Err_ass_glob);
fprintf("L'errore relativo globale vale: \n\t %.5e \n", Err_rel_glob);


% Se il metodo è su mesh non omogenea stampo anche il numero di nodi 
if strcmp(metodo, 'Gauss_Legendre_6_mesh_non_omogenea') || strcmp(metodo, 'Crank_Nicolson_mesh_non_omogenea')
    fprintf("Il numero di nodi della mesh e': \n\t %d \n", numel(t_vett));
end


% Se il metodo è adattivo stampo anche il numero di nodi ed il passo minimo e massimo:
if strcmp(metodo, 'Gauss_Legendre_2_6') || strcmp(metodo, 'Runge_Kutta_Fehlberg_4_5')

    fprintf("Il numero di nodi della mesh e': \n\t %d \n", numel(t_vett));
    
    h_min = realmax;  % valore massimo rappresentabile in matlab in precisione doppia
    h_max = 0;
    
    for n=2:numel(t_vett)

        h_n = t_vett(n) - t_vett(n-1);
        
        if h_n < h_min
            h_min = h_n;
        end


        if h_n > h_max
            h_max = h_n;
        end


    end  % end for

    fprintf("Il passo di discretizzazione minimo h_min della mesh e': \n\t %.5f \n", h_min);
    fprintf("Il passo di discretizzazione massimo h_max della mesh e': \n\t %.5f \n", h_max);

end % end if(metodo adattivo)


fprintf("\n\n");


% Ora determino gli errori sulle singole componenti y_1, y_2, y_15:
Err_ass_1 = max( abs(y_1 - u_1) );
Err_ass_2 = max( abs(y_2 - u_2) );
Err_ass_15 = max( abs(y_15 - u_15) );

Err_rel_1 = Err_ass_1 / max( abs(y_1) );
Err_rel_2 = Err_ass_2 / max( abs(y_2) );
Err_rel_15 = Err_ass_15 / max( abs(y_15) );



% *************************************************************************
% RAPPRESENTO I GRAFICI:


% Disegno y_1(t), per maggior intuitività solo nell'intervallo [0,10] per
% visualizzare il picco
str_err_ass_1 = sprintf('Err Ass y_1 = %.3e', Err_ass_1);
str_err_rel_1 = sprintf('Err Rel y_1 = %.3e', Err_rel_1);
figure(1), plot(t_vett, u_1, 'r*-'), hold on,
    % Oppure 'r*-' per enfatizzare i nodi degli adattivi o mesh non omogenee
plot(t_vett, y_1, 'b-'), legend('approssimata', 'esatta'), 
xlabel('t'), ylabel('y(t)'), xlim([t_0, 50]), 
title(sprintf('y_1(t)\n%s\n%s', str_err_ass_1, str_err_rel_1));

% Per matlab (e non Octave) invece si può usare anche (più elegante):
% title('y_1(t)')
% subtitle({str_err_ass_1, str_err_rel_1});





% Disegno y_2(t) in [t_0,5]
str_err_ass_2 = sprintf('Err Ass y_2 = %.3e', Err_ass_2);
str_err_rel_2 = sprintf('Err Rel y_2 = %.3e', Err_rel_2);
figure(2), plot(t_vett, u_2, 'r-'), hold on,
plot(t_vett, y_2, 'b-'), legend('approssimata', 'esatta'),
xlabel('t'), ylabel('y(t)'), xlim([t_0, 5]), 
title(sprintf('y_2(t)\n%s\n%s', str_err_ass_2, str_err_rel_2));



% Disegno y_15(t) in [t_0,4]
str_err_ass_15 = sprintf('Err Ass y_{15} = %.3e', Err_ass_15);
str_err_rel_15 = sprintf('Err Rel y_{15} = %.3e', Err_rel_15);
figure(3), plot(t_vett, u_15, 'r-'), hold on,
plot(t_vett, y_15, 'b-'), legend('approssimata', 'esatta'),
xlabel('t'), ylabel('y(t)'), xlim([t_0, 4])
title(sprintf('y_{15}(t)\n%s\n%s', str_err_ass_15, str_err_rel_15));


