

function U = Pareschi_Russo(M,y_0,t_vett)

% La presente funzione risolve un sistema autonomo di ODE della forma:
% y'(t) = M \cdot y(t) , con M matrice 15x15 a coefficienti costanti.
% La presente funzione restituisce la matrice U avente numel(t_vett)
% colonne e 15 righe. Ciascuna colonna individua un istante della mesh,
% mentre ciascuna delle 15 righe individua una componente y_i del vettore
% soluzione y \in R^15. 
% In tal senso, y_5 = U(5,:) è un vettore di dimensione numel(t_vett) che 
% contiene tutti i valori di y_5 sulla mesh t_vett.
% M è la matrice di dimensione 15x15 che descrive il sistema di ODE, mentre 
% y_0 è il valore iniziale del problema di Cauchy, t_vett la mesh ed h il 
% passo di discretizzazione. 

% Questo codice applica un metodo di Runge Kutta di Pareschi-Russo.
% L'ordine e la stabilità del metodo dipendono dal parametro R, e per tale
% discussione si rimanda ai commenti nella dedicata sezione sottostante.
% Questo metodo numerico è inoltre un metodo di Runge-Kutta diagonalmente
% implicito, ed il fatto che sia solo diagonalmente implicito è stato 
% sfruttato per ottimizzare il codice.

% NB: Questo codice è stato ottimizzato ad hoc per il problema in esame e 
% non è stato pensato per funzionare in generale.
% In particolare richiede che:
% -  Il sistema di ODE sia della forma:  y'(t) = M \cdot y(t)
% -  M abbia dimensione 15x15
% -  la mesh in t_vett sia omogenea (tutti i nodi distino esattamente h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inizializzo il passo di discretizzazione e la soluzione numerica:
    
    h = t_vett(2) - t_vett(1);

    U=zeros(15,numel(t_vett));
    U(:,1)=y_0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definisco le costanti del metodo di Pareschi-Russo
    
    R = 1 + sqrt(2)/2;
    % R = 1/2;
    % R = 1/8;
    
    c_1 = R;
    c_2 = 1-R;


    a_11 = R;
    a_12 = 0;

    a_21 = 1-2*R;
    a_22 = R;


    b_1 = 1/2;
    b_2 = 1/2;


% Il metodo di Pareschi-Russo ha ordine 2, per ogni valore del parametro R.
% Il metodo è A-stabile se e solo se R >= 1/4.
% Il metodo è L-stabile se e solo se R = 1 \pm sqrt(2)/2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Risolvo il sistema di ODE applicando il metodo numerico

    % Definisco la matrice identità:
    I_15 = eye(15);
    

    % Oss: I calcoli per ricavare le espressioni di k_1 e k_2 sono già
    % stati parzialmente svolti con carta e penna, arrivando alle
    % formule utilizzate in seguito.


    % Definisco le matrici necessarie per risolvere i due sistemi lineari.
    % Su tali matrici uso poi la fattorizazione LU al fine di ottimizzare 
    % il costo computazionale dal momento che le matrici A_1 ed A_2
    % non variano nel ciclo sottostante.
    % (dal momento che la variabile U è già utilizzata per altri scopi, 
    % in seguito chiamo V la matrice triagolare superiore)

    A_1 = I_15 - h*a_11*M;
    [L_1,V_1,P_1] = lu(A_1);

    A_2 = I_15 - h*a_22*M;
    [L_2,V_2,P_2] = lu(A_2);



    % Inizializzo prima del ciclo:

    % Piccola osservazione notazionale: il valore iniziale sarebbe forse 
    % stato meglio chiamarlo y_1, perchè è vero che comunque sarebbe 
    % y_1 = f(t_0, y(t_0)), ma d'altronde è anche vero che sarebbe il primo 
    % elemento della mesh, perchè in matlab i vettori non partono da 0
        u_n_meno_1 = y_0;


    for n=2:numel(t_vett)

        t = t_vett(n-1);    % Qui t=t_{n-1}, lo uso per ricavare 
                            % u_{n}, nuovo nodo

        % Oss: il termine noto del sistema A_1 * k_1 = b_1, anzichè b_1, 
        % lo chiamo z_1, dal momento che b_1 indica già uno dei pesi del 
        % metodo di RK
        z_1 = M * u_n_meno_1;

        % Risolvo il sistema k_1 = A_1 \ z_1 con la fattorizzazione LU:
        w_1 = L_1 \ (P_1 * z_1);
        k_1 = V_1 \ w_1;



        % Ora ricavo z_2, per il cui calcolo mi serve k_1
        z_2 = M * (u_n_meno_1 + h*a_21*k_1);

        % Risolvo il sistema k_2 = A_2 \ z_2 con la fattorizzazione LU:
        w_2 = L_2 \ (P_2 * z_2);
        k_2 = V_2 \ w_2;



        u_n = u_n_meno_1 + h*b_1*k_1 + h*b_2*k_2;


        U(:,n) = u_n;

    
        % Aggiorno per la prossima iterazione:
        u_n_meno_1 = u_n;


    end


return

