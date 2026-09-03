

function U = Eulero_Implicito(M,y_0,t_vett)

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

% Questo codice applica un metodo di Eulero implicito.
% Il presente metodo ha ordine 1
% Il presente metodo numerico è L-stabile.

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
    % Risolvo il sistema di ODE applicando il metodo numerico

    % Definisco la matrice identità:
    I_15 = eye(15);
    

    % Oss: I calcoli per ricavare l'espressione di u_n sono già
    % stati parzialmente svolti con carta e penna, arrivando alle
    % formule utilizzate in seguito.


    % Definisco la matrice necessaria per risolvere il sistema lineare.
    % Su tale matrice uso poi la fattorizazione LU al fine di ottimizzare 
    % il costo computazionale dal momento che la matrice A non variano nel 
    % ciclo sottostante.
    % (dal momento che la variabile U è già utilizzata per altri scopi, 
    % in seguito chiamo V la matrice triangolare superiore)

    A = I_15 - h*M;

    [L,V,P] = lu(A);



    % Inizializzo prima del ciclo:

    % Piccola osservazione notazionale: il valore iniziale sarebbe forse 
    % stato meglio chiamarlo y_1, perchè è vero che comunque sarebbe 
    % y_1 = f(t_0, y(t_0)), ma d'altronde è anche vero che sarebbe il primo 
    % elemento della mesh, perchè in matlab i vettori non partono da 0
        u_n_meno_1 = y_0;


    for n=2:numel(t_vett)

        t = t_vett(n-1);    % Qui t=t_{n-1}, lo uso per ricavare 
                            % u_{n}, nuovo nodo


        % Risolvo il sistema u_n = A \ u_n_meno_1 con la fattorizzazione LU
        w = L \ (P * u_n_meno_1);
        u_n = V \ w;


        U(:,n) = u_n;

    
        % Aggiorno per la prossima iterazione:
        u_n_meno_1 = u_n;


    end


return

