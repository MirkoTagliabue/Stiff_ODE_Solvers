

function U = Crank_Nicolson_mesh_non_omogenea(M,y_0,t_vett)

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

% Questo codice applica un metodo di Crank_Nicolson ed è applicabile anche
% su mesh non omogenea.
% Il presente metodo ha ordine 2
% Il presente metodo numerico è A-stabile.
% Questo metodo numerico è inoltre un metodo di Runge-Kutta diagonalmente
% implicito, ed il fatto che sia solo diagonalmente implicito è stato 
% sfruttato per ottimizzare il codice.

% NB: Questo codice è stato ottimizzato ad hoc per il problema in esame e 
% non è stato pensato per funzionare in generale.
% In particolare richiede che:
% -  Il sistema di ODE sia della forma:  y'(t) = M \cdot y(t)
% -  M abbia dimensione 15x15

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definisco ed inizializzo la soluzione numerica:
  
    U=zeros(15,numel(t_vett));
    U(:,1)=y_0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Risolvo il sistema di ODE applicando il metodo numerico

    % Definisco la matrice identità:
    I_15 = eye(15);


    % Inizializzo prima del ciclo:

    % Piccola osservazione notazionale: il valore iniziale sarebbe forse 
    % stato meglio chiamarlo y_1, perchè è vero che comunque sarebbe 
    % y_1 = f(t_0, y(t_0)), ma d'altronde è anche vero che sarebbe il primo 
    % elemento della mesh, perchè in matlab i vettori non partono da 0
        u_n_meno_1 = y_0;


    
    for n=2:numel(t_vett)

        t = t_vett(n-1);    % Qui t=t_{n-1}, lo uso per ricavare 
                            % u_{n}, nuovo nodo
        
        h = t_vett(n) - t_vett(n-1);


        % Oss: I calcoli per ricavare le espressioni di k_1 e k_2 sono già
        % stati parzialmente svolti con carta e penna, arrivando alle
        % seguenti formule:
        
        % Determino k_1: 
        k_1 = M * u_n_meno_1;

        % Determino k_2: 
        A = 2*I_15 - h*M;
        b = 2*k_1 + h*(M*k_1);
        k_2 = A\b;
        
        % determino y_n
        u_n  =  u_n_meno_1  +  h * 1/2 * k_1  +  h * 1/2 * k_2;


        U(:,n) = u_n;

    
        % Aggiorno per la prossima iterazione:
        u_n_meno_1 = u_n;


    end

return 

