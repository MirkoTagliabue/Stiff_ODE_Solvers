

function U = Gauss_Legendre_6_mesh_non_omogenea(M,y_0,t_vett)

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

% Questo codice applica un metodo di Runge Kutta basato su interpolazione 
% di Gauss-Legendre ed è applicabile anche su mesh non omogenea.
% Il presente metodo ha ordine 6, con p=2s, p=6, s=3. 
% Il presente metodo numerico è A-stabile.

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
    % Definisco le costanti del metodo di Gauss-Legendre, p=2s, p=6, s=3

    c_1 = 1/2 - sqrt(15)/10;
    c_2 = 1/2;
    c_3 = 1/2 + sqrt(15)/10;

    a_11 = 5/36;
    a_12 = 2/9 - sqrt(15)/15;
    a_13 = 5/36 - sqrt(15)/30;


    a_21 = 5/36 + sqrt(15)/24;
    a_22 = 2/9;
    a_23 = 5/36 - sqrt(15)/24;


    a_31 = 5/36 + sqrt(15)/30;
    a_32 = 2/9 + sqrt(15)/15;
    a_33 = 5/36;


    b_1 = 5/18;
    b_2 = 4/9;
    b_3 = 5/18;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Risolvo il sistema di ODE applicando il metodo numerico

    % Definisco la matrice identità:
    I_45 = eye(45);


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

        
        % Definisco le matrici a blocchi necessarie per risolvere il 
        % sistema lineare:
        B = zeros(45,45);

        B(1:15, 1:15) = h*a_11*M;
        B(1:15, 16:30) = h*a_12*M;
        B(1:15, 31:45) = h*a_13*M;

        B(16:30, 1:15) = h*a_21*M;
        B(16:30, 16:30) = h*a_22*M;
        B(16:30, 31:45) = h*a_23*M;

        B(31:45, 1:15) = h*a_31*M;
        B(31:45, 16:30) = h*a_32*M;
        B(31:45, 31:45) = h*a_33*M;


        A = I_45 - B;
        
        
        
        b = zeros(45,1);
        b(1:15) = M * u_n_meno_1;
        b(16:30) = b(1:15);
        b(31:45) = b(1:15);

        % Risolvo il sistema K = A\b 
        K = A\b;

        k_1 = K(1:15);
        k_2 = K(16:30);
        k_3 = K(31:45);


        u_n = u_n_meno_1 + h*b_1*k_1 + h*b_2*k_2 + h*b_3*k_3;


        U(:,n) = u_n;

    
        % Aggiorno per la prossima iterazione:
        u_n_meno_1 = u_n;


    end

return 

