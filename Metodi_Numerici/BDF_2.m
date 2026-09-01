

function U = BDF_2(M,y_0,t_vett)

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

% Questo codice applica un metodo multistep a 2 passi di tipo Backward 
% Differentiation Formulae, il metodo ha ordine 2. 
% Il presente metodo numerico è L-stabile 

% NB: Questo codice è stato ottimizzato ad hoc per il problema in esame e 
% non è stato pensato per funzionare in generale.
% In particolare richiede che:
% -  Il sistema di ODE sia della forma:  y'(t) = M \cdot y(t)
% -  M abbia dimensione 15x15
% -  la mesh in t_vett sia omogenea (tutti i nodi distino esattamente h)

% -  Oss: a differenza dei metodi di Runge Kutta, questo metodo numerico 
%    non è facilmente implementabile su mesh non omogenea senza prima 
%    derivare il nuovo metodo teoricamente con carta e penna, poichè 
%    differenti lunghezze di intervalli modificano necessariamente le 
%    differenze finite usate per approssimare le derivate.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inizializzo il passo di discretizzazione e la soluzione numerica:
    t_0 = t_vett(1);
    h = t_vett(2) - t_vett(1);

    U=zeros(15,numel(t_vett));
    U(:,1)=y_0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determino i passi iniziali

    % Per far partire il metodo mi servono i 2 nodi iniziali ricavati con 
    % un metodo numerico ad un passo di ordine almeno 2. 
    % Utilizzo a tal fine ad esempio il metodo di Radau_IIA_5 (L-stabile) 
    % di ordine 5
    % Oss: il primo di questi 2 punti è y_0 ed è comunque già noto

    mesh_aux = [t_0, t_0+h];
    U_sol_iniz = Radau_IIA_5(M,y_0,mesh_aux);

    u_n_meno_2 = y_0;  % = U_sol_iniz(:,1);
    u_n_meno_1 = U_sol_iniz(:,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                  IDEA DI RISOLUZIONE: 
    % Metodo BDF2:
    % u_{n} = 4/3*u_{n-1} - 1/3*u_{n-2} + 2/3*h*f_{n}
    % dove con f_n si intende:  f(t_n, u_n) = M * u_n
    % ma: f_{n} = M * u_{n}   ,  quindi otteniamo: 
    % (I-2/3*h*M) * u_{n}  =  4/3*u_{n-1} - 1/3*u_{n-2}
    % ottengo quindi un sistema lineare del tipo:  Ax=b, dove x=u_{n}.
    % Ciclo quindi su n, da n=3 fino ad numel(t_vett)

    % Inoltre, dal momento che dovremo risolvere molti sistemi della forma
    % Ax=b, con A matrice che non viene mai modificata all'interno del
    % ciclo, conviene utilizzare la fattorizzazione LU per ottimizzare il
    % costo computazionale. (Dal momento che la variabile U è già 
    % utilizzata per altri scopi, chiamo V la matrice triagolare superiore)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inizializzo la matrice A ed uso la fattorizazzione LU:

    I_15 = eye(15);
    A = I_15 - 2/3*h*M;
    [L,V,P] = lu(A);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Risolvo il sistema di ODE applicando il metodo numerico

    for n=3:numel(t_vett)

        t = t_vett(n-1);    % Qui t=t_{n-1}, lo uso per ricavare 
                            % u_{n}, nuovo nodo

        b = 4/3*u_n_meno_1 - 1/3*u_n_meno_2;


        % Risolvo il sistema u_n = A\b con la fattorizzazione LU:
        w = L\(P*b);
        u_n = V\w;


        U(:,n) = u_n;


        % Aggiorno per la prossima iterazione:
        u_n_meno_2 = u_n_meno_1;
        u_n_meno_1 = u_n;

    end

return

