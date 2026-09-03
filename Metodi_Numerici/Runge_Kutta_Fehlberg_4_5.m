

function [U, t_vett] = Runge_Kutta_Fehlberg_4_5(M, y_0, t_0, T, h_iniz, toll, r)

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

% Questo codice applica un metodo di Runge Kutta Fehlberg 4(5) esplicito 
% adattivo: i nodi della mesh ed il passo di discretizzazione sono 
% determinati ad ogni iterazione.
% Il presente metodo ha ordine 5.

% Il presente metodo numerico non è A-stabile, anzi, essendo esplicito è
% caratterizzato da una altissima instabilità. Il motivo per il quale 
% introduciamo questo metodo è solo quello di far notare come anche uno tra
% i migliori metodi di Runge Kutta non vada bene per risolvere un sistema
% di ODE caraterizzato da un elevato coefficiente di Stiffness.

% NB: Questo codice è stato ottimizzato ad hoc per il problema in esame e 
% non è stato pensato per funzionare in generale.
% In particolare richiede che:
% -  Il sistema di ODE sia della forma:  y'(t) = M \cdot y(t)
% -  M abbia dimensione 15x15

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definisco ed inizializzo la soluzione numerica:
  
    U=zeros(15,1);
    U(:,1)=y_0;

    t_vett = [t_0];

    % Tengo traccia per comodità del numero di nodi in una variabile:
    N = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definisco le costanti del metodo di Runge Kutta Fehlberg 4(5)

    p_min_piu_1 = 5;  % valore di 1 superiore all'ordine del metodo meno preciso

    % nodi di discretizzazione temporale
    c = [0, 1/4, 3/8, 12/13, 1, 1/2];
    
    % Matrice di Butcher
    A = zeros(6,6);

    A(2,1) = 1/4;
    
    A(3,1) = 3/32;
    A(3,2) = 9/32;
    
    A(4,1) = 1932/2197;
    A(4,2) = -7200/2197;
    A(4,3) = 7296/2197;
    
    A(5,1) = 439/216;
    A(5,2) = -8;
    A(5,3) = 3680/513;
    A(5,4) = -845/4104;
    
    A(6,1) = -8/27;
    A(6,2) = 2;
    A(6,3) = -3544/2565;
    A(6,4) = 1859/4104;
    A(6,5) = -11/40;
    
    % Pesi per la soluzione (Ordine 5 - per avanzare temporalmente)
    b = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];
    
    % Pesi per la stima dell'errore (Ordine 4 - per il controllo del passo)
    b_hat = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Risolvo il sistema di ODE applicando il metodo numerico


    % Inizializzo prima del ciclo:

    % Piccola osservazione notazionale: il valore iniziale sarebbe forse 
    % stato meglio chiamarlo y_1, perchè è vero che comunque sarebbe 
    % y_1 = f(t_0, y(t_0)), ma d'altronde è anche vero che sarebbe il primo 
    % elemento della mesh, perchè in matlab i vettori non partono da 0
        u_n_meno_1 = y_0;
        h = h_iniz;

    

    while t_vett(end) < T

        t = t_vett(end);    % Qui t=t_{n-1}, lo uso per ricavare 
                            % u_{n}, nuovo nodo
        
        
        % Calcolo i vari coefficienti k
        k = zeros(15,6);   % matrice dei coeffienti k come vettori colonna

        for i=1:6

            % Variabile ausiliaria per calcolare k_i = k(:,i)
            Sum_parz = u_n_meno_1;
               
            for j=1:(i-1)
                Sum_parz = Sum_parz  +  h * A(i,j) * k(:,j);  
            end

            k(:,i) = M * Sum_parz;

        end
        

        % Calcolo u_n
        u_n = u_n_meno_1; 
        for i=1:6
            u_n  =  u_n  +  h * b(i) * k(:,i);
        end


        % Calcolo u_n_hat
        u_n_hat = u_n_meno_1; 
        for i=1:6
            u_n_hat  =  u_n_hat  +  h * b_hat(i) * k(:,i);
        end



        Err_stima = norm( u_n - u_n_hat , 2);
        
        % Coinvolgo la precisione di macchina eps per evitare di dividere
        % per zero in caso di errore nullo
        q = (r * toll / (Err_stima + eps) ) ^ (1/p_min_piu_1);  % coefficiente moltiplicativo
        h_new = h * q;

        

        if Err_stima <= toll
            N=N+1;
            U(:,N) = u_n;
            t_vett(end+1) = t_vett(end) + h;
            h = min(h_new, T-t_vett(end));  % per non sforare T in caso di 
                                            % ultima iterazione
            u_n_meno_1 = u_n;   % aggiorno per la prossima iterazione

        else  % Rifaccio tutto ma con un nuovo h 
            h = h_new;
        end   % end if-else



    end   % end while

return   % end function

