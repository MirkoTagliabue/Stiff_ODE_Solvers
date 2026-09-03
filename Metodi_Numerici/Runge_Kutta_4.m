

function U = Runge_Kutta_4(M,y_0,t_vett)

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

% Questo codice applica un metodo di Runge Kutta a 4 stadi esplicito 
% Il presente metodo ha ordine 4. Il metodo in questione è il metodo di
% Runge Kutta a 4 stadi esplicito più famoso, sebbene non abbia un nome
% particolare è chiamato semplicemente "Runge-Kutta 4"

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
    % Inizializzo il passo di discretizzazione e la soluzione numerica:
    
    h = t_vett(2) - t_vett(1);

    U=zeros(15,numel(t_vett));
    U(:,1)=y_0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definisco le costanti del metodo di Runge Kutta 4


    % nodi di discretizzazione temporale
    c = [0, 1/2, 1/2, 1];
    

    % Matrice di Butcher
    A = zeros(4,4);

    A(2,1) = 1/2;
    A(3,2) = 1/2;
    A(4,3) = 1;
    

    % Pesi per la combinazione lineare:
    b = [1/6, 1/3, 1/3, 1/6];
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Risolvo il sistema di ODE applicando il metodo numerico


    % Inizializzo prima del ciclo:

    % Piccola osservazione notazionale: il valore iniziale sarebbe forse 
    % stato meglio chiamarlo y_1, perchè è vero che comunque sarebbe 
    % y_1 = f(t_0, y(t_0)), ma d'altronde è anche vero che sarebbe il primo 
    % elemento della mesh, perchè in matlab i vettori non partono da 0
        u_n_meno_1 = y_0;


    for n=2:numel(t_vett)

        t = t_vett(n-1);    % Qui t=t_{n-1}, lo uso per ricavare 
                            % u_{n}, nuovo nodo
        
        
        % Calcolo i vari coefficienti k
        k = zeros(15,4);   % matrice dei coeffienti k come vettori colonna

        for i=1:4

            % Variabile ausiliaria per calcolare k_i = k(:,i)
            Sum_parz = u_n_meno_1;
               
            for j=1:(i-1)
                Sum_parz = Sum_parz  +  h * A(i,j) * k(:,j);  
            end

            k(:,i) = M * Sum_parz;

        end
        


        % Calcolo u_n
        u_n = u_n_meno_1; 
        for i=1:4
            u_n  =  u_n  +  h * b(i) * k(:,i);
        end


        % Inserisco il valore appena trovato nella matrice soluzione
        U(:,n) = u_n;


        % aggiorno per la prossima iterazione
        u_n_meno_1 = u_n;


    end   % end for

return   % end function

