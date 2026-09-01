

function [M,y_0] = genera_dati_iniziali()


    % Crea la matrice M
    M = zeros(15,15);

    for j=1:14
        M(j,j) = -(j^2);
        M(j,j+1) = 150;
    end
    M(15,15)= -15^2;



    % Crea il vettore y_0
    y_0 = ones(15,1);
    y_0(15,1)=10;


return

