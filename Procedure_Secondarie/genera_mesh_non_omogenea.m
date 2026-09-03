

function t_vett = genera_mesh_non_omogenea(t_0,T,passo_1,passo_2,passo_3)

    mesh_1 = t_0 : passo_1 : 20;
    mesh_2 = (20+passo_2) : passo_2 : 100;
    mesh_3 = (100+passo_3) : passo_3 : T;

    
    if mesh_3(end) ~= T
        mesh_3(end+1) = T;
    end


    t_vett = [mesh_1, mesh_2, mesh_3];
    
return

