
% TP2 de Statistiques : fonctions a completer et rendre sur Moodle
% Nom :grethen
% Prénom : clémentine
% Groupe : 1SN-N

function varargout = fonctions_TP2_stat(varargin)

    [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});

end

% Fonction centrage_des_donnees (exercice_1.m) ----------------------------
function [x_G, y_G, x_donnees_bruitees_centrees, y_donnees_bruitees_centrees] = ...
                centrage_des_donnees(x_donnees_bruitees,y_donnees_bruitees)
     x_G=mean(x_donnees_bruitees);
     y_G=mean(y_donnees_bruitees);
     x_donnees_bruitees_centrees=x_donnees_bruitees-x_G;
     y_donnees_bruitees_centrees=y_donnees_bruitees-y_G;
end

% Fonction estimation_Dyx_MV (exercice_1.m) -------------------------------
function [a_Dyx,b_Dyx] = ...
           estimation_Dyx_MV(x_donnees_bruitees,y_donnees_bruitees,n_tests)
[x_G, y_G, x_donnees_bruitees_centrees, y_donnees_bruitees_centrees] =centrage_des_donnees(x_donnees_bruitees,y_donnees_bruitees)
psi= -(pi/2) + pi*rand(n_tests,1);
somme=sum((y_donnees_bruitees_centrees-tan(psi)*x_donnees_bruitees_centrees).^2,2);
[~,k]=min(somme);
a_Dyx= tan(psi(k));
b_Dyx = y_G-a_Dyx*x_G;
end

% Fonction estimation_Dyx_MC (exercice_2.m) -------------------------------
function [a_Dyx,b_Dyx] = ...
                   estimation_Dyx_MC(x_donnees_bruitees,y_donnees_bruitees)
[~,n]=size(x_donnees_bruitees);
%création de A
A=ones(n,2)
A(:,1) = x_donnees_bruitees;
%création de b 
B=transpose(y_donnees_bruitees);
pseudo_inverse_A=inv(transpose(A)*A)*transpose(A);
X=pseudo_inverse_A*B
a_Dyx=X(1)
b_Dyx=X(2)
end

% Fonction estimation_Dorth_MV (exercice_3.m) -----------------------------
function [theta_Dorth,rho_Dorth] = ...
         estimation_Dorth_MV(x_donnees_bruitees,y_donnees_bruitees,n_tests)

[x_G, y_G, x_donnees_bruitees_centrees, y_donnees_bruitees_centrees] = ...
                centrage_des_donnees(x_donnees_bruitees,y_donnees_bruitees)

 theta = pi*rand(n_tests,1);
 somme = sum((x_donnees_bruitees_centrees.*cos(theta)+y_donnees_bruitees_centrees.*sin(theta)).^2,2);
 [~,k] = min(somme);
 theta_Dorth= theta(k);
 rho_Dorth=x_G*cos(theta(k)) + y_G*sin(theta(k));
end

% Fonction estimation_Dorth_MC (exercice_4.m) -----------------------------
function [theta_Dorth,rho_Dorth] = ...
                 estimation_Dorth_MC(x_donnees_bruitees,y_donnees_bruitees)

[x_G, y_G, x_donnees_bruitees_centrees, y_donnees_bruitees_centrees] = ...
               centrage_des_donnees(x_donnees_bruitees,y_donnees_bruitees)


    C = transpose([x_donnees_bruitees_centrees;y_donnees_bruitees_centrees]);
    
    % Calcul de lambda
    [Vecteurs_p,Valeurs_p] = eig(transpose(C)*C);
    % on construit une matrice diagonale contenant les aleurs propres sur
    % la diagonale:
    Valeurs_p = diag(Valeurs_p);
    [~,k] = min(Valeurs_p);
    cos_theta = Vecteurs_p(k,1);
    sin_theta = Vecteurs_p(k,2);
    rho_Dorth = x_G*cos_theta + y_G*sin_theta;
    
    if cos_theta>0
        if sin_theta>0
            theta_Dorth = acos(cos_theta);
        else
            theta_Dorth = asin(sin_theta);
        end
    else
        if sin_theta>0
            theta_Dorth= acos(cos_theta);
        else
            theta_Dorth = -acos(cos_theta);
        end
    end
    end





