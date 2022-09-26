
% TP1 de Statistiques : fonctions a completer et rendre sur Moodle
% Nom :GRETHEN
% Prénom : Clémentine
% Groupe : 1SN-N

function varargout = fonctions_TP1_stat(varargin)

    [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});

end

% Fonction G_et_R_moyen (exercice_1.m) ----------------------------------
function [G, R_moyen, distances] = ...
         G_et_R_moyen(x_donnees_bruitees,y_donnees_bruitees)
G_x=mean(x_donnees_bruitees);
G_y=mean(y_donnees_bruitees);
G=[G_x G_y]
R_x = x_donnees_bruitees-G_x;
R_y = y_donnees_bruitees-G_y;
R_moyen = mean(sqrt(R_y.^2 + R_x.^2));
distances=sqrt(R_y.^2 + R_x.^2);
end

% Fonction estimation_C_uniforme (exercice_1.m) ---------------------------
function [C_estime, R_moyen] = ...
         estimation_C_uniforme(x_donnees_bruitees,y_donnees_bruitees,n_tests)
 
    [G,R_moyen,~]= G_et_R_moyen(x_donnees_bruitees,y_donnees_bruitees)
    G_x = mean(x_donnees_bruitees);
    G_y = mean(y_donnees_bruitees);
    X=R_moyen*(rand(n_tests,1)-0.5);
    X = X+G_x;
    Y = R_moyen*(rand(n_tests,1)-0.5);
    Y = Y+G_y;
    [~,N] = size(x_donnees_bruitees);
    X_c = repmat(X,1,N);
    X_db = repmat(x_donnees_bruitees,n_tests,1);
    Y_c = repmat(Y,1,N);
    Y_db = repmat(y_donnees_bruitees,n_tests,1);
    
    distances = (X_c-X_db).^2 + (Y_c-Y_db).^2;
    
    calcul_dis = (distances - R_moyen).^2;
    C = sum(transpose(calcul_dis));  
    [~,i] = min(C);
    C_estime = [X(i),Y(i)];
     

end

% Fonction estimation_C_et_R_uniforme (exercice_2.m) ----------------------
function [C_estime, R_estime] = ...
         estimation_C_et_R_uniforme(x_donnees_bruitees,y_donnees_bruitees,n_tests)

    [G,R_moyen,~]= G_et_R_moyen(x_donnees_bruitees,y_donnees_bruitees);
    G_x = mean(x_donnees_bruitees);
    G_y = mean(y_donnees_bruitees);
    X = rand(n_tests,1)+G_x;
    Y = rand(n_tests,1)+G_y;
    
    % Tirage aléatoire de n_test  de moyenne R_moyen
    R = rand(n_tests,1)+R_moyen;

    [~,N] = size(x_donnees_bruitees);
    X_c = repmat(X,1,N);
    X_db = repmat(x_donnees_bruitees,n_tests,1);
    Y_c = repmat(Y,1,N);
    Y_db = repmat(y_donnees_bruitees,n_tests,1);
    distances = (X_c-X_db).^2 + (Y_c-Y_db).^2;
    calcul = zeros(n_tests,N);

    for k = 1:N
        d = (distances(:,k) - R(k)).^2;
        calcul_dis(:,k) = d;
    end
  s = sum(transpose(calcul));   
    [~,z] = min(s);
    C_estime = [X(z),Y(z)];
    R_estime = R(z);
    
    
    


end

% Fonction occultation_donnees (donnees_occultees.m) ----------------------
function [x_donnees_bruitees, y_donnees_bruitees] = ...
         occultation_donnees(x_donnees_bruitees, y_donnees_bruitees, theta_donnees_bruitees)

theta_1=2*pi*rand(1,1);
theta_2=2*pi*rand(1,1);

if theta_1 <0 theta_2
    indice = find(theta_donnees_bruitees > theta_1 & theta_donnees_bruitees < theta_2);
else
    indice = find((theta_donnees_bruitees < theta_2) | (theta_donnees_bruitees > theta_1));
end
x_donnees_bruitees = x_donnees_bruitees(indice);
y_donnees_bruitees = y_donnees_bruitees(indice);

end

% Fonction estimation_C_et_R_normale (exercice_4.m, bonus) ----------------
function [C_estime, R_estime] = ...
         estimation_C_et_R_normale(x_donnees_bruitees,y_donnees_bruitees,n_tests)
 [G,R_moyen,~]= G_et_R_moyen(x_donnees_bruitees,y_donnees_bruitees);
    G_x = mean(x_donnees_bruitees);
    G_y = mean(y_donnees_bruitees);
    X = randn(n_tests,1)+G_x;
    Y = randn(n_tests,1)+G_y;
    
    % Tirage aléatoire de n_test  de moyenne R_moyen
    R = rand(n_tests,1)+R_moyen;

    [~,N] = size(x_donnees_bruitees);
    X_c = repmat(X,1,N);
    X_db = repmat(x_donnees_bruitees,n_tests,1);
    Y_c = repmat(Y,1,N);
    Y_db = repmat(y_donnees_bruitees,n_tests,1);
    distances = (X_c-X_db).^2 + (Y_c-Y_db).^2;
    calcul = zeros(n_tests,N);

    for k = 1:N
        d = (distances(:,k) - R(k)).^2;
        calcul_dis(:,k) = d;
    end
  s = sum(transpose(calcul));   
    [~,z] = min(s);
    C_estime = [X(z),Y(z)];
    R_estime = R(z);

end
