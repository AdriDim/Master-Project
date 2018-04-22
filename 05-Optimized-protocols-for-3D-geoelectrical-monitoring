
%% --------- Optimized-protocols-for-3D-geoelectrical-monitoring ----------
%
% Adrien Dimech - Master Project - 22/04/2018
%
% -------------------------------------------------------------------------
% Matlab codes to prepare - process - visualize and interpret 3D time-lapse
% geolelectrical monitoring of a waste rock pile.
% -------------------------------------------------------------------------
%
% This Matlab code was used to select optimized protocols for the
% geoelectrical monitoring of the experimental waste rock pile. The
% methodology follows the example of the 'COMPARE R' method developed by
% the British Geological Survey (Stummer et al., 2004; Wilkinson et al.,
% 2006b, 2012; Loke et al.,2014a, 2014b, 2015). Optimized protocols with
% 1000 configurations provide better coverage and senstivities than
% standard protocols of 4000 configurations. This Matlab code can be used
% for any 3D electrode distribution and the number of configurations can be
% selected by the user. 3D visualizations of the optimized protocol
% sensitivies are used to assess the quality of the selected protocols.
%
% Feel free to visit : https://www.researchgate.net/profile/Adrien_Dimech
% for more information about my research or contact me for more information
% and data files : adrien.dimech@gmail.com
%
%%% VERSION SIMPLIFIÉE DE LA METHODE 'COMPARE R' : IDENTIFICATION DES PROTOCOLES DE MESURE

%% I) CHARGEMENT DES DONNEES
actif=1;
for o=1
    if actif==1
        close all
        clear all
        figure
        p = get(gcf,'Position');
        set(0,'DefaultFigurePosition',p);
        close all
        load COORD_ELEC_SORT.mat
        load Maillage3DSurface.mat
        
        %         load COORD_ELEC_SORT_2.mat
        %         load COORD_ELEC_SORT_TEST.mat
        
        ELEC_FAUSSES=[7;16;27;41;46;160;191];
        ELEC_MAUVAISES=[4;10;12;14;20;22;37;53;62;71];
        WRONG_ELEC=[ELEC_FAUSSES;ELEC_MAUVAISES];
        
        % Affichage des électrodes
        ELEC=COORD_ELEC_SORT;
        figure('Color', [ 1 1 1])
        plot3(ELEC(:,1),ELEC(:,2),ELEC(:,3),'.k','MarkerSize',20)
        hold on
        for i=1:length(ELEC_FAUSSES(:,1))
            plot3(ELEC(ELEC_FAUSSES(i,1),1),ELEC(ELEC_FAUSSES(i,1),2),ELEC(ELEC_FAUSSES(i,1),3),'.r','MarkerSize',20)
        end
        for i=1:length(ELEC_MAUVAISES(:,1))
            plot3(ELEC(ELEC_MAUVAISES(i,1),1),ELEC(ELEC_MAUVAISES(i,1),2),ELEC(ELEC_MAUVAISES(i,1),3),'.m','MarkerSize',20)
        end
        plot3(Maillage3DSurface(:,1),Maillage3DSurface(:,2),Maillage3DSurface(:,3),'.-k','LineWidth',1)
        daspect([1 1 1])
    end
end

%% II) CREATION DU COMPREHENSIVE DATA SET
%--------------------------------------------------------------------------
actif_comprehensive_data_set=0;
for o=1
    if actif_comprehensive_data_set==1
        %% 2.1) Tous les alpha et beta arrays possibles pour une ligne de 24 electrodes
        % Creation le 10/04/2017
        % Modification le 11/04/2017
        for o=1
            % Configurations alpha          C1 -- P1 -- P2 -- C2
            clear ALPHA_ARRAYS
            ALPHA_ARRAYS(1,1:4)=[0,0,0,0];
            for C1=1:24
                for P1=C1+1:24
                    for P2=P1+1:24
                        for C2=P2+1:24
                            ALPHA_ARRAYS(end+1,1:4)=[C1,C2,P1,P2];
                        end
                    end
                end
            end
            ALPHA_ARRAYS(1,:)=[];
            
            % Configurations beta           C2 -- C1 -- P1 -- P2
            clear BETA_ARRAYS
            BETA_ARRAYS(1,1:4)=[0,0,0,0];
            for C2=1:24
                for C1=C2+1:24
                    for P1=C1+1:24
                        for P2=P1+1:24
                            BETA_ARRAYS(end+1,1:4)=[C1,C2,P1,P2];
                        end
                    end
                end
            end
            BETA_ARRAYS(1,:)=[];
            
            % Comprehensive data set avec 24 électrodes
            clear COMPREHENSIVE_24_ELEC
            for o=1:1
                COMPREHENSIVE_24_ELEC=[ALPHA_ARRAYS;BETA_ARRAYS];
                % Transformation de COMPREHENSIVE_24_ELEC pour etre coherent avec la
                % numerotation des electrodes sur le terrain
                for i=1:4
                    COMPREHENSIVE_24_ELEC(:,i)=(COMPREHENSIVE_24_ELEC(:,i)-1)*4+1;
                end
            end
            
            % Liberation de la memoire
            clear BETA_ARRAYS
            clear ALPHA_ARRAYS
        end
        
        %% 2.2.1) Extension a deux lignes (1 x base et 1 x surface) avec roll along vertical
        % Creation le 10/04/2017
        % Modification le 11/04/2017
        for o=1
            % Matrice de passage Roll_Along en Z
            clear ROLL_ALONG_Z
            ROLL_ALONG_Z=[1,1,1,1;1,2,1,1;1,2,1,2;1,2,2,2;2,2,2,2;1,1,2,2;2,2,1,1;2,1,1,1;2,1,2,1;2,1,2,2];
            
            % Comprehensive data set avec 48x2 électrodes (4 secondes)
            clear COMPREHENSIVE_48_ELEC
            COMPREHENSIVE_48_ELEC_AVANCEMENT=0;tic
            L_Z=length(ROLL_ALONG_Z(:,1));
            COMPREHENSIVE_48_ELEC=zeros(length(COMPREHENSIVE_24_ELEC(:,1))*L_Z,4);
            for i=1:length(COMPREHENSIVE_24_ELEC(:,1))
                clear Mat_24
                for j=1:4
                    Mat_24(:,j)=COMPREHENSIVE_24_ELEC(i,j)*ones(L_Z,1);
                end
                COMPREHENSIVE_48_ELEC((i-1)*L_Z+1:i*L_Z,1:4)=Mat_24+(ROLL_ALONG_Z-1)*96;
                COMPREHENSIVE_48_ELEC_AVANCEMENT=round(i/length(COMPREHENSIVE_24_ELEC(:,1))*100)
            end
            toc
            
            % Ajout des configurations dipole dipole equatorial
            % Configurations verticale
            %   C2 -- P2
            %   C1 -- P1
            clear DIP_DIP_EQUA_VERT
            DIP_DIP_EQUA_VERT(1,1:4)=[0,0,0,0];
            for C1=1:24
                for P1=C1+1:24
                    DIP_DIP_EQUA_VERT(end+1,1:4)=[C1,C1+96,P1,P1+96];
                end
            end
            DIP_DIP_EQUA_VERT(1,:)=[];
            
            % Configurations horizontale
            %   P1 -- P2
            %   C1 -- C2
            clear DIP_DIP_EQUA_HOR
            DIP_DIP_EQUA_HOR(1,1:4)=[0,0,0,0];
            for C1=1:24
                for C2=C1+1:24
                    DIP_DIP_EQUA_HOR(end+1,1:4)=[C1,C2,C1+96,C2+96];
                end
            end
            DIP_DIP_EQUA_HOR(1,:)=[];
            
            % Ajout des configurations Dipole Dipole Equatorial
            COMPREHENSIVE_48_ELEC=[COMPREHENSIVE_48_ELEC;DIP_DIP_EQUA_VERT;DIP_DIP_EQUA_HOR];
            
            % Liberation de la memoire
            clear DIP_DIP_EQUA_HOR
            clear DIP_DIP_EQUA_VERT
        end
        
        %% 2.2.2) Extension a huit lignes (2 x base et 2 x surface) avec roll along horizontal
        % Création le 11/04/2017
        for o=1
            % Matrice de passage Roll_Along en Y pour Nb_Lignes lignes
            clear ROLL_ALONG_Y
            Nb_Lignes=2;
            ROLL_ALONG_Y(1,1:4)=[0,0,0,0];
            for i=1:Nb_Lignes
                for j=1:Nb_Lignes
                    for k=1:Nb_Lignes
                        for l=1:Nb_Lignes
                            ROLL_ALONG_Y(end+1,1:4)=[i,j,k,l];
                        end
                    end
                end
            end
            ROLL_ALONG_Y(1,:)=[];
            
            % Comprehensive data set avec 48x2 électrodes (400 secondes)
            clear COMPREHENSIVE_96_ELEC
            COMPREHENSIVE_96_ELEC_AVANCEMENT=0;tic
            COMPREHENSIVE_96_ELEC=zeros(length(COMPREHENSIVE_48_ELEC(:,1))*length(ROLL_ALONG_Y(:,1)),4);
            L_Y=length(ROLL_ALONG_Y(:,1));
            for i=1:length(COMPREHENSIVE_48_ELEC(:,1))
                clear Mat_48
                for j=1:4
                    Mat_48(:,j)=COMPREHENSIVE_48_ELEC(i,j)*ones(L_Y,1);
                end
                COMPREHENSIVE_96_ELEC((i-1)*L_Y+1:i*L_Y,1:4)=Mat_48+(ROLL_ALONG_Y-1);
                COMPREHENSIVE_96_ELEC_AVANCEMENT=round(i/length(COMPREHENSIVE_48_ELEC(:,1))*100)
            end
            toc
            
            % Liberation de la memoire
            clear COMPREHENSIVE_24_ELEC
            clear COMPREHENSIVE_48_ELEC
            clear COMPREHENSIVE_48_ELEC_AVANCEMENT
            clear COMPREHENSIVE_96_ELEC_AVANCEMENT
            
            % Extension de COMPREHENSIVE_96_ELEC à toutes les lignes (192 elec)
            COMPREHENSIVE_192_ELEC=[COMPREHENSIVE_96_ELEC;COMPREHENSIVE_96_ELEC+1;COMPREHENSIVE_96_ELEC+2];
            COMPREHENSIVE_192_ELEC=unique(COMPREHENSIVE_192_ELEC,'rows','stable');
            
            clear COMPREHENSIVE_96_ELEC
        end
        
        %% 2.2.3) Ajout des configurations de surface entre les 4 lignes du haut puis du bas
        % Création le 18/04/2017
        % Modification le 22/04/2017
        for o=1
            % Alpha arrays entre les 4 lignes
            ALPHA_ARRAYS_4_LIGNES=0;
            for C1=1:4:13
                for P1=2:4:14
                    for P2=3:4:15
                        for C2=4:4:16
                            ALPHA_ARRAYS_4_LIGNES(end+1,1:4)=[C1,C2,P1,P2];
                        end
                    end
                end
            end
            ALPHA_ARRAYS_4_LIGNES(1,:)=[];
            
            % Beta arrays entre les 4 lignes
            BETA_ARRAYS_4_LIGNES=0;
            for C2=1:4:13
                for C1=2:4:14
                    for P1=3:4:15
                        for P2=4:4:16
                            BETA_ARRAYS_4_LIGNES(end+1,1:4)=[C1,C2,P1,P2];
                        end
                    end
                end
            end
            BETA_ARRAYS_4_LIGNES(1,:)=[];
            
            % Propagation aux 96 électrodes de la base de la halde
            ALPHA_ARRAYS_4_LIGNES_96=ALPHA_ARRAYS_4_LIGNES;
            BETA_ARRAYS_4_LIGNES_96=BETA_ARRAYS_4_LIGNES;
            for i=1:20
                ALPHA_ARRAYS_4_LIGNES_96=[ALPHA_ARRAYS_4_LIGNES_96;ALPHA_ARRAYS_4_LIGNES+i*4];
                BETA_ARRAYS_4_LIGNES_96=[BETA_ARRAYS_4_LIGNES_96;BETA_ARRAYS_4_LIGNES+i*4];
            end
            
            % Propagation aux 192 électrodes de la halde
            ALPHA_ARRAYS_4_LIGNES_192=[ALPHA_ARRAYS_4_LIGNES_96;ALPHA_ARRAYS_4_LIGNES_96+96];
            BETA_ARRAYS_4_LIGNES_192=[BETA_ARRAYS_4_LIGNES_96;BETA_ARRAYS_4_LIGNES_96+96];
            ALPHA_ARRAYS_4_LIGNES_192=unique(ALPHA_ARRAYS_4_LIGNES_192,'rows','stable');
            BETA_ARRAYS_4_LIGNES_192=unique(BETA_ARRAYS_4_LIGNES_192,'rows','stable');
            
            % Ajout des configurations de 4 lignes au COMPREHENSIVE_192_ELEC
            PRO_4_LIGNES=unique([ALPHA_ARRAYS_4_LIGNES_192;BETA_ARRAYS_4_LIGNES_192],'rows','stable');
            % ATTENTION A UNIQUE QUI TRIE LES CONFIGURATIONS!!!!
            COMPREHENSIVE_192_ELEC=unique([COMPREHENSIVE_192_ELEC;PRO_4_LIGNES],'rows','stable');
            
            % Liberation de la memoire
            clear ALPHA_ARRAYS_4_LIGNES
            clear BETA_ARRAYS_4_LIGNES
            clear ALPHA_ARRAYS_4_LIGNES_96
            clear BETA_ARRAYS_4_LIGNES_96
            clear ALPHA_ARRAYS_4_LIGNES_192
            clear BETA_ARRAYS_4_LIGNES_192
        end
        
        %% 2.2.4) Ajout des configurations standard pour les 192 electrodes (PROA, PROB, PROC)
        % Création le 21/04/2017
        % Modification le 22/04/2017
        for o=1
            % PROA
            load PROA.mat
            % PROB
            load PROB.mat
            % PROC
            load PROC.mat
            
            % Ajout des configurations standard au COMPREHENSIVE_192_ELEC
            %COMPREHENSIVE_192_ELEC=[PROA;PROB;PROC];
            PRO_STAND=[PROA;PROB;PROC];
            PRO_STAND=unique(PRO_STAND,'rows','stable');
            % ATTENTION A UNIQUE QUI TRIE LES CONFIGURATIONS!!!!
            COMPREHENSIVE_192_ELEC=unique([COMPREHENSIVE_192_ELEC;PRO_STAND],'rows','stable');
            
            
            
            % Liberation de la memoire
            clear ALPHA_ARRAYS_4_LIGNES
            clear BETA_ARRAYS_4_LIGNES
            clear ALPHA_ARRAYS_4_LIGNES_96
            clear BETA_ARRAYS_4_LIGNES_96
            clear ALPHA_ARRAYS_4_LIGNES_192
            clear BETA_ARRAYS_4_LIGNES_192
        end
        
        %% 2.3) Filtrage avec facteur géométrique
        % Création le 11/04/2017
        % Modification le 12/04/2017
        K_max=2000;
        for o=1
            % Chargement des coordonnées des électrodes
            clear ELEC
            ELEC=COORD_ELEC_SORT;
            N_Elec=length(ELEC(:,1));
            
            % Calcul des coordonnées des électrodes "miroir" au dessus de la surface
            % Calcul de la surface
            pente=(5.665-7.563)/(91.07-56.74); % cf. FIGURE INTITULÉE PENTE_SURFACE
            for i=1:192
                ELEC(192+i,1:2)=ELEC(i,1:2);
                ELEC(192+i,3)=ELEC(i,3)+2*(7.563+pente*(ELEC(i,1)-56.74)-ELEC(i,3));
            end
            
            % Visualisation de la surface, des électrodes et des électrodes "miroirs"
            for o=1:1
                figure('Color', [ 1 1 1])
                N_Elec=length(ELEC(:,1));
                plot3(ELEC(1:N_Elec/2,1),ELEC(1:N_Elec/2,2),ELEC(1:N_Elec/2,3),'.k','MarkerSize',20)
                hold on
                X_Pente=(30:1:100)';
                Z_Pente=7.563+pente*(X_Pente-56.74);
                for i=21:31
                    plot3(X_Pente(:,1),i*ones(length(X_Pente(:,1)),1),Z_Pente(:,1),'.r','MarkerSize',5)
                end
                plot3(ELEC(N_Elec/2+1:end,1),ELEC(N_Elec/2+1:end,2),ELEC(N_Elec/2+1:end,3),'.g','MarkerSize',20)
                xlim([30 95])
                ylim([21 31])
                daspect([1,1,1])
            end
            
            % Calcul de la matrice des distances
            N_Elec=length(ELEC(:,1));
            for i=1:N_Elec
                for j=1:N_Elec
                    DISTANCE(i,j)=pdist([ELEC(i,1),ELEC(i,2),ELEC(i,3);ELEC(j,1),ELEC(j,2),ELEC(j,3)],'euclidean');
                end
            end
            
            % Assemblage pour calcul des distances
            clear K_Mat_1;clear K_Mat_2;clear K_Mat_3;
            clear K_Mat_4; % matrice de distance
            K_Mat_1=[1,3;1,4;2,3;2,4;1,3;1,4;2,3;2,4];
            K_Mat_2=[0;0;0;0;1;1;0;0];
            K_Mat_3=[1;-1;-1;1;1;-1;-1;1];
            
            % Calcul des facteurs géométriques
            clear H; clear E1; clear E2;
            clear COORD_E1; clear COORD_E2
            H=0;
            for i=1:8
                E1=COMPREHENSIVE_192_ELEC(:,K_Mat_1(i,1))+192*K_Mat_2(i,1);
                E2=COMPREHENSIVE_192_ELEC(:,K_Mat_1(i,2))+192*K_Mat_2(i,1);
                COORD_E1(:,1:3)=[ELEC(E1,1),ELEC(E1,2),ELEC(E1,3)];
                COORD_E2(:,1:3)=[ELEC(E2,1),ELEC(E2,2),ELEC(E2,3)];
                clear somme
                somme=0;
                for j=1:3
                    somme=somme+(COORD_E1(:,j)-COORD_E2(:,j)).^2;
                end
                K_Mat_4(:,i)=sqrt(somme(:,1));
                H=H+K_Mat_3(i,1)*1./K_Mat_4(:,i);
            end
            K=4*pi./H;
            
            % Elimination des configurations avec un K trop important
            %K_max=10000; % Loke et al. 2014 (meme dimensions que mon modele)
            clear elimination_1
            elimination_1=find(abs(K(:,1))>K_max);
            COMPREHENSIVE_192_ELEC(:,5)=K;
            COMPREHENSIVE_192_ELEC_filtK=COMPREHENSIVE_192_ELEC;
            %clear COMPREHENSIVE_192_ELEC
            COMPREHENSIVE_192_ELEC_filtK(elimination_1,:)=[];
            
            H_filtK=H;
            H_filtK(elimination_1,:)=[];
            K_Mat_4_filtK=K_Mat_4;
            K_Mat_4_filtK(elimination_1,:)=[];
            COORD_E1(elimination_1,:)=[];
            COORD_E2(elimination_1,:)=[];
        end
        
        %% 2.4) Filtrage des configurations abberantes
        % Création le 15/04/2017
        for o=1
            % Suppression des NaN
            clear elimination_3
            elimination_3=find(COMPREHENSIVE_192_ELEC_filtK(:,5)==0);
            COMPREHENSIVE_192_ELEC_filtK(elimination_3,:)=[];
            H_filtK(elimination_3,:)=[];
            K_Mat_4_filtK(elimination_3,:)=[];
            COORD_E1(elimination_3,:)=[];
            COORD_E2(elimination_3,:)=[];
        end
        
        %% 2.5) Filtrage avec sensibilité du facteur géométrique aux erreurs de positions d'électrodes
        % Création le 12/04/2017
        % Modification le 13/04/2017
        K_err_max=5;
        for o=1
            % Mise en forme des coordonnees des electrodes
            clear COORD_ELEC
            for i=1:4
                Ei=COMPREHENSIVE_192_ELEC_filtK(:,i);
                COORD_ELEC(:,(i-1)*3+1:i*3)=ELEC(COMPREHENSIVE_192_ELEC_filtK(:,i),1:3);
            end
            
            % Calcul de la sensibilite de K (S_2) en fonction des erreurs de position
            clear K_Mat_5;
            K_Mat_5=[1,2;3,4;1,3;2,4];
            clear dK_dEi_2; % erreur de position au carré (col 1 = erreur % C1 ... )
            clear S_2 % sensibilité au carré
            S_2=0;
            for i=1:4 % erreur de position de C1, C2, P1, P2
                somme_1=0;
                for j=1:3 % pour chaque dimension (x, y et z)
                    somme_2=0;
                    for k=1:2
                        ind_Ci_Pi=K_Mat_5(i,k);
                        ind_E=K_Mat_1(ind_Ci_Pi,1:2);
                        col_COORD_ELEC=(ind_E-1)*3+j;
                        somme_2=somme_2+-1./((K_Mat_4_filtK(:,ind_Ci_Pi)).^2).*(-2*(COORD_ELEC(:,col_COORD_ELEC(1,2))-COORD_ELEC(:,col_COORD_ELEC(1,1))));
                    end
                    somme_1=somme_1+(-2*pi./(H_filtK.^2).*somme_2).^2;
                end
                dK_dEi_2(:,i)=somme_1;
                S_2=S_2+dK_dEi_2(:,i);
            end
            
            % Calcul de Re le Geometric Factor Relative Error
            clear Re
            Re = sqrt(S_2)./COMPREHENSIVE_192_ELEC_filtK(:,5);
            COMPREHENSIVE_192_ELEC_filtK(:,6)=Re;
            
            % Elimination des configurations avec un Re trop important
            % K_err_max=0.55; % Loke et al. 2014 (meme dimensions que mon modele)
            COMPREHENSIVE_192_ELEC_filtK_filtRe=COMPREHENSIVE_192_ELEC_filtK;
            %clear COMPREHENSIVE_96_ELEC_filtK
            clear elimination_2
            elimination_2=find(abs(Re(:,1))>K_err_max);
            COMPREHENSIVE_192_ELEC_filtK_filtRe(elimination_2,:)=[];
        end
        
        %% 2.5.1) Dépendance entre Re et K_max_observe
        % Création le 13/04/2017
        actif=0;
        for o=1
            if actif==1
                K_err_max_test=(0:1:50)';
                for i=1:length(K_err_max_test(:,1))
                    % Elimination des configurations avec un Re trop important
                    K_err_max=K_err_max_test(i,1);
                    COMPREHENSIVE_192_ELEC_filtK_filtRe=COMPREHENSIVE_192_ELEC_filtK_filtRe;
                    %clear COMPREHENSIVE_192_ELEC_filtK_filtRe
                    clear elimination_2
                    elimination_2=find(abs(Re(:,1))>K_err_max);
                    COMPREHENSIVE_192_ELEC_filtK_filtRe(elimination_2,:)=[];
                    
                    K_maxi_observe=max(abs(COMPREHENSIVE_192_ELEC_filtK_filtRe(:,5)));
                    K_err_max_test(i,2)=K_maxi_observe;
                end
                figure('Color', [ 1 1 1])
                plot(K_err_max_test(:,1),K_err_max_test(:,2),'-o')
                %ylim([-10 2500])
                grid
                figure('Color', [ 1 1 1])
                plot(K_err_max_test(:,1),K_err_max_test(:,2)./K_err_max_test(:,1),'-o')
                %ylim([-10 2500])
                grid
            end
        end
        
        %% 2.6) Visualisation des résultats
        % Création le 13/04/2017
        actif=0;
        for o=1
            if actif==1
                figure('Color', [ 1 1 1])
                plot(sqrt(S_2),'.','MarkerSize',1)
                yyaxis right
                plot(COMPREHENSIVE_192_ELEC_filtK(:,5),'.','MarkerSize',1)
                figure('Color', [ 1 1 1])
                plot(COMPREHENSIVE_192_ELEC_filtK(:,6),'.','MarkerSize',1)
                
                figure('Color', [ 1 1 1])
                plot(COMPREHENSIVE_192_ELEC_filtK_filtRe(:,6),'--.','MarkerSize',1)
            end
        end
        
    else
        load COMPREHENSIVE_192_ELEC_filtK_filtRe
    end
end

%% III) APPLICATION DE LA COMPARE R METHOD A UN PROTOCOLE
%--------------------------------------------------------------------------
actif_compare_r=1;
for o=1
    if actif_compare_r==1
        %% 3.1) Maillage de la halde pour identification des meilleures configurations
        % Creation le 02/04/2017
        % Modification le 15/04/2017
        actif=1;
        for o=1
            if actif==1
                dX = 1; %résolution spatiale en X
                dY = 1; %résolution spatiale en Y
                dZ = 0.5; %résolution spatiale en Z
                [X,Y,Z] = meshgrid(30:dX:95,21:dY:31,2:dZ:9);
                
                % Création points
                COORD=0;
                for i=1:length(X(1,:,1))
                    for j=1:length(X(:,1,1))
                        for k=1:length(X(1,1,:))
                            COORD(end+1,1)=X(1,i,1);
                            COORD(end,2)=Y(j,1,1);
                            COORD(end,3)=Z(1,1,k);
                        end
                    end
                end
                COORD(1,:)=[];
                
                pente=(5.665-7.563)/(91.07-56.74); % cf. FIGURE INTITULÉE PENTE_SURFACE
                
                
                COORD_2=0;
                for i=1:length(COORD(:,1))
                    % if COORD(i,3)<8.725-(COORD(i,1)-30)*0.05
                    if COORD(i,3)<7.563+pente*(COORD(i,1)-56.74)
                        COORD_2(end+1,1:3)=COORD(i,1:3);
                    end
                end
                
                figure('Color', [ 1 1 1])
                plot3(COORD_ELEC_SORT(:,1),COORD_ELEC_SORT(:,2),COORD_ELEC_SORT(:,3),'.k','MarkerSize',20)
                hold on
                plot3(COORD_2(:,1),COORD_2(:,2),COORD_2(:,3),'.r','MarkerSize',5)
                %                 xlim([30 95])
                %                 ylim([21 31])
                %                 zlim([1 9])
                plot3(Maillage3DSurface(:,1),Maillage3DSurface(:,2),Maillage3DSurface(:,3),'.-k')
                
                daspect([1,1,1])
            end
        end
        
        %% 3.2) Préparation du calcul de la sensibilité des configurations
        % Creation le 02/04/2017
        % Modification le 15/04/2017
        
        % Calcul des sensibilités pour chaque electrode
        ELEC=COORD_ELEC_SORT;
        MAILLAGE=COORD_2;
        for o=1
            % Matrice des duos possibles
            DUOS=zeros(192*191,2);
            count=0;
            tic
            for elec1=1:192
                for elec2=1:192
                    if elec1~=elec2
                        count=count+1;
                        %DUOS_AVANCEMENT=round(count/(192*191)*100)
                        DUOS(count,1:2)=[elec1,elec2];
                    end
                end
            end
            toc
            
            % (X-Xi)(X-Xj)+(Y-Yi)(Y-Yj)+(Z-Zi)(Z-Zj)
            clear Mat1
            Mat1=zeros(length(MAILLAGE(:,1)),length(DUOS(:,1)));
            tic
            for o=1:1
                for j=1:length(DUOS(:,1))
                    clear somme
                    somme=0;
                    for colonne=1:3
                        somme=somme+(MAILLAGE(:,colonne)-ELEC(DUOS(j,1),colonne)).*(MAILLAGE(:,colonne)-ELEC(DUOS(j,2),colonne));
                    end
                    Mat1(:,j)=somme;
                    %MAT1_AVANCEMENT=round(j/length(DUOS(:,1))*100)
                end
            end
            toc
            
            % ((X-Xi)^2+(Y-Yi)^2+(Z-Zi)^2)^1.5
            clear Mat2
            tic
            for j=1:length(ELEC(:,1))
                somme=0;
                for colonne=1:3
                    somme=somme+(MAILLAGE(:,colonne)-ELEC(j,colonne)).^2;
                end
                Mat2(:,j)=somme.^1.5;
            end
            toc
            
            % Matrice de sensibilité
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % REMARQUE OLA + ABDERREZAK : utiliser COMSOL avec modèle résistivités pour
            % prendre en compte hétérogénéité du modèle de la halde dans calcul de sensibilités
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tic
            for o=1:1
                clear SENSIBILITE
                SENSIBILITE=zeros(length(MAILLAGE(:,1)),length(DUOS(:,1)));
                for j=1:length(DUOS(:,1))
                    SENSIBILITE(:,j)=1/(4*pi^2)*Mat1(:,j)./(Mat2(:,DUOS(j,1)).*Mat2(:,DUOS(j,2)));
                    %SENSIBILITE_AVANCEMENT=round(j/length(DUOS(:,1))*100)
                end
            end
            toc
        end
        
        %% 3.3) Calcul des sensibilités pour chaque configuration
        % Creation le 02/04/2017
        % Modification le 15/04/2017
        for o=1
            % Sensibilités pour les quadripoles
            tic
            
            % Réecriture de COMPREHENSIVE_192_ELEC_filtK_filtRe
            COMPREHENSIVE_192_ELEC_filtK_filtRe_filtWRONG=COMPREHENSIVE_192_ELEC_filtK_filtRe;
            elimination_4=zeros(length(COMPREHENSIVE_192_ELEC_filtK_filtRe_filtWRONG(:,1)),1);
            count=1;
            for i=1:length(COMPREHENSIVE_192_ELEC_filtK_filtRe_filtWRONG(:,1))
                elecOK=1;
                for j=1:4
                    if isempty(find(COMPREHENSIVE_192_ELEC_filtK_filtRe_filtWRONG(i,j)==WRONG_ELEC(:,1)))==0
                        elecOK=0;
                    end
                end
                
                if elecOK==0
                    elimination_4(count,1)=i;
                    count=count+1;
                end
            end
            elimination_4(count:end,:)=[];
            COMPREHENSIVE_192_ELEC_filtK_filtRe_filtWRONG(elimination_4,:)=[];
            
            CONFIGURATIONS=COMPREHENSIVE_192_ELEC_filtK_filtRe_filtWRONG;
            %           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             SYNTHESE_VALUES=zeros(length(MAILLAGE(:,1)),5);
            %             SYNTHESE_INDICE=zeros(length(MAILLAGE(:,1)),5);
            %           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load SYNTHESE_VALUES.mat
            load SYNTHESE_INDICE.mat
            %                         SYNTHESE_VALUES=SYNTHESE_VALUES_655521;
            %                         SYNTHESE_INDICE=SYNTHESE_INDICE_655521;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=length(CONFIGURATIONS(:,1)):length(CONFIGURATIONS(:,1))
                clear SENSIBILITE_TOT
                C1_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS(i,1)&DUOS(:,2)==CONFIGURATIONS(i,3)));
                C2_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS(i,2)&DUOS(:,2)==CONFIGURATIONS(i,3)));
                C1_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS(i,1)&DUOS(:,2)==CONFIGURATIONS(i,4)));
                C2_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS(i,2)&DUOS(:,2)==CONFIGURATIONS(i,4)));
                SENSIBILITE_TOT(:,1)=abs(C1_P1-C2_P1-C1_P2+C2_P2);
                SENS_TOT_AVANCEMENT=round(i/length(CONFIGURATIONS(:,1))*1000000)/10000
                
                for j=1:length(MAILLAGE(:,1))
                    trouve=0;
                    for k=1:5
                        if trouve==0
                            if SENSIBILITE_TOT(j,1)>SYNTHESE_VALUES(j,k)
                                SYNTHESE_VALUES(j,k+1:end)=SYNTHESE_VALUES(j,k:end-1);
                                SYNTHESE_INDICE(j,k+1:end)=SYNTHESE_INDICE(j,k:end-1);
                                
                                SYNTHESE_VALUES(j,k)=SENSIBILITE_TOT(j,1);
                                SYNTHESE_INDICE(j,k)=i;
                                trouve=1;
                            end
                        end
                    end
                end
                
            end
            toc
        end
        
        %% 3.4) Identification des meilleures configurations
        % Création le 15/04/2017
        for o=1
            % Pondération de chaque cellule
            clear PONDERATION;
            for i=1:5
                %PONDERATION(:,i)=(SYNTHESE_VALUES(:,i)-SYNTHESE_VALUES(:,end))./(SYNTHESE_VALUES(:,1)-SYNTHESE_VALUES(:,end))*100*(6-i)*(6-i);
                %PONDERATION(:,i)=(SYNTHESE_VALUES(:,i)-SYNTHESE_VALUES(:,end))./(SYNTHESE_VALUES(:,1)-SYNTHESE_VALUES(:,end))*100*(6-i);
                %PONDERATION(:,i)=(SYNTHESE_VALUES(:,i)-SYNTHESE_VALUES(:,end))./(SYNTHESE_VALUES(:,1)-SYNTHESE_VALUES(:,end))*100;
                PONDERATION(:,i)=(SYNTHESE_VALUES(:,i)-SYNTHESE_VALUES(:,end))./(SYNTHESE_VALUES(:,1)-SYNTHESE_VALUES(:,end))*100+((6-i)*5);
                %PONDERATION(:,i)=zeros(length(MAILLAGE(:,1)),1)+((6-i)*5);
            end
            
            % Réécriture des matrices de pondération, de valeurs et d'indice
            % pour faciliter les calculs
            PONDERATION_NEW=PONDERATION(:,1);
            SYNTHESE_VALUES_NEW=SYNTHESE_VALUES(:,1);
            SYNTHESE_INDICE_NEW=SYNTHESE_INDICE(:,1);
            for i=2:5
                PONDERATION_NEW=[PONDERATION_NEW;PONDERATION(:,i)];
                SYNTHESE_VALUES_NEW=[SYNTHESE_VALUES_NEW;SYNTHESE_VALUES(:,i)];
                SYNTHESE_INDICE_NEW=[SYNTHESE_INDICE_NEW;SYNTHESE_INDICE(:,i)];
            end
            
            % Ranking de chaque configuration identifiée
            clear CONFIGURATION_UNIQUE
            CONFIGURATION_UNIQUE=unique(SYNTHESE_INDICE_NEW(:,1),'stable');
            CONFIGURATION_UNIQUE(:,2)=zeros(length(CONFIGURATION_UNIQUE(:,1)),1);
            for i=1:length(CONFIGURATION_UNIQUE(:,1))
                clear trouve
                trouve=find(SYNTHESE_INDICE_NEW(:,1)==CONFIGURATION_UNIQUE(i,1));
                for j=1:length(trouve(:,1))
                    CONFIGURATION_UNIQUE(i,2)=CONFIGURATION_UNIQUE(i,2)+PONDERATION_NEW(trouve(j,1),1);
                end
            end
            CONFIGURATION_UNIQUE=sortrows(CONFIGURATION_UNIQUE,-2);
            
            % Visualisation de la distribution de l'information
            figure('Color', [ 1 1 1])
            semilogy(CONFIGURATION_UNIQUE(:,end),'.k','MarkerSize',1)
            grid
            figure('Color', [ 1 1 1])
            plot(CONFIGURATION_UNIQUE(:,end),'.k','MarkerSize',1)
            grid
        end
        
    end
end

%% IV) CHOIX DU PROTOCOLE ET VISUALISATION DES SENSIBILITÉS
%--------------------------------------------------------------------------
actif=1;
for o=1
    if actif==1
        %% 4.1) Définition des protocoles de mesure
        % Creation le 16/04/2017
        % Modification le 20/04/2017
        count=0;
        auto=0;
        clear SYNTHESE_TEST
        Nb_CONFIG=1000;
        for o=1
            %for count1=1:0.25:4
            for count1=1
                %Nb_CONFIG=10^count1;
                %Nb_CONFIG=2000;
                actif=1;
                for o=1
                    if actif==1
                        % Protocole total du comprehensive data set
                        clear PROTOCOLE_TOTAL
                        
                        PROTOCOLE_TOTAL=CONFIGURATIONS(CONFIGURATION_UNIQUE(:,1),:);
                        PROTOCOLE_TOTAL(:,7)=CONFIGURATION_UNIQUE(:,2); % Colonne de ranking
                        % Matrice de sensibilité pour le protocole total
                        SENSIBILITE_TOTAL=[MAILLAGE,SYNTHESE_VALUES(:,1)];
                        
                        % Protocole avec nombre de configuration choisi
                        PROTOCOLE_LIMITE=CONFIGURATIONS(CONFIGURATION_UNIQUE(1:Nb_CONFIG,1),:);
                        PROTOCOLE_LIMITE(1:Nb_CONFIG,7)=CONFIGURATION_UNIQUE(1:Nb_CONFIG,2); % Colonne de ranking
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        % RENFORCEMENT DU PROTOCOLE AVEC CONFIGURATIONS POUR N_Channels CANAUX
                        Renforcement = menu('Voulez vous compléter votre protocole optimisé par des mesures supplémentaires (pas de temps de mesure supplémentaire!) ?','OUI','NON');
                        RenforcementD = [1 , 0];
                        Renforcement=RenforcementD(1,Renforcement);
                        
                        if Renforcement==1
                            
                            prompt = {sprintf('Combien de canaux utilisez vous ?')};
                            dlg_title = 'NbCanaux';
                            num_lines = 1;
                            defaultans = {'10'};
                            N_Channels = inputdlg(prompt,dlg_title,num_lines,defaultans);
                            N_Channels = str2num(N_Channels{:});
                            
                            %                         N_Channels=10;
                            
                            CONFIGURATION_UNIQUE_2=CONFIGURATION_UNIQUE;
                            for i=1:length(CONFIGURATION_UNIQUE_2(:,1))
                                CONFIGURATION_UNIQUE_2(i,1:6)=CONFIGURATIONS(CONFIGURATION_UNIQUE(i,1),1:6);
                                CONFIGURATION_UNIQUE_2(i,7)=CONFIGURATION_UNIQUE(i,2);
                            end
                            
                            % TX Matrix for Optimized selected protocol
                            PROTOCOLE_LIMITE_EXT=[0,0,0,0,0,0];
                            TX_LIM=unique(PROTOCOLE_LIMITE(:,1:2),'rows','stable');
                            for i=1:length(TX_LIM(:,1))
                                % Protocole optimisé
                                RX_LIM_rows=find(PROTOCOLE_LIMITE(:,1)==TX_LIM(i,1) & PROTOCOLE_LIMITE(:,2)==TX_LIM(i,2));
                                RX_LIM=PROTOCOLE_LIMITE(find(PROTOCOLE_LIMITE(:,1)==TX_LIM(i,1) & PROTOCOLE_LIMITE(:,2)==TX_LIM(i,2)),3:4);
                                for j=1:length(RX_LIM(:,1))
                                    PROTOCOLE_LIMITE_EXT(end+1,1:4)=[TX_LIM(i,1:2),RX_LIM(j,1:2)];
                                    PROTOCOLE_LIMITE_EXT(end,5:7)=CONFIGURATION_UNIQUE_2(find(PROTOCOLE_LIMITE_EXT(end,1)==CONFIGURATION_UNIQUE_2(:,1)&PROTOCOLE_LIMITE_EXT(end,2)==CONFIGURATION_UNIQUE_2(:,2)&PROTOCOLE_LIMITE_EXT(end,3)==CONFIGURATION_UNIQUE_2(:,3)&PROTOCOLE_LIMITE_EXT(end,4)==CONFIGURATION_UNIQUE_2(:,4)),5:7);
                                end
                                
                                % Protocole optimisé extended
                                RX_LIM_EXT=CONFIGURATION_UNIQUE_2(find(CONFIGURATION_UNIQUE_2(:,1)==TX_LIM(i,1) & CONFIGURATION_UNIQUE_2(:,2)==TX_LIM(i,2)),3:4);
                                elimination=0;
                                for j=1:length(RX_LIM_EXT(:,1))
                                    if isempty(find(RX_LIM_EXT(j,1)==RX_LIM(:,1)&RX_LIM_EXT(j,2)==RX_LIM(:,2)))==0
                                        elimination(end+1,1)=j;
                                    end
                                end
                                elimination(1,:)=[];
                                RX_LIM_EXT(elimination,:)=[];
                                
                                PROTOCOLE_LIMITE_EXT_possibles=[0,0,0,0,0,0];
                                for j=1:length(RX_LIM_EXT(:,1))
                                    PROTOCOLE_LIMITE_EXT_possibles(end+1,1:4)=[TX_LIM(i,1:2),RX_LIM_EXT(j,1:2)];
                                    PROTOCOLE_LIMITE_EXT_possibles(end,5:7)=CONFIGURATION_UNIQUE_2(find(PROTOCOLE_LIMITE_EXT_possibles(end,1)==CONFIGURATION_UNIQUE_2(:,1)&PROTOCOLE_LIMITE_EXT_possibles(end,2)==CONFIGURATION_UNIQUE_2(:,2)&PROTOCOLE_LIMITE_EXT_possibles(end,3)==CONFIGURATION_UNIQUE_2(:,3)&PROTOCOLE_LIMITE_EXT_possibles(end,4)==CONFIGURATION_UNIQUE_2(:,4)),5:7);
                                end
                                PROTOCOLE_LIMITE_EXT_possibles(1,:)=[];
                                if isempty(PROTOCOLE_LIMITE_EXT_possibles)==0
                                    PROTOCOLE_LIMITE_EXT_possibles=sortrows(PROTOCOLE_LIMITE_EXT_possibles,-7);
                                    
                                    count=0;
                                    while mod(length(RX_LIM(:,1))+count,N_Channels)~=0 & count<length(PROTOCOLE_LIMITE_EXT_possibles(:,1))
                                        count=count+1;
                                        PROTOCOLE_LIMITE_EXT(end+1,1:7)=PROTOCOLE_LIMITE_EXT_possibles(count,1:7);
                                    end
                                end
                            end
                            PROTOCOLE_LIMITE_EXT(1,:)=[];
                            PROTOCOLE_LIMITE=PROTOCOLE_LIMITE_EXT;
                        end
                        
                        
                        % Nouveau calcul des sensibilités
                        tic
                        
                        % Chargement des protocoles utiles
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Configurations standard 192 elec simplifié
                        load CONFIGURATION_STANDARD.mat
                        
                        
                        
                        % CHOIX DU PROTOCOLE
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        %                     PROTOCOLE_LIMITE=PROTOCOLE_LIMITE_EXT;
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        %CONFIGURATIONS_SELECT=CONFIGURATION_STANDARD;
                        CONFIGURATIONS_SELECT=PROTOCOLE_LIMITE;
                        %CONFIGURATIONS_SELECT=PROA;
                        %CONFIGURATIONS_SELECT=PROB;
                        %CONFIGURATIONS_SELECT=PROC;
                        
                        
                        SYNTHESE_VALUES_LIMITE=zeros(length(MAILLAGE(:,1)),10);
                        SYNTHESE_INDICE_LIMITE=zeros(length(MAILLAGE(:,1)),10);
                        for i=1:length(CONFIGURATIONS_SELECT(:,1))
                            clear SENSIBILITE_TOT
                            C1_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT(i,1)&DUOS(:,2)==CONFIGURATIONS_SELECT(i,3)));
                            C2_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT(i,2)&DUOS(:,2)==CONFIGURATIONS_SELECT(i,3)));
                            C1_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT(i,1)&DUOS(:,2)==CONFIGURATIONS_SELECT(i,4)));
                            C2_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT(i,2)&DUOS(:,2)==CONFIGURATIONS_SELECT(i,4)));
                            SENSIBILITE_TOT(:,1)=abs(C1_P1-C2_P1-C1_P2+C2_P2);
                            SENS_TOT_AVANCEMENT=round(i/length(CONFIGURATIONS_SELECT(:,1))*100,0)
                            
                            for j=1:length(MAILLAGE(:,1))
                                trouve=0;
                                for k=1:10
                                    if trouve==0
                                        if SENSIBILITE_TOT(j,1)>SYNTHESE_VALUES_LIMITE(j,k)
                                            SYNTHESE_VALUES_LIMITE(j,k)=SENSIBILITE_TOT(j,1);
                                            SYNTHESE_INDICE_LIMITE(j,k)=i;
                                            trouve=1;
                                        end
                                    end
                                end
                            end
                            
                        end
                        toc
                        
                        % Matrice de sensibilité pour le protocole limité
                        SENSIBILITE_LIMITE=MAILLAGE;
                        SENSIBILITE_LIMITE(:,4)=SYNTHESE_VALUES_LIMITE(:,1);
                        
                        % Calcul de l'efficacité du protocole de mesure
                        EFFICACITE=mean(100-(SENSIBILITE_TOTAL(:,4)-SENSIBILITE_LIMITE(:,4))./SENSIBILITE_TOTAL(:,4)*100);
                        % Comparaison graphique simple des deux sensibilités totales
                        if auto==0
                            figure('Color', [ 1 1 1])
                            plot(SENSIBILITE_TOTAL(:,4),SENSIBILITE_LIMITE(:,4),'.k','MarkerSize',5)
                            xlabel(['Sensibilité du comprehensive data set (',num2str(length(CONFIGURATIONS(:,1))),' configurations)'])
                            ylabel(['Sensibilité du protocole identifié (',num2str(length(CONFIGURATIONS_SELECT(:,1))),' configurations)'])
                            xlim([0 0.06])
                            ylim([0 0.06])
                            text(0.005,0.045,['Efficacité du protocole = ',num2str(round(EFFICACITE,1)),' %'],'FontSize',14)
                            grid
                            hold on
                            plot([0 0.06],[0 0.06],'-r')
                        end
                    end
                end
                count=count+1;
                SYNTHESE_TEST(count,1:2)=[Nb_CONFIG,EFFICACITE];
            end
            
            if auto==1
                % Représentation graphique de l'efficacité du protocole versus
                % taille du protocole
                figure('Color', [ 1 1 1])
                semilogx(SYNTHESE_TEST(:,1),SYNTHESE_TEST(:,2),'--ok','MarkerSize',5)
                xlabel('Nombre de configurations du protocole')
                ylabel('Efficactité du protocole en %')
                grid
                ylim([10 100])
            end
        end
        
        %% 4.2) Visualisation 3D des sensibilités des protocoles de mesure
        % Creation le 16/04/2017
        % Modification le 19/04/2017
        actif=1;
        for o=1
            if actif==1
                figure('color',[1 1 1])
                hold on
                
                SENSIBILITE_DATA=SENSIBILITE_LIMITE;
                SENSIBILITE_DATA(:,4)=log(SENSIBILITE_DATA(:,4));
                % Interpolation
                clearvars F;
                F = scatteredInterpolant(SENSIBILITE_DATA(:,1),SENSIBILITE_DATA(:,2),SENSIBILITE_DATA(:,3),SENSIBILITE_DATA(:,4),'natural');
                vq = F(X,Y,Z);
                ValAffichee = vq;
                
                % Suppression des données au dessus de la surface
                for xi=1:length(X(1,:,1))
                    for yi=1:length(X(:,1,1))
                        for zi=1:length(X(1,1,:))
                            if Z(yi,xi,zi)>7.563+pente*(X(yi,xi,zi)-56.74)+0.5
                                ValAffichee(yi,xi,zi)=NaN;
                            end
                        end
                    end
                end
                
                % Mise en forme de la figure
                daspect([1,1,1])
                axis tight
                ax = gca;
                ax.FontSize = 13;
                view(0,0)
                %view(-37.5,30)
                light('Position',[40 -10 50],'Style','local')
                %camzoom(1.8)
                camproj perspective
                colormap (flipud(hot(10)))
                cmap = colormap;
                
                % Plans d'interet
                %                 p2 = slice(X,Y,Z,ValAffichee,[],26,[]);
                %                 p2.FaceColor = 'interp';
                %                 p2.EdgeColor = 'none';
                
                %                 p4 = slice(X,Y,Z,ValAffichee,95,[],[]);
                %                 p4.FaceColor = 'interp';
                %                 p4.EdgeColor = 'none';
                
                %                 p3 = slice(X,Y,Z,ValAffichee,[],[],0);
                %                 p3.FaceColor = 'interp';
                %                 p3.EdgeColor = 'none';
                
                % Affichage des volumes (UNE seule valeur de coupure)
                Valeur=[-7;-5];
                V_Max=-5;
                V_Min=-10;
                Valeur_Indice = round((Valeur-V_Min)*10/(V_Max-V_Min));
                for i=1:length(Valeur_Indice(:,1))
                    p1=patch(isocaps(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
                        'FaceColor','interp','EdgeColor','none','FaceAlpha',0.7)
                    p1 = patch(isosurface(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
                        'FaceColor',cmap(Valeur_Indice(i,1),:),'EdgeColor','none','FaceAlpha',0.7)
                end
                %isonormals(X,Y,Z,ValAffichee,p1)
                %
                % % Affichage des volumes (DEUX valeurs de coupure)
                % % val = ValAffichee;
                % % val(val>2.09999)=nan;
                % % v4 = round((1.9-1.5)*100/0.8);
                % % p1=patch(isocaps(X,Y,Z,val,1.9,'enclose','above'),'FaceColor','interp','EdgeColor','none');
                % % p1 = patch(isosurface(X,Y,Z,val,1.9,'enclose','above'),'FaceColor',cmap(v4,:),'EdgeColor','none');
                % % isonormals(X,Y,Z,ValAffichee,p1)
                % %
                % val2 = ValAffichee;
                % val2(val2<-50)=nan;
                % % v5 = round((50-1.5)*100/0.8);
                % p2=patch(isocaps(X,Y,Z,val2,-50,'enclose','below'),'FaceColor','interp','EdgeColor','none');
                % %p2 = patch(isosurface(X,Y,Z,val2,2.1,'enclose','below'),'FaceColor',cmap(v5,:),'EdgeColor','none');
                
                
                % Légende et titre
                xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
                ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
                zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
                
                
                % Affichage des électrodes (approximation halde rectangulaire)
                for i=1:length(ELEC(:,1))
                    plot3(ELEC(i,1),ELEC(i,2),ELEC(i,3),'k.','MarkerSize',20)
                end
                for i=1:length(ELEC_FAUSSES(:,1))
                    plot3(ELEC(ELEC_FAUSSES(i,1),1),ELEC(ELEC_FAUSSES(i,1),2),ELEC(ELEC_FAUSSES(i,1),3),'.r','MarkerSize',20)
                end
                for i=1:length(ELEC_MAUVAISES(:,1))
                    plot3(ELEC(ELEC_MAUVAISES(i,1),1),ELEC(ELEC_MAUVAISES(i,1),2),ELEC(ELEC_MAUVAISES(i,1),3),'.m','MarkerSize',20)
                end
                
                % Limites de la halde
                %                 xlim([30 95]);
                %                 ylim([21 31]);
                %                 zlim([1 9]);
                plot3(Maillage3DSurface(:,1),Maillage3DSurface(:,2),Maillage3DSurface(:,3),'.-k','LineWidth',1)
                
                % Légende
                clear colorbar
                colorbar
                caxis([-10 -5])
                TicksM = [-10;-9;-8;-7;-6;-5];
                TicksL = round(round(exp(TicksM),4)*10000,0);
                
                c=colorbar;
                c.Label.String = 'Sensitivity cell value (x10000)';
                c.Label.FontWeight='bold';
                c.Label.FontSize=14;
                c.Ticks=TicksM
                c.TickLabels={num2str(TicksL)}
            end
        end % vue 1
        for o=1
            %             if actif==1
            %                 figure('color',[1 1 1])
            %                 hold on
            %
            %                 SENSIBILITE_DATA=SENSIBILITE_LIMITE;
            %                 SENSIBILITE_DATA(:,4)=log(SENSIBILITE_DATA(:,4));
            %                 % Interpolation
            %                 clearvars F;
            %                 F = scatteredInterpolant(SENSIBILITE_DATA(:,1),SENSIBILITE_DATA(:,2),SENSIBILITE_DATA(:,3),SENSIBILITE_DATA(:,4),'natural');
            %                 vq = F(X,Y,Z);
            %                 ValAffichee = vq;
            %
            %                 % Suppression des données au dessus de la surface
            %                 for xi=1:length(X(1,:,1))
            %                     for yi=1:length(X(:,1,1))
            %                         for zi=1:length(X(1,1,:))
            %                             if Z(yi,xi,zi)>7.563+pente*(X(yi,xi,zi)-56.74)+0.5
            %                                 ValAffichee(yi,xi,zi)=NaN;
            %                             end
            %                         end
            %                     end
            %                 end
            %
            %                 % Mise en forme de la figure
            %                 daspect([1,1,1])
            %                 axis tight
            %                 ax = gca;
            %                 ax.FontSize = 13;
            %                 %view(0,0)
            %                 view(-37.5,30)
            %                 light('Position',[40 -10 50],'Style','local')
            %                 %camzoom(1.8)
            %                 camproj perspective
            %                 colormap (flipud(hot(10)))
            %                 cmap = colormap;
            %
            %                 % Plans d'interet
            %                 %                 p2 = slice(X,Y,Z,ValAffichee,[],26,[]);
            %                 %                 p2.FaceColor = 'interp';
            %                 %                 p2.EdgeColor = 'none';
            %
            %                 %                 p4 = slice(X,Y,Z,ValAffichee,95,[],[]);
            %                 %                 p4.FaceColor = 'interp';
            %                 %                 p4.EdgeColor = 'none';
            %
            %                 %                 p3 = slice(X,Y,Z,ValAffichee,[],[],0);
            %                 %                 p3.FaceColor = 'interp';
            %                 %                 p3.EdgeColor = 'none';
            %
            %                 % Affichage des volumes (UNE seule valeur de coupure)
            %                 Valeur=[-7;-6.5;-6;-5.5;-5];
            %                 V_Max=-5;
            %                 V_Min=-10;
            %                 Valeur_Indice = round((Valeur-V_Min)*10/(V_Max-V_Min));
            %                 for i=1:length(Valeur_Indice(:,1))
            %                     p1=patch(isocaps(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
            %                         'FaceColor','interp','EdgeColor','none','FaceAlpha',0.7)
            %                     p1 = patch(isosurface(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
            %                         'FaceColor',cmap(Valeur_Indice(i,1),:),'EdgeColor','none','FaceAlpha',0.7)
            %
            %                     isonormals(X,Y,Z,ValAffichee,p1)
            %                 end
            %                 %isonormals(X,Y,Z,ValAffichee,p1)
            %                 %
            %                 % % Affichage des volumes (DEUX valeurs de coupure)
            %                 % % val = ValAffichee;
            %                 % % val(val>2.09999)=nan;
            %                 % % v4 = round((1.9-1.5)*100/0.8);
            %                 % % p1=patch(isocaps(X,Y,Z,val,1.9,'enclose','above'),'FaceColor','interp','EdgeColor','none');
            %                 % % p1 = patch(isosurface(X,Y,Z,val,1.9,'enclose','above'),'FaceColor',cmap(v4,:),'EdgeColor','none');
            %                 % % isonormals(X,Y,Z,ValAffichee,p1)
            %                 % %
            %                 % val2 = ValAffichee;
            %                 % val2(val2<-50)=nan;
            %                 % % v5 = round((50-1.5)*100/0.8);
            %                 % p2=patch(isocaps(X,Y,Z,val2,-50,'enclose','below'),'FaceColor','interp','EdgeColor','none');
            %                 % %p2 = patch(isosurface(X,Y,Z,val2,2.1,'enclose','below'),'FaceColor',cmap(v5,:),'EdgeColor','none');
            %
            %
            %                 % Légende et titre
            %                 xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
            %                 ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
            %                 zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
            %
            %
            %                 % Affichage des électrodes (approximation halde rectangulaire)
            %                 for i=1:length(ELEC(:,1))
            %                     plot3(ELEC(i,1),ELEC(i,2),ELEC(i,3),'k.','MarkerSize',20)
            %                 end
            %
            %                 % Limites de la halde
            % %                 xlim([30 95]);
            % %                 ylim([26 31]);
            % %                 zlim([1 9]);
            %                 plot3(Maillage3DSurface(:,1),Maillage3DSurface(:,2),Maillage3DSurface(:,3),'.-k')
            %
            %
            %                 % Légende
            %                 clear colorbar
            %                 colorbar
            %                 caxis([-10 -5])
            %                 TicksM = [-10;-9;-8;-7;-6;-5];
            %                 TicksL = round(round(exp(TicksM),4)*10000,0);
            %
            %                 c=colorbar;
            %                 c.Label.String = 'Sensitivity cell value (x10000)';
            %                 c.Label.FontWeight='bold';
            %                 c.Label.FontSize=14;
            %                 c.Ticks=TicksM
            %                 c.TickLabels={num2str(TicksL)}
            %             end
        end % vue 2
        for o=1
            %             if actif==1
            %                 figure('color',[1 1 1])
            %                 hold on
            %
            %                 SENSIBILITE_DATA=SENSIBILITE_LIMITE;
            %                 SENSIBILITE_DATA(:,4)=log(SENSIBILITE_DATA(:,4));
            %                 % Interpolation
            %                 clearvars F;
            %                 F = scatteredInterpolant(SENSIBILITE_DATA(:,1),SENSIBILITE_DATA(:,2),SENSIBILITE_DATA(:,3),SENSIBILITE_DATA(:,4),'natural');
            %                 vq = F(X,Y,Z);
            %                 ValAffichee = vq;
            %
            %                 % Suppression des données au dessus de la surface
            %                 for xi=1:length(X(1,:,1))
            %                     for yi=1:length(X(:,1,1))
            %                         for zi=1:length(X(1,1,:))
            %                             if Z(yi,xi,zi)>7.563+pente*(X(yi,xi,zi)-56.74)+0.5
            %                                 ValAffichee(yi,xi,zi)=NaN;
            %                             end
            %                         end
            %                     end
            %                 end
            %
            %                 % Mise en forme de la figure
            %                 daspect([1,1,1])
            %                 axis tight
            %                 ax = gca;
            %                 ax.FontSize = 13;
            %                 %view(0,0)
            %                 view(-37.5,30)
            %                 light('Position',[40 -10 50],'Style','local')
            %                 %camzoom(1.8)
            %                 camproj perspective
            %                 colormap (flipud(hot(10)))
            %                 cmap = colormap;
            %
            %                 % Plans d'interet
            %                 %                 p2 = slice(X,Y,Z,ValAffichee,[],26,[]);
            %                 %                 p2.FaceColor = 'interp';
            %                 %                 p2.EdgeColor = 'none';
            %
            %                 %                 p4 = slice(X,Y,Z,ValAffichee,95,[],[]);
            %                 %                 p4.FaceColor = 'interp';
            %                 %                 p4.EdgeColor = 'none';
            %
            %                 %                 p3 = slice(X,Y,Z,ValAffichee,[],[],0);
            %                 %                 p3.FaceColor = 'interp';
            %                 %                 p3.EdgeColor = 'none';
            %
            %                 % Affichage des volumes (UNE seule valeur de coupure)
            %                 Valeur=[-7;-6.5;-6;-5.5;-5];
            %                 V_Max=-5;
            %                 V_Min=-10;
            %                 Valeur_Indice = round((Valeur-V_Min)*10/(V_Max-V_Min));
            %                 for i=1:length(Valeur_Indice(:,1))
            %                     p1=patch(isocaps(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
            %                         'FaceColor','interp','EdgeColor','none','FaceAlpha',0.7)
            %                     p1 = patch(isosurface(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
            %                         'FaceColor',cmap(Valeur_Indice(i,1),:),'EdgeColor','none','FaceAlpha',0.7)
            %
            %                     isonormals(X,Y,Z,ValAffichee,p1)
            %                 end
            %                 %isonormals(X,Y,Z,ValAffichee,p1)
            %                 %
            %                 % % Affichage des volumes (DEUX valeurs de coupure)
            %                 % % val = ValAffichee;
            %                 % % val(val>2.09999)=nan;
            %                 % % v4 = round((1.9-1.5)*100/0.8);
            %                 % % p1=patch(isocaps(X,Y,Z,val,1.9,'enclose','above'),'FaceColor','interp','EdgeColor','none');
            %                 % % p1 = patch(isosurface(X,Y,Z,val,1.9,'enclose','above'),'FaceColor',cmap(v4,:),'EdgeColor','none');
            %                 % % isonormals(X,Y,Z,ValAffichee,p1)
            %                 % %
            %                 % val2 = ValAffichee;
            %                 % val2(val2<-50)=nan;
            %                 % % v5 = round((50-1.5)*100/0.8);
            %                 % p2=patch(isocaps(X,Y,Z,val2,-50,'enclose','below'),'FaceColor','interp','EdgeColor','none');
            %                 % %p2 = patch(isosurface(X,Y,Z,val2,2.1,'enclose','below'),'FaceColor',cmap(v5,:),'EdgeColor','none');
            %
            %
            %                 % Légende et titre
            %                 xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
            %                 ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
            %                 zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
            %
            %
            %                 % Affichage des électrodes (approximation halde rectangulaire)
            %                 for i=1:length(ELEC(:,1))
            %                     plot3(ELEC(i,1),ELEC(i,2),ELEC(i,3),'k.','MarkerSize',20)
            %                 end
            %
            %                 % Limites de la halde
            % %                 xlim([30 95]);
            % %                 ylim([21 31]);
            % %                 zlim([1 9]);
            %                 plot3(Maillage3DSurface(:,1),Maillage3DSurface(:,2),Maillage3DSurface(:,3),'.-k')
            %
            %
            %                 % Légende
            %                 clear colorbar
            %                 colorbar
            %                 caxis([-10 -5])
            %                 TicksM = [-10;-9;-8;-7;-6;-5];
            %                 TicksL = round(round(exp(TicksM),4)*10000,0);
            %
            %                 c=colorbar;
            %                 c.Label.String = 'Sensitivity cell value (x10000)';
            %                 c.Label.FontWeight='bold';
            %                 c.Label.FontSize=14;
            %                 c.Ticks=TicksM
            %                 c.TickLabels={num2str(TicksL)}
            %             end
        end % vue 3
    end
end
