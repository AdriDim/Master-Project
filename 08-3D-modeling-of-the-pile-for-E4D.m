%% ------------------- 3D-modeling-of-the-pile-for-E4D --------------------
%
% Adrien Dimech - Master Project - 22/04/2018
%
% -------------------------------------------------------------------------
% Matlab codes to prepare - process - visualize and interpret 3D time-lapse
% geolelectrical monitoring of a waste rock pile.
% -------------------------------------------------------------------------
%
% This Matlab code was used to load surveying data to create a complex 3D
% model for the inversion software E4D (Johnson, 2010). The 3D tetrahedron
% mesh is generated with Tetgen. This code provides 3D visualization of the
% different geometries and structure you wish to model with Tetgen and
% helps verifying the consistency of the model to find meshes errors for
% instance. Some useful tools are also provided to complement and
% interpolate surveying datasets. 
%
% Feel free to visit : https://www.researchgate.net/profile/Adrien_Dimech
% for more information about my research or contact me for more information
% and data files : adrien.dimech@gmail.com
%
%% CODE DE L'ARPENTAGE DE LA HALDE POUR E4D
% -------------------------------------------------------------------------
%

%
%                           Suivi du code
% -------------------------------------------------------------------------
% Creation          |       14/02/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       28/07/2017        |     Adrien Dimech
% -------------------------------------------------------------------------

%% 0) Chargement des donnees
% Creation le 14/02/2017
actif=1;
for o=1:1
    if actif==1
        close all
        clear all
        count = 0;
        indfig=0;
        FichierExcel_maillage = 'MAILLAGEV6.xlsx';
        FichierExcel_arpentage = 'ARPENTAGEV4.xlsx';
        
        figure(1)
        p = get(gcf,'Position');
        set(0,'DefaultFigurePosition',p);
        close all
        
        % Chargement des données
        
        [a_arpentage,b_arpentage]=xlsfinfo(FichierExcel_arpentage);
        NbPages_arpentage = length(b_arpentage(1,:));
        
        % Selection des donnees d'arpentage
        Choix = ['1 = ',b_arpentage{1}];
        for i=2:NbPages_arpentage
            Choix = [Choix, ' \n ', num2str(i), ' = ',b_arpentage{i}];
        end
        Choix = sprintf(Choix);
        
        prompt = {sprintf(['Quelles donnees d"arpentage pour les BORDS ?  \n \n  Pour selection multiple utilisez des espaces \n \n ',Choix, '\n'])};
        dlg_title = 'Donnees_BORD';
        num_lines = 1;
        defaultans = {'1'};
        Donnees_BORD = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Donnees_BORD = str2num(Donnees_BORD{:});
        
        prompt = {sprintf(['Quelles donnees d"arpentage pour la SURFACE ?  \n \n  Pour selection multiple utilisez des espaces \n \n ',Choix, '\n'])};
        dlg_title = 'Donnees_SURFACE';
        num_lines = 1;
        defaultans = {'2'};
        Donnees_SURFACE = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Donnees_SURFACE = str2num(Donnees_SURFACE{:});
        
        prompt = {sprintf(['Quelles donnees d"arpentage pour les POINTS INTERNES ?  \n \n  Pour selection multiple utilisez des espaces \n \n ',Choix, '\n'])};
        dlg_title = 'Donnees_POINTSINTERNES';
        num_lines = 1;
        defaultans = {'5 6 7'};
        Donnees_POINTSINTERNES = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Donnees_POINTSINTERNES = str2num(Donnees_POINTSINTERNES{:});
        
        % Chargement des données de maillage
        [a_maillage,b_maillage]=xlsfinfo(FichierExcel_maillage);
        NbPages_maillage = length(b_maillage(1,:));
        
        % Selection des donnees de maillage
        Choix = ['1 = ',b_maillage{1}];
        for i=2:NbPages_maillage
            Choix = [Choix, ' \n ', num2str(i), ' = ',b_maillage{i}];
        end
        Choix = sprintf(Choix);
        
        prompt = {sprintf(['Quelles donnees pour le maillage des ZONES ?  \n \n  Pour selection multiple utilisez des espaces \n \n ',Choix, '\n'])};
        dlg_title = 'Donnees_ZONES';
        num_lines = 1;
        defaultans = {''};
        Donnees_ZONES = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Donnees_ZONES = str2num(Donnees_ZONES{:});
    end
end

%% 0.1) Traitement des donnees d'arpentage automatique
% on supprime les points de l'arpentage automatique s'ils sont plus près de
% "dist_min" des aretes pour éviter de générer des artefacts de maillage
% (noter que aretes50cm a la priorité dans ce cas)
actif=0;
for o=1:1
    if actif==1
        dist_min=0.75;
        clear dist
        
        load reg50cm.mat
        load reg1m.mat
        load aretes50cm.mat
        reg=reg1m;
        reg_traite=0;
        for i=1:length(reg(:,1))
            round(i/length(reg(:,1))*100,2)
            for j=1:length(aretes50cm(:,1))
                dist(j,1)=pdist([reg(i,1:3);aretes50cm(j,1:3)]);
            end
            if min(dist(:,1))>dist_min
                reg_traite(end+1,1:3)=reg(i,1:3);
            end
        end
        reg_traite(1,:)=[];
        % Affichage des donnees traitees
        indfig = indfig +1 ;
        figure('Color', [ 1 1 1])
        hold on
        plot3(reg(:,1),reg(:,2),reg(:,3),'.r','MarkerSize',3)
        plot3(reg_traite(:,1),reg_traite(:,2),reg_traite(:,3),'.g','MarkerSize',5)
        plot3(aretes50cm(:,1),aretes50cm(:,2),aretes50cm(:,3),'.b','MarkerSize',6)
        legend('avant traitement','apres traitement','aretes aux 50 cm');
        xlim([0 110])
        ylim([0 50])
        zlim([0 10])
    end
end

%% GÉNÉRATION DE COUCHES RAFFINÉES
actif=0;
for o=1
    if actif==1
        load anortho.mat
        load sable.mat
        MAILLAGE_Raf=sable;
        X_min=min(MAILLAGE_Raf(:,1));
        X_max=max(MAILLAGE_Raf(:,1));
        Y_min=min(MAILLAGE_Raf(:,2));
        Y_max=max(MAILLAGE_Raf(:,2));
        
        pas_X=0.25; % en m => 20 cm
        pas_Y=0.25; % en m => 20 cm
        
        X_reg=(floor(X_min):pas_X:floor(X_max)+1)';
        L_X=length(X_reg);
        Y_reg=(floor(Y_min):pas_Y:floor(Y_max)+1)';
        L_Y=length(Y_reg);
        
        REG_Raf=zeros(L_X*L_Y,3);
        for i=1:L_X
            REG_Raf((i-1)*L_Y+1:i*L_Y,1)=ones(L_Y,1)*X_reg(i,1);
            REG_Raf((i-1)*L_Y+1:i*L_Y,2)=Y_reg;
        end
        
        REG_Raf_2=zeros(L_X*L_Y,3);
        for i=1:length(MAILLAGE_Raf(:,1))/3
            lignes=inpolygon(REG_Raf(:,1),REG_Raf(:,2),MAILLAGE_Raf((i-1)*3+1:i*3,1),MAILLAGE_Raf((i-1)*3+1:i*3,2));
            count=1;
            clear lignes_2
            for j=1:L_X*L_Y
                if lignes(j,1)==1
                    lignes_2(count,1)=j;
                    count=count+1;
                end
            end
            
            REG_Raf_2(lignes_2,1:2)=REG_Raf(lignes_2,1:2);
            % Reconstruction de la hauteur par interpolation
            A=MAILLAGE_Raf((i-1)*3+1,:)';
            B=MAILLAGE_Raf((i-1)*3+2,:)';
            C=MAILLAGE_Raf((i-1)*3+3,:)';
            AB=B-A;
            AC=C-A;
            n=cross(AB,AC);
            X_2=REG_Raf_2(lignes_2,1)-A(1);
            Y_2=REG_Raf_2(lignes_2,2)-A(2);
            REG_Raf_2(lignes_2,3)=-(n(1)*X_2+n(2)*Y_2)/n(3)+A(3);
        end
        REG_Raf_2(REG_Raf_2(:,1)==0,:)=[];
        figure('Color', [ 1 1 1])
        plot3(MAILLAGE_Raf(:,1),MAILLAGE_Raf(:,2),MAILLAGE_Raf(:,3),'r.','markersize',10)
        hold on
        plot3(REG_Raf_2(:,1),REG_Raf_2(:,2),REG_Raf_2(:,3),'k.','markersize',5)
        
        
    end
end

%% 1) Affichage des donnees d'arpentage et du maillage des zones
% Creation le 14/02/2017
% Modification le 18/02/2017
actif=1;
raffinage=0;
for o=1:1
    if actif==1
        prec=3;
        dimraf=0.05; % cubes de 20 cm de côté
        maillage_true=1;
        
        indfig = indfig +1 ;
        figure('Color', [ 1 1 1])
        hold on
        Xname = '\fontsize{18}\bf Coord X';
        Yname = '\fontsize{18}\bf Coord Y';
        % Récupération des données de Excel
        for i=1:length(Donnees_BORD(1,:))
            DATA{i,1} = round(xlsread(FichierExcel_arpentage,Donnees_BORD(1,i)),prec);
            DATA{i,2}=2;
        end
        for i=1:length(Donnees_SURFACE(1,:))
            DATA{end+1,1} = round(xlsread(FichierExcel_arpentage,Donnees_SURFACE(1,i)),prec);
            DATA{end,2}=1;
        end
        if isempty(Donnees_POINTSINTERNES)==0
            DATA_interne=[0,0,0];
            for i=1:length(Donnees_POINTSINTERNES(1,:))
                DATA_interne=[DATA_interne;round(xlsread(FichierExcel_arpentage,Donnees_POINTSINTERNES(1,i)),3)];
                DATA{end+1,1} = round(xlsread(FichierExcel_arpentage,Donnees_POINTSINTERNES(1,i)),3);
                DATA{end,2}=0;
            end
            DATA_interne(1,:)=[];
            
            
            % Ajout de points sur cube <3 pour raffiner maillage interne
            % Modification le 18/02/2017
            if raffinage==1
                DATA_interne_raf=0;
                for i=1:length(DATA_interne(:,1))
                    for j=1:2
                        for k=1:2
                            for l=1:2
                                DATA_interne_raf(end+1,1)=DATA_interne(i,1)-dimraf+(j==2)*dimraf*2;
                                DATA_interne_raf(end,2)=DATA_interne(i,2)-dimraf+(k==2)*dimraf*2;
                                DATA_interne_raf(end,3)=DATA_interne(i,3)-dimraf+(l==2)*dimraf*2;
                            end
                        end
                    end
                end
                DATA_interne_raf(1,:)=[];
                DATA_interne_raf=unique(DATA_interne_raf,'rows');
                DATA{end+1,1} = DATA_interne_raf;
                DATA{end,2}=0;
            end
        end
        
        DATA2(:,:) = DATA{1,1}(:,:);
        mat=ones(length(DATA{1,1}(:,1)));
        DATA2(:,4)=DATA{1,2}*mat(:,1);
        if length(DATA)>1
            for i=2:length(DATA)
                clear data
                data(:,:)=DATA{i,1}(:,:);
                mat=ones(length(DATA{i,1}(:,1)));
                data(:,4)=DATA{i,2}*mat(:,1);
                DATA2 = [DATA2 ; data];
            end
        end
        DATA2=unique(DATA2,'rows');
        % Affichage des données d'arpentage
        title('\fontsize{20}\bf Arpentage de la halde');
        xlabel(Xname);
        ylabel(Yname);
        plot3(DATA2(find(DATA2(:,4)==0),1),DATA2(find(DATA2(:,4)==0),2),DATA2(find(DATA2(:,4)==0),3),'.r','MarkerSize',10)
        plot3(DATA2(find(DATA2(:,4)==1),1),DATA2(find(DATA2(:,4)==1),2),DATA2(find(DATA2(:,4)==1),3),'.b','MarkerSize',10)
        plot3(DATA2(find(DATA2(:,4)==2),1),DATA2(find(DATA2(:,4)==2),2),DATA2(find(DATA2(:,4)==2),3),'.g','MarkerSize',10)
        legend('internal points','surface topo points','external boudaries');
        xlim([0 110])
        ylim([0 50])
        zlim([0 10])
        
        % Récupération des données de maillage Excel
        if isempty(Donnees_ZONES)==0
            for i=1:length(Donnees_ZONES(1,:))
                DATA_MAILLAGE{i,1} = round(xlsread(FichierExcel_maillage,Donnees_ZONES(1,i)),10);
            end
            
            if maillage_true==1
                % Affichage des données de maillage
                for i=1:length(Donnees_ZONES(1,:))
                    plot3(DATA_MAILLAGE{i,1}(:,1),DATA_MAILLAGE{i,1}(:,2),DATA_MAILLAGE{i,1}(:,3),'k.','MarkerSize',8)
                    for j=1:length(DATA_MAILLAGE{i,1}(:,1))/3
                        line(DATA_MAILLAGE{i,1}((j-1)*3+1:(j-1)*3+2,1),DATA_MAILLAGE{i,1}((j-1)*3+1:(j-1)*3+2,2),DATA_MAILLAGE{i,1}((j-1)*3+1:(j-1)*3+2,3),'Color','k','LineWidth',2)
                        line(DATA_MAILLAGE{i,1}((j-1)*3+2:(j-1)*3+3,1),DATA_MAILLAGE{i,1}((j-1)*3+2:(j-1)*3+3,2),DATA_MAILLAGE{i,1}((j-1)*3+2:(j-1)*3+3,3),'Color','k','LineWidth',2)
                        line([DATA_MAILLAGE{i,1}((j-1)*3+3,1),DATA_MAILLAGE{i,1}((j-1)*3+1,1)],[DATA_MAILLAGE{i,1}((j-1)*3+3,2),DATA_MAILLAGE{i,1}((j-1)*3+1,2)],[DATA_MAILLAGE{i,1}((j-1)*3+3,3),DATA_MAILLAGE{i,1}((j-1)*3+1,3)],'Color','k','LineWidth',2)
                        % Remplissage des triangles du maillage
                        X = DATA_MAILLAGE{i,1}((j-1)*3+1:(j-1)*3+3,1);
                        Y = DATA_MAILLAGE{i,1}((j-1)*3+1:(j-1)*3+3,2);
                        Z = DATA_MAILLAGE{i,1}((j-1)*3+1:(j-1)*3+3,3);
                        C = [0.5 0.5 0.5 ];
                        fill3(X,Y,Z,C,'FaceAlpha',0.5);
                    end
                end
            end
        end
        view(-16.42,37.12)
    end
end

%% 2) Generation automatique du .cfg pour E4D
% Creation le 14/02/2016
actif=1;
regularisation=0;
for o=1:1
    if actif==1
        % CONTROL POINTS
        % points de controle des données d'arpentage
        CONTROL_POINTS=DATA2;
        if isempty(Donnees_ZONES)==0
            for i=1:length(DATA_MAILLAGE)
                % points de controle des defs de zones
                clear c_p
                c_p=DATA_MAILLAGE{i};
                c_p(:,4)=zeros(length(c_p(:,1)),1);
                c_p=unique(c_p,'rows');
                CONTROL_POINTS=[CONTROL_POINTS ; c_p];
            end
            CONTROL_POINTS=round(CONTROL_POINTS,prec);
            CONTROL_POINTS=unique(CONTROL_POINTS,'rows');
            
            Nb_CONTROL_POINTS=length(CONTROL_POINTS(:,1));
            
            % DEF_PLAN
            for i=1:length(DATA_MAILLAGE)
                Maill = round(DATA_MAILLAGE{i,1}(:,:),prec);
                for j=1:length(Maill(:,1))/3
                    for k=1:3
                        ind=(j-1)*3+k;
                        Def_Plan{i,1}(j,k)=find(Maill(ind,1)==CONTROL_POINTS(:,1)& Maill(ind,2)==CONTROL_POINTS(:,2) & Maill(ind,3)==CONTROL_POINTS(:,3));
                    end
                end
            end
        end
        
        % RÉDACTION DU DOCUMENT
        clear DOC
        ligne=1;
        Nb_CONTROL_POINTS=length(CONTROL_POINTS(:,1));
        DOC(ligne,1)=Nb_CONTROL_POINTS;
        DOC(ligne,2:5)=NaN;
        ligne=ligne+1;
        for i=1:Nb_CONTROL_POINTS
            DOC(ligne,1)=i;
            DOC(ligne,2:5)=CONTROL_POINTS(i,1:4);
            ligne=ligne+1;
        end
        
        Nb_PLANS=0;
        % SI DEFINITION DE ZONES
        Ligne_Nb_PLANS=0;
        if isempty(Donnees_ZONES)==0 % Si def de zones
            for i=1:length(Def_Plan)
                Nb_PLANS=Nb_PLANS+length(Def_Plan{i}(:,1));
            end
            % Nombre de plans pour l'ensemble du doc
            DOC(ligne,1)=Nb_PLANS;
            % SAUVEGARDE DE LA LIGNE POUR MODIFICATION POSSIBLE
            Ligne_Nb_PLANS=ligne;
            DOC(ligne,2:5)=NaN;
            ligne=ligne+1;
            
            for i=1:length(Def_Plan)
                ind_ZONE=10+i;
                for j=1:length(Def_Plan{i}(:,1))
                    % Nb points du plan et indice de la zone
                    DOC(ligne,1)=3;
                    DOC(ligne,2)=ind_ZONE;
                    DOC(ligne,3:5)=NaN;
                    ligne=ligne+1;
                    % Points définissant le plan
                    DOC(ligne,1:3)=Def_Plan{i}(j,1:3);
                    DOC(ligne,4:5)=NaN;
                    ligne=ligne+1;
                end
            end
        end
        
        % SI REGULARISATION
        if regularisation==1
            % mise en forme DATA GS3 fichier de données de régularisation
            
            load DATA_GS3.mat
            
            % SPECIAL POUR TROUVER COORD OK GS3
            DATA_GS3(:,2:4)=round(xlsread(FichierExcel_arpentage,6),prec);
            
            DATA_GS3(isnan(DATA_GS3(:,5)),:)=[];
            
            NB_GS3_reg=length(DATA_GS3(:,1));
            % modification du nb de plans
            if Ligne_Nb_PLANS==0
                Nb_PLANS=NB_GS3_reg*6;
                % Nombre de plans pour l'ensemble du doc
                DOC(ligne,1)=Nb_PLANS;
                DOC(ligne,2:5)=NaN;
                ligne=ligne+1;
            else
                Nb_PLANS=Nb_PLANS+NB_GS3_reg*6;
                % Nombre de plans pour l'ensemble du doc
                DOC(Ligne_Nb_PLANS,1)=Nb_PLANS;
            end
            
            for i=1:NB_GS3_reg
                % Coord de la GS3
                X_GS3=DATA_GS3(i,2);
                Y_GS3=DATA_GS3(i,3);
                Z_GS3=DATA_GS3(i,4);
                
                % trouver l'indice des nodes des coins correspondant
                Mat_CUBE=[ones(8,1)*X_GS3,ones(8,1)*Y_GS3,ones(8,1)*Z_GS3];
                dimraf=0.05;
                ligne_cube=1;
                for j=1:2
                    for k=1:2
                        for l=1:2
                            Mat_CUBE(ligne_cube,1)=Mat_CUBE(ligne_cube,1)-dimraf+(j==2)*dimraf*2;
                            Mat_CUBE(ligne_cube,2)=Mat_CUBE(ligne_cube,2)-dimraf+(k==2)*dimraf*2;
                            Mat_CUBE(ligne_cube,3)=Mat_CUBE(ligne_cube,3)-dimraf+(l==2)*dimraf*2;
                            ligne_cube=ligne_cube+1;
                        end
                    end
                end
                Mat_CUBE=round(Mat_CUBE,3);
                
                Mat_CUBE_Node=ones(8,1);
                for j=1:8
                    Mat_CUBE_Node(j,1)=DOC(find(DOC(:,2)==Mat_CUBE(j,1)&DOC(:,3)==Mat_CUBE(j,2)&DOC(:,4)==Mat_CUBE(j,3)),1);
                end
                
                % Matrice des faces des cubes
                MAT_FACES=[1,2,6,5;5,6,8,7;5,1,3,7;1,2,4,3;4,8,7,3;2,4,8,6];
                MAT_FACES_Node=zeros(6,4);
                for j=1:4
                    MAT_FACES_Node(:,j)=Mat_CUBE_Node(MAT_FACES(:,j));
                end
                
                % Rédaction du .cfg
                ind_ZONE=10000+i;
                for j=1:6
                    % Nb points du plan et indice de la zone
                    DOC(ligne,1)=4;
                    DOC(ligne,2)=ind_ZONE;
                    DOC(ligne,3:5)=NaN;
                    ligne=ligne+1;
                    % Points définissant le plan
                    DOC(ligne,1:4)=MAT_FACES_Node(j,:);
                    DOC(ligne,5)=NaN;
                    ligne=ligne+1;
                end
            end
            
            % REDACTION DOC INFO ZONES
            DOC_ZONES=zeros(NB_GS3_reg+2+Donnees_ZONES,6);
            % Corps de la halde
            DOC_ZONES(1,1:6)=[NB_GS3_reg+1+Donnees_ZONES,NaN,NaN,NaN,NaN,NaN];
            DOC_ZONES(2,1:6)=[1,33,26,0.75,1000000,0.001];
            if Donnees_ZONES==1 %LYSI
                DOC_ZONES(3,1:6)=[2,48.75,25.0,1.78,1000000,0.00000001];
            end
            for j=1:NB_GS3_reg
                DOC_ZONES(j+2+Donnees_ZONES,1)=j+1+Donnees_ZONES;
                DOC_ZONES(j+Donnees_ZONES+2,2:4)=round(DATA_GS3(j,2:4),3);
                DOC_ZONES(j+Donnees_ZONES+2,5)=1000000;
                DOC_ZONES(j+Donnees_ZONES+2,6)=round(1/DATA_GS3(j,5),4);
            end
            
        end
        
        % SI AUCUNE DEF DE ZONE NI REGULARISATION
        if isempty(Donnees_ZONES)==1 & regularisation==0
            DOC(ligne,1)=0;
            DOC(ligne,2:5)=NaN;
            ligne=ligne+1;
        end
        
        
        % nb de holes in the mesh
        DOC(ligne,1)=0;
        DOC(ligne,2:5)=NaN;
        ligne=ligne+1;
    end
end

DOC(:,6)=ones(length(DOC(:,1)),1)*NaN;
DOC_TOT=[DOC;DOC_ZONES];


%% 2.1) Generation automatique du .inv pour E4D
% Creation le 02/08/2017
actif=1;
for o=1:1
    if actif==1
        clear INV_ZONES
        
        % Intégration des lysimètres ?
        Nb_Zones=0;
        if Donnees_ZONES==1
            Nb_Zones=1;
        end
        % Nombre total de régularisations de zones
        Nb_Reg=1+Nb_Zones+NB_GS3_reg;
        INV_ZONES=[Nb_Reg,NaN,NaN,NaN];
        
        % LISSAGE SPATIAL DANS HALDE
        INV_ZONES(end+1,:)=[1,NaN,NaN,NaN];
        INV_ZONES(end+1,:)=[2,1,1,1];
        INV_ZONES(end+1,:)=[1,2,0.001,NaN];
        INV_ZONES(end+1,:)=[0,NaN,NaN,NaN];
        INV_ZONES(end+1,:)=[0,NaN,NaN,NaN];
        INV_ZONES(end+1,:)=[1,NaN,NaN,NaN];
        
        % MIN VALUE DANS LYSI
        INV_ZONES(end+1,:)=[2,NaN,NaN,NaN];
        INV_ZONES(end+1,:)=[3,1,1,1];
        INV_ZONES(end+1,:)=[1,0,0.05,NaN];
        INV_ZONES(end+1,:)=[0,NaN,NaN,NaN];
        INV_ZONES(end+1,:)=[log(0.00001),NaN,NaN,NaN];
        INV_ZONES(end+1,:)=[1,NaN,NaN,NaN];
        
        % KNOWN VALUE DANS LES GS3
        for i=1:NB_GS3_reg
            % ligne blanche
            INV_ZONES(end+1,:)=[NaN,NaN,NaN,NaN];
            
            INV_ZONES(end+1,:)=[i+2,NaN,NaN,NaN];
            INV_ZONES(end+1,:)=[3,1,1,1];
            INV_ZONES(end+1,:)=[3,0,0.01,NaN];
            INV_ZONES(end+1,:)=[0,NaN,NaN,NaN];
            INV_ZONES(end+1,:)=[log(1/DATA_GS3(i,5)),NaN,NaN,NaN];
            INV_ZONES(end+1,:)=[1,NaN,NaN,NaN];
        end
        
        complement_LISS(1)=NB_GS3_reg;
        for i=1:NB_GS3_reg
            complement_LISS(1,end+1)=i+1;
        end
    end
end

%% 2) Generation automatique de donnees d'arpentage
% Creation le 15/02/2017
actif=0;
for o=1:1
    if actif==1
        
        %Interpolation = menu('Creer des données par interpolation ?','OUI','NON');
        Interpolation = 1;
        
        
        if Interpolation==1
            fig = figure(indfig);
            dcm_obj = datacursormode(fig);
            % set(dcm_obj,'DisplayStyle','datatip',...
            %     'SnapToDataVertex','off','Enable','on')
            
            c_info = getCursorInfo(dcm_obj);
            for i=1:length(c_info)
                Coord(i,1:3) = c_info(i).Position;
            end
            
            %Nombre de points a aujouter
            %NbPoints = menu('Combien de points a ajouter ?','2 points','5 points','7 points','10 points','15 points','20 points');
            NbPoints = 20;
            
            
            for i=1:3
                dCoord(1,i) = max(Coord(:,i))- min(Coord(:,i));
                Pas(1,i) = dCoord(1,i)/(NbPoints+1);
                dCoord2(1,i) = Coord(2,i)-Coord(1,i);
                if dCoord2(1,i)>0
                    cCoord(1,i)=1;
                else
                    cCoord(1,i)=-1;
                end
            end
            
            for i=1:NbPoints
                for j=1:3
                    Interp2(i,j)=Coord(1,j)+i*cCoord(1,j)*Pas(1,j);
                end
            end
            
            if count < 1
                Interp = Interp2;
                count =1;
            end
            Interp=[Interp;Interp2];
            %Interp=[Interp;Interp2;Coord];
            plot3(Interp(:,1),Interp(:,2),Interp(:,3),'r.','MarkerSize',10)
        end
        
        % xlim([0 110])
        % ylim([0 50])
        % zlim([0 10])
    end
end

%% 3) Generation automatique du maillage triangulaire
% Creation le 15/02/2017
actif=1;
for o=1:1
    if actif==1
        fig = figure(1);
        dcm_obj = datacursormode(fig);
        
        c_info = getCursorInfo(dcm_obj);
        for i=1:length(c_info)
            Coord(i,1:3) = c_info(i).Position;
        end
        
        Maill2=Coord;
        if count < 1
            Maill = Maill2;
            count =1;
        else
            Maill=[Maill;Maill2];
        end
        plot3(Maill2(:,1),Maill2(:,2),Maill2(:,3),'r.','MarkerSize',20)
        line(Maill2(1:2,1),Maill2(1:2,2),Maill2(1:2,3),'Color','r','LineWidth',4)
        line(Maill2(2:3,1),Maill2(2:3,2),Maill2(2:3,3),'Color','r','LineWidth',4)
        line([Maill2(3,1),Maill2(1,1)],[Maill2(3,2),Maill2(1,2)],[Maill2(3,3),Maill2(1,3)],'Color','r','LineWidth',4)
    end
end

%% 4) Création maillage automatique pour E4D
% Creation le 07/02/2017
actif=0;
for o=1:1
    if actif==1
        Pas=0.1;
        load Maillage3DSurface.mat
        DATA2=Maillage3DSurface;
        % Maillage régulier de la halde
        [X,Y]=meshgrid(0:Pas:110,0:Pas:50);
        
        Interp_Zmemoire=zeros(length(X(:,1)),length(X(1,:)));
        INmemoire=zeros(length(X(:,1)),length(X(1,:)));
        for i=1:length(DATA2(:,1))/3 % test pour chaque polygone constituant le maillage
            ii=(i-1)*3+1;
            IN = inpolygon(X(:,:),Y(:,:),DATA2(ii:ii+2,1),DATA2(ii:ii+2,2)); % test pour chaque point du maillage régulier
            IN=IN-INmemoire;
            IN(IN(:,:)==-1)=0;
            INmemoire=INmemoire+IN;
            % Equation du plan
            A=DATA2(ii,1:3);
            B=DATA2(ii+1,1:3);
            C=DATA2(ii+2,1:3);
            AB=B-A;
            AC=C-A;
            n=cross(AB,AC);
            Interp_Z=A(3)-(n(1)*(X-A(1))+n(2)*(Y-A(2)))/n(3);
            Interp_Z=Interp_Z.*IN;
            Interp_Zmemoire=Interp_Zmemoire+Interp_Z;
        end
        clear ARP_NEW
        ARP_NEW=0;
        for i=1:length(X(:,1))
            for j=1:length(X(1,:))
                ARP_NEW(end+1,1)=X(i,j);
                ARP_NEW(end,2)=Y(i,j);
                ARP_NEW(end,3)=Interp_Zmemoire(i,j);
            end
        end
        ARP_NEW(1,:)=[];
        ARP_NEW(ARP_NEW(:,1)==0,:)=[];
        ARP_NEW(ARP_NEW(:,2)==0,:)=[];
        ARP_NEW(ARP_NEW(:,1)==110,:)=[];
        ARP_NEW(ARP_NEW(:,2)==50,:)=[];
        ARP_NEW(ARP_NEW(:,3)==0,:)=[];
        
        figure('Color', [ 1 1 1])
        hold on
        plot3(ARP_NEW(:,1),ARP_NEW(:,2),ARP_NEW(:,3),'.')
        daspect([1,1,0.5])
        view(101.34,20.48)
        %camzoom(2)
    end
end


