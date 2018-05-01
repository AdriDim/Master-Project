## --------------- 3D-modeling-of-the-waste-rock-pile ---------------------

 Adrien Dimech - Master Project - 21/04/2018
 -------------------------------------------------------------------------
 Matlab codes to prepare - process - visualize and interpret 3D time-lapse 
 geolelectrical monitoring of a waste rock pile.
 -------------------------------------------------------------------------

 This Matlab code was designed to use surveying of the pile to create a
 complex 3D model of the pile with both external topography and internal
 structure and instrumentation positions. The 3D model is used to provide
 better inversion results for the 3D time-lapse geoelectrical monitoring.
 of the experimental waste rock pile.
 
 Feel free to visit : https://www.researchgate.net/profile/Adrien_Dimech
 for more information about my research or contact me for more information
 and data files : adrien.dimech@gmail.com
 
                           Suivi du code
 -------------------------------------------------------------------------
 Creation          |       01/10/2016        |     Adrien Dimech
 -------------------------------------------------------------------------
 Modification      |       09/11/2016        |     Adrien Dimech
 -------------------------------------------------------------------------
 Modification      |       17/11/2016        |     Adrien Dimech
 -------------------------------------------------------------------------
 Modification      |       25/01/2017        |     Adrien Dimech
 -------------------------------------------------------------------------
 Modification      |       31/01/2017        |     Adrien Dimech
 -------------------------------------------------------------------------
 Modification      |       03/02/2017        |     Adrien Dimech
 -------------------------------------------------------------------------
 Modification      |       07/02/2017        |     Adrien Dimech
 -------------------------------------------------------------------------
 Modification      |       13/02/2017        |     Adrien Dimech
 -------------------------------------------------------------------------

## 0.1) Rotation des donneés brutes
 Creation le 17/11/2016
actif=0;
for o=1:1 % non utilisé : Rotation des données
    if actif==1
        close all
        clear all
        FichierExcel = 'ARPENTAGEV2.xlsx';
        [a,b]=xlsfinfo(FichierExcel);
        NbPages = length(b(1,:));
        
        Ymin = 2.435360*10^6;
        Xmin = -5.06059*10^6;
        Zmin = 213;
        
        X0 = 3.824469562857140e+05;
        Y0 = 5.603052734142860e+06;
        
        X1 = 3.824732754761910e+05;
        Y1 = 5.602986912476190e+06;
        
        Theta = -acosd((X1-X0)/sqrt((X1-X0)^2+(Y1-Y0)^2));
        
        for i=1:NbPages
            DATA{i} = xlsread(FichierExcel,i);
            DATA2 = DATA{i};
            DATA3=DATA2;
            DATA3(:,1) = DATA2(:,1).*cosd(Theta)+DATA2(:,2).*sind(Theta)-Xmin;
            DATA3(:,2) = DATA2(:,2).*cosd(Theta)-DATA2(:,1).*sind(Theta)-Ymin;
            DATA3(:,3) = DATA2(:,3)-Zmin;
            DATA{i}(:,:) = DATA3(:,:);
        end
    end
end

%% 0) Chargement des donnees
% Creation le 01/10/2016
actif=1;
for o=1:1
    if actif==1
        close all
        clear all
        indfig = 0;
        count = 0;
        
        figure(1)
        p = get(gcf,'Position');
        set(0,'DefaultFigurePosition',p);
        close all
        
        Normalisation = 2;
        for o=1:1 % Non utilisé : Normalisation des données
            % % Elements de normalisation
            % Xo = 3.824431970000000e+05;
            % Yo = 5.602979262000000e+06;
            % Zo = 2.176130000000000e+02-3;
            %
            % Normalisation = menu('Normaliser les données ?','OUI','NON');
            %
            % XE = [6.40;8.31;10.13];
            % YE = [59.85;60.77;61.61];
            % ZE = [6.614;6.629;6.636];
            %
            % if Normalisation==1
            %     XE=XE+Xo;
            %     YE=YE+Yo;
            %     ZE=ZE+Zo;
            % end
        end
        
        % Chargement des données
        FichierExcel = 'ARPENTAGEV6.xlsx';
        [a,b]=xlsfinfo(FichierExcel);
        NbPages = length(b(1,:));
        
        % Selection des donnees d'arpentage
        Choix = ['1 = ',b{1}];
        for i=2:NbPages
            Choix = [Choix, ' \n ', num2str(i), ' = ',b{i}];
        end
        Choix = sprintf(Choix);
        
        prompt = {sprintf(['Quelles donnees d"arpentage ?  \n \n  Pour selection multiple utilisez des espaces \n \n ',Choix, '\n'])};
        dlg_title = 'Donnees';
        num_lines = 1;
        defaultans = {'1 2 3 4'};
        Donnees = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Donnees = str2num(Donnees{:});
    end
end

%% 1) Affichage des donnees d'arpentage
% Creation le 09/11/2016
actif=1;
for o=1:1
    if actif==1
        indfig = indfig +1 ;
        figure('Color', [ 1 1 1])
        hold on
        
        for i=1:length(Donnees(1,:))
            DATA{i} = xlsread(FichierExcel,Donnees(1,i));
            if Normalisation==1
                DATA{i}(:,1) = DATA{i}(:,1) - Xo;
                DATA{i}(:,2) = DATA{i}(:,2) - Yo;
                DATA{i}(:,3) = DATA{i}(:,3) - Zo;
                Xname = '\fontsize{18}\bf Coord X normalisees (m)';
                Yname = '\fontsize{18}\bf Coord Y normalisees (m)';
            else
                Xname = '\fontsize{18}\bf Coord X';
                Yname = '\fontsize{18}\bf Coord Y';
            end
        end
        
        DATA2(:,:) = DATA{1}(:,:);
        if length(Donnees(1,:))>1
            for i=2:length(Donnees(1,:))
                DATA2 = [DATA2 ; DATA{i}(:,:)];
            end
        end
        
        % % Rotation des données
        % rotation = menu('Effectuer une rotation des données ?','OUI','NON');
        %
        % if rotation ==1
        %     X0 = 3.824469562857140e+05;
        %     Y0 = 5.603052734142860e+06;
        %
        %     X1 = 3.824732754761910e+05;
        %     Y1 = 5.602986912476190e+06;
        %
        %     Theta = -acosd((X1-X0)/sqrt((X1-X0)^2+(Y1-Y0)^2));
        %
        %     DATA3=DATA2;
        %     DATA3(:,1) = DATA2(:,1).*cosd(Theta)+DATA2(:,2).*sind(Theta);
        %     DATA3(:,2) = DATA2(:,2).*cosd(Theta)-DATA2(:,1).*sind(Theta);
        %     DATA2 = DATA3;
        % end
        
        title('\fontsize{20}\bf Arpentage de la halde');
        xlabel(Xname);
        ylabel(Yname);
        plot3(DATA2(:,1),DATA2(:,2),DATA2(:,3),'.','MarkerSize',5)
        
        xlim([0 110])
        ylim([0 50])
        zlim([0 10])
    end
end

%% 2) Generation automatique de donnees d'arpentage avec nb de points
% Creation le 09/11/2016
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
            NbP = [1,5 ,7, 10,15,20];
            NbPoints=NbP(1,NbPoints);
            
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
% Creation le 18/11/2016
% Modification le 13/02/2017
actif=1;
%count=0
for o=1:1
    if actif==1
        fig = figure(indfig);
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

%% 4) Visualisation of the 3D meshing
% Creation le 18/11/2016
actif=1;
%figure('Color', [ 1 1 1])
        hold on
load lysi_adate.mat
Data_Maillage = [lysi_adate] ;
%Data_Maillage = [Maillage3DSable_Anortho] ;
for o=1:1
    if actif==1
        for i=1:length(Data_Maillage(:,1))/3
            X = Data_Maillage((i-1)*3+1:(i-1)*3+3,1);
            Y = Data_Maillage((i-1)*3+1:(i-1)*3+3,2);
            Z = Data_Maillage((i-1)*3+1:(i-1)*3+3,3);
            C = [0.5000 0.5 0.5 ];
            fill3(X,Y,Z,C,'FaceAlpha',0.1);
            daspect([1 1 1])
        end
        view(0,0)
        xlim([0 110])
        ylim([0 50])
        zlim([0 20])
        
        ax = gca;
        ax.FontSize = 13;
    end
end

%% 5) Détermination des coordonnées des GS3 et MPS en surface et à la base
% Creation le 25/01/2017
% Modification le 26/01/2017
actif=0;
for o=1:1
    if actif==1
        
        for o=1:1   % Détermination des centres des lysimètres à la surface de la halde
            % clear Coord
            % figure(1)
            % dcm_obj = datacursormode(fig);
            %
            % c_info = getCursorInfo(dcm_obj);
            % if length(c_info)>0
            %     for j=1:length(c_info)
            %         Coord(j,1:3) = c_info(j).Position;
            %     end
            %     for i=1:3
            %         CENTRE_SURFACE(lysi,i)=mean(Coord(:,i));
            %     end
            %     hold on
            %     plot3(CENTRE_SURFACE(lysi,1),CENTRE_SURFACE(lysi,2),CENTRE_SURFACE(lysi,3),'r.','MarkerSize',20)
            % lysi=lysi+1;
            % end
        end
        
        for o=1:1   % Détermination des coordonnées des sondes à la base de la halde
            % clear Coord
            % figure(1)
            % dcm_obj = datacursormode(fig);
            %
            % c_info = getCursorInfo(dcm_obj);
            % if length(c_info)>0
            %     for j=1:length(c_info)
            %         Coord(j,1:3) = c_info(j).Position;
            %     end
            %
            %     for j=1:length(c_info)/2
            %         for i=1:3
            %             COORD_BASE(j,i)=mean(Coord((j-1)*2+1:j*2,i));
            %         end
            %     end
            %
            %     hold on
            %     plot3(COORD_BASE(:,1),COORD_BASE(:,2),COORD_BASE(:,3),'r.','MarkerSize',20)
            % end
        end
        
        load CENTRE_SURFACE.mat
        load COORD_BASE.mat
        
        % SONDES GS3
        % Coordonnées du premier niveau d'instrumentation GS3
        GS3{1}(:,:)=CENTRE_SURFACE;
        GS3{1}(:,3)=GS3{1}(:,3)-0.25;
        
        % Coordonnées du deuxième niveau d'instrumentation GS3
        GS3{2}(:,:)=CENTRE_SURFACE;
        GS3{2}(:,3)=GS3{2}(:,3)-0.55;
        
        % Coordonnées du troisième niveau d'instrumentation GS3
        GS3{3}(:,:)=CENTRE_SURFACE;
        GS3{3}(:,3)=GS3{3}(:,3)-0.90;
        
        % Coordonnées du quatrième niveau d'instrumentation GS3
        GS3{4}(:,:)=CENTRE_SURFACE;
        GS3{4}(:,3)=GS3{4}(:,3)-1.2;
        
        % Coordonnées du cinquième niveau d'instrumentation GS3 (base de la halde)
        GS3{5}(:,:)=COORD_BASE;
        
        
        
        % SONDES MPS
        MPS{1}=GS3{2};
        MPS{2}=GS3{3};
        MPS{3}=GS3{5};
        % Coordonnées du premier niveau d'instrumentation MPS
        MPS{1}(:,1)=MPS{1}(:,1)+0.2;
        
        % Coordonnées du deuxième niveau d'instrumentation MPS
        MPS{2}(:,1)=MPS{2}(:,1)+0.2;
        
        % Coordonnées du troisième niveau d'instrumentation MPS (base de la halde)
        MPS{3}(:,2)=MPS{3}(:,2)+0.2;
        
        
        % Affichage des GS3 et MPS
        for i=1:5
            plot3(GS3{i}(:,1),GS3{i}(:,2),GS3{i}(:,3),'g.','MarkerSize',16);
        end
        for i=1:3
            plot3(MPS{i}(:,1),MPS{i}(:,2),MPS{i}(:,3),'r.','MarkerSize',16);
        end
    end
end

%% 6) Calcul de l'aire de la halde
% Creation le 07/02/2017
actif=0;
for o=1:1
    if actif==1
        clear Coord
        figure(1)
        dcm_obj = datacursormode(fig);
        
        c_info = getCursorInfo(dcm_obj);
        if length(c_info)>0
            for j=1:length(c_info)
                Coord(j,1:3) = c_info(j).Position;
            end
            
            hold on
            plot3(Coord(:,1),Coord(:,2),Coord(:,3),'r.','MarkerSize',20)
            for i=2:length(Coord(:,1))
                line(Coord(i-1:i,1),Coord(i-1:i,2),Coord(i-1:i,3),'Color','r','LineWidth',4)
            end
            line([Coord(1,1),Coord(end,1)],[Coord(1,2),Coord(end,2)],[Coord(1,3),Coord(end,3)],'Color','r','LineWidth',4)
            AIRE=polyarea(Coord(:,1),Coord(:,2))/(cos(atan(5/100)))
        end
    end
end

%% 8) Generation automatique de donnees d'arpentage avec pas constant
% Creation le 09/11/2016
actif=0;
Pas=0.5;
for o=1:1
    if actif==1
        clear Coord
        fig = figure(indfig);
        dcm_obj = datacursormode(fig);
        
        c_info = getCursorInfo(dcm_obj);
        for i=1:length(c_info)
            Coord(i,1:3) = c_info(i).Position;
        end
        
        %Nombre de points a aujouter
        dist=pdist(Coord(:,:));
        NbPoints=round(dist/Pas,0);
        
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
        Interp=[Interp;Interp2;Coord];
        plot3(Interp(:,1),Interp(:,2),Interp(:,3),'r.','MarkerSize',10)
        daspect([1,1,0.5])
        %view(101.34,20.48)
    end
end



