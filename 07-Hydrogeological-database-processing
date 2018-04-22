%% ------------------- Hydrogeological-database-processing ----------------
%
% Adrien Dimech - Master Project - 22/04/2018
%
% -------------------------------------------------------------------------
% Matlab codes to prepare - process - visualize and interpret 3D time-lapse
% geolelectrical monitoring of a waste rock pile.
% -------------------------------------------------------------------------
%
% This Matlab code was used to load and process the hydrogeological
% database measured in the experimental waste rock pile. Both GS3 and MPS
% sensors are used (Decagon Devices) to monitor water content, temperature,
% resistivity and suction in the waste rock. Precipitation data are also
% processed. This code provide useful tools to handle complexe database and
% to visualize different datasets. Hydrogeological data are then used as a
% validation of geoelectrical results.
%
% Feel free to visit : https://www.researchgate.net/profile/Adrien_Dimech
% for more information about my research or contact me for more information
% and data files : adrien.dimech@gmail.com
%
%% CODE DE PRISE EN MAIN DES DONNÉEES HYDROGEOLOGIQUES
% -------------------------------------------------------------------------
%
%
%                           Suivi du code
% -------------------------------------------------------------------------
% Creation          |       20/01/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       22/01/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       24/01/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       30/01/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       28/02/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       02/03/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       24/03/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       29/03/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       27/07/2017        |     Adrien Dimech

%% 0.1) Chargement des donnees
% Creation le 20/01/2017
% Modification le 31/01/2017
for o=1:1
    % Chargement des données
    figure(1)
    p = get(gcf,'Position');
    set(0,'DefaultFigurePosition',p);
    close all
    clear all
    
    load INFO
    
    T_Debut_Statique=4.263248888888889e+04-1;
    T_Fin_Statique=42635.5881944444-1;
    
    T_Debut_Terrain=42626;
    T_Fin_Terrain=42643;
    
    T1 = 42636.41667;
    T2 = 42640.54167;
    
    formatOut1 = 'dd/mm';
    formatOut2 = 'dd/mm/yy';
    formatOut3 = 'mm/yy';
    
    FichierExcel = 'DATA_HYDROGEO.xlsx';
    [~,b]=xlsfinfo(FichierExcel);
    NbPages = length(b(1,:));
    
    for o=1:1 % Code d'écriture de la matrice DATAH
        % Enregistrement des donnees dans DATAH
        %         for i=1:NbPages
        %             DATAH{i} = xlsread(FichierExcel,i);
        %             round(i/NbPages*100)
        %         end
        load 'DATAH.mat'
    end
end

%% 0.2) Chargement des donnees METEO 2015-2016
% Creation le 30/01/2017
% Modification le 27/07/2017
for o=1:1
    % DATA_METEO{1} = 24_heures
    %   DATA_METEO{1}(:,1) = Heures du 02/01/2015 au 31/10/2016
    %   DATA_METEO{1}(:,2) = T moy °C
    %   DATA_METEO{1}(:,3) = T min °C
    %   DATA_METEO{1}(:,4) = T max °C
    %   DATA_METEO{1}(:,5) = Humidité
    %   DATA_METEO{1}(:,6) = Rafales
    %   DATA_METEO{1}(:,7) = Précipitation sommée en mm
    
    % DATA_METEO{1} = xlsread('METEO_2015_2016.xlsx',1);
    
    % DATA_METEO{2} = 60_minutes
    %   DATA_METEO{2}(:,1) = Heures du 02/01/2015 au 31/10/2016
    %   DATA_METEO{2}(:,2) = T moy °C
    %   DATA_METEO{2}(:,3) = Humidité
    %   DATA_METEO{2}(:,4) = T max °C
    %   DATA_METEO{2}(:,5) = T min °C
    %   DATA_METEO{2}(:,6) = Rafales
    %   DATA_METEO{2}(:,7) = Humidité
    %   DATA_METEO{2}(:,8) = Précipitation sommée en mm
    
    % DATA_METEO{2} = xlsread('METEO_2015_2016.xlsx',2);
    
    
    DATA_METEO = xlsread('DATA_METEO.xlsx',1);
    
    % Création de DATA_METEO_2 : données météo utiles
    DATA_METEO_2(:,1)=DATA_METEO(:,1);
    DATA_METEO_2(:,2)=DATA_METEO(:,8);
end

%% 1.1) Concaténation des données
% Creation le 20/01/2017
% Modification le 30/01/2017
for o=1:1 % Code d'écriture de la matrice DATAH2 : concaténation temporelle
    for o=1:1
        %             DATAH2 = DATAH{1,1};
        %             for i=2:NbPages
        %                 dates=DATAH2(:,1);
        %                 L_DATAH2=length(DATAH2(1,:));
        %                 L_DATAHi=length(DATAH{1,i}(1,:));
        %                 for j=1:length(DATAH{1,i}(:,1))
        %                     if isempty(find(DATAH{1,i}(j,1)==dates(:,1)))==0
        %                         DATAH2(find(DATAH{1,i}(j,1)==dates(:,1)),L_DATAH2+1:L_DATAH2+L_DATAHi-1)=DATAH{1,i}(j,2:L_DATAHi);
        %                     else
        %                         DATAH2(end+1,1)=DATAH{1,i}(j,1);
        %                         DATAH2(end,L_DATAH2+1:L_DATAH2+L_DATAHi-1)=DATAH{1,i}(j,2:L_DATAHi);
        %                     end
        %                 end
        %                 round(i/NbPages*100)
        %             end
        %
        %             DATAH2=sortrows(DATAH2,1);
        %             DATAH2(DATAH2==0)=NaN;
    end
    
    load 'DATAH2.mat'
end

%% 1.1.1) Couverture des données
% Creation le 27/01/2017
% Modification le 27/07/2017
actif=0;
pas_de_temps=30;
format=formatOut2;
for o=1:1
    if actif==1
        close all
        
        LYSI_INFO{1}=['Lysimètre 1 : '];
        LYSI_INFO{2}=['Lysimètre 2 : '];
        LYSI_INFO{3}=['Lysimètre 3 : '];
        LYSI_INFO{4}=['Lysimètre 4 : '];
        LYSI_INFO{5}=['Lysimètre 5 : '];
        LYSI_INFO{6}=['Lysimètre 6 : '];
        
        INFO_SONDE{1}=['GS3 - Anorthosite Haut, Z = 0.15 m'];
        INFO_SONDE{2}=['GS3 - Sable Haut, Z = 0.6 m'];
        INFO_SONDE{3}=['GS3 - Sable Bas, Z = 0.9 m'];
        INFO_SONDE{4}=['GS3 - Ilmenite Haut, Z = 1.2 m'];
        INFO_SONDE{5}=['GS3 - Ilmenite Bas, Z = 1.5 m'];
        INFO_SONDE{6}=['MPS - Sable Bas, Z = 0.9 m'];
        INFO_SONDE{7}=['MPS - Ilmenite Haut, Z = 1.2 m'];
        INFO_SONDE{8}=['GS3 - Sable Base N1, Z = 7 m'];
        INFO_SONDE{9}=['GS3 - Sable Base N2, Z = 7 m'];
        INFO_SONDE{10}=['MPS - Sable Base N1, Z = 7 m'];
        INFO_SONDE{11}=['MPS - Sable Base N2, Z = 7 m'];
        INFO_SONDE{12}=['Debitmetre'];
        
        % ordre des données pour un Lysi dans la matrice concaténée DATAH2
        Column_DATA=[13,10,7,4,1,18,16,23,20,28,26,30];
        for j=1:6
            figure('Color', [ 1 1 1])
            line([T_Debut_Terrain T_Debut_Terrain],[0 12],'LineWidth',2,'Color','k','DisplayName','Début du terrain sept 2016')
            line([T_Fin_Terrain T_Fin_Terrain],[0 12],'LineWidth',2,'Color','k','DisplayName','Fin du terrain sept 2016')
            
            T_Debut_Terrain_N1=4.251492708333262e+04;
            T_Fin_Terrain_N1=4.252897916666666e+04;
            
            line([T_Debut_Terrain_N1 T_Debut_Terrain_N1],[0 12],'LineWidth',2,'Color','r','DisplayName','Début du terrain mai 2016')
            line([T_Fin_Terrain_N1 T_Fin_Terrain_N1],[0 12],'LineWidth',2,'Color','r','DisplayName','Fin du terrain mai 2016')
            
            line([T1 T1],[0 12],'LineWidth',2,'Color','b','DisplayName','Essai Infiltration')
            line([T2 T2],[0 12],'LineWidth',2,'Color','b','DisplayName','Essai Infiltration')
            
            hold on
            for i=1:length(Column_DATA)
                clear ISNAN
                ISNAN(:,1)=double(isnan(DATAH2(:,1+27*(j-1)+Column_DATA(i))));
                ISNAN(:,1)=(ISNAN(:,1)-1)*(-1);
                plot(DATAH2(:,1),ISNAN(:,1)*i,'o','DisplayName',[LYSI_INFO{j},INFO_SONDE{i}])
                hold on
            end
            legend show
            grid
            ax = gca;
            ax.XTick=[pas_de_temps*unique(round(DATAH2(:,1)/pas_de_temps))];
            ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATAH2(:,1)/pas_de_temps)),format)];
            
            % TOUTE LA PERIODE DE DONNEES
            %xlim([4.219238541666666e+04-1 4.267141666666666e+04+1])
            
            %             % ZOOM SUR DEBITMETRE
            %             xlim([4.253145833333334e+04-1 4.264250000000000e+04+1])
            
            ylim([0.5 11.5])
        end
    end
end

%% 1.2) Mise en forme des données
% Creation le 24/03/2017
% Modification le 27/07/2017
for o=1:1
    %     for o=1:1
    %             % INTERPOLATION DES DONNEES
    %             DATAH2_intra=DATAH2;
    %             for col=2:length(DATAH2_intra(1,:))
    %                 for i=1:length(DATAH2_intra(:,1))-1
    %                     if isnan(DATAH2_intra(i,col))
    %                         ISNAN=isnan(DATAH2_intra(i:end,col));
    %                         clear trouve
    %                         trouve=find(ISNAN==0);
    %                         if isempty(trouve)==0
    %                             VALUE_APRES=DATAH2_intra(trouve(1)+i-1,col);
    %                             if i>1
    %                                 VALUE_AVANT=DATAH2_intra(i-1,col);
    %                             else
    %                                 VALUE_AVANT=VALUE_APRES;
    %                             end
    %                             TAILLE=trouve(1)-1;
    %                             for j=1:TAILLE
    %                                 DATAH2_intra(i+j-1,col)=VALUE_AVANT+(VALUE_APRES-VALUE_AVANT)/(TAILLE+1);
    %                             end
    %                         end
    %                     end
    %                 end
    %             end
    %
    %             % SUPPRESSION DES DONNEES EN TROP
    %             limite_ecart=2/24;
    %             DATAH2_intra_simpl=DATAH2_intra(1,:);
    %             for i=2:length(DATAH2_intra(:,1))
    %                 if DATAH2_intra(i,1)-DATAH2_intra_simpl(end,1)>limite_ecart
    %                     DATAH2_intra_simpl(end+1,:)=DATAH2_intra(i,:);
    %                 end
    %             end
    %
    %             DATAH3=DATAH2_intra_simpl;
    %     end
    %     DATAH3(:,1)=DATAH3(:,1)-1.21;
    load DATAH3.mat
end

%% 1.3) Mise en forme des données DATAH3
% Creation le 24/03/2017
% Modification le 37/03/2017
% Modification le 27/07/2017
for o=1:1
    % Matrices d'informations
    LYSI_INFO{1}=['Lysimètre 1'];
    LYSI_INFO{2}=['Lysimètre 2'];
    LYSI_INFO{3}=['Lysimètre 3'];
    LYSI_INFO{4}=['Lysimètre 4'];
    LYSI_INFO{5}=['Lysimètre 5'];
    LYSI_INFO{6}=['Lysimètre 6'];
    
    NOM_SONDE{1}=['Anorthosite Haut, Z = 0.15 m'];
    NOM_SONDE{2}=['Sable Haut, Z = 0.6 m'];
    NOM_SONDE{3}=['Sable Bas, Z = 0.9 m'];
    NOM_SONDE{4}=['Anorthosite Haut, Z = 1.2 m'];
    NOM_SONDE{5}=['Ilmenite Bas, Z = 1.5 m'];
    NOM_SONDE{6}=['Sable Base N1, Z = 7 m'];
    NOM_SONDE{7}=['Sable Base N2, Z = 7 m'];
    NOM_SONDE{8}=['Sable Bas, Z = 0.9 m'];
    NOM_SONDE{9}=['Ilmenite Haut, Z = 1.2 m'];
    NOM_SONDE{10}=['Sable Base N1, Z = 7 m'];
    NOM_SONDE{11}=['Sable Base N2, Z = 7 m'];
    
    % Donnees de conductivite
    % position_conduct=[12;9;6;3;19;22];
    clear CONDUCTIVITE
    position_conduct=[15;12;9;6;3;22;25];
    for lysi=1:6
        for indZ=1:length(position_conduct)
            CONDUCTIVITE{indZ,lysi}(:,1)=DATAH3(:,1);
            CONDUCTIVITE{indZ,lysi}(:,2)=DATAH3(:,1+(lysi-1)*30+position_conduct(indZ));
        end
    end
    
    % Donnees de résistivité
    clear RESISTIVITE
    for lysi=1:6
        for indZ=1:length(position_conduct)
            RESISTIVITE{indZ,lysi}(:,1)=CONDUCTIVITE{indZ,lysi}(:,1);
            RESISTIVITE{indZ,lysi}(:,2)=(1./(CONDUCTIVITE{indZ,lysi}(:,2)/10));
        end
    end
    
    % Donnees de température
    % position_temperature=[11;8;5;2;18;21;16;14;24;26];
    position_temperature=[14;11;8;5;2;21;23;19;17;27;29];
    clear TEMPERATURE
    for lysi=1:6
        for indZ=1:length(position_temperature)
            TEMPERATURE{indZ,lysi}(:,1)=DATAH3(:,1);
            TEMPERATURE{indZ,lysi}(:,2)=DATAH3(:,1+(lysi-1)*30+position_temperature(indZ));
        end
    end
    
    %     % Donnees de succion
    %     % position_succion=[0;15;13;0;23;25];
    %     position_succion=[0;15;13;0;23;25];
    %     for lysi=1:6
    %         for indZ=1:6
    %             if position_succion(indZ)==0
    %                 SUCCION{indZ,lysi}(1,1)=NaN;
    %             else
    %             SUCCION{indZ,lysi}(:,1)=DATAH3(:,1);
    %             SUCCION{indZ,lysi}(:,2)=DATAH3(:,1+(lysi-1)*27+position_succion(indZ));
    %             end
    %         end
    %     end
    
    % Donnees de Débitmetre
    clear DEBITMETRE
    for lysi=1:6
        for indZ=1:1
            DEBITMETRE{indZ,lysi}(:,1)=DATAH3(:,1);
            DEBITMETRE{indZ,lysi}(:,2)=DATAH3(:,1+(lysi-1)*30+30);
        end
    end
    
end

%% 1.2.1) Affichage des valeurs de conductivité lors de l'imagerie statique
% Création le 28/02/2017
for o=1
    STATIQUE_ESSAI=[25,28,30,33,46,52,27,29,32,47,48,49,50,51,34,37,38,39,40,41,42,44,45];
    for i=1:length(STATIQUE_ESSAI(1,:))
        DATES_ESSAI(i,1)=STATIQUE_ESSAI(1,i);
        DATES_ESSAI(i,2)=INFO(find(INFO(:,1)==STATIQUE_ESSAI(1,i)),4)-1;
        DATES_ESSAI(i,3)=INFO(find(INFO(:,1)==STATIQUE_ESSAI(1,i)),6)+DATES_ESSAI(i,2);
        DATES_ESSAI(i,4)=(DATES_ESSAI(i,2)+DATES_ESSAI(i,3))/2;
    end
end

%% 1.2.2) Affichage des donnees pour chaque sonde = FIGURE SAGEEP
% Creation le 20/01/2017
% Modification le 24/01/2017
% Modification le 28/02/2017
% Modification le 24/03/2017
% Modification le 27/07/2017
for o=1
    DATA_TO_PLOT=RESISTIVITE;
    meteo=0;
    essai_infilt=0;
    
    pas_de_temps=1;
    format=formatOut1;
    delay=5;
    
    actif=1;
    close all
    for o=1:1
        if actif==1
            VALEURS_UNIQUES_COND=0;
            Z=[0.15,0.6,0.9,1.2,1.5];
            for lysi=1:6
                for sonde=2
                    figure('Color', [ 1 1 1])
                    hold on
                    plot(DATA_TO_PLOT{sonde,lysi}(:,1),DATA_TO_PLOT{sonde,lysi}(:,2),'.-','LineWidth',4,'MarkerSize',10)
                    
                    % CHOIX DE LA PERIODE D'INTERET
                    T_Debut=T_Debut_Terrain-delay;
                    T_Fin=T_Fin_Terrain+delay;
                    % xlim([T_Debut T_Fin])
                    %xlim([42510 42560])
                    % PRÉCIPITATION MONTRÉE DANS SAGEEP
                    %xlim([42602 42609])
                    % TERRAIN N3
                    xlim([42890 42901])
                    
                    
                    xlabel({'Elapsed time (dd/mm/yy)'},'FontSize',14,'FontWeight','bold')
                    
                    if essai_infilt==1
                        clear MIN
                        clear MAX
                        MIN=min(DATA_TO_PLOT{sonde,lysi}(find(DATA_TO_PLOT{sonde,lysi}(:,1)>T_Debut&DATA_TO_PLOT{sonde,lysi}(:,1)<T_Fin),2));
                        MAX=max(DATA_TO_PLOT{sonde,lysi}(find(DATA_TO_PLOT{sonde,lysi}(:,1)>T_Debut&DATA_TO_PLOT{sonde,lysi}(:,1)<T_Fin),2));
                        
                        line([T1 T1],[MIN MAX],'LineWidth',2,'Color','b','DisplayName','Essai Infiltration N1')
                        line([T2 T2],[MIN MAX],'LineWidth',2,'Color','b','DisplayName','Essai Infiltration N2')
                    end
                    
                    title([LYSI_INFO{lysi},' : ',NOM_SONDE{sonde}],'FontSize',12,'FontWeight','bold' )
                    
                    if DATA_TO_PLOT{1}==RESISTIVITE{1}
                        ylabel({'Resistiviy';'(Ohm.m)'},'FontSize',14,'FontWeight','bold')
                    else
                        ylabel({'Conductivite';'(mS/cm)'},'FontSize',14,'FontWeight','bold')
                    end
                    
                    if meteo==1
                        % Affichage des précipitations
                        yyaxis right
                        plot(DATA_METEO_2(:,1),DATA_METEO_2(:,2),'.--','LineWidth',2,'MarkerSize',20)
                        ylabel({'Precipitation every hour';'(mm)'},'FontSize',14,'FontWeight','bold')
                    end
                    
                    grid
                    ax = gca;
                    ax.XTick=[pas_de_temps*unique(round(DATA_TO_PLOT{1,1}(:,1)/pas_de_temps))];
                    ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATA_TO_PLOT{1,1}(:,1)/pas_de_temps)),format)];
                    
                    
                end
            end
        end
    end
end

%% 1.2.3) Création de la matrice des données de résistivité GS3 pour comparaison E4D
% Creation le 27/07/2017
for o=1
    load GS3_location.mat
    GS3_DATA_mes=GS3_location;
    
    T_MIN_2017=42890;
    T_MAX_2017=42901;
    Ligne_DATAH3=find(DATAH3(:,1)>T_MIN_2017&DATAH3(:,1)<T_MAX_2017);
    T_GS3=DATAH3(Ligne_DATAH3,1)';
    indCol=[5,4,3,2,1,6,7];
    for lysi=1:6
        for sonde=1:7
            ind_ligne=(lysi-1)*7+indCol(sonde);
            GS3_DATA_mes(ind_ligne,5:length(T_GS3)+4)=RESISTIVITE{sonde,lysi}(Ligne_DATAH3,2)';
        end
    end
end

%% 1.3) Affichage des donnees de température pour chaque sonde
% Creation le 30/01/2017
actif=0;
%Name_Probe=NOM_SONDE;
for o=1:1
    if actif==1
        figure('Color', [ 1 1 1])
        hold on
        for i=1:length(Name_Probe)
            if mod(i,6)==1
                plot(DATA_TEMP_GS3(:,1),DATA_TEMP_GS3(:,i+1),'.','DisplayName',Name_Probe{i})
            end
        end
        title(['Évolution de la température au cours du temps'])
        xlabel('Date de l"evenement en j')
        ylabel('Température en °C')
        
        legend show
        grid
        ax = gca;
        pas_de_temps=1;
        ax.XTick=[pas_de_temps*unique(round(DATA_TEMP_GS3(:,1)/pas_de_temps))];
        ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATA_TEMP_GS3(:,1)/pas_de_temps)),formatOut2)];
    end
end

%% 2.1) Identification des évenements d'infiltration
% Creation le 23/01/2017
actif=0;
for o=1:1   % code de récupération des dates des évenements d'infiltration
    if actif==1
        for i=1:36
            clear Coord
            fig = figure(i);
            figure(i)
            dcm_obj = datacursormode(fig);
            
            c_info = getCursorInfo(dcm_obj);
            if length(c_info)>0
                for j=1:length(c_info)
                    Coord(j,1:2) = c_info(j).Position;
                end
                Evenement{i}=Coord;
            else
                Evenement{i}=0;
            end
        end
    end
end

%% 2.2.1) Affichage des précipitations
% Creation le 23/01/2017
actif=0;
for o=1:1
    if actif==1
        figure('Color', [ 1 1 1])
        hold on
        plot(DATA_METEO_2(:,1),DATA_METEO_2(:,2),'o-','DisplayName','Précipitations cumulées aux 60 min (en mm)')
        legend('show')
        grid
        ax = gca;
        ax.XTick=[unique(round(DATA_METEO_2(:,1)))];
        ax.XTickLabel =[datestr(unique(round(DATA_METEO_2(:,1))),formatOut2)];
        title(['Précipitations sur le site cumulées aux 60 min en mm'])
        xlabel('Date de l"evenement de précipitations en j')
        ylabel('Précipitation en mm')
        xlim([T_Debut_Terrain T_Fin_Terrain])
        hold off
    end
end

%% 2.2.2) Donnees METEO 2015-2016 : Calcul de l'amplitude de précipitations particulières en mm
% Creation le 30/01/2017
actif=0; % Calcul de l'amplitude de la précipitation en mm
for o=1:1
    if actif==1
        clear x y Precipitation
        [x,y] = ginput(2);
        [~,t1]=min(abs(DATA_METEO_2(:,1)-x(1,1)));
        [~,t2]=min(abs(DATA_METEO_2(:,1)-x(2,1)));
        Precipitation = sum(DATA_METEO_2(t1:t2,2));
    end
end

%% 2.2.3) Superposition donnees meteo + donnees GS3
% Creation le 30/01/2017
actif=0;

for o=1:1
    if actif==1
        DATA_TO_PLOT=DATA_RESIST_simpl;
        figure('Color', [ 1 1 1])
        hold on
        title(['Précipitations sur le site cumulées aux 60 min en mm et resistivités en Ohms.m'])
        xlabel('Date de l"evenement de précipitations en j')
        
        % Affichage des résistivités
        for i=1:6
            for j=4
                Nprobe=(i-1)*6+j;
                plot(DATA_TO_PLOT(:,1),DATA_TO_PLOT(:,Nprobe+1),'.-','DisplayName',Name_Probe{Nprobe})
            end
        end
        ylabel('Résistivité en Ohms.m')
        
        % Affichage des précipitations
        yyaxis right
        plot(DATA_METEO_2(:,1),DATA_METEO_2(:,2),'.-k','DisplayName','Précipitations cumulées aux 60 min (en mm)')
        ylabel('Précipitation en mm')
        
        
        
        legend('show')
        grid
        ax = gca;
        ax.XTick=[unique(round(DATA_METEO_2(:,1)))];
        ax.XTickLabel =[datestr(unique(round(DATA_METEO_2(:,1))),formatOut2)];
        
        xlim([T_Debut_Terrain T_Fin_Terrain])
        hold off
    end
end

%% 2.3) Affichage des temps d'infiltration
% Creation le 23/01/2017
actifGlob=0;
for o=1:1
    if actifGlob==1
        actif=0;
        for o=1:1 % Affichage des sables haut et bas pour tous les lysimètres
            if actif==1
                figure('Color', [ 1 1 1])
                hold on
                count=0;
                for i=1:6
                    for j=1:2
                        count=count+1;
                        Nprobe=(i-1)*6+j+2;
                        MAT=ones(length(Evenement{Nprobe}(:,1)));
                        plot(Evenement{Nprobe}(:,1),count*MAT(:,1),'o','DisplayName',Name_Probe{Nprobe})
                    end
                end
            end
        end
        
        actif=0;
        for o=1:1 % Affichage des sables et ilmenite haut et bas pour chaque lysimètre
            if actif==1
                for i=1:6
                    figure('Color', [ 1 1 1])
                    hold on
                    count=0;
                    for j=1:4
                        count=count+1;
                        Nprobe=(i-1)*6+j;
                        MAT=ones(length(Evenement{Nprobe}(:,1)));
                        plot(Evenement{Nprobe}(:,1),count*MAT(:,1),'o','DisplayName',Name_Probe{Nprobe})
                    end
                    legend('show')
                    for j=1:22
                        line([EVENT(j,1) EVENT(j,1)],[0 4],'Color','r');
                    end
                end
            end
        end
        
        actif=0;
        if actif==1
            for i=1:1000
                [x,y] = ginput(1);
                EVENT(i,1)=x;
                line([x x],[0 1]);
                ylim([0 1])
            end
        end
    end
end

%% 2.4) Affichage des temps d'infiltration V2 entre sable surface et base
% Creation le 23/01/2017
actif=0;

meteo=1;
for o=1:1
    if actif==1
        DATA_TO_PLOT=DATA_RESIST_simpl;
        for i=1:6
            figure('Color', [ 1 1 1])
            hold on
            count=0;
            for j=1:6
                count=count+1;
                Nprobe=(i-1)*6+j;
                % Normalisation
                clear g
                MIN=min(DATA_TO_PLOT(:,Nprobe+1));
                g=DATA_TO_PLOT(:,Nprobe+1)-MIN;
                MAX=max(g);
                g=g/MAX;
                plot(DATA_TO_PLOT(:,1),g,'o-','DisplayName',Name_Probe{Nprobe})
            end
            if meteo==1
                % Affichage des précipitations
                yyaxis right
                plot(DATA_METEO_2(:,1),DATA_METEO_2(:,2),'.-','DisplayName','Précipitations cumulées aux 60 min (en mm)')
                ylabel('Précipitation en mm')
            end
            legend('show')
            
            grid
            ax = gca;
            pas_de_temps=1;
            ax.XTick=[pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps))];
            ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps)),formatOut2)];
        end
    end
end

%% 2.5) Normalisation locale des évenements d'infiltration
% Creation le 24/01/2017
for o=1
    actif =0;
    
    meteo=1;
    for o=1:1
        if actif==1
            DATA_TO_PLOT=DATA_RESIST_simpl;
            LYSI = 2;
            
            clear EVENT
            clear data
            for i=1:2
                [x,y] = ginput(1);
                EVENT(i,1)=x;
                line([x x],[0 1]);
                ylim([0 1])
            end
            figure('Color', [ 1 1 1])
            data = DATA_TO_PLOT(find(DATA_TO_PLOT(:,1)>min(EVENT(:,1)) & DATA_TO_PLOT(:,1)<max(EVENT(:,1))),:);
            for i=1:6
                Nprobe=(LYSI-1)*6+i;
                g = data(:,Nprobe+1);
                g=g-min(g);
                g=g/max(g);
                plot(data(:,1),g,'.-','DisplayName',Name_Probe{Nprobe})
                hold on
            end
            
            
            grid
            ax = gca;
            pas_de_temps=1;
            ax.XTick=[pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps))];
            ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps)),formatOut2)];
            
            title(['Evenement d"infiltration particulier et réponse des capteurs au LYSI n°',num2str(LYSI)])
            xlabel('Date de l"evenement d"infiltration en j')
            ylabel('Evolution relative des valeurs de résistivités')
            xlim([min(EVENT(:,1)) max(EVENT(:,1))])
            %ylim([0 1])
            
            if meteo==1
                % Affichage des précipitations
                yyaxis right
                plot(DATA_METEO_2(:,1),DATA_METEO_2(:,2),'.-','DisplayName','Précipitations cumulées aux 60 min (en mm)')
                ylabel('Précipitation en mm')
            end
            legend show
            hold off
        end
    end
end

%% 3) Traitement des donnees d'infiltration
for o=1
    % Creation le 22/01/2017
    load DATA_INFILT1.mat
end

%% 5) Test pour retrouver conductivité électrique de l'eau dans les pores
for o=1
    p = get(gcf,'Position');
    set(0,'DefaultFigurePosition',p);
    actif=1;
    resistivite=1;
    meteo=1;
    load ARROSAGE.mat
    for o=1:1
        if actif==1
            % Donnees de conductivite Sbulk (mS/cm)
            clear Sbulk
            position_conduct=[15;12;9;6;3;22;25];
            for lysi=1:6
                for indZ=1:length(position_conduct)
                    Sbulk{indZ,lysi}(:,1)=DATAH3(:,1);
                    Sbulk{indZ,lysi}(:,2)=DATAH3(:,1+(lysi-1)*30+position_conduct(indZ));
                end
            end
            
            % Donnees de Volumetric Water Content VWC (m3/m3)
            clear VWC
            position_conduct=[15;12;9;6;3;22;25];
            for lysi=1:6
                for indZ=1:length(position_conduct)
                    VWC{indZ,lysi}(:,1)=DATAH3(:,1);
                    VWC{indZ,lysi}(:,2)=DATAH3(:,1+(lysi-1)*30+position_conduct(indZ)-2);
                end
            end
            
            % Donnees de Temperature du sol Tsoil (°C)
            clear Tsoil
            position_conduct=[15;12;9;6;3;22;25];
            for lysi=1:6
                for indZ=1:length(position_conduct)
                    Tsoil{indZ,lysi}(:,1)=DATAH3(:,1);
                    Tsoil{indZ,lysi}(:,2)=DATAH3(:,1+(lysi-1)*30+position_conduct(indZ)-1);
                end
            end
            
            
            % Calcul de conductivite de l'eau dans les pores Spore (mS/cm)
            load GS3_location.mat
            GS3_DATA_EAU_mes_raw=GS3_location;
            T_GS3_DATA_EAU_mes_raw=Sbulk{1,1}(:,1);
            L_GS3_DATA_EAU_mes_raw=length(T_GS3_DATA_EAU_mes_raw)+4;
            position=[5,4,3,2,1,6,7];
            for lysi=1:6
                for indZ=1:length(position_conduct)
                    ligne_GS3_DATA_EAU_mes=(lysi-1)*7+position(indZ);
                    Spore{indZ,lysi}(:,1)=Sbulk{indZ,lysi}(:,1);
                    Spore{indZ,lysi}(:,2)=((80.3-0.37*(Tsoil{indZ,lysi}(:,2)-20)).*Sbulk{indZ,lysi}(:,2))./(((VWC{indZ,lysi}(:,2)+0.117)/0.118).^2-4.1);
                    GS3_DATA_EAU_mes_raw(ligne_GS3_DATA_EAU_mes,5:L_GS3_DATA_EAU_mes_raw)=Spore{indZ,lysi}(:,2)';
                end
            end
            
            % Affichage
            for o=1
                DATA_TO_PLOT=Spore;
                pas_de_temps=1;
                format=formatOut1;
                
                close all
                VALEURS_UNIQUES_COND=0;
                Z=[0.15,0.6,0.9,1.2,1.5];
                
                NOM_SONDE{4}=['GS3 - Anorthosite Haut, Z = 1.2 m'];
                NOM_SONDE{5}=['GS3 - Anorthosite Bas, Z = 1.5 m'];
                
                for sonde=1:7
                    figure('Color', [ 1 1 1])
                    hold on
                    for lysi=1:6
                        if resistivite==1
                            plot(DATA_TO_PLOT{sonde,lysi}(:,1),10./DATA_TO_PLOT{sonde,lysi}(:,2),'.-','LineWidth',4,'MarkerSize',10,'DisplayName',[LYSI_INFO{lysi},' : ',NOM_SONDE{sonde}])
                        else
                            plot(DATA_TO_PLOT{sonde,lysi}(:,1),1000*DATA_TO_PLOT{sonde,lysi}(:,2),'.-','LineWidth',4,'MarkerSize',10,'DisplayName',[LYSI_INFO{lysi},' : ',NOM_SONDE{sonde}])
                        end
                        
                        % CHOIX DE LA PERIODE D'INTERET
                        % TERRAIN N3
                        xlim([42891 42904])
                        % PRÉCIPITATION MONTRÉE DANS SAGEEP
                        % xlim([42602 42609])
                        
                        xlabel({'Elapsed time (dd/mm/yy)'},'FontSize',12,'FontWeight','bold')
                        
                        
                    end
                    legend show
                    %title([LYSI_INFO{lysi},' : ',NOM_SONDE{sonde}],'FontSize',12,'FontWeight','bold' )
                    title('Résistivité électrique du fluide intersticiel dans le sable (Z = 7 m)','FontSize',14,'FontWeight','bold' )
                    %ylim([2 6])
                    %ylim([2.5 4.5])
                    if resistivite==1
                        ylabel({'Résistivité du fluide intersticiel';'(Ohms/m)'},'FontSize',14,'FontWeight','bold')
                    else
                        ylabel({'Conductivite du fluide intersticiel';'(\muS/cm)'},'FontSize',14,'FontWeight','bold')
                    end
                    
                    grid
                    ax = gca;
                    ax.XTick=[pas_de_temps*unique(round(DATA_TO_PLOT{1,1}(:,1)/pas_de_temps))];
                    ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATA_TO_PLOT{1,1}(:,1)/pas_de_temps)),format)];
                    
                    %                     yyaxis right
                    %                     plot(ARROSAGE(:,1),ARROSAGE(:,2),'.--b','LineWidth',2,'MarkerSize',10)
                    %                     ylabel({'Arrosage de la halde';'(Nb de passage par heure)'},'FontSize',14,'FontWeight','bold','Color','b')
                    %                     yyaxis left
                end
            end
        end
    end
end

%% DONNEES DE TEMPÉRATURE
for o=1
    p = get(gcf,'Position');
    set(0,'DefaultFigurePosition',p);
    % Affichage
    for o=1
        DATA_TO_PLOT=TEMPERATURE;
        pas_de_temps=1;
        format=formatOut1;
        
        close all
        VALEURS_UNIQUES_COND=0;
        Z=[0.15,0.6,0.9,1.2,1.5];
        
        NOM_SONDE{4}=['GS3 - Anorthosite Haut, Z = 1.2 m'];
        NOM_SONDE{5}=['GS3 - Anorthosite Bas, Z = 1.5 m'];
        
        for sonde=1:7
            figure('Color', [ 1 1 1])
            hold on
            for lysi=1:6
                plot(DATA_TO_PLOT{sonde,lysi}(:,1),DATA_TO_PLOT{sonde,lysi}(:,2),'.-','LineWidth',4,'MarkerSize',10,'DisplayName',[LYSI_INFO{lysi},' : ',NOM_SONDE{sonde}])
                
                
                % CHOIX DE LA PERIODE D'INTERET
                % TERRAIN N3
                %xlim([42891 42904])
                % PRÉCIPITATION MONTRÉE DANS SAGEEP
                % xlim([42602 42609])
                
                xlabel({'Elapsed time (dd/mm/yy)'},'FontSize',12,'FontWeight','bold')
                
                
            end
            legend show
            %title([LYSI_INFO{lysi},' : ',NOM_SONDE{sonde}],'FontSize',12,'FontWeight','bold' )
            %title('Température du fluide intersticiel dans le sable (Z = 7 m)','FontSize',14,'FontWeight','bold' )
            %ylim([2 6])
            %ylim([2.5 4.5])
            if resistivite==1
                ylabel({'Résistivité du fluide intersticiel';'(Ohms/m)'},'FontSize',14,'FontWeight','bold')
            else
                ylabel({'Conductivite du fluide intersticiel';'(\muS/cm)'},'FontSize',14,'FontWeight','bold')
            end
            
            grid
            ax = gca;
            ax.XTick=[pas_de_temps*unique(round(DATA_TO_PLOT{1,1}(:,1)/pas_de_temps))];
            ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATA_TO_PLOT{1,1}(:,1)/pas_de_temps)),format)];
            
            %                     yyaxis right
            %                     plot(ARROSAGE(:,1),ARROSAGE(:,2),'.--b','LineWidth',2,'MarkerSize',10)
            %                     ylabel({'Arrosage de la halde';'(Nb de passage par heure)'},'FontSize',14,'FontWeight','bold','Color','b')
            %                     yyaxis left
        end
    end
end

%% 6) Création de la matrice des données de résistivité DE L'EAU GS3 pour comparaison E4D
% Creation le 12/10/2017
for o=1
    load GS3_location.mat
    GS3_DATA_EAU_mes=GS3_DATA_EAU_mes_raw(:,1:4);
    
    T_MIN_2017=42890;
    T_MAX_2017=42901;
    Colonnes_GS3_DATA_EAU_mes_raw=find(T_GS3_DATA_EAU_mes_raw(:,1)>T_MIN_2017&T_GS3_DATA_EAU_mes_raw(:,1)<T_MAX_2017);
    T_GS3_RES_EAU=T_GS3_DATA_EAU_mes_raw(Colonnes_GS3_DATA_EAU_mes_raw,1)';
    
    GS3_DATA_EAU_mes(:,5:length(T_GS3_RES_EAU)+4)=10./GS3_DATA_EAU_mes_raw(:,Colonnes_GS3_DATA_EAU_mes_raw+4);
end

%% 7) Figure pour visualisation données de résistivité en sortie des lysimètres (Bissé 2018)
p = get(gcf,'Position');
set(0,'DefaultFigurePosition',p);
actif=1;
for o=1:1
    if actif==1
        close all
        % Chargement des données
        for i=1:6
            DATA_BISSE{i} = xlsread('DATA_BISSE.xlsx',i);
        end
        load ARROSAGE.mat
        pas_de_temps=1;
        format=formatOut1;
        
        % Affichage données de pH
        figure('Color', [ 1 1 1])
        hold on
        for i=1:6
            DATA_TO_PLOT=[DATA_BISSE{i}(:,1),DATA_BISSE{i}(:,2)];
            plot(DATA_TO_PLOT(:,1),DATA_TO_PLOT(:,2),'.-','LineWidth',3,'MarkerSize',20,'DisplayName',['Lysimètre ',num2str(i)])
        end
        xlim([42891 42901])
        grid
        ax = gca;
        ax.XTick=[pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps))];
        ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps)),format)];
        
        ylabel({'pH mesuré à la sortie des lysimètres';'(-)'},'FontSize',14,'FontWeight','bold')
        legend show
        title('Évolution du pH de l''eau en sortie des lysimètres au cours du temps','FontSize',14,'FontWeight','bold' )
        
        % Affichage données de conductivité
        figure('Color', [ 1 1 1])
        hold on
        for i=1:6
            DATA_TO_PLOT=[DATA_BISSE{i}(:,1),DATA_BISSE{i}(:,3)];
            plot(DATA_TO_PLOT(:,1),DATA_TO_PLOT(:,2),'.-','LineWidth',3,'MarkerSize',20,'DisplayName',['Lysimètre ',num2str(i)])
        end
        xlim([42891 42901])
        grid
        ax = gca;
        ax.XTick=[pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps))];
        ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps)),format)];
        
        ylabel({'Conductivité électrique mesurée à la sortie des lysimètres';'(\muS/cm)'},'FontSize',14,'FontWeight','bold')
        legend show
        title('Évolution de la conductivité électrique de l''eau en sortie des lysimètres au cours du temps','FontSize',14,'FontWeight','bold' )
        
        % Affichage données de résistivité
        figure('Color', [ 1 1 1])
        hold on
        for i=1:6
            DATA_TO_PLOT=[DATA_BISSE{i}(:,1),10000./DATA_BISSE{i}(:,3)];
            plot(DATA_TO_PLOT(:,1),DATA_TO_PLOT(:,2),'.-','LineWidth',3,'MarkerSize',20,'DisplayName',['Lysimètre ',num2str(i)])
        end
        xlim([42891 42901])
        grid
        ax = gca;
        ax.XTick=[pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps))];
        ax.XTickLabel =[datestr(pas_de_temps*unique(round(DATA_TO_PLOT(:,1)/pas_de_temps)),format)];
        
        ylabel({'Résistivité électrique mesurée à la sortie des lysimètres';'(\Omega.m)'},'FontSize',14,'FontWeight','bold')
        legend show
        title('Évolution de la résistivité électrique de l''eau en sortie des lysimètres au cours du temps','FontSize',14,'FontWeight','bold' )
        
        
    end
end
