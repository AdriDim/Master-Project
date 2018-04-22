%% ------------- Optimized-protocols-for-ABEM-Terrameter-LS ---------------
%
% Adrien Dimech - Master Project - 22/04/2018
%
% -------------------------------------------------------------------------
% Matlab codes to prepare - process - visualize and interpret 3D time-lapse
% geolelectrical monitoring of a waste rock pile.
% -------------------------------------------------------------------------
%
% This Matlab code was used to select optimized protocols for any
% Electrical Resistivity Tomography measurements with surface electrodes in
% 1D or 2D grids (regular or not). The optimized protocols are identified
% following the methodology of the 'COMPARE R' method developed by the 
% British Geological Survey (Stummer et al., 2004; Wilkinson et al., 2006b, 
% 2012; Loke et al.,2014a, 2014b, 2015). Optimized protocols provide better
% resolution and coverage compared to standard protocols and they need less
% acquisition time. The Matlab code automatically generates .xml files to
% upload the optimized protocols on the Terramter LS (ABEM) which can be
% completed by one or more ES1064C electrode selecter (ABEM).
%
% Feel free to visit : https://www.researchgate.net/profile/Adrien_Dimech
% for more information about my research or contact me for more information
% and data files : adrien.dimech@gmail.com
%
%%% 3D VERTICAL SCAN : V1 => Maximum resolution 3D imaging
% -------------------------------------------------------------------------
% Automatic generation of optimized 3D ERT protocols with maximal coverage
% and resolution and minimal execution time : integrated tool from the 
% office survey to the field data acquisition.
% Designed for ABEM Terrameter with ES1064C.

%                           Code Activity
% -------------------------------------------------------------------------
% Creation          |       09/05/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       18/05/2017        |     Adrien Dimech

%% I) PARAMETERS DEFINITION
%--------------------------------------------------------------------------
actif_parameters_definition=1;
for o=1
    if actif_parameters_definition==1
        % Creation : 09/05/2017
        actif=1;
        for o=1
            if actif==1
                close all
                clear all
                figure
                p = get(gcf,'Position');
                set(0,'DefaultFigurePosition',p);
                clear p
                close all
                
                % PARAMETERS
                
                % 1) GRID SIZE
                
                % X : number of electrodes (#)
                prompt = {'How many electrodes in the X direction ?'};
                dlg_title = 'Number of electrodes';
                num_lines = 1;
                defaultans = {'0'};
                X_Nb_E = inputdlg(prompt,dlg_title,num_lines,defaultans);
                X_Nb_ELEC = str2num(X_Nb_E{:});
                %X_Nb_ELEC=21;
                % X : grid spacing (m)
                X_GRID_SPACING=2;
                % Y : number of electrodes (#)
                prompt = {'How many electrodes in the Y direction ?'};
                dlg_title = 'Number of electrodes';
                num_lines = 1;
                defaultans = {'0'};
                Y_Nb_E = inputdlg(prompt,dlg_title,num_lines,defaultans);
                Y_Nb_ELEC = str2num(Y_Nb_E{:});
                %Y_Nb_ELEC=1;
                % Y : grid spacing (m)
                Y_GRID_SPACING=1;
                % Position of the first electrode (m)
                X_Elec_1=1;
                Y_Elec_1=1;
                
                % 2) SENSITIVITY CALCULATION GRID SIZE
                
                % X : grid spacing for sensitivity calculation (m)
                % Suggered value : X_grid spacing / 2
                X_SENSI_GRID_SPACING=X_GRID_SPACING/2;
                % Y : grid spacing for sensitivity calculation (m)
                % Suggered value : Y_grid spacing / 2
                Y_SENSI_GRID_SPACING=Y_GRID_SPACING/2;
                % Z : grid spacing for sensitivity calculation (m)
                % Suggered value : min(X_GRID_SPACING,Y_GRID_SPACING) / 4
                Z_SENSI_GRID_SPACING=min(X_GRID_SPACING,Y_GRID_SPACING)/4;
                % Z : Maximum depth for sensitivity calculation (m)
                % Suggered value :
                % max(X_GRID_SPACING*X_Nb_ELEC,Y_GRID_SPACING*Y_Nb_ELEC)/6
                Z_MAX=max(X_GRID_SPACING*X_Nb_ELEC,Y_GRID_SPACING*Y_Nb_ELEC)/8;
                
                % Standard measurement protocols
                standard_exist=1;
                load WENNER_84x1.mat
                load DIPDIP_84x1.mat
                load MULTIGRAD_84x1.mat
                load PRO.mat
                PRO_STAND_84x1=unique([WENNER_84x1;DIPDIP_84x1],'rows');
                
                load STANDARD_21x3.mat
                load STANDARD_21x1.mat
                load STANDARD_9x9.mat
                load STANDARD_21x4.mat
                load STANDARD_84x1.mat
                load STANDARD_42x1.mat
                load STANDARD_21x2.mat
                load STANDARD_8x5.mat
                PRO_STANDARD=PRO;
            end
        end
        clear defaultans dlg_title num_lines prompt X_Nb_E Y_Nb_E 
    end
end

%% II)CONCEPTION OF THE COMPREHENSIVE DATA SET
%--------------------------------------------------------------------------
actif_comprehensive_data_set=1;
for o=1
    if actif_comprehensive_data_set==1
        %% 2.1) Electrode grid definition (numero, coordinates...)
        % Creation : 09/05/2017
        actif=1;
        for o=1
            if actif==1
                ELEC=zeros(X_Nb_ELEC*Y_Nb_ELEC,3);
                ELEC_MATRIX=zeros(X_Nb_ELEC,Y_Nb_ELEC);
                count=1;
                for i=1:X_Nb_ELEC
                    for j=1:Y_Nb_ELEC
                        ELEC_MATRIX(i,j)=count;
                        ELEC(count,1)=(i-1)*X_GRID_SPACING+X_Elec_1;
                        ELEC(count,2)=(j-1)*Y_GRID_SPACING+Y_Elec_1;
                        ELEC(count,3)=0;
                        count=count+1;
                    end
                end
                clear count i j o actif
                
                o_2_1_Electrode_grid_definition=1
                clear o_2_1_Electrode_grid_definition
                
            end
        end
       
        %% 2.2) Alpha and Beta arrays with one line in X direction
        % Creation : 09/05/2017
        actif=1;
        for o=1
            if actif==1 & X_Nb_ELEC>1    
                % Alpha configuration          C1 -- P1 -- P2 -- C2
                clear ALPHA_ARRAYS
                ALPHA_ARRAYS=zeros(X_Nb_ELEC*X_Nb_ELEC*X_Nb_ELEC*X_Nb_ELEC,4);
                count=1;
                for C1=1:X_Nb_ELEC
                    for P1=C1+1:X_Nb_ELEC
                        for P2=P1+1:X_Nb_ELEC
                            for C2=P2+1:X_Nb_ELEC
                                ALPHA_ARRAYS(count,1:4)=[C1,C2,P1,P2];
                                count=count+1;
                            end
                        end
                    end
                end
                ALPHA_ARRAYS(count:end,:)=[];
            
               	% Beta configuration          C2 -- C1 -- P1 -- P2
                clear BETA_ARRAYS
                BETA_ARRAYS=zeros(X_Nb_ELEC*X_Nb_ELEC*X_Nb_ELEC*X_Nb_ELEC,4);
                count=1;
                for C2=1:X_Nb_ELEC
                    for C1=C2+1:X_Nb_ELEC
                        for P1=C1+1:X_Nb_ELEC
                            for P2=P1+1:X_Nb_ELEC
                                BETA_ARRAYS(count,1:4)=[C1,C2,P1,P2];
                                count=count+1;
                            end
                        end
                    end
                end
                BETA_ARRAYS(count:end,:)=[];
                
                % Comprehensive data set with one X line
                clear COMPREHENSIVE_X_1_line
                COMPREHENSIVE_X_1_line=[ALPHA_ARRAYS;BETA_ARRAYS];
                COMPREHENSIVE_X_1_line=unique(COMPREHENSIVE_X_1_line,'rows','stable');
                
                clear C1 C2 P1 P2 ALPHA_ARRAYS BETA_ARRAYS actif o
                
                o_2_2_Alpha_and_Beta_arrays_with_one_line_in_X_direction=1
                clear o_2_2_Alpha_and_Beta_arrays_with_one_line_in_X_direction count
            end
        end
        
        %% 2.3) Horizontal Roll Along in Y direction with all the lines
        % Creation : 09/05/2017
        actif=1;
        for o=1
            if actif==1 & X_Nb_ELEC>1           
                % Roll Along Y Matrix
                clear ROLL_ALONG_Y
                 ROLL_ALONG_Y=zeros(Y_Nb_ELEC*Y_Nb_ELEC*Y_Nb_ELEC*Y_Nb_ELEC,4);
                count=1;
                for i=1:Y_Nb_ELEC
                    for j=1:Y_Nb_ELEC
                        for k=1:Y_Nb_ELEC
                            for l=1:Y_Nb_ELEC
                                ROLL_ALONG_Y(count,1:4)=[i,j,k,l];
                                count=count+1;
                            end
                        end
                    end
                end
                ROLL_ALONG_Y(count:end,:)=[];
                
                % Comprehensive data set with all the X lines
                clear COMPREHENSIVE_X_all_lines
                COMPREHENSIVE_X_all_lines=zeros(length(COMPREHENSIVE_X_1_line(:,1))*length(ROLL_ALONG_Y(:,1)),4);
                L_Y=length(ROLL_ALONG_Y(:,1));
                for i=1:length(COMPREHENSIVE_X_1_line(:,1))
                    clear Mat_X_all_lines
                    for j=1:4
                        Mat_X_all_lines(:,j)=COMPREHENSIVE_X_1_line(i,j)*ones(L_Y,1);
                    end
                    COMPREHENSIVE_X_all_lines((i-1)*L_Y+1:i*L_Y,1:4)=Mat_X_all_lines+(ROLL_ALONG_Y-1)*X_Nb_ELEC;
                end
                
                COMPREHENSIVE_X_all_lines=unique(COMPREHENSIVE_X_all_lines,'rows','stable');
                
                clear COMPREHENSIVE_X_1_line actif actif_comprehensive_data_set i j k l L_Y Mat_X_all_lines o ROLL_ALONG_Y
              
                o_2_3_Horizontal_Roll_Along_in_Y_direction_with_all_the_lines=1
                clear o_2_3_Horizontal_Roll_Along_in_Y_direction_with_all_the_lines count
            end
        end
        
        %% 2.4) Alpha and Beta arrays with one line in Y direction
        % Creation : 09/05/2017
        actif=1;
        for o=1
            if actif==1 & Y_Nb_ELEC>4
                % Alpha configuration          C1 -- P1 -- P2 -- C2
                clear ALPHA_ARRAYS
                ALPHA_ARRAYS=zeros(Y_Nb_ELEC*Y_Nb_ELEC*Y_Nb_ELEC*Y_Nb_ELEC,4);
                count=1;
                for C1=1:Y_Nb_ELEC
                    for P1=C1+1:Y_Nb_ELEC
                        for P2=P1+1:Y_Nb_ELEC
                            for C2=P2+1:Y_Nb_ELEC
                                ALPHA_ARRAYS(count,1:4)=[(C1-1)*X_Nb_ELEC+1,(C2-1)*X_Nb_ELEC+1,(P1-1)*X_Nb_ELEC+1,(P2-1)*X_Nb_ELEC+1];
                                count=count+1;
                            end
                        end
                    end
                end
                 ALPHA_ARRAYS(count:end,:)=[];
            
               	% Beta configuration          C2 -- C1 -- P1 -- P2
                clear BETA_ARRAYS
                BETA_ARRAYS=zeros(Y_Nb_ELEC*Y_Nb_ELEC*Y_Nb_ELEC*Y_Nb_ELEC,4);
                count=1;
                for C2=1:Y_Nb_ELEC
                    for C1=C2+1:Y_Nb_ELEC
                        for P1=C1+1:Y_Nb_ELEC
                            for P2=P1+1:Y_Nb_ELEC
                                BETA_ARRAYS(count,1:4)=[(C1-1)*X_Nb_ELEC+1,(C2-1)*X_Nb_ELEC+1,(P1-1)*X_Nb_ELEC+1,(P2-1)*X_Nb_ELEC+1];
                                count=count+1;
                            end
                        end
                    end
                end
                BETA_ARRAYS(count:end,:)=[];
                
                % Comprehensive data set with one Y line
                clear COMPREHENSIVE_Y_1_line
                COMPREHENSIVE_Y_1_line=[ALPHA_ARRAYS;BETA_ARRAYS];
                COMPREHENSIVE_Y_1_line=unique(COMPREHENSIVE_Y_1_line,'rows','stable');
                
                clear C1 C2 P1 P2 ALPHA_ARRAYS BETA_ARRAYS actif o
                
                o_2_4_Alpha_and_Beta_arrays_with_one_line_in_Y_direction=1
                clear o_2_4_Alpha_and_Beta_arrays_with_one_line_in_Y_direction count actif o
            end
        end
        
        %% 2.5) Horizontal Roll Along in X direction with all the lines
        % Creation : 09/05/2017
        actif=1;
        for o=1
            if actif==1 & Y_Nb_ELEC>4
                % Roll Along X Matrix
                clear ROLL_ALONG_X
                ROLL_ALONG_X=zeros(X_Nb_ELEC*X_Nb_ELEC*X_Nb_ELEC*X_Nb_ELEC,4);
                count=1;
                for i=1:X_Nb_ELEC
                    for j=1:X_Nb_ELEC
                        for k=1:X_Nb_ELEC
                            for l=1:X_Nb_ELEC
                                ROLL_ALONG_X(count,1:4)=[i,j,k,l];
                                count=count+1;
                            end
                        end
                    end
                end
                ROLL_ALONG_X(count:end,:)=[];

                % Comprehensive data set with all the Y lines
                clear COMPREHENSIVE_Y_all_lines
                COMPREHENSIVE_Y_all_lines=zeros(length(COMPREHENSIVE_Y_1_line(:,1))*length(ROLL_ALONG_X(:,1)),4);
                L_X=length(ROLL_ALONG_X(:,1));
                for i=1:length(COMPREHENSIVE_Y_1_line(:,1))
                    clear Mat_Y_all_lines
                    for j=1:4
                        Mat_Y_all_lines(:,j)=COMPREHENSIVE_Y_1_line(i,j)*ones(L_X,1);
                    end
                    COMPREHENSIVE_Y_all_lines((i-1)*L_X+1:i*L_X,1:4)=Mat_Y_all_lines+(ROLL_ALONG_X-1);
                end
                
                COMPREHENSIVE_Y_all_lines=unique(COMPREHENSIVE_Y_all_lines,'rows','stable');
                
                clear COMPREHENSIVE_Y_1_line actif actif_comprehensive_data_set i j k l L_X Mat_Y_all_lines o ROLL_ALONG_X
                
                o_2_5_Horizontal_Roll_Along_in_X_direction_with_all_the_lines=1
                clear o_2_5_Horizontal_Roll_Along_in_X_direction_with_all_the_lines actif o
            end
        end
        
        %% 2.6) Cross diagonal dipole dipole with all the lines in X direction
        % Creation : 11/05/2017
        actif=1;
        for o=1
            if actif==1 & X_Nb_ELEC>1 & Y_Nb_ELEC>1
                % Cross diagonal configuration          C1 -- C2 
                %                                       P1 -- P2                            
                clear CROSS_DIAG_ARRAYS_X_2_lines
                CROSS_DIAG_ARRAYS_X_2_lines=zeros(X_Nb_ELEC*X_Nb_ELEC*X_Nb_ELEC*X_Nb_ELEC,4);
                count=1;
                for C1=1:X_Nb_ELEC
                    for C2=C1+1:X_Nb_ELEC
                        for P1=X_Nb_ELEC+1:2*X_Nb_ELEC
                            for P2=P1+1:2*X_Nb_ELEC
                                CROSS_DIAG_ARRAYS_X_2_lines(count,1:4)=[C1,C2,P1,P2];
                                count=count+1;
                            end
                        end
                    end
                end
                CROSS_DIAG_ARRAYS_X_2_lines(count:end,:)=[];
                
                % Propagation to all of the X lines
                % Matrix from two to all the X-lines
                clear ROLL_ALONG_Y
                ROLL_ALONG_Y(1,1:4)=[0,0,0,0];;
                for i=1:Y_Nb_ELEC-1
                    for j=i:Y_Nb_ELEC-1
                        ROLL_ALONG_Y(end+1,1:4)=[i i j j];
                    end
                end
                ROLL_ALONG_Y(1,:)=[];
                
                % Cross diag data set with all the Y lines
                clear CROSS_DIAG_ARRAYS_X_all_lines
                CROSS_DIAG_ARRAYS_X_all_lines=zeros(length(CROSS_DIAG_ARRAYS_X_2_lines(:,1))*length(ROLL_ALONG_Y(:,1)),4);
                L_Y=length(ROLL_ALONG_Y(:,1));
                for i=1:length(CROSS_DIAG_ARRAYS_X_2_lines(:,1))
                    clear Mat_X_all_lines
                    for j=1:4
                        Mat_X_all_lines(:,j)=CROSS_DIAG_ARRAYS_X_2_lines(i,j)*ones(L_Y,1);
                    end
                    CROSS_DIAG_ARRAYS_X_all_lines((i-1)*L_Y+1:i*L_Y,1:4)=Mat_X_all_lines+(ROLL_ALONG_Y-1)*X_Nb_ELEC;
                end
                
                CROSS_DIAG_ARRAYS_X_all_lines=unique(CROSS_DIAG_ARRAYS_X_all_lines,'rows','stable');
                
                clear C1 C2 P1 P2 ROLL_ALONG_Y CROSS_DIAG_ARRAYS_X_2_lines actif i j L_Y Mat_X_all_lines o
                
                o_2_6_Cross_diagonal_dipole_dipole_with_all_the_lines_in_X_direction=1
                clear o_2_6_Cross_diagonal_dipole_dipole_with_all_the_lines_in_X_direction
            end
        end
        
        %% 2.7) Cross diagonal dipole dipole with all the lines in Y direction
        % Creation : 11/05/2017
        actif=1;
        for o=1
            if actif==1 & X_Nb_ELEC>1 & Y_Nb_ELEC>1
                % Cross diagonal configuration          C1 -- P1 
                %                                       C2 -- P2                            
                clear CROSS_DIAG_ARRAYS_Y_2_lines
                CROSS_DIAG_ARRAYS_Y_2_lines(1,1:4)=[0,0,0,0];
                for C1=1:X_Nb_ELEC:length(ELEC(:,1))
                    for C2=C1+X_Nb_ELEC:X_Nb_ELEC:length(ELEC(:,1))
                        for P1=C1+1:X_Nb_ELEC:length(ELEC(:,1))
                            for P2=P1+X_Nb_ELEC:X_Nb_ELEC:length(ELEC(:,1))
                                CROSS_DIAG_ARRAYS_Y_2_lines(end+1,1:4)=[C1,C2,P1,P2];
                            end
                        end
                    end
                end
                CROSS_DIAG_ARRAYS_Y_2_lines(1,:)=[];
                
                % Propagation to all of the Y lines
                % Matrix from two to all the Y-lines
                clear ROLL_ALONG_X
                ROLL_ALONG_X(1,1:4)=[0,0,0,0];;
                for i=1:X_Nb_ELEC-1
                    for j=i:X_Nb_ELEC-1
                        ROLL_ALONG_X(end+1,1:4)=[i i j j];
                    end
                end
                ROLL_ALONG_X(1,:)=[];
                
                % Cross diag data set with all the X lines
                clear CROSS_DIAG_ARRAYS_Y_all_lines
                CROSS_DIAG_ARRAYS_Y_all_lines=zeros(length(CROSS_DIAG_ARRAYS_Y_2_lines(:,1))*length(ROLL_ALONG_X(:,1)),4);
                L_X=length(ROLL_ALONG_X(:,1));
                for i=1:length(CROSS_DIAG_ARRAYS_Y_2_lines(:,1))
                    clear Mat_Y_all_lines
                    for j=1:4
                        Mat_Y_all_lines(:,j)=CROSS_DIAG_ARRAYS_Y_2_lines(i,j)*ones(L_X,1);
                    end
                    CROSS_DIAG_ARRAYS_Y_all_lines((i-1)*L_X+1:i*L_X,1:4)=Mat_Y_all_lines+(ROLL_ALONG_X-1);
                end
                
                CROSS_DIAG_ARRAYS_Y_all_lines=unique(CROSS_DIAG_ARRAYS_Y_all_lines,'rows','stable');
                
                clear C1 C2 P1 P2 ROLL_ALONG_X CROSS_DIAG_ARRAYS_Y_2_lines actif i j L_X Mat_Y_all_lines o
            
                o_2_7_Cross_diagonal_dipole_dipole_with_all_the_lines_in_Y_direction=1
                clear o_2_7Cross_diagonal_dipole_dipole_with_all_the_lines_in_Y_direction
            end
        end
       
        %% 2.8) Total Comprehensive Data Set
        % Creation : 09/05/2017
        actif=1;
        for o=1
            if actif==1
                
                COMPREHENSIVE(1,1:4)=[0,0,0,0];
                if X_Nb_ELEC>4
                    COMPREHENSIVE=[COMPREHENSIVE;COMPREHENSIVE_X_all_lines];
                end
                if Y_Nb_ELEC>4
                    COMPREHENSIVE=[COMPREHENSIVE;COMPREHENSIVE_Y_all_lines];
                end
                if X_Nb_ELEC>1 & Y_Nb_ELEC>1
                    COMPREHENSIVE=[COMPREHENSIVE;CROSS_DIAG_ARRAYS_X_all_lines;CROSS_DIAG_ARRAYS_Y_all_lines];
                end
                COMPREHENSIVE(1,:)=[];
                
                
                COMPREHENSIVE=unique(COMPREHENSIVE,'rows','stable');
                
                clear COMPREHENSIVE_X_all_lines COMPREHENSIVE_Y_all_lines actif CROSS_DIAG_ARRAYS_X_all_lines CROSS_DIAG_ARRAYS_Y_all_lines
                
                o_2_8_Total_Comprehensive_Data_Set=1
                clear o_2_8_Total_Comprehensive_Data_Set
            end
        end
        
        %% 2.9) Geometric factor filtering
        % Creation : 09/05/2017
        K_max=1000;
        actif=1;
        for o=1
            if actif==1
                % Number of electrodes
                N_Elec=length(ELEC(:,1));

                % Electrode visualization
                for o=1:1
                    figure('Color', [ 1 1 1])
                    plot3(ELEC(:,1),ELEC(:,2),ELEC(:,3),'.k','MarkerSize',20)
                    hold on
                    xlim([0 X_Nb_ELEC*X_GRID_SPACING+2*X_Elec_1])
                    ylim([0 Y_Nb_ELEC*Y_GRID_SPACING+2*Y_Elec_1])
                    daspect([1,1,1])
                end
                
                % Distances matrix calculation
                for i=1:N_Elec
                    for j=1:N_Elec
                        DISTANCE(i,j)=pdist([ELEC(i,1),ELEC(i,2),ELEC(i,3);ELEC(j,1),ELEC(j,2),ELEC(j,3)],'euclidean');
                    end
                end
                
                % Distance matrix
                clear K_Mat_1 K_Mat_2 K_Mat_3 K_Mat_4
                K_Mat_1=[1,3;1,4;2,3;2,4;1,3;1,4;2,3;2,4];
                K_Mat_2=[0;0;0;0;0;0;0;0];
                K_Mat_3=[1;-1;-1;1;1;-1;-1;1];
                
                % Geometric factor calculation
                clear H; clear E1; clear E2;
                clear COORD_E1; clear COORD_E2
                H=0;
                for i=1:8
                    E1=COMPREHENSIVE(:,K_Mat_1(i,1))+N_Elec*K_Mat_2(i,1);
                    E2=COMPREHENSIVE(:,K_Mat_1(i,2))+N_Elec*K_Mat_2(i,1);
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
                
                % K filtering : suppr of configurations with too high K values
                clear elimination_1
                elimination_1=find(abs(K(:,1))>K_max);
                COMPREHENSIVE(:,5)=K;
                COMPREHENSIVE_filtK=COMPREHENSIVE;
                COMPREHENSIVE_filtK(elimination_1,:)=[];
                
                H_filtK=H;
                H_filtK(elimination_1,:)=[];
                K_Mat_4_filtK=K_Mat_4;
                K_Mat_4_filtK(elimination_1,:)=[];
                COORD_E1(elimination_1,:)=[];
                COORD_E2(elimination_1,:)=[];
                
                clear actif E1 E2 elimination_1 i j o 
                
                o_2_9_Geometric_factor_filtering=1
                clear o_2_9_Geometric_factor_filtering
            end
        end
        
        %% 2.10) Abberant configurations filtering
        % Creation : 10/05/2017
        actif=1;
        for o=1
            if actif==1
                % Suppression of NaN values (if existing)
                clear elimination_3
                elimination_3=find(COMPREHENSIVE_filtK(:,5)==0);
                COMPREHENSIVE_filtK(elimination_3,:)=[];
                H_filtK(elimination_3,:)=[];
                K_Mat_4_filtK(elimination_3,:)=[];
                COORD_E1(elimination_3,:)=[];
                COORD_E2(elimination_3,:)=[];
                
                clear o actif
                
                o_2_10_Abberant_configurations_filtering=1
                clear o_2_10_Geometric_factor_filtering
            end
        end
        
        %% 2.11) Geometric Factor Relative Error filtering
        % Creation : 10/05/2017
        K_err_max=5;
        actif=1;
        for o=1
            if actif==1
                % Pre processing of electrodes coordinates
                clear COORD_ELEC
                for i=1:4
                    Ei=COMPREHENSIVE_filtK(:,i);
                    COORD_ELEC(:,(i-1)*3+1:i*3)=ELEC(COMPREHENSIVE_filtK(:,i),1:3);
                end
                
                % K sensitivity calculation (S_2) relative to position errors
                clear K_Mat_5;
                K_Mat_5=[1,2;3,4;1,3;2,4];
                clear dK_dEi_2; % square of the error prosition (col 1 = erreur % C1 ... )
                clear S_2 % square of the sensitivities
                S_2=0;
                for i=1:4 % error position of C1, C2, P1 and P2
                    somme_1=0;
                    for j=1:3 % for each dimension (x, y and z)
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
                
                % Re calculation : Geometric Factor Relative Error
                clear Re
                Re = sqrt(S_2)./COMPREHENSIVE_filtK(:,5);
                COMPREHENSIVE_filtK(:,6)=Re;
                
                % Re filtering : suppr of configurations with too high re values
                COMPREHENSIVE_filtK_filtRe=COMPREHENSIVE_filtK;
                clear elimination_2
                elimination_2=find(abs(Re(:,1))>K_err_max);
                COMPREHENSIVE_filtK_filtRe(elimination_2,:)=[];
                
                clear col_COORD_ELEC COORD_E1 COORD_E2 COORD_ELEC DISTANCE dK_dEi_2 Ei elimination_2 elimination_3
                clear H H_filtK i ind_Ci_Pi ind_E j k K K_Mat_1 K_Mat_2 K_Mat_3 K_Mat_4 K_Mat_4_filtK
                clear K_Mat_5 o Re S_2 somme somme_1 somme_2 actif
                
                o_2_11_Geometric_Factor_Relative_Error_filtering=1
                clear o_2_11_Geometric_Factor_Relative_Error_filtering
            end
        end
        
      	%% 2.12) Abberant configurations filtering
        % Creation : 10/05/2017
        actif=1;
        for o=1
            if actif==1
                % Suppression of NaN values (if existing)
                COMPREHENSIVE_filtK_filtRe(isnan(COMPREHENSIVE_filtK_filtRe(:,5)),:)=[];
                
                clear o actif COMPREHENSIVE COMPREHENSIVE_filtK K_err_max K_max 
                
                o_2_12_Abberant_configurations_filtering=1
                clear o_2_12_Abberant_configurations_filtering
            end
        end
        
    end
end

%% III) COMPARE R METHOD
%--------------------------------------------------------------------------
actif_compare_r=1;
for o=1
    if actif_compare_r==1
        %% 3.1) Sensitivity calculation grid definition
        % Creation : 10/05/2017
        actif=1;
        for o=1
            if actif==1
                [X,Y,Z] = meshgrid(0:X_SENSI_GRID_SPACING:X_Nb_ELEC*X_GRID_SPACING+X_Elec_1,0:Y_SENSI_GRID_SPACING:Y_Nb_ELEC*Y_GRID_SPACING+Y_Elec_1,0:-Z_SENSI_GRID_SPACING:-Z_MAX);
                
                % Point grid matrix
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
                COORD_2=COORD;
                              
                % Sensitivity grid calculation visualization
                figure('Color', [ 1 1 1])
                plot3(ELEC(:,1),ELEC(:,2),ELEC(:,3),'.k','MarkerSize',20)
                hold on
                plot3(COORD_2(:,1),COORD_2(:,2),COORD_2(:,3),'.r','MarkerSize',5)
                xlim([0 X_Nb_ELEC*X_GRID_SPACING+X_Elec_1])
                ylim([0 Y_Nb_ELEC*Y_GRID_SPACING+Y_Elec_1])
                zlim([-Z_MAX 0])
                daspect([1,1,1])
                 
                clear actif actif_compare_r COORD i j k o 
            end
        end
        
        %% 3.2) Compare R preparation
        % Creation : 10/05/2017
        actif=1;
        for o=1
            if actif==1
                % Sensitivity calculation for each electrode (homogeneous 3D model)
                MAILLAGE=COORD_2;
                
                % Matrix of all possible DUOS
                DUOS=zeros(N_Elec*(N_Elec-1),2);
                count=0;
                tic
                for elec1=1:N_Elec
                    for elec2=1:N_Elec
                        if elec1~=elec2
                            count=count+1;
                            DUOS_ADVANCE=round(count/(192*191)*100,0)
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
                        MAT1_ADVANCE=round(j/length(DUOS(:,1))*100,0)
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
                
                % Sensitivity Matrix
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
                        SENSIBILITE_ADVANCE=round(j/length(DUOS(:,1))*100,0)
                    end
                end
                toc
                
                clear actif colonne count DUOS_ADVANCE elec1 elec2 j MAT1_ADVANCE somme o SENSIBILITE_ADVANCE
            end
        end
        
        %% 3.3) Sensitivity calculation for each configuration
        % Creation : 10/05/2017
        actif=1;
        for o=1
            if actif==1
                % Sensitivity calculation for each quadripole
                tic
                N_COLUMNS=5;
                
                CONFIGURATIONS=COMPREHENSIVE_filtK_filtRe;
                SYNTHESE_VALUES=zeros(length(MAILLAGE(:,1)),N_COLUMNS);
                SYNTHESE_INDICE=zeros(length(MAILLAGE(:,1)),N_COLUMNS);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                load SYNTHESE_VALUES_84x1.mat
                load SYNTHESE_INDICE_84x1.mat
                clear SYNTHESE_INDICE
                clear SYNTHESE_VALUES
                SYNTHESE_INDICE=SYNTHESE_INDICE_84x1;
                SYNTHESE_VALUES=SYNTHESE_VALUES_84x1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=length(CONFIGURATIONS(:,1)):length(CONFIGURATIONS(:,1))
                    clear SENSIBILITE_TOT
                    C1_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS(i,1)&DUOS(:,2)==CONFIGURATIONS(i,3)));
                    C2_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS(i,2)&DUOS(:,2)==CONFIGURATIONS(i,3)));
                    C1_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS(i,1)&DUOS(:,2)==CONFIGURATIONS(i,4)));
                    C2_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS(i,2)&DUOS(:,2)==CONFIGURATIONS(i,4)));
                    SENSIBILITE_TOT(:,1)=abs(C1_P1-C2_P1-C1_P2+C2_P2);
                    SENS_TOT_ADVANCE=round(i/length(CONFIGURATIONS(:,1))*100,3)
                    
                    for j=1:length(MAILLAGE(:,1))
                        trouve=0;
                        for k=1:N_COLUMNS
                            if trouve==0
                                if SENSIBILITE_TOT(j,1)>SYNTHESE_VALUES(j,k)
                                    SYNTHESE_VALUES(j,k+1:end)=SYNTHESE_VALUES(j,k:end-1);
                                    SYNTHESE_VALUES(j,k)=SENSIBILITE_TOT(j,1);
                                    SYNTHESE_INDICE(j,k+1:end)=SYNTHESE_INDICE(j,k:end-1);
                                    SYNTHESE_INDICE(j,k)=i;
                                    trouve=1;
                                end
                            end
                        end
                    end
                    
                end
                toc
                TIME_TAKEN(1)=length(COMPREHENSIVE_filtK_filtRe(:,1));
                TIME_TAKEN(2)=toc
                
                clear C1_P1 C1_P2 C2_P1 C2_P2 i j k Mat1 Mat2 o trouve
            end
        end
        
        %% 3.4) Ranking of the configurations
        % Creation : 10/05/2017
        actif=1;
        for o=1
            if actif==1
                % Ponderation of each cell of the SYNTHESE_VALUES
                clear PONDERATION;
                for i=1:N_COLUMNS
                    %PONDERATION(:,i)=(SYNTHESE_VALUES(:,i)-SYNTHESE_VALUES(:,N_COLUMNS))./(SYNTHESE_VALUES(:,1)-SYNTHESE_VALUES(:,N_COLUMNS))*100*(N_COLUMNS+1-i)*(N_COLUMNS+1-i);
                    %PONDERATION(:,i)=(SYNTHESE_VALUES(:,i)-SYNTHESE_VALUES(:,N_COLUMNS))./(SYNTHESE_VALUES(:,1)-SYNTHESE_VALUES(:,N_COLUMNS))*100*(N_COLUMNS+1-i);
                    %PONDERATION(:,i)=(SYNTHESE_VALUES(:,i)-SYNTHESE_VALUES(:,N_COLUMNS))./(SYNTHESE_VALUES(:,1)-SYNTHESE_VALUES(:,N_COLUMNS))*100;
                    PONDERATION(:,i)=(SYNTHESE_VALUES(:,i)-SYNTHESE_VALUES(:,N_COLUMNS))./(SYNTHESE_VALUES(:,1)-SYNTHESE_VALUES(:,N_COLUMNS))*100+((N_COLUMNS+1-i)*N_COLUMNS);
                    %PONDERATION(:,i)=zeros(length(MAILLAGE(:,1)),1)+((N_COLUMNS+1-i)*N_COLUMNS);
                end
                
                % Matrix (VALUES, INDICE and PONDERATION) re-writing for calculations
                PONDERATION_NEW=PONDERATION(:,1);
                SYNTHESE_VALUES_NEW=SYNTHESE_VALUES(:,1);
                SYNTHESE_INDICE_NEW=SYNTHESE_INDICE(:,1);
                for i=2:N_COLUMNS
                    PONDERATION_NEW=[PONDERATION_NEW;PONDERATION(:,i)];
                    SYNTHESE_VALUES_NEW=[SYNTHESE_VALUES_NEW;SYNTHESE_VALUES(:,i)];
                    SYNTHESE_INDICE_NEW=[SYNTHESE_INDICE_NEW;SYNTHESE_INDICE(:,i)];
                end
                
                % Ranking of each configuration
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
                CONFIGURATION_UNIQUE=sortrows(CONFIGURATION_UNIQUE,-2); % decreasing order
                CONFIGURATION_UNIQUE(isnan(CONFIGURATION_UNIQUE(:,2)),:)=[];
                % Information distribution visualization
                figure('Color', [ 1 1 1])
                semilogy(CONFIGURATION_UNIQUE(:,end),'.k','MarkerSize',1)
                grid
                figure('Color', [ 1 1 1])
                plot(CONFIGURATION_UNIQUE(:,end),'.k','MarkerSize',1)
                grid
                
                clear actif i j o trouve 
            end
        end
        
    end
end

%% IV) Identification of the optimized protocol number of configurations 
%--------------------------------------------------------------------------
actif_optimized_protocol_identification=1;
for o=1
    if actif_optimized_protocol_identification==1
        close all
        %% 4.1) Calculation of efficiency VS number of configuration graph
        % Creation : 10/05/2017
        actif=1;
        for o=1
            if actif==1
                count=0;
                clear SYNTHESE_TEST
                for count1=1:0.25:4
                  	Nb_CONFIG=10^count1;
                    if Nb_CONFIG<length(CONFIGURATION_UNIQUE(:,1))
                        
                        % Maximum sensitivity values : total protocol
                        clear PROTOCOLE_TOTAL
                        
                        PROTOCOLE_TOTAL=CONFIGURATIONS(CONFIGURATION_UNIQUE(:,1),:);
                        PROTOCOLE_TOTAL(:,7)=CONFIGURATION_UNIQUE(:,2); % Ranking column
                        % Sensitivity matrix for the total protocol
                        SENSIBILITE_TOTAL=[MAILLAGE,SYNTHESE_VALUES(:,1)];
                        
                        % Protocol with a limited number of configurations
                        PROTOCOLE_LIMITE=CONFIGURATIONS(CONFIGURATION_UNIQUE(1:Nb_CONFIG,1),:);
                        PROTOCOLE_LIMITE(1:Nb_CONFIG,7)=CONFIGURATION_UNIQUE(1:Nb_CONFIG,2); % Colonne de ranking
                        
                        % New calculation of the sensitivities for the
                        % limited protocol
                        tic
                        
                        CONFIGURATIONS_SELECT=PROTOCOLE_LIMITE;
                        clear SENSIBILITE_EACH
                        SENSIBILITE_EACH=zeros(length(MAILLAGE(:,1)),length(CONFIGURATIONS_SELECT(:,1)));
                        
                        SYNTHESE_VALUES_LIMITE=zeros(length(MAILLAGE(:,1)),N_COLUMNS);
                        SYNTHESE_INDICE_LIMITE=zeros(length(MAILLAGE(:,1)),N_COLUMNS);
                        for i=1:length(CONFIGURATIONS_SELECT(:,1))
                            clear SENSIBILITE_TOT
                            C1_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT(i,1)&DUOS(:,2)==CONFIGURATIONS_SELECT(i,3)));
                            C2_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT(i,2)&DUOS(:,2)==CONFIGURATIONS_SELECT(i,3)));
                            C1_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT(i,1)&DUOS(:,2)==CONFIGURATIONS_SELECT(i,4)));
                            C2_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT(i,2)&DUOS(:,2)==CONFIGURATIONS_SELECT(i,4)));
                            SENSIBILITE_TOT(:,1)=abs(C1_P1-C2_P1-C1_P2+C2_P2);
                            SENS_TOT_AVANCEMENT=round(i/length(CONFIGURATIONS_SELECT(:,1))*100,0)
                            
                            SENSIBILITE_EACH(:,i)=SENSIBILITE_TOT(:,1);
                            
                            for j=1:length(MAILLAGE(:,1))
                                trouve=0;
                                for k=1:N_COLUMNS
                                    if trouve==0
                                        if SENSIBILITE_TOT(j,1)>SYNTHESE_VALUES_LIMITE(j,k)
                                            SYNTHESE_VALUES_LIMITE(j,k+1:end)=SYNTHESE_VALUES_LIMITE(j,k:end-1);
                                            SYNTHESE_VALUES_LIMITE(j,k)=SENSIBILITE_TOT(j,1);
                                            SYNTHESE_INDICE_LIMITE(j,k+1:end)=SYNTHESE_INDICE_LIMITE(j,k:end-1);
                                            SYNTHESE_INDICE_LIMITE(j,k)=i;
                                            trouve=1;
                                        end
                                    end
                                end
                            end
                        end
                        toc
                        
                        % Sensibility matrix for the limited protocol
                        SENSIBILITE_LIMITE=MAILLAGE;
                        SENSIBILITE_LIMITE(:,4)=SYNTHESE_VALUES_LIMITE(:,1);
                        
                        % Efficiency calculation for the limited protocol
                        EFFICACITE=mean(100-(SENSIBILITE_TOTAL(:,4)-SENSIBILITE_LIMITE(:,4))./SENSIBILITE_TOTAL(:,4)*100);
                        
                    end
                    count=count+1;
                    SYNTHESE_TEST(count,1:2)=[Nb_CONFIG,EFFICACITE];
                end
                
                close all
                
                % Graphic representation of the efficiency of the protocol
                % versus the number of configurations
                figure('Color', [ 1 1 1])
                semilogx(SYNTHESE_TEST(:,1),SYNTHESE_TEST(:,2),'--ok','MarkerSize',5)
                xlabel('Number of configurations of the protocol','FontSize',12)
                ylabel('Efficiency of the protocole in %','FontSize',12)
                grid
                hold on
                ylim([10 100])
                title(['Efficiency graph for the ',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),' optimized protocol'],'FontSize',14)
                
                clear actif C1_P1 C1_P2 C2_P1 C2_P2 count count1 EFFICACITE i j k o trouve
            end
        end
        
        %% 4.2) User selection of the number of configurations
        % Creation : 11/05/2017
        actif=1;
        for o=1
            if actif==1
                % Dialog box to choose the number of configurations
                p = get(gcf,'Position');
                for o=1
                    d = dialog('Position',[p(1)+800,p(2)+300,400,400],'Name','3D Vertical Scan');
                    
                    txt = uicontrol('Parent',d,...
                        'Style','text',...
                        'Position',[100 200 210 80],...
                        'String',['Select the number of configurations needed by clicking on the dashed line. Less than ',num2str(length(CONFIGURATION_UNIQUE(:,1))),' configurations !'],...
                        'Fontsize',10);
                    
                    btn = uicontrol('Parent',d,...
                        'Position',[170 20 70 25],...
                        'String','Close',...
                        'Callback','delete(gcf)');
                end % Dialog box 
                
                % Selection of the number of configurations wanted
                N_Config=ginput(1);
                N_Config(1)=round(N_Config(1)/10,0)*10;
                N_Config(2)=round(N_Config(2),2);
                
                 % Dialog box to inform the user about his choice
                 p = get(gcf,'Position');
                 for o=1
%                     d = dialog('Position',[p(1)+800,p(2)+300,400,400],'Name','3D Vertical Scan');
%                     
%                     txt = uicontrol('Parent',d,...
%                         'Style','text',...
%                         'Position',[100 200 210 90],...
%                         'String',['You have selected an optimized protocol with : ',num2str(N_Config(1)),' configurations which corresponds to an efficiency of ',num2str(N_Config(2)),' %.'],...
%                         'Fontsize',10);
%                     
%                     btn = uicontrol('Parent',d,...
%                         'Position',[170 20 70 25],...
%                         'String','Close',...
%                         'Callback','delete(gcf)');
                 end
                
                % Calculation of the sensitivities for the user selected optimized protocol
                Nb_CONFIG=N_Config(1);
                
                % Protocol with a limited number of configurations
                PROTOCOLE_LIMITE=CONFIGURATIONS(CONFIGURATION_UNIQUE(1:Nb_CONFIG,1),:);
                PROTOCOLE_LIMITE(1:Nb_CONFIG,7)=CONFIGURATION_UNIQUE(1:Nb_CONFIG,2); % Colonne de ranking
                
                % New calculation of the sensitivities for the limited protocol
                tic
                
                CONFIGURATIONS_SELECT_OPT=PROTOCOLE_LIMITE;
                clear SENSIBILITE_EACH
                SENSIBILITE_EACH=zeros(length(MAILLAGE(:,1)),length(CONFIGURATIONS_SELECT_OPT(:,1)));
                
                SYNTHESE_VALUES_LIMITE=zeros(length(MAILLAGE(:,1)),N_COLUMNS);
                SYNTHESE_INDICE_LIMITE=zeros(length(MAILLAGE(:,1)),N_COLUMNS);
                erreur_calcul_1=0;
                for i=1:length(CONFIGURATIONS_SELECT_OPT(:,1))
                    clear SENSIBILITE_TOT C1_P1 C1_P2 C2_P1 C2_P2
                    C1_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT_OPT(i,1)&DUOS(:,2)==CONFIGURATIONS_SELECT_OPT(i,3)));
                    C2_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT_OPT(i,2)&DUOS(:,2)==CONFIGURATIONS_SELECT_OPT(i,3)));
                    C1_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT_OPT(i,1)&DUOS(:,2)==CONFIGURATIONS_SELECT_OPT(i,4)));
                    C2_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT_OPT(i,2)&DUOS(:,2)==CONFIGURATIONS_SELECT_OPT(i,4)));
                    if isempty(C1_P1) | isempty(C2_P1) | isempty(C1_P2) | isempty(C2_P2)
                        SENSIBILITE_TOT(:,1)=zeros(length(MAILLAGE(:,1)),1);
                        erreur_calcul_1=erreur_calcul_1+1;
                    else
                        SENSIBILITE_TOT(:,1)=abs(C1_P1-C2_P1-C1_P2+C2_P2);
                    end
                    SENS_TOT_AVANCEMENT=round(i/length(CONFIGURATIONS_SELECT_OPT(:,1))*100,0)
                    
                    SENSIBILITE_EACH(:,i)=SENSIBILITE_TOT(:,1);
                    
                    for j=1:length(MAILLAGE(:,1))
                        trouve=0;
                        for k=1:N_COLUMNS
                            if trouve==0
                                if SENSIBILITE_TOT(j,1)>SYNTHESE_VALUES_LIMITE(j,k)
                                    SYNTHESE_VALUES_LIMITE(j,k+1:end)=SYNTHESE_VALUES_LIMITE(j,k:end-1);
                                    SYNTHESE_VALUES_LIMITE(j,k)=SENSIBILITE_TOT(j,1);
                                    SYNTHESE_INDICE_LIMITE(j,k+1:end)=SYNTHESE_INDICE_LIMITE(j,k:end-1);
                                    SYNTHESE_INDICE_LIMITE(j,k)=i;
                                    trouve=1;
                                end
                            end
                        end
                    end
                end
                erreur_calcul_1
                toc
                
                % Sensibility matrix for the limited protocol
                SENSIBILITE_LIMITE=MAILLAGE;
                SENSIBILITE_LIMITE(:,4)=SYNTHESE_VALUES_LIMITE(:,1);
                
                % Efficiency calculation for the limited protocol
                EFFICACITE_OPT=mean(100-(SENSIBILITE_TOTAL(:,4)-SENSIBILITE_LIMITE(:,4))./SENSIBILITE_TOTAL(:,4)*100);
                
                % Update of the figure
                plot(N_Config(1),N_Config(2),'.b','MarkerSize',30)
                text(N_Config(1)+5,N_Config(2)-2.5,['Efficiency = ',num2str(round(EFFICACITE_OPT,1)),' %'],'Color','blue','FontSize',14)
                text(N_Config(1)+5,N_Config(2)-5.5,[num2str(N_Config(1)),' configurations'],'Color','blue','FontSize',14)
                
               
                % Simple representation of the efficiency of the selected optimized protocol
                if standard_exist==0
                    figure('Color', [ 1 1 1])
                    % Optimized protocol sensitivities in blue
                    loglog(SENSIBILITE_TOTAL(:,4),SENSIBILITE_LIMITE(:,4),'.b','MarkerSize',65)
                    hold on
                    % labels
                    xlabel(['Sensitivity values for the comprehensive data set (',num2str(length(CONFIGURATIONS(:,1))),' configurations)'],'FontSize',12)
                    ylabel(['Sensitivity values for the optimized protocol in blue'],'FontSize',12)
                    xlim([1*10^(-5) 0.5])
                    ylim([1*10^(-5) 0.5])
                    title(['Efficiency of the ',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),' optimized protocol compared to the standard protocol'],'FontSize',14)
                    text(4*10^(-5),0.03,['Efficiency of the optimized protocol = ',num2str(round(EFFICACITE_OPT,1)),' % (with ',num2str(length(CONFIGURATIONS_SELECT_OPT(:,1))),' configurations)'],'Color','blue','FontSize',14)
                    grid
                    hold on
                    plot([10^(-10) 10],[10^(-10) 10],'-r')
                end
                
                % Calculation of the sensitivity and efficiency for the standard protocol (if existing)
                for o=1
                    if standard_exist==1
                        tic
                        
                        CONFIGURATIONS_SELECT_STAND=PRO_STANDARD;
                        clear SENSIBILITE_EACH
                        SENSIBILITE_EACH=zeros(length(MAILLAGE(:,1)),length(CONFIGURATIONS_SELECT_STAND(:,1)));
                        
                        SYNTHESE_VALUES_STANDARD=zeros(length(MAILLAGE(:,1)),N_COLUMNS);
                        SYNTHESE_INDICE_STANDARD=zeros(length(MAILLAGE(:,1)),N_COLUMNS);
                        erreur_calcul_2=0;
                        for i=1:length(CONFIGURATIONS_SELECT_STAND(:,1))
                            clear SENSIBILITE_TOT
                            C1_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT_STAND(i,1)&DUOS(:,2)==CONFIGURATIONS_SELECT_STAND(i,3)));
                            C2_P1=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT_STAND(i,2)&DUOS(:,2)==CONFIGURATIONS_SELECT_STAND(i,3)));
                            C1_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT_STAND(i,1)&DUOS(:,2)==CONFIGURATIONS_SELECT_STAND(i,4)));
                            C2_P2=SENSIBILITE(:,find(DUOS(:,1)==CONFIGURATIONS_SELECT_STAND(i,2)&DUOS(:,2)==CONFIGURATIONS_SELECT_STAND(i,4)));
                            if isempty(C1_P1) | isempty(C2_P1) | isempty(C1_P2) | isempty(C2_P2)
                                SENSIBILITE_TOT(:,1)=zeros(length(MAILLAGE(:,1)),1);
                                erreur_calcul_2=erreur_calcul_2+1;
                            else
                                SENSIBILITE_TOT(:,1)=abs(C1_P1-C2_P1-C1_P2+C2_P2);
                            end
                            SENSIBILITE_TOT(:,1)=abs(C1_P1-C2_P1-C1_P2+C2_P2);
                            SENS_TOT_AVANCEMENT=round(i/length(CONFIGURATIONS_SELECT_STAND(:,1))*100,0)
                            
                            SENSIBILITE_EACH(:,i)=SENSIBILITE_TOT(:,1);
                            
                            for j=1:length(MAILLAGE(:,1))
                                trouve=0;
                                for k=1:N_COLUMNS
                                    if trouve==0
                                        if SENSIBILITE_TOT(j,1)>SYNTHESE_VALUES_STANDARD(j,k)
                                            SYNTHESE_VALUES_STANDARD(j,k+1:end)=SYNTHESE_VALUES_STANDARD(j,k:end-1);
                                            SYNTHESE_VALUES_STANDARD(j,k)=SENSIBILITE_TOT(j,1);
                                            SYNTHESE_INDICE_STANDARD(j,k+1:end)=SYNTHESE_INDICE_STANDARD(j,k:end-1);
                                            SYNTHESE_INDICE_STANDARD(j,k)=i;
                                            trouve=1;
                                        end
                                    end
                                end
                            end
                        end
                        erreur_calcul_2
                        toc
                        
                        % Sensibility matrix for the limited protocol
                        SENSIBILITE_STANDARD=MAILLAGE;
                        SENSIBILITE_STANDARD(:,4)=SYNTHESE_VALUES_STANDARD(:,1);
                        
                        % Efficiency calculation for the limited protocol
                        EFFICACITE_STAND=mean(100-(SENSIBILITE_TOTAL(:,4)-SENSIBILITE_STANDARD(:,4))./SENSIBILITE_TOTAL(:,4)*100);
                        
                        % Simple representation of the efficiency of the selected standard protocol
                        figure('Color', [ 1 1 1])
                        % Standard protocol sensitivities in green
                        loglog(SENSIBILITE_TOTAL(:,4),SENSIBILITE_STANDARD(:,4),'.g','MarkerSize',6)
                        hold on
                        % Optimized protocol sensitivities in blue
                        loglog(SENSIBILITE_TOTAL(:,4),SENSIBILITE_LIMITE(:,4),'.b','MarkerSize',6)
                        % labels
                        xlabel(['Sensitivity values for the comprehensive data set (',num2str(length(CONFIGURATIONS(:,1))),' configurations)'],'FontSize',12)
                        ylabel(['Sensitivity values for the protocols (standard in green, optimized in blue)'],'FontSize',12)
                        title(['Efficiency of the ',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),' optimized protocol compared to the standard protocol'],'FontSize',14)
                        xlim([1*10^(-5) 0.5])
                        ylim([1*10^(-5) 0.5])
                        text(4*10^(-5),0.03-0.01,['Efficiency of the standard protocol = ',num2str(round(EFFICACITE_STAND,1)),' % (with ',num2str(length(CONFIGURATIONS_SELECT_STAND(:,1))),' configurations)'],'Color','green','FontSize',14)
                        text(4*10^(-5),0.03,['Efficiency of the optimized protocol = ',num2str(round(EFFICACITE_OPT,1)),' % (with ',num2str(length(CONFIGURATIONS_SELECT_OPT(:,1))),' configurations)'],'Color','blue','FontSize',14)
                        grid
                        hold on
                        plot([10^(-10) 10],[10^(-10) 10],'-r')
                    end
                end
                clear actif btn C1_P1 C1_P2 C2_P1 C2_P2 d i j k o p trouve txt 
            end
        end
        
        %% 4.2) 3D Visualization of the sensivity distribution for the selected protocols
        % Creation : 11/05/2017
        actif=1;
        for o=1
            if actif==1
                if standard_exist==1
                    figure('color',[1 1 1])
                    hold on
                    
                    SENSIBILITE_DATA=SENSIBILITE_STANDARD;
                    SENSIBILITE_DATA(:,4)=log(SENSIBILITE_DATA(:,4));
                    % Interpolation
                    clearvars F;
                    F = scatteredInterpolant(SENSIBILITE_DATA(:,1),SENSIBILITE_DATA(:,2),SENSIBILITE_DATA(:,3),SENSIBILITE_DATA(:,4),'natural');
                    vq = F(X,Y,Z);
                    ValAffichee = vq;
                    
                    % Figure editing
                    daspect([1,1,1])
                    axis tight
                    ax = gca;
                    ax.FontSize = 13;
                    %view(0,0)
                    view(-37.5,30)
                    light('Position',[-1 -1 0.75],'Style','infinite')
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
                    
                    % Volume visualization (ONE cutting value)
                    Valeur=[-7;-5];
                    V_Max=-5;
                    V_Min=-10;
                    Valeur_Indice = round((Valeur-V_Min)*10/(V_Max-V_Min));
                    for i=1:length(Valeur_Indice(:,1))
                        p1=patch(isocaps(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
                            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.7);
                        p1 = patch(isosurface(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
                            'FaceColor',cmap(Valeur_Indice(i,1),:),'EdgeColor','none','FaceAlpha',0.7);
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
                    
                    
                    % Labels
                    xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
                    ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
                    zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
                    
                    title(['3D sensitivity distribution for the ',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),' standard protocol'],'FontSize',14)
                        
                    % Electrode visualization
                    for i=1:length(ELEC(:,1))
                        plot3(ELEC(i,1),ELEC(i,2),ELEC(i,3),'k.','MarkerSize',20)
                    end
                    
                    % Limits of the model
                    xlim([0 X_Nb_ELEC*X_GRID_SPACING+X_Elec_1])
                    ylim([0 Y_Nb_ELEC*Y_GRID_SPACING+Y_Elec_1])
                    zlim([-Z_MAX 0])
                    
                    % Legend
                    clear colorbar
                    colorbar
                    caxis([-10 -5]);
                    TicksM = [-10;-9;-8;-7;-6;-5];
                    TicksL = round(round(exp(TicksM),4)*10000,0);
                    
                    c=colorbar;
                    c.Label.String = 'Sensitivity cell value (x10000)';
                    c.Label.FontWeight='bold';
                    c.Label.FontSize=14;
                    c.Ticks=TicksM;
                    c.TickLabels={num2str(TicksL)};
                end
            end
        end % vue 2 with STANDARD PROTOCOL
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
                
                % Figure editing
                daspect([1,1,1])
                axis tight
                ax = gca;
                ax.FontSize = 13;
                %view(0,0)
                view(-37.5,30)
                light('Position',[-1 -1 0.75],'Style','infinite')
                %camzoom(1.8)
                camproj perspective
                colormap (flipud(hot(10)));
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
                
                % Volume visualization (ONE cutting value)
                Valeur=[-7;-5];
                V_Max=-5;
                V_Min=-10;
                Valeur_Indice = round((Valeur-V_Min)*10/(V_Max-V_Min));
                for i=1:length(Valeur_Indice(:,1))
                    p1=patch(isocaps(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
                        'FaceColor','interp','EdgeColor','none','FaceAlpha',0.7);
                    p1 = patch(isosurface(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
                        'FaceColor',cmap(Valeur_Indice(i,1),:),'EdgeColor','none','FaceAlpha',0.7);
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
                
                
                % Labels
                xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
                ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
                zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
                
                title(['3D sensitivity distribution for the ',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),' optimized protocol (',num2str(round(EFFICACITE_OPT,1)),'%)'],'FontSize',14)
                
                
                % Electrode visualization
                for i=1:length(ELEC(:,1))
                    plot3(ELEC(i,1),ELEC(i,2),ELEC(i,3),'k.','MarkerSize',20)
                end
                
                % Limits of the model
                xlim([0 X_Nb_ELEC*X_GRID_SPACING+X_Elec_1])
                ylim([0 Y_Nb_ELEC*Y_GRID_SPACING+Y_Elec_1])
                zlim([-Z_MAX 0])
                
                % Legend
                clear colorbar
                colorbar
                caxis([-10 -5]);
                TicksM = [-10;-9;-8;-7;-6;-5];
                TicksL = round(round(exp(TicksM),4)*10000,0);
                
                c=colorbar;
                c.Label.String = 'Sensitivity cell value (x10000)';
                c.Label.FontWeight='bold';
                c.Label.FontSize=14;
                c.Ticks=TicksM;
                c.TickLabels={num2str(TicksL)};
            end
        end % vue 2 with OPTIMIZED PROTOCOL
        
        %% 4.3) 3D Visualization of the difference between standard and optimized sensitivities
        % Creation : 16/05/2017
        actif=1;
        for o=1
            if actif==1
                if standard_exist==1
                    figure('color',[1 1 1])
                    hold on
                    
                    SENSIBILITE_DATA=SENSIBILITE_LIMITE;
                    SENSIBILITE_DATA(:,4)=(SENSIBILITE_LIMITE(:,4)-SENSIBILITE_STANDARD(:,4))./SENSIBILITE_STANDARD(:,4)*100;
                    %SENSIBILITE_DATA(:,4)=log(SENSIBILITE_DATA(:,4));
                    % Interpolation
                    clearvars F;
                    F = scatteredInterpolant(SENSIBILITE_DATA(:,1),SENSIBILITE_DATA(:,2),SENSIBILITE_DATA(:,3),SENSIBILITE_DATA(:,4),'natural');
                    vq = F(X,Y,Z);
                    ValAffichee = vq;
                    
                    % Figure editing
                    daspect([1,1,1])
                    axis tight
                    ax = gca;
                    ax.FontSize = 13;
                    %view(0,0)
                    view(-37.5,30)
                    light('Position',[-1 -1 0.75],'Style','infinite')
                    %camzoom(1.8)
                    camproj perspective
                    colormap ((jet(10)));
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
                    
                    % Volume visualization (ONE cutting value)
                    Valeur=[0;100;200];
                    V_Max=max(SENSIBILITE_DATA(:,4));
                    V_Min=min(SENSIBILITE_DATA(:,4));
                    Valeur_Indice = round((Valeur-V_Min)*10/(V_Max-V_Min))+1;
                    for i=1:length(Valeur_Indice(:,1))
                        p1=patch(isocaps(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
                            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.7);
                        p1 = patch(isosurface(X,Y,Z,ValAffichee,Valeur(i,1),'enclose','above'),...
                            'FaceColor',cmap(Valeur_Indice(i,1),:),'EdgeColor','none','FaceAlpha',0.7);
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
          
                    % Labels
                    xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
                    ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
                    zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
                    
                    title(['3D sensitivity improvement between the ',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),' optimized (',num2str(round(EFFICACITE_OPT,1)),'%) and standard protocol'],'FontSize',14)
                    
                    % Electrode visualization
                    for i=1:length(ELEC(:,1))
                        plot3(ELEC(i,1),ELEC(i,2),ELEC(i,3),'k.','MarkerSize',20)
                    end
                    
                    % Limits of the model
                    xlim([0 X_Nb_ELEC*X_GRID_SPACING+X_Elec_1])
                    ylim([0 Y_Nb_ELEC*Y_GRID_SPACING+Y_Elec_1])
                    zlim([-Z_MAX 0])
                    
                    
                    % Legend
                    clear colorbar
                    colorbar
                    caxis([0 300]);
                    TicksM = [0;50;100;200;300];
                    TicksL = TicksM;
                    
                    c=colorbar;
                    c.Label.String = 'Sensitivity improvement (%)';
                    c.Label.FontWeight='bold';
                    c.Label.FontSize=14;
                    c.Ticks=TicksM;
                    c.TickLabels={num2str(TicksL)};
                end
            end
            
            clear actif actif_optimized_protocol_identification ax c cmap F i Mat1 Mat2 o o_2_10_Abberant_configurations_filtering o_2_6_Cross_diagonal_dipole_dipole_with_all_the_lines_in_X_dire o_2_7_Cross_diagonal_dipole_dipole_with_all_the_lines_in_Y_dire
            clear p1 TicksL TicksM V_Max V_Min ValAffichee Valeur Valeur_Indice vq 
        end % vue 2 of the difference of resistivities
        
    end
end

%% V) Automated redaction of the SPREAD and PROTOCOL files (direct and reciprocal) 
%--------------------------------------------------------------------------
actif_Automated_redaction_of_the_SPREAD_and_PROTOCOL_files=1;
for o=1
    if actif_Automated_redaction_of_the_SPREAD_and_PROTOCOL_files==1
        %% 5.1) Selection of the electrode setup on the field
        % Creation : 17/05/2017
        actif=1;
        for o=1
            if actif==1
                CABLE_MATRIX_ind_cable=zeros(X_Nb_ELEC,Y_Nb_ELEC);
                CABLE_MATRIX_ind_electrode=zeros(X_Nb_ELEC,Y_Nb_ELEC);
                ELEC_MATRIX_FIELD=zeros(X_Nb_ELEC,Y_Nb_ELEC);
                % Selection of the number of cables used
                Nb_Cables = menu([num2str(N_Elec),' electrodes are beeing used, how many cables do you want to use ?'],'1 cable : capacity = [ 1 - 32 ]','2 cables : capacity = [ 1 - 64 ]','3 cables : capacity = [ 1 - 96 ]','4 cables : capacity = [ 1 - 128 ]','5 cables : capacity = [ 1 - 164 ]','6 cables : capacity = [ 1 - 192 ]');
                Nb_C = [1,2,3,4,5,6];
                Nb_Cables=Nb_C(1,Nb_Cables);
                
                if Nb_Cables>2
                    % Warning Message
                     % Selection of the number of cables used
                     Warning = menu(['You chose to use ',num2str(Nb_Cables),' cables which means you need to use ',num2str(floor((Nb_Cables+1)/2)-1),' ES1064C'],'OK');
                     Warning = [1];
                end
                    
                % Electrode visualization
                for o=1:1
                    figure('Color', [ 1 1 1])
                    plot(ELEC(:,1),ELEC(:,2),'.k','MarkerSize',20)
                    hold on
                    xlim([0 X_Nb_ELEC*X_GRID_SPACING+X_Elec_1])
                    ylim([0 Y_Nb_ELEC*Y_GRID_SPACING+Y_Elec_1+5])
                    xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
                    ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
                    daspect([1,1,1])
                    title(['Disposition of electrodes on the field for the protocol ',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC)],'FontSize',14)
                end
                
                color_sel{1}='b';
                color_sel{2}='g';
                color_sel{3}='m';
                color_sel{4}='c';
                color_sel{5}='y';
                color_sel{6}='r';
                
                % Selection of the electrodes organization for each cable
                ELEC_MATRIX_FIELD_ind=1;
                for i=1:Nb_Cables
                    % Warning Message
                    prompt = {sprintf(['You are going to choose the electrode distribution for the cable N°%d \n \n How many electrodes do you want to select for the cable N°%d ?  \n \n '],i,i)};
                    dlg_title = 'Number of electrodes';
                    num_lines = 1;
                    defaultans = {'0'};
                    Nb_E = inputdlg(prompt,dlg_title,num_lines,defaultans);
                    Nb_Electrodes(i,1) = str2num(Nb_E{:});
                    
                    Warning = menu(sprintf(['Are you sure that you want to use %d electrodes on the cable N°%d ?  \n \n '],Nb_Electrodes(i,1) ,i),'YES','NO');
                   
                    if Warning==2
                        % Warning Message
                    prompt = {sprintf(['You are going to choose the electrode distribution for the cable N°%d \n \n How many electrodes do you want to select for the cable N°%d ?  \n \n '],i,i)};
                    dlg_title = 'Number of electrodes';
                    num_lines = 1;
                    defaultans = {'0'};
                    Nb_E = inputdlg(prompt,dlg_title,num_lines,defaultans);
                    Nb_Electrodes(i,1) = str2num(Nb_E{:});
                    
                    end
                    
                    % Electrode selection
                    count=0;
                    for j=1:Nb_Electrodes(i,1)
                        % Warning Message
                        if count==0
                            % Selection of the jth electrode
                            Warning = menu(sprintf(['Please click on the electrode N°%d of the cable N°%d ?  \n \n '],j,i),'OK');
                            Warning = [1];
                        end
                        
                        ELEC_SELECTED=ginput(1);
                        
                        % Identification of the closest electrode
                        clear Dist
                        for k=1:N_Elec
                            Dist(k)=pdist([ELEC_SELECTED;ELEC(k,1:2)],'euclidean');
                        end
                        [~,ind_elec]=min(Dist);
                        
                        [ind_row,ind_column]=find(ELEC_MATRIX(:,:)==ind_elec);
                        
                        % Updating of the CABLE_MATRIX matrix
                        CABLE_MATRIX_ind_cable(ind_row,ind_column)=i;
                        CABLE_MATRIX_ind_electrode(ind_row,ind_column)=j;
                        ELEC_MATRIX_FIELD(ind_row,ind_column)=ELEC_MATRIX_FIELD_ind;
                        ELEC_MATRIX_FIELD_ind=ELEC_MATRIX_FIELD_ind+1;
                        
                        plot(ELEC(ind_elec,1),ELEC(ind_elec,2),['.',color_sel{i}],'MarkerSize',25)
                        text(ELEC(ind_elec,1)+0.1,ELEC(ind_elec,2)+1,[num2str(i),' - ',num2str(j)],'Color','black','FontSize',12,'Rotation',90)
                        
                        if count>0
                            line([X_line_mem ELEC(ind_elec,1)],[Y_line_mem ELEC(ind_elec,2)],'LineWidth',0.5,'Color',color_sel{i})
                            X_line_mem=ELEC(ind_elec,1);
                            Y_line_mem=ELEC(ind_elec,2);
                        else
                            X_line_mem=ELEC(ind_elec,1);
                            Y_line_mem=ELEC(ind_elec,2);
                            count=1;
                        end
                    end
                end
                clear actif color_sel 
            end
        end

        %% 5.2) Preparation for the redaction of spread and protocol files
        % Creation : 17/05/2017
        actif=1;
        for o=1
            if actif==1
                % Matrix to re-write protocols according to the setup on
                % the field
                for i=1:N_Elec
                    ELEC_PASSAGE_MATRIX(i,1)=i;
                    [xi,yi]=find(ELEC_MATRIX(:,:)==i);
                    ELEC_PASSAGE_MATRIX(i,2)=ELEC_MATRIX_FIELD(xi,yi);
                end
                
                % New writing of the identified optimized protocol
                PROTOCOLE_LIMITE_FIELD=PROTOCOLE_LIMITE;
                for i=1:length(PROTOCOLE_LIMITE(:,1)) 
                    for j=1:4
                        PROTOCOLE_LIMITE_FIELD(i,j)=ELEC_PASSAGE_MATRIX(PROTOCOLE_LIMITE(i,j),2);
                    end
                end
                
                % New matrix of electrode coordinates on the field
                ELEC_FIELD=ELEC;
                for i=1:N_Elec
                    for j=1:3
                        ELEC_FIELD(i,j)=ELEC(ELEC_PASSAGE_MATRIX(i,2),j);
                    end
                end 
            end
        end 
        
        %% 5.3) Redaction of the SPREAD file
        % Creation : 18/05/2017
        actif=1;
        for o=1
            if actif==1
                % Name of the spread file
                Name_Spread=['ADT3_S_',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC)];
                
                % Description of the spread file
                Description_Spread=['AD Lac TIO : Terrain N3 : Spread avec ',num2str(N_Elec),' electrodes. Mesure pour grille de ',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC)];
                
                % Number of Switch used
                Nb_SWITCH=floor(Nb_Cables/2);
                
                % Number of Cables used
                Nb_CABLES=Nb_Cables;
                
                % Tab
                tab='     ';
                
                % Creation of the .txt file
                fileID = fopen([Name_Spread,'.xml'],'w');
                
                % <?xml version="1.0" encoding="UTF-8" ?>
                fprintf(fileID,'%6s \r\n','<?xml version="1.0" encoding="UTF-8" ?>');
                
                % New Line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Spread>
                fprintf(fileID,'%6s \r\n','<Spread>');
                
                % <Name>
                fprintf(fileID,'%6s \r\n',[tab,'<Name> ',Name_Spread,' </Name>']);
                
                % New Line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Description>
                fprintf(fileID,'%6s \r\n',[tab,'<Description> ',Description_Spread,' </Description>']);
                
                % New Line
                fprintf(fileID,'%6s \r\n',' ');
                
                ind_cable=1;
                ind_elec=1;
                elec_previous_cables=1;
                ind_elec_SwitchAddress=1;
                while ind_cable < Nb_CABLES+1
                    if mod(ind_cable+1,2)==1 & ind_elec_SwitchAddress<33
                        ind_elec_SwitchAddress=33;
                    end
                    if mod(ind_cable+1,2)==0 & ind_elec_SwitchAddress>32
                        ind_elec_SwitchAddress=1;
                    end
                    
                    % <Cable>
                    fprintf(fileID,'%6s \r\n',[tab,'<Cable>']);
                    
                    % <Name>
                    fprintf(fileID,'%6s \r\n',[tab,tab,'<Name> ',num2str(ind_cable),' </Name> ']);
                    
                    % New Line
                    fprintf(fileID,'%6s \r\n',' ');
                    
                    while ind_elec<elec_previous_cables+Nb_Electrodes(ind_cable) && ind_elec<N_Elec+1
                        % <Electrode>
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Electrode>']);
                        
                        % <Id>
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Id> ',num2str(ind_elec),' </Id>']);
                        
                        % <X>
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<X> ',num2str(ind_elec),' </X>']);
                        
                        % <Name>
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Name> ',num2str(ind_elec),' </Name>']);
                        
                        % <SwitchId>
                        % if floor((ind_cable-1)/2)>0
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<SwitchId> ',num2str(floor((ind_cable-1)/2)+1),' </SwitchId>']);
                        %end
                        
                        % <SwitchAddress>
                        % from 1 to 64 for each IdSwitch
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<SwitchAddress> ',num2str(ind_elec_SwitchAddress),' </SwitchAddress>']);
                        
                        % <Electrode>
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'</Electrode>']);
                        
                        % New line
                        fprintf(fileID,'%6s \r\n',' ');
                        
                        ind_elec=ind_elec+1;
                        ind_elec_SwitchAddress=ind_elec_SwitchAddress+1;
                        
                        
                    end
                    elec_previous_cables=elec_previous_cables+Nb_Electrodes(ind_cable);
                    
                    % </Cable>
                    fprintf(fileID,'%6s \r\n',[tab,'</Cable>']);
                    
                    % New Line
                    fprintf(fileID,'%6s \r\n',' ');
                    
                    ind_cable=ind_cable+1;
                end
                
                % </Spread>
                fprintf(fileID,'%6s \r\n','</Spread>');
                
                fclose(fileID);
            end
        end

        %% 5.4) Redaction of the DIRECT PROTOCOL file
        % Creation : 18/05/2017
        actif=1;
        for o=1
            if actif==1
                PROTOCOLE=PROTOCOLE_LIMITE_FIELD;
                % Injections : TX
                TX=unique(PROTOCOLE(:,1:2),'rows');
                
                % Name of the protocol file
                Name_Protocole_D=['ADT3_P_D_',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),'_',num2str(length(PROTOCOLE(:,1)))];
                
                % Description of the protocol file
                Description_Protocole=['AD Lac TIO : Terrain N3 : Protocole optimise : ',num2str(N_Elec),' elecs et ',num2str(length(PROTOCOLE(:,1))),' configs'];
                
                % tab
                tab='     ';
                
                % Creation of the file .txt
                fileID = fopen([Name_Protocole_D,'.xml'],'w');
                
                % <?xml version="1.0" encoding="UTF-8" ?>
                fprintf(fileID,'%6s \r\n','<?xml version="1.0" encoding="UTF-8" ?>');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Protocol>
                fprintf(fileID,'%6s \r\n','<Protocol>');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Name>
                fprintf(fileID,'%6s \r\n',[tab,'<Name> ',Name_Protocole_D,' </Name>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Description>
                fprintf(fileID,'%6s \r\n',[tab,'<Description> ',Description_Protocole,' </Description>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Arraycode>
                fprintf(fileID,'%6s \r\n',[tab,'<Arraycode> 11 </Arraycode>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <SpreadFile>
                fprintf(fileID,'%6s \r\n',[tab,'<SpreadFile> ',Name_Spread,'.xml',' </SpreadFile>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Sequence>
                fprintf(fileID,'%6s \r\n',[tab,'<Sequence>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                for i=1:length(TX(:,1))
                    clear RX
                    RX=PROTOCOLE(find(PROTOCOLE(:,1)==TX(i,1)&PROTOCOLE(:,2)==TX(i,2)),3:4);
                    
                    % <Measure>
                    fprintf(fileID,'%6s \r\n',[tab,tab,'<Measure>']);
                    
                    % <Tx>
                    fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Tx> ',num2str(TX(i,1)),' ',num2str(TX(i,2)),' </Tx>']);
                    
                    for j=1:length(RX(:,1))
                        % <Rx>
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Rx> ',num2str(RX(j,1)),' ',num2str(RX(j,2)),' </Rx>']);
                    end
                    
                    % </Measure>
                    fprintf(fileID,'%6s \r\n',[tab,tab,'</Measure>']);
                    
                    % New line
                    fprintf(fileID,'%6s \r\n',' ');
                end
                
                % </Sequence>
                fprintf(fileID,'%6s \r\n',[tab,'</Sequence>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % </Protocol>
                fprintf(fileID,'%6s \r\n',['</Protocol>']);
            end
        end
        
        %% 5.5) Redaction of the RECIPROQUAL PROTOCOL file
        % Creation : 18/05/2017
        actif=1;
        for o=1
            if actif==1
                PROTOCOLE=PROTOCOLE_LIMITE_FIELD;
                % Injections : TX (reciproqual)
                TX=unique(PROTOCOLE(:,3:4),'rows');
                
                % Name of the protocol file
                Name_Protocole_R=['ADT3_P_R_',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),'_',num2str(length(PROTOCOLE(:,1)))];
                
                % Description of the protocol file
                Description_Protocole=['AD Lac TIO : Terrain N3 : Protocole optimise reciproqual: ',num2str(N_Elec),' elecs et ',num2str(length(PROTOCOLE(:,1))),' configs'];
                
                % tab
                tab='     ';
                
                % Creation of the file .txt
                fileID = fopen([Name_Protocole_R,'.xml'],'w');
                
                % <?xml version="1.0" encoding="UTF-8" ?>
                fprintf(fileID,'%6s \r\n','<?xml version="1.0" encoding="UTF-8" ?>');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Protocol>
                fprintf(fileID,'%6s \r\n','<Protocol>');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Name>
                fprintf(fileID,'%6s \r\n',[tab,'<Name> ',Name_Protocole_R,' </Name>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Description>
                fprintf(fileID,'%6s \r\n',[tab,'<Description> ',Description_Protocole,' </Description>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Arraycode>
                fprintf(fileID,'%6s \r\n',[tab,'<Arraycode> 11 </Arraycode>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <SpreadFile>
                fprintf(fileID,'%6s \r\n',[tab,'<SpreadFile> ',Name_Spread,'.xml',' </SpreadFile>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Sequence>
                fprintf(fileID,'%6s \r\n',[tab,'<Sequence>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                for i=1:length(TX(:,1))
                    clear RX %(reciproqual)
                    RX=PROTOCOLE(find(PROTOCOLE(:,3)==TX(i,1)&PROTOCOLE(:,4)==TX(i,2)),1:2);
                    
                    % <Measure>
                    fprintf(fileID,'%6s \r\n',[tab,tab,'<Measure>']);
                    
                    % <Tx>
                    fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Tx> ',num2str(TX(i,1)),' ',num2str(TX(i,2)),' </Tx>']);
                    
                    for j=1:length(RX(:,1))
                        % <Rx>
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Rx> ',num2str(RX(j,1)),' ',num2str(RX(j,2)),' </Rx>']);
                    end
                    
                    % </Measure>
                    fprintf(fileID,'%6s \r\n',[tab,tab,'</Measure>']);
                    
                    % New line
                    fprintf(fileID,'%6s \r\n',' ');
                end
                
                % </Sequence>
                fprintf(fileID,'%6s \r\n',[tab,'</Sequence>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % </Protocol>
                fprintf(fileID,'%6s \r\n',['</Protocol>']);
            end
            
            
            % Warning Message
            % End of the 3D VERTICAL SCAN programm
            Warning = menu(sprintf('Dear user, we are pleased to inform you that 3D VERTICAL SCAN executed without errors : you are ready to transfer the SPREAD and PROTOCOL files to the Terrameter LS device ! \r \n  \n Your files are : \r \n \n - spread file : %6s \n \r - direct protocol file : %6s \n \r - reverse protocol file : %6s \r\n \n Have a nice field work ! \r\n \r\n \n \n Adrien Dimech ',Name_Spread,Name_Protocole_D,Name_Protocole_R),'Thanks');
            Warning = [1];
            
        end
        
        %% 5.6) Preparation for the redaction of spread and standard protocol files
        % Creation : 02/06/2017
        actif=1;
        for o=1
            if actif==1 & standard_exist==1
                % New writing of the identified optimized protocol
                PRO_STANDARD_FIELD=PRO_STANDARD;
                for i=1:length(PRO_STANDARD(:,1)) 
                    for j=1:4
                        PRO_STANDARD_FIELD(i,j)=ELEC_PASSAGE_MATRIX(PRO_STANDARD(i,j),2);
                    end
                end
            end
        end 
        
        %% 5.7) Redaction of the STANDARD DIRECT PROTOCOL file
        % Creation : 02/06/2017
        actif=1;
        for o=1
            if actif==1 & standard_exist==1
                PROTOCOLE=PRO_STANDARD_FIELD;
                % Injections : TX
                TX=unique(PROTOCOLE(:,1:2),'rows');
                
                % Name of the protocol file
                Name_Protocole_D=['ADT3_P_stand_D_',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),'_',num2str(length(PROTOCOLE(:,1)))];
                
                % Description of the protocol file
                Description_Protocole=['AD Lac TIO : Terrain N3 : Protocole standard : ',num2str(N_Elec),' elecs et ',num2str(length(PROTOCOLE(:,1))),' configs'];
                
                % tab
                tab='     ';
                
                % Creation of the file .txt
                fileID = fopen([Name_Protocole_D,'.xml'],'w');
                
                % <?xml version="1.0" encoding="UTF-8" ?>
                fprintf(fileID,'%6s \r\n','<?xml version="1.0" encoding="UTF-8" ?>');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Protocol>
                fprintf(fileID,'%6s \r\n','<Protocol>');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Name>
                fprintf(fileID,'%6s \r\n',[tab,'<Name> ',Name_Protocole_D,' </Name>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Description>
                fprintf(fileID,'%6s \r\n',[tab,'<Description> ',Description_Protocole,' </Description>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Arraycode>
                fprintf(fileID,'%6s \r\n',[tab,'<Arraycode> 11 </Arraycode>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <SpreadFile>
                fprintf(fileID,'%6s \r\n',[tab,'<SpreadFile> ',Name_Spread,'.xml',' </SpreadFile>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Sequence>
                fprintf(fileID,'%6s \r\n',[tab,'<Sequence>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                for i=1:length(TX(:,1))
                    clear RX
                    RX=PROTOCOLE(find(PROTOCOLE(:,1)==TX(i,1)&PROTOCOLE(:,2)==TX(i,2)),3:4);
                    
                    % <Measure>
                    fprintf(fileID,'%6s \r\n',[tab,tab,'<Measure>']);
                    
                    % <Tx>
                    fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Tx> ',num2str(TX(i,1)),' ',num2str(TX(i,2)),' </Tx>']);
                    
                    for j=1:length(RX(:,1))
                        % <Rx>
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Rx> ',num2str(RX(j,1)),' ',num2str(RX(j,2)),' </Rx>']);
                    end
                    
                    % </Measure>
                    fprintf(fileID,'%6s \r\n',[tab,tab,'</Measure>']);
                    
                    % New line
                    fprintf(fileID,'%6s \r\n',' ');
                end
                
                % </Sequence>
                fprintf(fileID,'%6s \r\n',[tab,'</Sequence>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % </Protocol>
                fprintf(fileID,'%6s \r\n',['</Protocol>']);
            end
        end
        
        %% 5.5) Redaction of the STANDARD RECIPROQUAL PROTOCOL file
        % Creation : 02/06/2017
        actif=1;
        for o=1
            if actif==1 & standard_exist==1
                PROTOCOLE=PRO_STANDARD_FIELD;
                % Injections : TX (reciproqual)
                TX=unique(PROTOCOLE(:,3:4),'rows');
                
                % Name of the protocol file
                Name_Protocole_R=['ADT3_P_stand_R_',num2str(X_Nb_ELEC),'x',num2str(Y_Nb_ELEC),'_',num2str(length(PROTOCOLE(:,1)))];
                
                % Description of the protocol file
                Description_Protocole=['AD Lac TIO : Terrain N3 : Protocole standard reciproqual: ',num2str(N_Elec),' elecs et ',num2str(length(PROTOCOLE(:,1))),' configs'];
                
                % tab
                tab='     ';
                
                % Creation of the file .txt
                fileID = fopen([Name_Protocole_R,'.xml'],'w');
                
                % <?xml version="1.0" encoding="UTF-8" ?>
                fprintf(fileID,'%6s \r\n','<?xml version="1.0" encoding="UTF-8" ?>');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Protocol>
                fprintf(fileID,'%6s \r\n','<Protocol>');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Name>
                fprintf(fileID,'%6s \r\n',[tab,'<Name> ',Name_Protocole_R,' </Name>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Description>
                fprintf(fileID,'%6s \r\n',[tab,'<Description> ',Description_Protocole,' </Description>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Arraycode>
                fprintf(fileID,'%6s \r\n',[tab,'<Arraycode> 11 </Arraycode>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <SpreadFile>
                fprintf(fileID,'%6s \r\n',[tab,'<SpreadFile> ',Name_Spread,'.xml',' </SpreadFile>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % <Sequence>
                fprintf(fileID,'%6s \r\n',[tab,'<Sequence>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                for i=1:length(TX(:,1))
                    clear RX %(reciproqual)
                    RX=PROTOCOLE(find(PROTOCOLE(:,3)==TX(i,1)&PROTOCOLE(:,4)==TX(i,2)),1:2);
                    
                    % <Measure>
                    fprintf(fileID,'%6s \r\n',[tab,tab,'<Measure>']);
                    
                    % <Tx>
                    fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Tx> ',num2str(TX(i,1)),' ',num2str(TX(i,2)),' </Tx>']);
                    
                    for j=1:length(RX(:,1))
                        % <Rx>
                        fprintf(fileID,'%6s \r\n',[tab,tab,tab,'<Rx> ',num2str(RX(j,1)),' ',num2str(RX(j,2)),' </Rx>']);
                    end
                    
                    % </Measure>
                    fprintf(fileID,'%6s \r\n',[tab,tab,'</Measure>']);
                    
                    % New line
                    fprintf(fileID,'%6s \r\n',' ');
                end
                
                % </Sequence>
                fprintf(fileID,'%6s \r\n',[tab,'</Sequence>']);
                
                % New line
                fprintf(fileID,'%6s \r\n',' ');
                
                % </Protocol>
                fprintf(fileID,'%6s \r\n',['</Protocol>']);
            end
            
            
            % Warning Message
            % End of the 3D VERTICAL SCAN programm
            Warning = menu(sprintf('Dear user, we are pleased to inform you that 3D VERTICAL SCAN executed without errors : you are ready to transfer the SPREAD and PROTOCOL files to the Terrameter LS device ! \r \n  \n Your files are : \r \n \n - spread file : %6s \n \r - direct protocol file : %6s \n \r - reverse protocol file : %6s \r\n \n Have a nice field work ! \r\n \r\n \n \n Adrien Dimech ',Name_Spread,Name_Protocole_D,Name_Protocole_R),'Thanks');
            Warning = [1];
            
        end
        
        
    end
end
