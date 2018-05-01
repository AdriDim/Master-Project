
%% --------------- 3D-modeling-of-the-pile-for-COMSOL ---------------------
%
% Adrien Dimech - Master Project - 21/04/2018
%
% -------------------------------------------------------------------------
% Matlab codes to prepare - process - visualize and interpret 3D time-lapse
% geolelectrical monitoring of a waste rock pile.
% -------------------------------------------------------------------------
%
% This Matlab code was designed to create a 3D model of the waste rock pile
% for the COMSOL Multiphysics modelisation software. External topography
% and internal structure was used to calculate geometric factors for both
% standard and optimized protocols. The approach developped to model the
% pile can be applied to any complex structure to generate a complex 3D
% high-resolution model for COMSOL Multiphysics with LiveLink for Matlab.
%
% Feel free to visit : https://www.researchgate.net/profile/Adrien_Dimech
% for more information about my research or contact me for more information
% and data files : adrien.dimech@gmail.com
%
%%% MODELISATION DE LA HALDE AVEC COMSOL
% -------------------------------------------------------------------------
%
%
%                           Suivi du code
% -------------------------------------------------------------------------
% Creation          |       01/11/2016        |     Eric T.K. Chou
% -------------------------------------------------------------------------
% Modification      |       09/11/2016        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       13/12/2016        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       12/01/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       18/01/2017        |     Adrien Dimech
% -------------------------------------------------------------------------
% Modification      |       03/02/2017        |     Adrien Dimech
% -------------------------------------------------------------------------

%% 0) Création du modèle COMSOL
% Creation le 01/11/2016
for o=1
    clear all
    figure
    p = get(gcf,'Position');
    set(0,'DefaultFigurePosition',p);
    close all
    clc
    
    import com.comsol.model.*
    import com.comsol.model.util.*
    model = ModelUtil.create('Model');
    
    % Connect the Model in the COMSOL Desktop
    ModelUtil.showProgress(true);
    ModelUtil.setServerBusyHandler(ServerBusyHandler(2));
end

%% 1.1) Création de la géométrie de la halde
% Creation le 01/11/2016
% Modification le 09/11/2016
% Modification le 13/12/2016
% Modification le 12/01/2017
% Modification le 03/02/2017
for o=1
    model.modelNode.create('comp1');
    
    model.geom.create('geom1', 3);
    model.mesh.create('mesh1', 'geom1');
    model.view('view1').set('transparency', 'on');
    
    for o=1:1   % Définition de la géométrie de la halde = 'dif1'
        % Definition du bloc principal
        model.geom('geom1').create('blk1', 'Block');
        model.geom('geom1').feature('blk1').set('size', {'110' '50' '25'});
        model.geom('geom1').feature('blk1').set('pos', {'0' '0' '-5'});
        model.geom('geom1').run('blk1');
        
        % Chargement de la géométrie de la halde
        load 'Maillage3D_V2.mat'
        Maillage3D=Maillage3D_V2;
        
        ModelUtil.showProgress(true);
        L = length(Maillage3D(:,1))/3;
        
        % Maillage tétraédrique de la halde
        TETRA=0;
        count=0;
        % Maillage primaire
        for i=1:L
            Name = ['tet',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            count=count+1;
            for j=1:3
                for k=1:3
                    % Modification le 13/12/2016
                    TETRA((count-1)*4+j,k)=Maillage3D((i-1)*3+j,k);
                    Coordstr = num2str(Maillage3D((i-1)*3+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            for j=1:2
                Coordstr = num2str(Maillage3D((i-1)*3+1,j));
                % Modification le 13/12/2016
                TETRA((count-1)*4+4,j)=Maillage3D((i-1)*3+1,j);
                model.geom('geom1').feature(Name).setIndex('p', Coordstr , j-1, 3);
            end
            % Modification le 13/12/2016
            TETRA((count-1)*4+4,3)=20;
            model.geom('geom1').feature(Name).setIndex('p', '20' , 2, 3);
            
            clear Coordnum
            model.geom('geom1').run(Name);
        end
        % Maillage secondaire
        TETRACOMPL1=TETRA;
        TETRACOMPL2=TETRA;
        for i=1:L
            TETRACOMPL1((i-1)*4+1,3)=20;
            TETRACOMPL1((i-1)*4+4,:)=TETRA((i-1)*4+3,:);
            TETRACOMPL1((i-1)*4+4,3)=20;
            
            TETRACOMPL2((i-1)*4+1,3)=20;
            TETRACOMPL2((i-1)*4+2,3)=20;
            TETRACOMPL2((i-1)*4+3,3)=20;
            TETRACOMPL2((i-1)*4+4,:)=TETRA((i-1)*4+2,:);
        end
        
        for i=1:L
            Name = ['tetcomp1_',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL1((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        for i=1:L
            Name = ['tetcomp2_',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL2((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        % Création du volume par différence
        selection1 = {'tet1'};
        selection2 = {'tetcomp1_1'};
        selection3 = {'tetcomp2_1'};
        for i=2:L
            selection1{1,i} = ['tet',num2str(i)];
            selection2{1,i} = ['tetcomp1_',num2str(i)];
            selection3{1,i} = ['tetcomp2_',num2str(i)];
        end
        selectionTET1=[selection1,selection2,selection3];
        
        model.geom('geom1').create('dif1', 'Difference');
        model.geom('geom1').feature('dif1').selection('input').set({'blk1'});
        model.geom('geom1').feature('dif1').selection('input2').set(selectionTET1);
        model.geom('geom1').feature('dif1').set('keep', 'on');
        model.geom('geom1').run('dif1');
    end
    
    for o=1:1   % Définition de la couche de sable + anorthosite = 'dif2'
        
        % Chargement de la géométrie de la halde
        load 'Maillage3DSteriles_Sable2.mat'
        Maillage3DSab_Ste=Maillage3DSteriles_Sable2;
        
        LSab_Ste = length(Maillage3DSab_Ste(:,1))/3;
        
        % Maillage tétraédrique de la halde
        TETRA=0;
        count=0;
        % Maillage primaire
        for i=1:LSab_Ste
            Name = ['tetSab_Ste',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            count=count+1;
            for j=1:3
                for k=1:3
                    TETRA((count-1)*4+j,k)=Maillage3DSab_Ste((i-1)*3+j,k);
                    Coordstr = num2str(Maillage3DSab_Ste((i-1)*3+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            for j=1:2
                Coordstr = num2str(Maillage3DSab_Ste((i-1)*3+1,j));
                TETRA((count-1)*4+4,j)=Maillage3DSab_Ste((i-1)*3+1,j);
                model.geom('geom1').feature(Name).setIndex('p', Coordstr , j-1, 3);
            end
            TETRA((count-1)*4+4,3)=-5;
            model.geom('geom1').feature(Name).setIndex('p', '-5' , 2, 3);
            
            clear Coordnum
            model.geom('geom1').run(Name);
        end
        % Maillage secondaire
        TETRACOMPL1=TETRA;
        TETRACOMPL2=TETRA;
        for i=1:LSab_Ste
            TETRACOMPL1((i-1)*4+1,3)=-5;
            TETRACOMPL1((i-1)*4+4,:)=TETRA((i-1)*4+3,:);
            TETRACOMPL1((i-1)*4+4,3)=-5;
            
            TETRACOMPL2((i-1)*4+1,3)=-5;
            TETRACOMPL2((i-1)*4+2,3)=-5;
            TETRACOMPL2((i-1)*4+3,3)=-5;
            TETRACOMPL2((i-1)*4+4,:)=TETRA((i-1)*4+2,:);
        end
        
        for i=1:LSab_Ste
            Name = ['tetcomp1_Sab_Ste',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL1((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        for i=1:LSab_Ste
            Name = ['tetcomp2_Sab_Ste',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL2((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        % Création du volume par différence
        clear selection
        selection1 = {'tetSab_Ste1'};
        selection2 = {'tetcomp1_Sab_Ste1'};
        selection3 = {'tetcomp2_Sab_Ste1'};
        for i=2:LSab_Ste
            selection1{1,i} = ['tetSab_Ste',num2str(i)];
            selection2{1,i} = ['tetcomp1_Sab_Ste',num2str(i)];
            selection3{1,i} = ['tetcomp2_Sab_Ste',num2str(i)];
        end
        selectionTET2=[selection1,selection2,selection3];
        
        model.geom('geom1').create('dif2', 'Difference');
        model.geom('geom1').feature('dif2').selection('input').set({'dif1'});
        model.geom('geom1').feature('dif2').selection('input2').set(selectionTET2);
        model.geom('geom1').feature('dif2').set('keep', 'on');
        model.geom('geom1').run('dif2');
        
    end
    
    for o=1:1   % Définition de la couche d'anorthosite = 'Anorthosite'
        % Chargement de la géométrie de la halde
        load 'Maillage3DSable1_Anortho.mat'
        Maillage3DS_A=Maillage3DSable1_Anortho;
        
        LS_A = length(Maillage3DS_A(:,1))/3;
        
        % Maillage tétraédrique de la halde
        TETRA=0;
        count=0;
        % Maillage primaire
        for i=1:LS_A
            Name = ['tetS_A',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            count=count+1;
            for j=1:3
                for k=1:3
                    TETRA((count-1)*4+j,k)=Maillage3DS_A((i-1)*3+j,k);
                    Coordstr = num2str(Maillage3DS_A((i-1)*3+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            for j=1:2
                Coordstr = num2str(Maillage3DS_A((i-1)*3+1,j));
                TETRA((count-1)*4+4,j)=Maillage3DS_A((i-1)*3+1,j);
                model.geom('geom1').feature(Name).setIndex('p', Coordstr , j-1, 3);
            end
            TETRA((count-1)*4+4,3)=-5;
            model.geom('geom1').feature(Name).setIndex('p', '-5' , 2, 3);
            
            clear Coordnum
            model.geom('geom1').run(Name);
        end
        % Maillage secondaire
        TETRACOMPL1=TETRA;
        TETRACOMPL2=TETRA;
        for i=1:LS_A
            TETRACOMPL1((i-1)*4+1,3)=-5;
            TETRACOMPL1((i-1)*4+4,:)=TETRA((i-1)*4+3,:);
            TETRACOMPL1((i-1)*4+4,3)=-5;
            
            TETRACOMPL2((i-1)*4+1,3)=-5;
            TETRACOMPL2((i-1)*4+2,3)=-5;
            TETRACOMPL2((i-1)*4+3,3)=-5;
            TETRACOMPL2((i-1)*4+4,:)=TETRA((i-1)*4+2,:);
        end
        
        for i=1:LS_A
            Name = ['tetcomp1_S_A',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL1((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        for i=1:LS_A
            Name = ['tetcomp2_S_A',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL2((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        % Création du volume par différence
        selection1 = {'tetS_A1'};
        selection2 = {'tetcomp1_S_A1'};
        selection3 = {'tetcomp2_S_A1'};
        for i=2:LS_A
            selection1{1,i} = ['tetS_A',num2str(i)];
            selection2{1,i} = ['tetcomp1_S_A',num2str(i)];
            selection3{1,i} = ['tetcomp2_S_A',num2str(i)];
        end
        selectionTET3=[selection1,selection2,selection3];
        
        model.geom('geom1').create('Anorthosite', 'Difference');
        model.geom('geom1').feature('Anorthosite').selection('input').set({'dif1'});
        model.geom('geom1').feature('Anorthosite').selection('input2').set(selectionTET3);
        model.geom('geom1').feature('Anorthosite').set('keep', 'on');
        model.geom('geom1').feature('Anorthosite').label('ANORTHOSITE');
        model.geom('geom1').run('Anorthosite');
        
    end
    
    for o=1:1   % Définition de la couche de sable du haut = 'SableH'
        model.geom('geom1').create('SableH', 'Difference');
        model.geom('geom1').feature('SableH').selection('input').set({'dif2'});
        model.geom('geom1').feature('SableH').selection('input2').set({'Anorthosite'});
        model.geom('geom1').feature('SableH').set('keep', 'on');
        model.geom('geom1').feature('SableH').label('SABLE_HAUT');
        model.geom('geom1').run('SableH');
    end
    
    for o=1:1   % Définition de la géométrie de la base de la halde = 'Base'
        
        % Chargement de la géométrie de la base de la halde
        load 'Maillage3DBase.mat'
        
        L_B = length(Maillage3DBase(:,1))/3;
        
        % Maillage tétraédrique de la halde
        TETRA=0;
        count=0;
        % Maillage primaire
        for i=1:L_B
            Name = ['tet_B',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            count=count+1;
            for j=1:3
                for k=1:3
                    % Modification le 13/12/2016
                    TETRA((count-1)*4+j,k)=Maillage3DBase((i-1)*3+j,k);
                    Coordstr = num2str(Maillage3DBase((i-1)*3+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            for j=1:2
                Coordstr = num2str(Maillage3DBase((i-1)*3+1,j));
                % Modification le 13/12/2016
                TETRA((count-1)*4+4,j)=Maillage3DBase((i-1)*3+1,j);
                model.geom('geom1').feature(Name).setIndex('p', Coordstr , j-1, 3);
            end
            % Modification le 13/12/2016
            TETRA((count-1)*4+4,3)=20;
            model.geom('geom1').feature(Name).setIndex('p', '20' , 2, 3);
            
            clear Coordnum
            model.geom('geom1').run(Name);
        end
        % Maillage secondaire
        TETRACOMPL1=TETRA;
        TETRACOMPL2=TETRA;
        for i=1:L_B
            TETRACOMPL1((i-1)*4+1,3)=20;
            TETRACOMPL1((i-1)*4+4,:)=TETRA((i-1)*4+3,:);
            TETRACOMPL1((i-1)*4+4,3)=20;
            
            TETRACOMPL2((i-1)*4+1,3)=20;
            TETRACOMPL2((i-1)*4+2,3)=20;
            TETRACOMPL2((i-1)*4+3,3)=20;
            TETRACOMPL2((i-1)*4+4,:)=TETRA((i-1)*4+2,:);
        end
        
        for i=1:L_B
            Name = ['tetcomp1__B',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL1((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        for i=1:L_B
            Name = ['tetcomp2__B',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL2((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        % Création du volume par différence
        selection1 = {'tet_B1'};
        selection2 = {'tetcomp1__B1'};
        selection3 = {'tetcomp2__B1'};
        for i=2:L_B
            selection1{1,i} = ['tet_B',num2str(i)];
            selection2{1,i} = ['tetcomp1__B',num2str(i)];
            selection3{1,i} = ['tetcomp2__B',num2str(i)];
        end
        selectionTET4=[selection1,selection2,selection3];
        
        model.geom('geom1').create('Base', 'Difference');
        model.geom('geom1').feature('Base').selection('input').set({'blk1'});
        model.geom('geom1').feature('Base').selection('input2').set(selectionTET4);
        model.geom('geom1').feature('Base').set('keep', 'on');
        model.geom('geom1').run('Base');
    end
    
    for o=1:1   % Définition de la couche de sable du bas = 'SableB'
        % Chargement de la géométrie de la halde
        load 'Maillage3DSable_Base.mat'
        Maillage3DS_Base=Maillage3DSable_Base;
        
        LS_Base = length(Maillage3DS_Base(:,1))/3;
        
        % Maillage tétraédrique de la halde
        TETRA=0;
        count=0;
        % Maillage primaire
        for i=1:LS_Base
            Name = ['tetS_Base',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            count=count+1;
            for j=1:3
                for k=1:3
                    TETRA((count-1)*4+j,k)=Maillage3DS_Base((i-1)*3+j,k);
                    Coordstr = num2str(Maillage3DS_Base((i-1)*3+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            for j=1:2
                Coordstr = num2str(Maillage3DS_Base((i-1)*3+1,j));
                TETRA((count-1)*4+4,j)=Maillage3DS_Base((i-1)*3+1,j);
                model.geom('geom1').feature(Name).setIndex('p', Coordstr , j-1, 3);
            end
            TETRA((count-1)*4+4,3)=-5;
            model.geom('geom1').feature(Name).setIndex('p', '-5' , 2, 3);
            
            clear Coordnum
            model.geom('geom1').run(Name);
        end
        % Maillage secondaire
        TETRACOMPL1=TETRA;
        TETRACOMPL2=TETRA;
        for i=1:LS_Base
            TETRACOMPL1((i-1)*4+1,3)=-5;
            TETRACOMPL1((i-1)*4+4,:)=TETRA((i-1)*4+3,:);
            TETRACOMPL1((i-1)*4+4,3)=-5;
            
            TETRACOMPL2((i-1)*4+1,3)=-5;
            TETRACOMPL2((i-1)*4+2,3)=-5;
            TETRACOMPL2((i-1)*4+3,3)=-5;
            TETRACOMPL2((i-1)*4+4,:)=TETRA((i-1)*4+2,:);
        end
        
        for i=1:LS_Base
            Name = ['tetcomp1_S_Base',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL1((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        for i=1:LS_Base
            Name = ['tetcomp2_S_Base',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL2((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        % Création du volume par différence
        selection1 = {'tetS_Base1'};
        selection2 = {'tetcomp1_S_Base1'};
        selection3 = {'tetcomp2_S_Base1'};
        for i=2:LS_Base
            selection1{1,i} = ['tetS_Base',num2str(i)];
            selection2{1,i} = ['tetcomp1_S_Base',num2str(i)];
            selection3{1,i} = ['tetcomp2_S_Base',num2str(i)];
        end
        selectionTET5=[selection1,selection2,selection3];
        
        model.geom('geom1').create('SableB', 'Difference');
        model.geom('geom1').feature('SableB').selection('input').set(selectionTET5);
        model.geom('geom1').feature('SableB').selection('input2').set('Base');
        model.geom('geom1').feature('SableB').set('keep', 'on');
        model.geom('geom1').feature('SableB').set('intbnd', 'off');
        model.geom('geom1').feature('SableB').label('SABLE_BAS');
        model.geom('geom1').run('SableB');
        
    end
    
    for o=1:1   % Définition du corps de la halde = 'Halde'
        model.geom('geom1').create('Halde', 'Difference');
        model.geom('geom1').feature('Halde').selection('input').set({'dif1'});
        model.geom('geom1').feature('Halde').selection('input2').set({'Anorthosite' 'SableH' 'SableB' 'Base'});
        model.geom('geom1').feature('Halde').set('keep', 'on');
        model.geom('geom1').feature('Halde').label('HALDE');
        model.geom('geom1').run('Halde');
    end
    
    for o=1:1   % Définition de la couche d'anorthosite = 'Anorthosite_Gross'
        % Chargement de la géométrie de la halde
        load 'Maillage3DAnortho_Gross.mat'
        
        LAnortho_Gross = length(Maillage3DAnortho_Gross(:,1))/3;
        
        % Maillage tétraédrique de la halde
        TETRA=0;
        count=0;
        % Maillage primaire
        for i=1:LAnortho_Gross
            Name = ['tetAnortho_Gross',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            count=count+1;
            for j=1:3
                for k=1:3
                    TETRA((count-1)*4+j,k)=Maillage3DAnortho_Gross((i-1)*3+j,k);
                    Coordstr = num2str(Maillage3DAnortho_Gross((i-1)*3+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            for j=1:2
                Coordstr = num2str(Maillage3DAnortho_Gross((i-1)*3+1,j));
                TETRA((count-1)*4+4,j)=Maillage3DAnortho_Gross((i-1)*3+1,j);
                model.geom('geom1').feature(Name).setIndex('p', Coordstr , j-1, 3);
            end
            TETRA((count-1)*4+4,3)=20;
            model.geom('geom1').feature(Name).setIndex('p', '20' , 2, 3);
            
            clear Coordnum
            model.geom('geom1').run(Name);
        end
        % Maillage secondaire
        TETRACOMPL1=TETRA;
        TETRACOMPL2=TETRA;
        for i=1:LAnortho_Gross
            TETRACOMPL1((i-1)*4+1,3)=20;
            TETRACOMPL1((i-1)*4+4,:)=TETRA((i-1)*4+3,:);
            TETRACOMPL1((i-1)*4+4,3)=20;
            
            TETRACOMPL2((i-1)*4+1,3)=20;
            TETRACOMPL2((i-1)*4+2,3)=20;
            TETRACOMPL2((i-1)*4+3,3)=20;
            TETRACOMPL2((i-1)*4+4,:)=TETRA((i-1)*4+2,:);
        end
        
        for i=1:LAnortho_Gross
            Name = ['tetcomp1_Anortho_Gross',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL1((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        for i=1:LAnortho_Gross
            Name = ['tetcomp2_Anortho_Gross',num2str(i)];
            model.geom('geom1').create(Name, 'Tetrahedron');
            model.geom('geom1').feature(Name).set('type', 'solid');
            for j=1:4
                for k=1:3
                    Coordstr = num2str(TETRACOMPL2((i-1)*4+j,k));
                    model.geom('geom1').feature(Name).setIndex('p', Coordstr , k-1, j-1);
                end
            end
            
            model.geom('geom1').run(Name);
        end
        
        % Création du volume par différence
        selection1 = {'tetAnortho_Gross1'};
        selection2 = {'tetcomp1_Anortho_Gross1'};
        selection3 = {'tetcomp2_Anortho_Gross1'};
        for i=2:LAnortho_Gross
            selection1{1,i} = ['tetAnortho_Gross',num2str(i)];
            selection2{1,i} = ['tetcomp1_Anortho_Gross',num2str(i)];
            selection3{1,i} = ['tetcomp2_Anortho_Gross',num2str(i)];
        end
        selectionTET6=[selection1,selection2,selection3];
        
        model.geom('geom1').create('dif3', 'Difference');
        model.geom('geom1').feature('dif3').selection('input').set(selectionTET6);
        model.geom('geom1').feature('dif3').selection('input2').set('Halde');
        model.geom('geom1').feature('dif3').set('keep', 'on');
        model.geom('geom1').feature('dif3').set('intbnd', 'off');
        model.geom('geom1').feature('dif3').label('dif3');
        model.geom('geom1').run('dif3');
        
        model.geom('geom1').create('Anorthosite_Gross', 'Difference');
        model.geom('geom1').feature('Anorthosite_Gross').selection('input').set(selectionTET6);
        model.geom('geom1').feature('Anorthosite_Gross').selection('input2').set('dif3');
        model.geom('geom1').feature('Anorthosite_Gross').set('keep', 'on');
        model.geom('geom1').feature('Anorthosite_Gross').set('intbnd', 'off');
        model.geom('geom1').feature('Anorthosite_Gross').label('ANORTHOSITE_GROSSIERE');
        model.geom('geom1').run('Anorthosite_Gross');
        
    end
    
    for o=1:1   % Définition du corps de la halde en ilmenite = 'Ilmenite'
        model.geom('geom1').create('Ilmenite', 'Difference');
        model.geom('geom1').feature('Ilmenite').selection('input').set('Halde');
        model.geom('geom1').feature('Ilmenite').selection('input2').set('Anorthosite_Gross');
        model.geom('geom1').feature('Ilmenite').set('keep', 'on');
        model.geom('geom1').feature('Ilmenite').label('ILMENITE');
        model.geom('geom1').run('Ilmenite');
    end
    
    % Suppression des tous les éléments de construction
    selectionTETTOT = [selectionTET1,selectionTET2,selectionTET3,selectionTET4,selectionTET5,selectionTET6,'dif1','dif2','blk1','dif3'];
    model.geom('geom1').create('del1', 'Delete');
    model.geom('geom1').feature('del1').selection('input').init;
    model.geom('geom1').feature('del1').selection('input').set(selectionTETTOT)
    
    model.geom('geom1').run('del1');
end

%% 2.1) Affichage de la géométrie
% Creation le 13/12/2016
for o=1
    for o=1:1 % Création d'une vue du modèle 3D (automatique)
        figure('color',[1 1 1])
        mphgeom(model);
        camzoom(1.8)
        xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
        ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
        zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
        xlim([0 110]);
        ylim([0 50]);
        zlim([-5 10]);
    end
end

%% 2.2) Affichage de la géométrie en vue de coupe
% Creation le 18/01/2017
for o=1
    for o=1:1 % Création d'une vue de coupe (automatique)
        model.geom('geom1').create('blk2', 'Block');
        model.geom('geom1').feature('blk2').set('pos', {'0' '26' '-5'});
        model.geom('geom1').feature('blk2').set('size', {'110' '1' '15'});
        model.geom('geom1').run('blk2');
        
        model.geom('geom1').create('uni1', 'Union');
        model.geom('geom1').feature('uni1').selection('input').set({'Anorthosite' 'Base' 'SableB' 'SableH' 'Anorthosite_Gross' 'Ilmenite'});
        model.geom('geom1').run('uni1');
        
        model.geom('geom1').create('Int1', 'Difference');
        model.geom('geom1').feature('Int1').selection('input').set('blk2');
        model.geom('geom1').feature('Int1').selection('input2').set('uni1');
        model.geom('geom1').feature('Int1').set('keep', 'on');
        model.geom('geom1').feature('Int1').label('Int1');
        model.geom('geom1').run('Int1');
        
        model.geom('geom1').create('Int2', 'Difference');
        model.geom('geom1').feature('Int2').selection('input').set('blk2');
        model.geom('geom1').feature('Int2').selection('input2').set('Int1');
        model.geom('geom1').feature('Int2').set('keep', 'off');
        model.geom('geom1').feature('Int2').label('Int2');
        model.geom('geom1').run('Int2');
        
        figure('color',[1 1 1])
        mphgeom(model);
        camzoom(1.8)
        xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
        ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
        zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
        xlim([0 110]);
        ylim([0 50]);
        zlim([-5 10]);
        
        model.geom('geom1').feature('blk2').active(false);
        model.geom('geom1').feature('uni1').active(false);
        model.geom('geom1').feature('Int1').active(false);
        model.geom('geom1').feature('Int2').active(false);
        model.geom('geom1').runPre('fin');
    end
end

%% 2.3) Insertion des électrodes au modele
% Creation le 13/12/2016
% Modification le 12/01/2017
for o=1
    for o=1:1       % Tri des electrodes => COORD_ELEC_SORT.mat
        % load COORD_ELEC.txt
        %
        % % Organisation des électrodes
        % Esort = sortrows(COORD_ELEC,3);
        % L1=Esort(1:96,:);
        % L1=sortrows(L1,1);
        % for i=1:24
        %     L1((i-1)*4+1:i*4,:)=sortrows(L1((i-1)*4+1:i*4,:),2);
        % end
        % L2=Esort(97:192,:);
        % L2=sortrows(L2,1);
        % for i=1:24
        %     L2((i-1)*4+1:i*4,:)=sortrows(L2((i-1)*4+1:i*4,:),2);
        % end
        % COORD_ELEC_SORT=[L1;L2];
    end
    load COORD_ELEC_SORT.mat
    
    % Création des electrodes
    for i=1:192
        Name = ['elec_',num2str(i)];
        model.geom('geom1').feature.create(Name,'Point');
        for j=1:3
            Coordstr = num2str(COORD_ELEC_SORT(i,j));
            model.geom('geom1').feature(Name).setIndex('p', Coordstr, j-1);
        end
        model.geom('geom1').run(Name);
    end
end

%% 2.4) Insertion des sondes MPS et GS3
% Creation le 26/01/2017
% Modification le 03/02/2017
for o=1
    GS3_MPS=1;
    point=1;
    
    for o=1:1
        if GS3_MPS==1
            load GS3.mat
            load MPS.mat
            
            % Création des sondes GS3
            for i=1:length(GS3)
                for j=1:length(GS3{1,i})
                    Name = ['GS3_niv',num2str(i),'_lysi',num2str(j)];
                    if point==1
                        model.geom('geom1').feature.create(Name,'Point');
                        for k=1:3
                            Coordstr = num2str(GS3{1,i}(j,k)-0.05*(k==3));
                            model.geom('geom1').feature(Name).setIndex('p', Coordstr, k-1);
                        end
                    else
                        model.geom('geom1').create(Name, 'Cylinder');
                        model.geom('geom1').feature(Name).set('r', '0.03');
                        model.geom('geom1').feature(Name).set('h', '0.1');
                        for k=1:3
                            Coordstr = num2str(GS3{1,i}(j,k)-0.05*(k==3));
                            model.geom('geom1').feature(Name).setIndex('pos', Coordstr, k-1);
                        end
                    end
                    model.geom('geom1').run(Name);
                end
            end
            % Création des sondes MPS
            for i=1:length(MPS)
                for j=1:length(MPS{1,i})
                    Name = ['MPS_niv',num2str(i),'_lysi',num2str(j)];
                    if point==1
                        model.geom('geom1').feature.create(Name,'Point');
                        for k=1:3
                            Coordstr = num2str(MPS{1,i}(j,k)-0.05*(k==3)-0.05*(k==3)*(j==3));
                            model.geom('geom1').feature(Name).setIndex('p', Coordstr, k-1);
                        end
                    else
                        model.geom('geom1').create(Name, 'Cylinder');
                        model.geom('geom1').feature(Name).set('r', '0.03');
                        model.geom('geom1').feature(Name).set('h', '0.1');
                        for k=1:3
                            Coordstr = num2str(MPS{1,i}(j,k)-0.05*(k==3)-0.05*(k==3)*(j==3));
                            model.geom('geom1').feature(Name).setIndex('pos', Coordstr, k-1);
                        end
                    end
                    model.geom('geom1').run(Name);
                end
            end
        end
    end
end

%% 2.5) Insertion des câbles DTS au modèle
% Creation le 03/02/2017
for o=1
    Cables_DTS=1;
    
    for o=1:1
        if Cables_DTS==1
            load('DTS.mat')
            for i=1:3
                % Création du câble DTS N°i
                Name = ['DTS_niv_',num2str(i)];
                model.geom('geom1').create(Name, 'Polygon');
                model.geom('geom1').feature(Name).set('x', DTS{i}(:,1));
                model.geom('geom1').feature(Name).set('y',  DTS{i}(:,2));
                model.geom('geom1').feature(Name).set('z',  DTS{i}(:,3));
                model.geom('geom1').run(Name);
            end
        end
    end
end

%% 2.5) Attribution des matériaux aux éléments
% Création le 12/01/2017
% Modification le 18/01/2017
% Modification le 26/01/2017
for o=1
    for o=1:1
        % Sable
        model.material.create('SableMat', 'Common', 'comp1');
        model.material('SableMat').label('SABLE');
        % prop physiques
        model.material('SableMat').propertyGroup('def').set('porosity', '0.3');
        model.material('SableMat').propertyGroup('def').set('resistivity', {'3000'});
        model.material('SableMat').propertyGroup('def').set('electricconductivity', {'1/3000'});
        model.material('SableMat').propertyGroup('def').set('relpermittivity', {'1'});
        model.material('SableMat').set('family', 'plastic');
        % selection
        model.material('SableMat').selection.set([3 5]);
        % couleur
        for o=1:1
            model.material('SableMat').set('family', 'custom');
            model.material('SableMat').set('specular', 'yellow');
            model.material('SableMat').set('customdiffuse', [1 0.8196078538894653 0.3607843220233917]);
            model.material('SableMat').set('customambient', [1 0.8196078538894653 0.3607843220233917]);
            model.material('SableMat').set('noise', 'off');
        end
        
        % Ilmenite
        model.material.create('IlmeniteMat', 'Common', 'comp1');
        model.material('IlmeniteMat').label('ILMENITE');
        model.material('IlmeniteMat').propertyGroup('def').set('porosity', '0.1');
        model.material('IlmeniteMat').propertyGroup('def').set('resistivity', {'600'});
        model.material('IlmeniteMat').propertyGroup('def').set('electricconductivity', {'1/600'});
        model.material('IlmeniteMat').propertyGroup('def').set('relpermittivity', {'1'});
        model.material('IlmeniteMat').set('family', 'plastic');
        % selection
        model.material('IlmeniteMat').selection.set([1 2]);
        % couleur
        for o=1:1
            model.material('IlmeniteMat').set('family', 'custom');
            model.material('IlmeniteMat').set('customspecular', [0.4117647111415863 0.4117647111415863 0.4117647111415863]);
            model.material('IlmeniteMat').set('customdiffuse', [0.4117647111415863 0.4117647111415863 0.4117647111415863]);
            model.material('IlmeniteMat').set('customambient', [0.4117647111415863 0.4117647111415863 0.4117647111415863]);
            model.material('IlmeniteMat').set('noise', 'off');
        end
        
        % Anorthosite
        model.material.create('AnorthositeMat', 'Common', 'comp1');
        model.material('AnorthositeMat').label('ANORTHOSITE');
        model.material('AnorthositeMat').propertyGroup('def').set('porosity', '0.1');
        model.material('AnorthositeMat').propertyGroup('def').set('resistivity', {'8000'});
        model.material('AnorthositeMat').propertyGroup('def').set('electricconductivity', {'1/8000'});
        model.material('AnorthositeMat').propertyGroup('def').set('relpermittivity', {'1'});
        model.material('AnorthositeMat').set('family', 'plastic');
        % selection
        model.material('AnorthositeMat').selection.set([4 6]);
        % couleur
        for o=1:1
            model.material('AnorthositeMat').set('family', 'custom');
            model.material('AnorthositeMat').set('noise', 'off');
            model.material('AnorthositeMat').set('customspecular', [0.7529411911964417 0.7529411911964417 0.7529411911964417]);
            model.material('AnorthositeMat').set('customdiffuse', [0.7529411911964417 0.7529411911964417 0.7529411911964417]);
            model.material('AnorthositeMat').set('customambient', [0.7529411911964417 0.7529411911964417 0.7529411911964417]);
            
        end
        
        if GS3_MPS==1 & point==0
            % Sondes GS3
            model.material.create('GS3Mat', 'Common', 'comp1');
            model.material('GS3Mat').label('GS3');
            model.material('GS3Mat').set('family', 'custom');
            model.material('GS3Mat').set('specular', 'red');
            model.material('GS3Mat').set('lighting', 'simple');
            model.material('GS3Mat').set('customdiffuse', [0.6980392336845398 0.13333334028720856 0.13333334028720856]);
            model.material('GS3Mat').set('customambient', [0.6980392336845398 0.13333334028720856 0.13333334028720856]);
            model.material('GS3Mat').set('noise', 'off');
            
            model.material('GS3Mat').propertyGroup('def').set('porosity', '0.1');
            model.material('GS3Mat').propertyGroup('def').set('resistivity', {'600'});
            model.material('GS3Mat').propertyGroup('def').set('electricconductivity', {'1/600'});
            model.material('GS3Mat').propertyGroup('def').set('relpermittivity', {'1'});
            
            % Sondes MPS
            model.material.create('MPSMat', 'Common', 'comp1');
            model.material('MPSMat').label('MPS');
            model.material('MPSMat').set('family', 'custom');
            model.material('MPSMat').set('specular', 'green');
            model.material('MPSMat').set('lighting', 'simple');
            model.material('MPSMat').set('customdiffuse', [0.03529411926865578 0.4627451002597809 0.03529411926865578]);
            model.material('MPSMat').set('customambient', [0.03529411926865578 0.4627451002597809 0.03529411926865578]);
            model.material('MPSMat').set('noise', 'off');
            
            model.material('MPSMat').propertyGroup('def').set('porosity', '0.1');
            model.material('MPSMat').propertyGroup('def').set('resistivity', {'600'});
            model.material('MPSMat').propertyGroup('def').set('electricconductivity', {'1/600'});
            model.material('MPSMat').propertyGroup('def').set('relpermittivity', {'1'});
        end
    end
end

%% 3) Maillage de la géométrie
% Creation le 13/12/2016
% Modification le 12/01/2017
for o=1
    automatique=0;
    
    for o=1:1 % Maillage avec les paramètres souhaités et affichage sur matlab
        if automatique ==1
            Finesse = menu('Quelle finesse du maillage COMSOL ?','extremement fin (51s / itération)','extra fin (25s / itération)','plus fin (14s / itération)','fin (11 / itération)','normal (10s / itération)','grossier (8s / itération)');
            FinesseChoix = [1,2,3,4,5,6];
            Finesse=FinesseChoix(1,Finesse);
            model.mesh('mesh1').autoMeshSize(Finesse);
            model.mesh('mesh1').run;
        else
            model.mesh('mesh1').automatic(false);
            model.mesh('mesh1').feature('size').set('table', 'default');
            model.mesh('mesh1').feature('size').set('custom', 'on');
            model.mesh('mesh1').feature('size').set('hmax', '20');
            model.mesh('mesh1').feature('size').set('hmin', '0.3');
            model.mesh('mesh1').feature('size').set('hgrad', '2.5');
            model.mesh('mesh1').feature('size').set('hcurve', '1');
            model.mesh('mesh1').feature('size').set('hnarrow', '2');
            model.mesh('mesh1').run;
        end
        figure('color',[1 1 1])
        mphmesh(model);
        camzoom(1.8)
        xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
        ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
        zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
        xlim([0 110]);
        ylim([0 50]);
        zlim([-5 10]);
    end
end

%% 4) Chargement des configurations utilisées
% Creation le 13/12/2016
for o=1
    for o=1:1       % Tri des configurations utilisées => CONFIG_TOT.mat
        %     load DATA2.mat
        %     CONFIG_TOT = DATA2{1}(:,1:4);
        %     CONFIG_TOT(:,5) = DATA2{1}(:,12);           % Facteur géometrique K
        %     for i=2:length(DATA2)
        %         CONFIG_TOTlocal = DATA2{i}(:,1:4);
        %         if length(DATA2{i}(1,:))==10
        %             Z=zeros(length(DATA2{i}(:,1)));
        %             CONFIG_TOTlocal(:,5)=Z(:,1);
        %         else
        %             CONFIG_TOTlocal(:,5) = DATA2{i}(:,12);
        %         end
        %         CONFIG_TOT=[CONFIG_TOT;CONFIG_TOTlocal];
        %         clear CONFIG_TOTlocal
        %     end
        %     CONFIG_TOT=unique(CONFIG_TOT,'rows');
    end
    load CONFIG_TOT.mat
    INJECTION=unique(CONFIG_TOT(:,1:2),'rows');
end

%% 5) Calcul du facteur géométrique réel
% Creation le 13/12/2016
% Modification le 12/01/2017
for o=1
    % Ajout de la physique ACDC a COMSOL
    model.physics.create('ec', 'ConductiveMedia', 'geom1');
    
    % Gestion des conditions aux frontières
    model.physics('ec').feature.create('pot1', 'ElectricPotential', 2);
    model.physics('ec').feature('pot1').selection.set([1 2 3 22 450]);
    
    % Création de l'étude
    model.study.create('std1');
    model.study('std1').create('stat', 'Stationary');
    model.study('std1').feature('stat').activate('ec', true);
    model.study('std1').feature('stat').set('showdistribute', false);
    model.study('std1').feature('stat').set('notlistsolnum', 1);
    model.study('std1').feature('stat').set('notsolnum', '1');
    model.study('std1').feature('stat').set('listsolnum', 1);
    model.study('std1').feature('stat').set('solnum', '1');
    
    I=1;     % courant injecté en A
    
    for i=1:1   %:length(COORD_ELEC_SORT(:,1))
        tic
        C1=i;
        
        % Injection de +I a C1
        if i==1
            model.geom('geom1').create('selC1', 'ExplicitSelection');
        end
        model.geom('geom1').feature('selC1').selection('selection').init(0);
        model.geom('geom1').feature('selC1').selection('selection').init;
        SelC1={['elec_',num2str(C1)]};
        model.geom('geom1').feature('selC1').selection('selection').set(SelC1);
        if i==1
            model.physics('ec').feature.create('elec_C1', 'PointCurrentSource', 0);
        end
        model.physics('ec').feature('elec_C1').selection.named('geom1_selC1_pnt');
        model.physics('ec').feature('elec_C1').set('Qjp', num2str(I));
        
        % Création de la solution
        if i==1
            model.sol.create('sol1');
            model.sol('sol1').study('std1');
            model.sol('sol1').create('st1', 'StudyStep');
            model.sol('sol1').feature('st1').set('study', 'std1');
            model.sol('sol1').feature('st1').set('studystep', 'stat');
            model.sol('sol1').create('v1', 'Variables');
            model.sol('sol1').feature('v1').set('control', 'stat');
            model.sol('sol1').create('s1', 'Stationary');
            model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
            model.sol('sol1').feature('s1').create('i1', 'Iterative');
            model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'cg');
            model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'amg');
            model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'i1');
            model.sol('sol1').feature('s1').feature.remove('fcDef');
            model.sol('sol1').attach('std1');
            model.sol('sol1').runFromTo('st1', 's1');
        end
        model.sol('sol1').runAll;
        
        % Affichage des résultats avec COMSOL
        if i==1
            model.result.create('pg1', 'PlotGroup3D');
            model.result('pg1').feature.create('mslc1', 'Multislice');
        end
        model.result('pg1').feature('mslc1').set('data', 'dset1');
        model.result('pg1').feature('mslc1').set('rangecoloractive', 'on');
        model.result('pg1').feature('mslc1').set('rangecolormin', '0');
        model.result('pg1').feature('mslc1').set('rangecolormax', '5');
        model.result('pg1').run;
        model.view('view1').set('transparency', 'off');
        model.result('pg1').feature('mslc1').set('ynumber', '0');
        model.result('pg1').feature('mslc1').set('zcoord', '1');
        model.result('pg1').run;
        model.result('pg1').feature('mslc1').set('xnumber', '15');
        model.result('pg1').run;
        
        % Affichage des résultats avec MATLAB
        figure('color',[1 1 1])
        mphplot(model,'pg1','rangenum',1)
        camzoom(1.8)
        xlabel('X (m)', 'FontSize', 12,'FontWeight','bold')
        ylabel('Y (m)', 'FontSize', 12,'FontWeight','bold')
        zlabel('Z (m)', 'FontSize', 12,'FontWeight','bold')
        xlim([0 110]);
        ylim([0 50]);
        zlim([-5 10]);
        
        % Mesure des potentiels aux autres electrodes
        for j=1:length(COORD_ELEC_SORT(:,1))
            if i==j
                POTENTIEL(j,i)=NaN;
            else
                POTENTIEL(j,i)=mphinterp(model,{'V'},'coord',[COORD_ELEC_SORT(j,1)' COORD_ELEC_SORT(j,2)' COORD_ELEC_SORT(j,3)']');
            end
        end
        i
        toc
    end
    
    % Calcul des facteurs géométriques
    load POTENTIEL.mat
    for i=1:length(CONFIG_TOT(:,1))
        C1=CONFIG_TOT(i,1);
        C2=CONFIG_TOT(i,2);
        P1=CONFIG_TOT(i,3);
        P2=CONFIG_TOT(i,4);
        
        C1P1=POTENTIEL(P1,C1);
        C1P2=POTENTIEL(P2,C1);
        C2P1=POTENTIEL(P1,C2);
        C2P2=POTENTIEL(P2,C2);
        DV=(C1P2-C1P1)-(C2P2-C2P1);
        K=-rho0/DV*I;
        CONFIG_TOT(i,6)=K;
        CONFIG_TOT(i,7)=abs((CONFIG_TOT(i,6)-CONFIG_TOT(i,5))/CONFIG_TOT(i,5))*100;
    end
    
    figure('color',[1 1 1])
    plot(CONFIG_TOT(1:100,5),CONFIG_TOT(1:100,6),'k*')
    xlim([0 10000])
    ylim([0 10000])
    
    figure('color',[1 1 1])
    plot(CONFIG_TOT(:,5),CONFIG_TOT(:,6),'k*')
    xlim([0 10000])
    ylim([0 10000])
    
end
