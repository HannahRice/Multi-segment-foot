% Hannah Rice August 2025

% Function that uses plantar pressure data to obtain the relative force acting under three segments of the foot 
% during running stance and the location of the centre of pressure of each segment in a QTM (Qualisys) coordinate system

% The function effectively treats the pressure plate as three force plates that represent the
% rearfoot, midfoot and forefoot for each trial, based on the relative location of
% each participant's foot markers during a static trial.

% Function is specific to a (Tekscan Mat 7101E) pressure plate with
% the following dimensions: 
% Matrix Height: 447.0 mm
% Matrix Width: 487.7 mm
% Thickness: 0.102 mm



%% File TYPE INFORMATION 
% Input pressure files are .csv files from Tekscan software and represent a matrix of force values for each pressure sensor where
% each row the location along the anterior-posterior plate axis in the direction of running; each column
% represents the location along the medial-lateral plate axis in the direction of running. Each data frame is represented in a new matrix
% separated by empty rows 

% Input static file provides the marker locations from the static trial
% from OpenSim, and presented in the QTM Coordinate system, saved as a .m file. 
% Marker names are called in 80 - 86 of the function. 


%% OUTPUT is a structure which contains for each of the three foot segments an n by 3 matrix, where n represents the number of frames 
% during stance and the columns are the force % under that foot region, the
% ML and the AP coordinate of the centre of pressure of that segment in the
% QTM coordinate system


function [ID_input] = MultSeg_Force_CoP(ST, DT) 
    % ST = static trial
    % DT = dynamic trial
    OrData = xlsread(DT); %Original Data
    ss=load(ST); %load Static File


    % change order of xyz coordinates from OpenSim to Qualisys 
    % clearvars -except OrData STATIC_MARKER
    FootSegment = {'RF','MF','FF'};
    
    %% Pressure plate dimensions
    PP_ML_width = 447;
    PP_AP_length = 487.7; %These dimensions are from the manufacturer information
    PP_Origin = mean(ss.STATIC_MARKER.PPTR); %marker representing the bottom left corner of the pressure plate sensors in the directon of running
    % Obtain plate surface height above the ground
    Plate_Mean = mean([mean(ss.STATIC_MARKER.PPBL);mean(ss.STATIC_MARKER.PPBR);mean(ss.STATIC_MARKER.PPTL);mean(ss.STATIC_MARKER.PPTR)]);
    Av_Z_coord = Plate_Mean(3); 
    
    
    %% Find frames separated by NaNs in file
    % Each data frame in the pressure .csv file is represented in a new matrix
    % separated by empty rows which are read as NaNs in OrData
    
    [a b] = find(isnan(OrData(:,1)));
    
        for i = 1:length(a)-1
            FF(i,1) = a(i+1)-a(i);
        end
        [c d] = find(FF>1);
        for j=1:length(c)
            FrameStartRow(j,1) = a(c(j))+1;
            FrameEndRow(j,1) = a(c(j)+1)-1;
        end
    
    s=size(OrData(FrameStartRow(1):FrameEndRow(1),:));


    %% Define polygons of the three foot zones from marker positions during static trial
    % 88 x 96 sensors in the medial-lateral and anterior-posterior directions respectively 
    % Calculate conversion factors for the plate size and number of sensors
    CF_ML = PP_ML_width/s(1); 
    CF_AP = PP_AP_length/s(2); 
    
    % Identify relevant markers used to create the polygons
    DistCalc = mean(ss.STATIC_MARKER.RDC);
    Navic = mean(ss.STATIC_MARKER.RNAV);
    M5base = mean(ss.STATIC_MARKER.RMB5); 
    M1head = mean(ss.STATIC_MARKER.RMH1);
    M5head = mean(ss.STATIC_MARKER.RMH5);
    Toe2 = mean(ss.STATIC_MARKER.RD2);
    Hallux = mean(ss.STATIC_MARKER.RD1);
    
    
    %% Find longest length of the whole foot, choosing between either the length from Distal Calc to Hallux or Distal Calc to second toe, depending on which is longest
    Len_to_T2 = sqrt((Toe2(1)-DistCalc(1))^2 + (Toe2(2) - DistCalc(2))^2);
    Len_to_Hx = sqrt((Hallux(1)-DistCalc(1))^2 + (Hallux(2) - DistCalc(2))^2);
    
    
    if Len_to_Hx > Len_to_T2
        FootLengthLine = Len_to_Hx;
        FootEnd = Hallux;
    else
        FootLengthLine = Len_to_T2;
        FootEnd = Toe2;
    end
    
    
    %% Calculate lengths of each polygon converted into the lab length units (i.e. true)
    Rf_left_edge_pp = (abs(DistCalc(2)-Navic(2)))/CF_AP;
    Rf_right_edge_pp = (abs(DistCalc(2)-M5base(2)))/CF_AP;
    Mf_left_edge_pp = (abs(M1head(2)-Navic(2)))/CF_AP;
    Mf_right_edge_pp = (abs(M5head(2)-M5base(2)))/CF_AP;
    
    %% Find the sum of the force under the whole foot for each frame in order to determine 
    % touchdown and takeoff

    for f = 1:length(FrameStartRow)
        PlateForce(f,1:s(1),1:s(2)) = OrData(FrameStartRow(f):FrameEndRow(f),:);
        SumPF(f) = sum(sum(PlateForce(f,:,:)));
    end
    TD  = find(SumPF>0,1,'first');
    TO = find(SumPF>0,1,'last');
    
    %% Consider by position rather than time to segment the foot
    % Find sum of forces per sensor to segment the foot   
    for t=1:s(1)
        SumSensors(t,:) = sum(PlateForce(:,t,:));
    end
    
    %% Find the sensors that mark the outer ML and AP edges of the foot contact based on force detection 

    [~, APpos_prox] = find(SumSensors(:,:)>10,1,'first');
    [~, APpos_dist] = find(SumSensors(:,:)>10,1,'last'); APpos_dist = APpos_dist+1;
    
    [~, MLpos_Left] = find(SumSensors(:,:)'>10,1,'first'); 
    if MLpos_Left>1 
        MLpos_Left-1;
    end
    [~, MLpos_Right] = find(SumSensors(:,:)'>10,1,'last');
    
    
          
    
    %% Create polygons that represent each region
    % Segment start and end points, in the AP direction, relative to the
    % pressure plate
    % p,d,m,l refer to proximal, distal, medial, lateral
    
    
    RF_p_m = [APpos_prox,MLpos_Left];
    RF_p_l = [APpos_prox,MLpos_Right];
    
    RF_d_m = [APpos_prox + round(Rf_left_edge_pp),MLpos_Left];
    RF_d_l = [APpos_prox + round(Rf_right_edge_pp),MLpos_Right];
    
    MF_d_m = [RF_d_m(1) + round(Mf_left_edge_pp),MLpos_Left];
    MF_d_l = [RF_d_m(1) + round(Mf_right_edge_pp),MLpos_Right];
    
    FF_d_m = [APpos_dist,MLpos_Left];
    FF_d_l = [APpos_dist,MLpos_Right];
    
    % To ensure the polygons do not exceed the pressure plate boundaries:
    if FF_d_m(1) > s(1)
        FF_d_m(1) = APpos_dist -1;
    end
    
   
   
    %% Figure showing division of foot as a check
    figure; imagesc(SumSensors)
    %vertical lines
    line([RF_p_m(1),RF_p_l(1)],[RF_p_m(2),RF_p_l(2)],'Color','w')
    line([RF_d_m(1),RF_d_l(1)],[RF_d_m(2),RF_d_l(2)],'Color','w')
    line([MF_d_m(1),MF_d_l(1)],[MF_d_m(2),MF_d_l(2)],'Color','w')
    line([FF_d_m(1),FF_d_l(1)],[FF_d_m(2),FF_d_l(2)],'Color','w')
    %horizontal lines
    line([RF_p_m(1),FF_d_m(1)],[RF_d_m(2),RF_d_m(2)],'Color','w')
    line([RF_p_l(1),FF_d_l(1)],[RF_d_l(2),RF_d_l(2)],'Color','w')

    
    % for each zone line divider find the angle relative to the vertical 
    RF_d_theta = atand(((RF_d_l(2)-RF_d_m(2))/((RF_d_l(1)-RF_d_m(1)))));
    MF_d_theta = atand(((MF_d_l(2)-MF_d_m(2))/((MF_d_l(1)-MF_d_m(1)))));
    if isinf(RF_d_theta)
        RF_d_theta = 90;
    end
    if isinf(MF_d_theta)
        MF_d_theta = 90;
    end

    % Now for each row in that polygon, find the x coordinate of the line
    % X-coordinate must be rounded because it will read the force value
    % from that cell
    for ii = TD:TO
        ForceData = OrData(FrameStartRow(ii):FrameEndRow(ii),:);
        %Create a new matrix that reads in the force values, starting with a
        %rectangular filled with NaNs so that CoP can be calculated from a
        %rectangular shape - effectively overlaid with the polygon that includes
        %the values 
        
        %% NaNs the size of the whole foot box
    
        Zone.RF(:,:,ii) = nan(RF_p_l(2)-RF_p_m(2),FF_d_l(1)-RF_p_l(1));
        Zone.MF(:,:,ii) = nan(RF_p_l(2)-RF_p_m(2),FF_d_l(1)-RF_p_l(1));
        Zone.FF(:,:,ii) = nan(RF_p_l(2)-RF_p_m(2),FF_d_l(1)-RF_p_l(1));
    
    
        %% Calculate Force in each zone to later obtain the necessary relative percentages 
        % Summing horizontally
         for h = RF_p_m(2):RF_p_l(2) % each sensor going down in ML direction 
            x_RF_d(h-RF_p_m(2)+1) = round(((h-RF_p_m(2))/tand(RF_d_theta)) + RF_d_m(1));
            x_MF_d(h-RF_p_m(2)+1) = round(((h-RF_p_m(2))/tand(MF_d_theta)) + MF_d_m(1));
         end
         
         
         for v = 1:length(RF_p_m(2):RF_p_l(2))-1      
             Start_row_MF = x_RF_d(v)-min(RF_d_m(1),RF_d_l(1));
             Start_row_FF = x_MF_d(v)-min(MF_d_m(1),MF_d_l(1));
             Zone.RF(v,1:length(RF_p_m(1):x_RF_d(v)),ii) = ForceData(RF_p_m(2)+v-1,RF_p_m(1):x_RF_d(v));
             Zone.MF(v,Start_row_MF+1:Start_row_MF+length(x_RF_d(v):x_MF_d(v)),ii) = ForceData(RF_p_m(2)+v-1,x_RF_d(v):x_MF_d(v));
             Zone.FF(v,Start_row_FF+1:Start_row_FF+length(x_MF_d(v):FF_d_m(1)),ii) = ForceData(RF_p_m(2)+v-1,x_MF_d(v):FF_d_m(1));
         end
    
         for fs=1:length(FootSegment)
    
            Force.(FootSegment{fs})(ii) = sum(sum(Zone.(FootSegment{fs})(:,:,ii),'omitnan'),'omitnan');
         end
        TotalForce(ii) = sum([Force.RF(ii),Force.MF(ii),Force.FF(ii)]);
        
    end
    
    %% Obtain the force under each foot segment as a % of the total foot force
    TotalForce(1:TD-1) = [];
    for fs = 1: length(FootSegment) 
        Force.(FootSegment{fs})(1:TD-1) = [];
        Force_pc.(FootSegment{fs}) = 100*(Force.(FootSegment{fs}))./TotalForce; 
    end
    

    % % 
    figure; hold on
    plot(Force.RF); plot(Force.MF), plot(Force.FF); plot(TotalForce,'k')
    legend('RF','MF','FF','TotalForce')

    jj=TD+5; %5 frames after touchdown, just as a check 
    figure; heatmap(OrData(FrameStartRow(jj):FrameEndRow(jj),:));
    figure; heatmap(Zone.RF(:,:,jj));
    figure; heatmap(Zone.MF(:,:,jj));
    figure; heatmap(Zone.FF(:,:,jj));
    
    

    %% For each timeframe, and for each foot zone (RF,MF,FF) find the sum of the force in each row and column 
    s_array = size(Zone.RF);
    for ii = TD:TO
        for fs=1:length(FootSegment)
        
           
            ForceSum_horiz.(FootSegment{fs}) = sum(Zone.(FootSegment{fs})(:,:,ii)','omitnan')';
            ForceSum_vert.(FootSegment{fs}) = sum(Zone.(FootSegment{fs})(:,:,ii),'omitnan')';
    
            clearvars h v
    
            for h = 1:s_array(1)
                Weighting_horiz.(FootSegment{fs})(h) = ForceSum_horiz.(FootSegment{fs})(h)*(h-1);
            end
            if sum(ForceSum_horiz.(FootSegment{fs}))>0
                X_coord_labCS.(FootSegment{fs})(ii,1) = sum(Weighting_horiz.(FootSegment{fs})')/sum(ForceSum_horiz.(FootSegment{fs}));
            else
                X_coord_labCS.(FootSegment{fs})(ii,1) = NaN;
            end
    
            for v = 1:FF_d_m(1)-RF_p_m(1)
                Weighting_vert.(FootSegment{fs})(v) = ForceSum_vert.(FootSegment{fs})(v)*(v-1);
            end
            if sum(ForceSum_vert.(FootSegment{fs}))>0
                Y_coord_labCS.(FootSegment{fs})(ii,1) = sum(Weighting_vert.(FootSegment{fs})')/sum(ForceSum_vert.(FootSegment{fs}));
            else
                Y_coord_labCS.(FootSegment{fs})(ii,1) = NaN;
            end
            %Remove distorted values due to only one pressure cell registering
            %pressure 
            [a b] = find(Weighting_horiz.(FootSegment{fs})'>0);
            [c d] = find(Weighting_vert.(FootSegment{fs})'>0);
            if length(a) && length(c) < 2
                X_coord_labCS.(FootSegment{fs})(ii,1)=X_coord_labCS.(FootSegment{fs})(ii-1,1);
                Y_coord_labCS.(FootSegment{fs})(ii,1)=Y_coord_labCS.(FootSegment{fs})(ii-1,1);
            end
    
            clearvars a b c d Weighting_vert Weighting_horiz
            Z_coord.(FootSegment{fs})(ii,1) = Av_Z_coord;
        end
    end
    
    for fs=1:length(FootSegment)
        X_coord_labCS.(FootSegment{fs})(1:TD-1) = [];
        Y_coord_labCS.(FootSegment{fs})(1:TD-1) = [];
        Z_coord.(FootSegment{fs})(1:TD-1) = [];
    end
    
    % Need to add the coordinates of X and Y for each segment to the
    % "origin of that region (i.e. its proximal medial line coordinates)
    
        X_coord_labCS.RF = X_coord_labCS.RF + RF_p_m(2);
        Y_coord_labCS.RF = Y_coord_labCS.RF + RF_p_m(1);
        X_coord_labCS.MF = X_coord_labCS.MF + RF_d_m(2);
        Y_coord_labCS.MF = Y_coord_labCS.MF + min(RF_d_m(1),RF_d_l(1));
        X_coord_labCS.FF = X_coord_labCS.FF + MF_d_m(2);
        Y_coord_labCS.FF = Y_coord_labCS.FF + min(MF_d_m(1),MF_d_l(1));
        % % 
    for fs=1:length(FootSegment)
        CoP.(FootSegment{fs})=[X_coord_labCS.(FootSegment{fs}),Y_coord_labCS.(FootSegment{fs}),Z_coord.(FootSegment{fs})];
    end
    
    
    
    
    %% Translate these into Qualisys/lab coordinate system
    
    figure; hold on
    for fs=1:length(FootSegment)
        for m = 1:length(CoP.(FootSegment{fs}))
            if CoP.(FootSegment{fs})(m,1) > 0
                CoP_final.(FootSegment{fs})(m,1) = PP_Origin(1)+ (-CoP.(FootSegment{fs})(m,1)*CF_ML);
                CoP_final.(FootSegment{fs})(m,2) = PP_Origin(2)+ (-CoP.(FootSegment{fs})(m,2)*CF_AP);
    
            else
                 CoP_final.(FootSegment{fs})(m,:) = NaN;  
            end
        end
        scatter(CoP_final.(FootSegment{fs})(:,2),CoP_final.(FootSegment{fs})(:,1))
        set (gca,'ydir','reverse')
        axis equal
        legend('RF','MF','FF')
    end
    
    %% Replace NaN with zero for use in OpenSim
    for fs=1:length(FootSegment)
        for c = 1
            [a b]=find(isnan(CoP_final.(FootSegment{fs})(:,c)));
            CoP_final.(FootSegment{fs})(a,:)=0;
            clearvars a b
            ID_input.(FootSegment{fs}) =[Force_pc.(FootSegment{fs})',CoP_final.(FootSegment{fs})];
        end 
    end
    
    
    % save('ScaledForce_S24_PostRun_0003_S95.mat','ID_input')
