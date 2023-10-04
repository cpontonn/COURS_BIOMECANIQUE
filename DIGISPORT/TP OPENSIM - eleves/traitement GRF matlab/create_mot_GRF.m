function []=create_mot_GRF(filename,q1,q2,q3,f_cutoff,fout)
% filename le fichier C3D a ouvrir sans l'exntesion 'myfile'
%f cutoff fréquence passe bas
%fout nom du fichier de sortie
% q1 q2 q3 les rotations sur x y z entre le repère mocap et opensim en
% degrés
%% calcul du CoP dans les repères plateformes
h = btkReadAcquisition([filename '.c3d']);
[Forceplates, ForceplatesInfo] = btkGetForcePlatforms(h);
f_mocap=btkGetPointFrequency(h);
for i=1:numel(Forceplates)
    for j=1:6
        FieldsName = fieldnames(Forceplates(i).channels);
        Data(i).FullData(:,j)=resample(eval(['Forceplates(i).channels.' FieldsName{j}]),f_mocap,ForceplatesInfo(i).frequency);
        Data(i).RawData(:,j) = Data(i).FullData(:,j);
    end
end
% changing N.mm to N.m
for i=1:numel(Data)
    Data(i).RawData(:,4:6) = Data(i).RawData(:,4:6)/1000;
end

%filtering

    for i=1:numel(Data)
        Data(i).Data = filt_data(Data(i).RawData,f_cutoff,f_mocap);
    end


%% passage des Fx Fy Fz et CoP dans le repère capture

%% tranformation dans le repère opensim
Ropcap=transfo_cap_to_opensim(q1,q2,q3) % 

for i=1:numel(Forceplates)
        Origin = mean(Forceplates(i).corners,2)/1000;
        x = ((Forceplates(i).corners(:,4) + Forceplates(i).corners(:,1))/2) - ((Forceplates(i).corners(:,2) + Forceplates(i).corners(:,3))/2); x=x/norm(x);
        y = ((Forceplates(i).corners(:,1) + Forceplates(i).corners(:,2))/2) - ((Forceplates(i).corners(:,3) + Forceplates(i).corners(:,4))/2); y=y/norm(y);
        z = cross(x,y); z = z/norm(z);
        y = cross(z,x);
        Rplatform = [x y z]'; 
for j=1:numel(Data(i).Data(:,1))
    F=Rplatform*Data(i).Data(j,1:3)';% effort dans la base mocap
     Mp=Rplatform*Data(i).Data(j,4:6)';%effort à l'origine du repère plateforme dans la base mocap
     M0=Mp + cross(Origin,F);%moment à l'origine du repère mocap, dans la base mocap

    CoP = cross(Data(i).Data(j,1:3)',Data(i).Data(j,4:6)')/(norm(Data(i).Data(j,1:3)')^2);%dans le repère plateforme
    CoP = CoP - (CoP(3)/Data(i).Data(j,3))*Data(i).Data(j,1:3)'; % point sur l'axe avec z=0;
    CoP0=Origin+Rplatform*CoP;
    F=Ropcap*F;
    M0=Ropcap*M0;
    CoP0=Ropcap*CoP0;
    Data(i).FCoPM(j,:) = [F' CoP0' M0'];   % enrichissement de la variable "external_forces"
end
end
%% écriture du fichier .mot

% Order:  rGRF(xyz), rCOP(xyz), lGRF(xyz), lCOP(xyz), rT(xyz), lT(xyz)
label{1} = 'ground_force_vx';
label{2} = 'ground_force_vy';
label{3} = 'ground_force_vz';
label{4} = 'ground_force_px';
label{5} = 'ground_force_py';
label{6} = 'ground_force_pz';
label{7} = 'ground_force_vx';
label{8} = 'ground_force_vy';
label{9} = 'ground_force_vz';
label{10} = 'ground_force_px';
label{11} = 'ground_force_py';
label{12} = 'ground_force_pz';
label{13} = 'ground_torque_x';
label{14} = 'ground_torque_y';
label{15} = 'ground_torque_z';
label{16} = 'ground_torque_x';
label{17} = 'ground_torque_y';
label{18} = 'ground_torque_z';
forceIndex = length(label);

Data(i).FCoPM(j,:)    
% Initialize 'motion file data matrix' for writing data of interest.
nRows = numel(Data(1).FCoPM(:,1));  
nCols = length(label)+1;   % plus time
motData = zeros(nRows, nCols);

% Write time array to data matrix.
time = [0:1/f_mocap:((nRows-1)/f_mocap)]'; 
motData(:, 1) = time;

% Write force data to data matrix.
% NOTE:  each field of mCS.forces has xyz components.
forceData = [Data(1).FCoPM(:,1) Data(1).FCoPM(:,2) Data(1).FCoPM(:,3) ...
             Data(1).FCoPM(:,4) Data(1).FCoPM(:,5) Data(1).FCoPM(:,6) ... 
             Data(2).FCoPM(:,1) Data(2).FCoPM(:,2) Data(2).FCoPM(:,3) ...
             Data(2).FCoPM(:,4) Data(2).FCoPM(:,5) Data(2).FCoPM(:,6) ...
             Data(1).FCoPM(:,7) Data(1).FCoPM(:,8) Data(1).FCoPM(:,9) ...
             Data(2).FCoPM(:,7) Data(2).FCoPM(:,8) Data(2).FCoPM(:,9)];

motData(:, 2:end) = forceData;          

% Open file for writing.
fid = fopen(fout, 'w');
if fid == -1
    error(['unable to open ', fout])
end

% Write header.
fprintf(fid, 'name %s\n', fout);
fprintf(fid, 'datacolumns %d\n', nCols);
fprintf(fid, 'datarows %d\n', nRows);
fprintf(fid, 'range %d %d\n', time(1), time(nRows));
fprintf(fid, 'endheader\n\n');

% Write column labels.
fprintf(fid, '%20s\t', 'time');
for i = 1:nCols-1
	fprintf(fid, '%20s\t', label{i});
end

% Write data.
for i = 1:nRows
    fprintf(fid, '\n'); 
	for j = 1:nCols
        fprintf(fid, '%20.8f\t', motData(i, j));
    end
end

fclose(fid);

end