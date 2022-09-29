%% dynamo2relion_particlesInMicrographs

% see comments in dynamo2relion_particles.m for running instructions



%% inputs
inputTable = ''; %Dynamo table
voltage = <>;
sphericalAberration = <>;
amplitudeContrast = <>;
detectorPixelSize = <>;
particleBaseName = ''; %the path to the particles in Relion, starting from the project folder 
particleExtension = '000001.mrc';
docFileName = ''; %the Dynamo doc file from the crop folder
calibratedPixelSize = <>;

headerFile = 'relionHeader.txt';
starFileName = ''; % output star file name
micrographStarFileName = 'micrographs.star';


%% get Particles

%shift particles into the center of their boxes;
table=dread(inputTable);
% table(:,24)=table(:,4)+table(:,24);
% table(:,25)=table(:,5)+table(:,25);
% table(:,26)=table(:,6)+table(:,26);
% table(:,4:6)=0;

%calculate the Eulers in Relion format (ZYZ)
tdrot=table(:,7);
tilt=table(:,8);
narot=table(:,9);
rot_xmipp=-tdrot+90;
tilt_xmipp=+tilt;
psi_xmipp=-narot+270;
rot_relion=-psi_xmipp;
tilt_relion=-tilt_xmipp;
psi_relion=-rot_xmipp;

%calculate actual magnification (NOT the nominal)
magnification=(detectorPixelSize*10000/calibratedPixelSize);

%read in the .doc file
delimiter = ' ';
formatSpec = '%f%s%[^\n\r]';
fileID = fopen(docFileName,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
dataArray(1) = cellfun(@(x) num2cell(x), dataArray(1), 'UniformOutput', false);
dataArray(2) = cellfun(@(x) mat2cell(x, ones(length(x), 1)), dataArray(2), 'UniformOutput', false);
bin1tomograms = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

tomoIndices=cell2mat(bin1tomograms(:,1));


%% Assemble the Star file

for i=1:length(table)
    tomogramIndex = table(i,20);
    tomoStringIndex = find(tomoIndices(:)==tomogramIndex);
    tomogramName = bin1tomograms{tomoStringIndex,2};
    starCellArray{i,1} = tomogramName;
    starCellArray{i,2} = voltage;
    starCellArray{i,3} = sphericalAberration;
    starCellArray{i,4} = amplitudeContrast;
    starCellArray{i,5} = magnification;
    starCellArray{i,6} = detectorPixelSize;
    starCellArray{i,7} = ([particleBaseName,sprintf('%06d',table(i,1)),particleExtension]);
    starCellArray{i,8} = round(table(i,24));
    starCellArray{i,9} = round(table(i,25));
    starCellArray{i,10} = round(table(i,26));
    starCellArray{i,11} = rot_relion(i);
    starCellArray{i,12} = tilt_relion(i);
    starCellArray{i,13} = psi_relion(i);
    starCellArray{i,14} = 0;
    starCellArray{i,15} = 0;
    starCellArray{i,16} = 0;
end


%% write StarFile

C = starCellArray.';
fid = fopen('temp.txt', 'wt');
fprintf(fid, '%s %5f %5f %6f %g %6f %s %5f %5f %5f %3f %3f %3f %3f %3f %3f\n', C{:});
fclose(fid);

system(['cat ',headerFile,' temp.txt > ',starFileName]);
!rm temp.txt

%% write Micrograph Star File

for i=1:size(bin1tomograms,1)
micrographArray{i,1} = bin1tomograms{i,2};
micrographArray{i,2} = voltage;
micrographArray{i,3} = sphericalAberration;
micrographArray{i,4} = amplitudeContrast;
micrographArray{i,5} = magnification;
micrographArray{i,6} = detectorPixelSize;
end

C = micrographArray.';
fid = fopen('temp.txt', 'wt');
fprintf(fid, '%s %5f %5f %6f %g %6f\n', C{:});
fclose(fid);

system(['cat ',headerFile,' temp.txt > ',micrographStarFileName]);


%% clean up

!rm temp.txt
clear
close all