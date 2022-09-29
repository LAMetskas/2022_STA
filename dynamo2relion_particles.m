%% dynamo2relion_particles

% dynamo2relion_particle.m is a script to generate a Relion star file from a
% Dynamo project.  This output is meant to point Relion toward a directory
% of pre-extracted subtomograms, for use in bin1 cases where the tomograms
% are too large to work with directly.

% To use:
% If working with bin1 tomograms that are unwieldly, crop mrc motls first
% before putting into relion.  Then use the particles star file to run a
% re-extraction on the subtomograms in Relion.  

% Next, find the Relion extraction path and input it into the 
% particlesInMicrographs version of this script to generate a star file with
% the paths to the subtomograms but linking them back to their parent tomogram
% (to preserve noise estimates). Use the particlesInMicrographs star file output
% for the 3D Classification run.


% Note this script assumes that the Dynamo tomogram index is the same as the 
% tomogram numerical string.

%% inputs
inputTable = '';
voltage = <>;
sphericalAberration = <>;
amplitudeContrast = <>;
detectorPixelSize = <>;
boxSize = <>;
particleBaseName = 'particle_';
particleExtension = '.mrc';
particlePath = ''; *path to Dynamo crop folder
calibratedPixelSize = <>;

headerFile = 'relionHeader.txt';
starFileName = ''; *chosen star file name for output
micrographStarFileName = 'micrographs.star';

shiftParticleCentersForRecropping=0; %Boolean toggle

%% get Particles

%shift particles into the center of their boxes;
table=dread(inputTable);
if shiftParticleCentersForRecropping==1
    table(:,24)=table(:,4)+table(:,24);
    table(:,25)=table(:,5)+table(:,25);
    table(:,26)=table(:,6)+table(:,26);
    table(:,4:6)=0;
end

%convert Eulers
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

for i=1:length(table)
    tomogramIndex = table(i,20);
    if tomogramIndex < 10
        tomoString = ['0',num2str(tomogramIndex)];
    else
        tomoString = num2str(tomogramIndex);
    end
    starCellArray{i,1} = ([particlePath,particleBaseName,sprintf('%06d',table(i,1)),particleExtension]);
    starCellArray{i,2} = voltage;
    starCellArray{i,3} = sphericalAberration;
    starCellArray{i,4} = amplitudeContrast;
    starCellArray{i,5} = magnification;
    starCellArray{i,6} = detectorPixelSize;
    starCellArray{i,7} = ([particleBaseName,sprintf('%06d',table(i,1)),particleExtension]);
    starCellArray{i,8} = (boxSize/2) + table(i,4);
    starCellArray{i,9} = (boxSize/2) + table(i,5);
    starCellArray{i,10} = (boxSize/2) + table(i,6);
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


%% clean up

!rm temp.txt
clear
close all