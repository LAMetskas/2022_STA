%% relion2dynamo_classes.m

% relion2dynamo_classes.m is a script that converts a Relion 3D classification
% star file into a Dynamo table.

% Inputs include 1-2 Selection star files containing the classes of
% interest. Only particle identities are transferred. Input is star file from 
% a Relion Select run.

% Written by Lauren Ann Metskas on 15 May 2020.  Please credit if using. This
% work is licensed under the Creative Commons Attribution-ShareAlike 4.0
% International License.  To view a copy of this license, visit
% http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative
% Commons, PO Box 1866, Mountain View, CA 94042, USA.



%% Inputs

mainDirectory = ''; %path to Dynamo project directory
starFilepath1 = ''; %path to folder containing Relion star file 1
%starFilepath2 = ''; %path to folder containing Relion star file 2. Uncomment to use.

dynamoTableName = '<*.tbl>';
newTableName = '<*.tbl>';

%% Make a list of particle ID numbers from Relion star files

cd(starFilepath1)
!grep -o "particle_.............mrc" particles.star > temp.txt
!cat temp.txt | cut -d '_' -f2 | cut -d '.' -f1 > temp2.txt
!sed 's/......$//' temp2.txt > particleIDs.txt
!rm temp.txt 
!rm temp2.txt

if exist('starFilepath2','var')~=0
    cd(starFilepath2)
    !grep -o "particle_.......mrc" particles.star > temp.txt
    !cat temp.txt | cut -d '_' -f2 | cut -d '.' -f1 > particleIDs.txt
    !rm temp.txt
    cd(mainDirectory)
    system(['cat ',starFilepath1,'particleIDs.txt ',starFilepath2,'particleIDs.txt > allSelectedParticles.txt']);
else
    system(['mv ',starFilepath1,'particleIDs.txt ',mainDirectory,'allSelectedParticles.txt']);
    cd(mainDirectory)
end

clear starFilepath1 starFilepath2



%% Reassemble the Dynamo table containing only the Relion selections

particles = dlmread('allSelectedParticles.txt');
dynamoTable = dread(dynamoTableName);

% populate the new particles table with the Euler angles and other 
% parameters from the Dynamo run

for i=1:length(particles)
    [r,c]=find(dynamoTable(:,1)==particles(i,1));
    particles(i,1:(size(dynamoTable,2))) = dynamoTable(r,1:(size(dynamoTable,2)));
end

dwrite(particles,newTableName);



%% Clean up

clear
close all