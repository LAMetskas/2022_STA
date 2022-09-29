%% BulkInputOfImodModels.m

% BulkInputOfImodModels is a script that can be run to convert a
% number of xyz text files to grids of points, and finally to Dynamo model
% files.

% This script generates a grid of xyz coordinates based on a polygon with 
% user-defined coordinates. The polygon can have any number of sides in 
% xyz space.

% To use: in imod, click anywhere on the lowest z slice you wish to include.
% Next, click a series of points to define the shape you wish, at any
% location in z.  Finally, click anywhere on the top slice you want to
% include. Save as an IMOD model, then convert to a text file using
% model2points.  If you were using a different binning than the tomogram
% to be cropped, note this in the binsize to convert.

% If you wish to input a more defined shape in 3D, in Slicer shift the x
% value to 90 and draw an outline of the widest part of the item. In the
% model2point, save as CB_xz(k).txt or a similar name (changing in the
% code).

% Tomogram numbers appear in column 20; different objects appear in column
% 21.

% Written by Lauren Ann Metskas on 20 May 2020.  Please credit if using. This
% work is licensed under the Creative Commons Attribution-ShareAlike 4.0
% International License.  To view a copy of this license, visit
% http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative
% Commons, PO Box 1866, Mountain View, CA 94042, USA.



%% inputs
%paths, tomos to include
mainDirectory = ''; %full path to tomogram directory
subtomoDirectory = ''; %full path to Dynamo project directory
catalogueName = ''; %Dynamo catalog name where tomograms are located
xyPlaneName = ''; %root name only
xzPlaneName = ''; %optional, root name only
listMin = ##; %the first tomogram index to include
listMax = ##; %the last tomogram index to include (continuous from first)

%grid choices
binsize = <>; %tomogram binning
check = <0/1>; %Boolean toggle, 0 = treat all pixel sizes as absolute, 1 = unbin
xzPlane = <0/1>; %Boolean toggle, 0 = xy file only, 1 = xy and xz files
spacing = <# pixels>; %(bin1 pixels if using check=1, otherwise binned pixels)
characters = 2; %assumes folders are in a ## naming scheme
structureName = ''; %name for the model files
trimEdges = 0; %Boolean toggle for xy polygon only. 1 = trimming away corners, 0=not
cornerLength = <#>; %number of unit cells to delete in and down




%% Main Code

cd(mainDirectory);
tomoList = dlmread('tomo_list.txt');
if listMin > 1
    tomoList(1:(listMin-1))=[];
end
tomoList(((listMax-listMin)+2):end)=[];
numTomos = length(tomoList);

for i=1:numTomos
    tomoString = sprintf(['%0',num2str(characters),'d'],(tomoList(i)));
    cd([mainDirectory,tomoString]);
    [status,j] = system(['ls ',xyPlaneName,'*.txt | wc -l']);
    if length(j)>3
        j=0;
    else
        j=str2num(j);
    end
    clear status
    k=1;
    while j>=1 && k<=j
        % get coordinates from imod output
        cd([mainDirectory,tomoString]);
        coordsxy=dlmread([xyPlaneName,num2str(k),'.txt']);
        zbottom=coordsxy(1,3);
        ztop=coordsxy(end,3);
        coordsxy(end,:)=[];
        coordsxy(1,:)=[];
        
        % convert to unbinned coordinates if desired
        if check == 1
            binned=coordsxy;
            coordsxy=binsize*coordsxy-(binsize/2);
            coordsxz=binsize*coordsxz-(binsize/2);
            zbottom=binsize*zbottom-(binsize/2);
            ztop=ztop*binsize-(binsize/2);
        end
        
        %generate a square point lattice that includes chosen area
        maxx=(max(coordsxy(:,1))-min(coordsxy(:,1)));
        maxy=(max(coordsxy(:,2))-min(coordsxy(:,2)));
        sidelength=ceil(max(maxx, maxy));
        xqv=min(coordsxy(:,1)):spacing:(min(coordsxy(:,1))+sidelength);
        yqv=min(coordsxy(:,2)):spacing:(min(coordsxy(:,2))+sidelength);
        [xq,yq]=meshgrid(xqv,yqv);
        
        %find points inside of chosen area
        xv=coordsxy(:,1);
        yv=coordsxy(:,2);
        in=inpolygon(xq,yq,xv,yv);
        
        %generate xyz coordinates of chosen points
        %try to come up with a way to preallocate the 'points' space
        zplanes=zbottom:spacing:ztop;
        points=[0 0 0];
        for itez=1:length(zplanes)
            currz=zplanes(itez);
            for itey=1:size(in,1)
                for itex=1:size(in,2)
                    if in(itex,itey)==1
                        points(end+1,:)=[xq(itex,itey) yq(itex,itey) currz];
                    end
                end
            end
        end
        points(1,:)=[];
        
        
        %trim corners or use an xz plane dataset to define the shape better
        %useful for shapes that get narrower near the top
        
        %identify points lying outside a model in the xz plane
        if xzPlane == 1
            coordsxz=dlmread([xzPlaneName,num2str(k),'.txt']);
            
            %bin/unbin if indicated
            if check == 1
                coordsxz=binsize*coordsxz-(binsize/2);
            end
            
            clear xv yv in
            xv=coordsxz(:,1);
            yv=coordsxz(:,3);
            xq=points(:,1);
            yq=points(:,3);
            in=inpolygon(xq,yq,xv,yv);
            newset=points(in,:);
            
        elseif trimEdges == 1
            %initialize loop values and params
            newset = points;
            l = cornerLength;
            m = 0;
            while l>0
                %define boundaries
                xlow=min(points(:,1))+(spacing*(l-1));
                xhigh=max(points(:,1))-(spacing*(l-1));
                ylow=min(points(:,2))+(spacing*(l-1));
                yhigh=max(points(:,2))-(spacing*(l-1));
                zlow=min(points(:,3))+(spacing*m);
                zhigh=max(points(:,3))-(spacing*m);
                
                %use logical indexing to locate points outside boundaries and
                %remove them
                cuts = ((newset(:,1)<=xlow | newset(:,1)>=xhigh) | (newset(:,2)<=ylow | newset(:,2)>=yhigh)) & (newset(:,3)==zlow | newset(:,3)==zhigh);
                newset(cuts,:,:)=[];
                l=l-1;
                m=m+1;
            end
            points=newset;
        end
        
        %write a file named points.txt with the output coordinates, plus files
        %recording the parameters used for the point generation for record-keeping
        dlmwrite(['points',num2str(k),'.txt'],points);
        dlmwrite(['vertices',num2str(k),'.txt'],coordsxy);
        params=[zbottom ztop spacing binsize check];
        dlmwrite(['pointparams',num2str(k),'.txt'],params);
        
        %make a Dynamo model
        cd(subtomoDirectory);
        m=dmodels.general();
        m.addPoint(points);
        m.updateCrop();
        dcm('-c',catalogueName,'-index',num2str(tomoString),'-add_model',m,'-modelname',[structureName,num2str(k),'.omd']);
%         dcm('-c',catalogueName,'-index',num2str(i),'-add_model',m,'-modelname',[structureName,num2str(k),'.omd']);
        
        %clean up and advance the loop
        clear currz in itex itey itez m maxx maxy params points sidelength xq xqv xv yq yqv yv zbottom zplanes ztop
        k=k+1;
    end
end


%% clean up
clear
close all