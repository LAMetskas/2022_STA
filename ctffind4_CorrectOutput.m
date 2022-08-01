%% ctffind4_CorrectOutput.m

% This is a script to go through the output of ctffind4 and find the tilts
% where it has performed badly (user's choice of automated or manual
% input).  Bad defoci are replaced with a near-neighbor average and
% visually checked to allow manual adjustment.

% The automated version of the script assumes that bad defocus estimations
% will always have a large max resolution number as well.  If this is not
% the case for your data, use the manual version of the script and mark bad
% estimates by hand.

% Dependencies: output files from ctffind4, assumes name and location
% conventions of ## directory names, ##.mrc tilt stack names, and ctffind4
% default output names. Also assumes the presence of "tomo_list.txt" in the
% root directory, which is a text file with one tomogram prefix (##) per
% line.

% Written by Lauren Ann Metskas in May 2018 and updated through December
% 2018. This work is licensed under the Creative Commons Attribution-
% NonCommercial-ShareAlike 4.0 International License. To view a copy of 
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% Inputs
clear
mode = 0; %0 for automatic localization, 1 for manual oversight
rootdir='<path>';
ctfsuffix='_output';
boxsize=0; %0 for 512, 1 for 1024

% For mode 0
threshold=700; %adjust based on last column of output. ~1000 for box size 1024, ~700 for box size 512


%% Initial params
list='tomo_list.txt';
list=dlmread(list);
digits=2;
headerlines=5; %number of header lines in output .txt file
ntomo = size(list,1);
numcols=ceil(ntomo/5);


%% Find bad defocus estimations
for i = 1:ntomo
    clear baddefs
    display(['Starting on tomogram ',num2str(i)]);
    tomostr = sprintf(['%0',num2str(digits),'d'],(list(i)));
    cd([rootdir,tomostr])

    % Skip if this tomogram has already been done
    if exist('ctfdone.txt','file')==2
        continue
    end    
    
    ctf=dlmread([tomostr,ctfsuffix,'.txt'],'',headerlines,0);
    ctf(:,8)=(ctf(:,2)+ctf(:,3))./2; %note this uses the average because the power spectrum is 1 axis. If get two-axis power spectra, replace this.
    numtilts=size(ctf,1);
    
    % Find the bad indices
    check(1:numtilts,1)=0;
    if mode==0
        for j=1:size(ctf,1)
            if ctf(j,7)>threshold
                check(j)=1;
            end
        end
    elseif mode==1
        system(['3dmod ',tomostr,ctfsuffix,'.mrc']);
        wait = 0;
        while wait == 0
            assess_string = input('Please list bad estimates, space delimiter (Start counting at one) \n','s');
            if isempty(assess_string)
                wait = 1;
            else
                baddefs = str2num(assess_string);
                wait = 1;
            end
        end
        if exist('baddefs','var')==0
            !touch ctfdone.txt
            continue
        end
        check(:,2)=1:length(check);
        for z=1:length(check)
            if ismember(check(z,2),baddefs)==1
                check(z,1)=1;
            end
        end
        check(:,2)=[];
    end
    
    %Fix the defocus estimate
    save('checkpoint.mat');
    clear j z
    
    % Get power spectrum and frequencies for Plotting
    tiltstack=dlmread([tomostr,'_output_avrot.txt'],' ',headerlines,0);
    params=dlmread([rootdir,'ctfparams.txt']);
    pixsize=params(1);
    volts=params(2)*1000;
    Cs=params(3)*1e-3;
    ampc=params(4);
    clear params
    lambda=((12.25e-10)/sqrt(volts))*1/(sqrt(1+(volts*1.6e-19)/(2*9.1e-31*(3e8)^2)));
    freq=tiltstack(1,:)*1e10;
    trim=freq(1,:)==0;
    trim(:,1)=0;
    freq(:,trim)=[];
    freqplot=freq*1e-10;
    for j=3:6:size(tiltstack,1)
        pspec(((j-3)/6)+1,:)=tiltstack(j,:);
    end
    
    if boxsize==1
        plotend=563;
        normstart=282;
    elseif boxsize==0
        plotend=274;
        normstart=165;
    end
    
    
    pspec(:,trim)=[];
    for k=1:size(pspec,1)
        pspec(k,1:plotend)=rescale(pspec(k,1:plotend),-3,3);
        normfactor=mean(pspec(k,normstart:plotend));
        pspec(k,:)=pspec(k,:)-normfactor+0.5;
    end
    clear j k
    
    for j=1:size(ctf,1)
        if check(j)==1
            if sum(check)==size(ctf,1)
            ctf(j,2)=3;
            ctf(j,3)=3;
            elseif j==1
                if check(j+1)==0
                    ctf(j,2)=ctf(j+1,2);
                    ctf(j,3)=ctf(j+1,3);
                elseif check(j+1)==1
                    test=0;
                    k=1;
                    trigger=0;
                    while trigger==0
                        test=sum(check(1:k));
                        if k~=test
                            trigger=1;
                        end
                        k=k+1;
                    end
                    ctf(j,2)=ctf(k-1,2);
                    ctf(j,3)=ctf(k-1,3);
                end
                
            elseif sum(check(1:j))==j
                if check(j+1)==0
                    ctf(j,2)=ctf(j+1,2);
                    ctf(j,3)=ctf(j+1,3);
                elseif check(j+1)==1
                    test=0;
                    k=1;
                    trigger=0;
                    while trigger==0
                        test=sum(check(1:k));
                        if k~=test
                            trigger=1;
                        end
                        k=k+1;
                    end
                    ctf(j,2)=ctf(k-1,2);
                    ctf(j,3)=ctf(k-1,3);
                end
                
            elseif j==length(check) || sum(check(j:end))==0
                if check(j-1)==0
                    ctf(j,2)=ctf(j-1,2);
                    ctf(j,3)=ctf(j-1,3);
                elseif check(j-1)==1
                    test=0;
                    k=j;
                    trigger=0;
                    while trigger==0
                        test=sum(check(k:j));
                        if j-k+1~=test
                            trigger=1;
                        end
                        k=k-1;
                    end
                    ctf(j,2)=ctf(k+1,2);
                    ctf(j,3)=ctf(k+1,3);
                end
                
            else
                if sum(check(j:end))==length(check)-j+1
                    if check(j-1)==0
                        ctf(j,2)=ctf(j-1,2);
                        ctf(j,3)=ctf(j-1,3);
                    elseif check(j-1)==1
                        test=0;
                        k=j;
                        trigger=0;
                        while trigger==0
                            test=sum(check(k:j));
                            if j-k+1~=test
                                trigger=1;
                            end
                            k=k-1;
                        end
                        ctf(j,2)=ctf(k+1,2);
                        ctf(j,3)=ctf(k+1,3);
                    end
                else
                    test=0;
                    k=j;
                    trigger=0;
                    while trigger==0
                        test=sum(check(j:k));
                        if abs(k-j+1)~=test
                            trigger=1;
                        end
                        k=k+1;
                    end
                    
                    l=j;
                    test=0;
                    trigger=0;
                    while trigger==0
                        test=sum(check(l:j));
                        if abs(j-l)+1~=test
                            trigger=1;
                        end
                        l=l-1;
                    end
                    
                    lowdist=j-(l+1);
                    highdist=abs(j-(k-1));
                    lowweight=highdist/(lowdist+highdist);
                    highweight=lowdist/(lowdist+highdist);
                    ctf(j,2)=(ctf(l+1,2)*lowweight)+(ctf(k-1,2)*highweight);
                    ctf(j,3)=(ctf(l+1,3)*lowweight)+(ctf(k-1,3)*highweight);
                end
            end
            
            % Plot the new ctf and fix it if needed
            fitctf=((ctf(j,2)+ctf(j,3))/2)*1e-10;
            set1=pi*lambda*freq(1,:).^2;
            set2=fitctf-(0.5*(lambda^2)*(freq(1,:).^2)*Cs);
            fit(:)=abs(sin((set1.*set2)+atan(ampc/sqrt(1-(ampc^2)))));
            
            happy=0;
            while happy==0
                figure(j)
                hold on
                plot(freqplot,pspec(j,:),'k')
                plot(freqplot,fit,'m')
                xlim([0 0.1])
                movegui('southeast')
                judge=input('Are you happy with this fit? Type y or n. ','s');
                if judge=='y'
                    ctf(j,8)=(ctf(j,2)+ctf(j,3))/2;
                    hold off
                    happy=1;
                    close all
                elseif isempty(judge)
                    ctf(j,8)=(ctf(j,2)+ctf(j,3))/2;
                    hold off
                    happy=1;
                    close all
                elseif judge=='n'
                    foo=0;
                    while foo==0
                        givevalue=['The current defocus fit is ',num2str(fitctf)];
                        disp(givevalue)
                        def=input('Please type a defocus value in m to test.');
                        set1=pi*lambda*freq(1,:).^2;
                        set2=def-(0.5*(lambda^2)*(freq(1,:).^2)*Cs);
                        ctfcalc=sin((set1.*set2)+atan(ampc/sqrt(1-(ampc^2))));
                        ctfcalc(:)=abs(ctfcalc(:));
                        figure(j)
                        hold on
                        if exist('h','var')
                            delete(h)
                        end
                        h=plot(freqplot,ctfcalc,'c','LineWidth',2);
                        hold off
                        ctfcheck=input('Are you happy now? Type y or n.   ','s');
                        if ctfcheck=='y'
                            happy=1;
                            foo=1;
                            ctf(j,8)=def*1e10;
                            ctf(j,2)=def*1e10;
                            ctf(j,3)=def*1e10;
                            close all
                        else
                            foo=0;
                            happy=0;
                        end
                    end
                else
                    disp('Invalid input, please try again.')
                    happy=0;
                end
            end
            clear h ctfcheck foo ctfcalc
            check(j)=0;
        end
    end
    
    %write in ctffind4 format with astigmatism info
    dlmwrite('temp1.txt',ctf(:,1:7),'delimiter',' ','precision','%.6f');
    system(['head -n 5 ',tomostr,'_output.txt > temp2.txt']);
    system('cat temp2.txt temp1.txt > defocus_corrected.txt');
    !rm temp*.txt
    
    %write in IMOD ctfphaseflip format without astigmatism info
    csvwrite('ctffinddefocus.csv',ctf(:,2:4));
    tilts=dlmread([tomostr,'.rawtlt']);
    defmean=ctf(:,8)./10;
    numtilts=size(defmean,1);
    setfocus=[0:numtilts-1;0:numtilts-1]';
    setfocus=cat(2,setfocus,tilts,tilts,defmean);
    dlmwrite('setfocus.txt',setfocus,'delimiter',' ','precision',4);
    
    clear check tilts defmean setfocus numtilts ctfoutput ctf pspec set1 set2 tiltstack trim freq freqplot fit defmean
    ! rm checkpoint.mat
    ! touch ctfdone.txt
    
end

%write list of average defocus per stack

%% Plot final outputs for review
cd(rootdir);

figure
xlim([0 80])
alldone=0;

while alldone==0
    for i=1:length(list)
        tomo_str = sprintf(['%0',num2str(digits),'d'],(list(i)));
        %     if exist  (['TS_', tomo_str,'/alldefocus.txt'],'file');
        defocusfilectffind=csvread([tomo_str,'/ctffinddefocus.csv']);
        defocusfilematlab=dlmread([tomo_str,'/setfocus.txt']);
        meandef(i)=mean(defocusfilematlab(:,5))/1000;
        
        
        defocusfilectffind(:,1:2)=round(defocusfilectffind(:,1:2)/10);
        ymax=max(defocusfilematlab(:,5))+200;
        ymin=min(defocusfilematlab(:,5))-200;
        subplot(5,numcols,i)
        plot(defocusfilematlab(:,3),(defocusfilematlab(:,5)).','k*')
        hold all
        subplot(5,numcols,i)
        
        plot(defocusfilematlab(:,3),defocusfilectffind(:,1),'o')
        subplot(5,numcols,i)
        
        plot(defocusfilematlab(:,3),defocusfilectffind(:,2),'o')
        title(tomo_str);
        ylim([ymin ymax])
        xlim([-70 70])
    end
    dlmwrite('allfoci.txt',meandef,'precision',3);
    clicktofinish=input('Press enter to exit','s');
    if isempty(clicktofinish)
        alldone=1;
    end
end



%% Clean up

clear
close all

