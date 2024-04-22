function one_view(inputfolder)
%% skypipe.m
%---------------------------------------------------------------
%       USAGE: one_view(inputfolder) %
%
%      AUTHOR: James Foster                 DATE: 2019 10 17
%    MODIFIED: James Foster                 DATE: 2024 04 22
%    
% DESCRIPTION: A one-stop pipeline for processing a set of single view images. 
%              Adapted from skypipe.m in https://github.com/JJFosterLab/light-pollution.
%
%      INPUTS:  1. Filename
%
%     OUTPUTS:  .pdf figures   
%
%  REFERENCES:  Foster et al., 2021 https://doi.org/10.1016/j.cub.2021.06.038
%---------------------------------------------------------------
%
% DON'T FORGET: EVERY FOLDER NEEDS A brackets.info FILE (see below)
%
% TO DO:
% - 
%
%---    JOCHEN SAYS...  ---%
%
% Filename must include a cell array im_filt_reproj with four elements. By
% default, these are assumed to be filtered at 2/4/8/16 degrees, but this
% only matters for labels. The second input variable 'rotate' should be
% true if the Milky Way is vertical in the image, as the algorithm assumes
% a roughly horizontal orientation.
%---    --- --- --- 	---%%%
%% add all the functions needed
filePath = matlab.desktop.editor.getActiveFilename; %path to this file
addpath(fullfile(fileparts(filePath),'jochensmolka-elf-461d23388273'));
addpath(fullfile(fileparts(filePath),'ELF4LP'));
addpath(fullfile(fileparts(filePath),'kakearney-cptcmap-pkg-845bf83'));
addpath(fullfile(fileparts(filePath),'altmany-export_fig-3175417'));
elf_paths;
%% select folder to process
%EVERY folder should contain a file called "brackets.info" with the number of
%the first and the last image in the bracket (I think in chronological order)
% so if the first backet contained 5 images, and the second one was three
% brackets.info would have the contents
% 1 5
% 6 9
if ~exist('inputfolder','var')
    %if no inputfolder has been specified, get the user to select one
    para = elf_para('%');
else %if the user has specificed an inputfolder, make that the root directory
    para = elf_para(inputfolder);
end
% this step never did anything, let the user select the folder
    % if ~exist('skydir','var')
    %     skydir = fullfile(getenv('HOME'),'/Documents/Starry Project/Light Pollution Imaging/');
    % end
    % if ~exist('filename','var')
    %     inputfolder = fullfile(skydir,'LPexper-20181127-LowMWay-Thornwood');
    % end
%don't be surprised if you get asked for the top folder
% para = elf_para('20181127StarryThornwood_1928_1940');
%find datasets in the folder
[~, ~, datasets] = elf_checkdata(para);
nohor = NaN(1, 4);%for now let's keep the horizon.
%In the past I used 70 to include only 90--20° (i.e. exclude 0--20°)

%% convert to HDR scenes
warning('off','imageio:tiffmexutils:libtiffWarning')%I don't think this is relevant
for i = 1:length(datasets)%[3 6]%
    dataset = datasets{i};
%     night_HDRscenes_quickfixJJF(dataset, nohor(i));%difficulty getting calibration right
    elf_main1_HdrAndInt(dataset, '*dng',false, 90)
    night_filter(dataset); close all
end
warning('on','imageio:tiffmexutils:libtiffWarning')%turn it back on for other steps
%% Calculate irradiance and display (can be skipped)
for i = 1:length(datasets)%[3 6]%
    % allfiles    = dir(fullfile(reprojfolder, '*.mat'));
    allfiles    = dir(fullfile(para.paths.root,datasets{i}, para.paths.filtfolder, '*filt.mat'));
    pdffolder   = fullfile(para.paths.root,datasets{i});%'E:\James data\04Jochen_pdfoutput';
    para.paths.datapath = fullfile(para.paths.root,datasets{i});
    
    % Calculate total irradiance
    for j = 1:length(allfiles)
        irrad = night_irradiance(strcat(allfiles(j).folder, '/',allfiles(j).name)); % in photons / s / m2 / nm
        disp('Irradiance (photons/s/cm2/nm)');
        disp(irrad/10000)
    end

    % Means of Johnsen et al., 2006 data for starlight:
    % 400-500 (34:66)   :    2.6e6
    % 500-600 (67:100)  :    4.9e6
    % 600-700 (101:133) :    5.1e6
    % 400-700: 4.21

end
%% Reproject (only needs to be done once)
for i = 1:length(datasets)%[3 6]%
    allfiles    = dir(fullfile(para.paths.root,datasets{i}, para.paths.filtfolder, '*filt.mat'));
    pdffolder   = fullfile(para.paths.root,datasets{i});%'E:\James data\04Jochen_pdfoutput';
    para.paths.datapath = fullfile(para.paths.root,datasets{i});
    rotation = zeros( 1,length(allfiles));%no rotation right now
    for j = 1:length(allfiles)
                        elf_support_logmsg('         Reprojecting Scene %03d\n', j);
        ims             = elf_readwrite(para, 'loadfilt_mat', ['scene' sprintf('%03d', j)]);
        im_filt_reproj  = cell(size(ims));
        roti = repmat(rotation(j), 1, length(ims));
        for ii = 1:length(ims)
            im_filt_reproj{ii} = elf_project_reproject2fisheye_simple(ims{ii}, [], [], [1000 1000], roti(ii)); % im, azi, ele, imsize
        end
        [~, f] = fileparts(allfiles(j).name);
        save(fullfile(allfiles(j).folder, [f '_reproj.mat']), 'im_filt_reproj');
    end
end

%%  Plot and save relative to an appropriate maximum
for i = 1:length(datasets)%[3 6]%
    allfiles    = dir(fullfile(para.paths.root,datasets{i}, para.paths.filtfolder, '*filt.mat'));
    pdffolder   = fullfile(para.paths.root,datasets{i});%'E:\James data\04Jochen_pdfoutput';
    para.paths.datapath = fullfile(para.paths.root,datasets{i});
    rotation = zeros( 1,length(allfiles));%no rotation right now
    % nmax is the maximum value in photons / s / m2 / nm
    % as used by Foster et al., 2018a: https://www.doi.org/10.1098/rspb.2017.2322
    % 3*10^10 is typical of clear moonless nights
    % 1*10^11 is for moderate light pollution or aurora borealis (norrsken)
    % 1*10^12 is for strong light pollution (e.g. Johannesburg)
    % 3*10^12 is good for full moon
    % 3*10^14 is good for security lighting (e.g. Mushrooms)
    % nmax =   3*10^10; 3*10^12; 1*10^11; 1*10^12;  %VALUE TO MAXIMISE TO
    nmoptions = [3*10^10; 3*10^12; 1*10^11; 1*10^12; 3*10^14];

    for j = 1:length(allfiles)
        irrad = night_irradiance(strcat(allfiles(j).folder, '/',allfiles(j).name)); % in photons / s / m2 / nm
        [mini, indi] = min(abs(log10(mean(irrad(4,:))) - log10(nmoptions)));
        if mini >1
            disp 'uh oh, no good display value found';
        end
        nmax = nmoptions(indi);
    %now plot them
        [~, f] = fileparts(allfiles(j).name);
        temp    = load(fullfile(allfiles(j).folder,[f,'_reproj']), 'im_filt_reproj');
        ims     = temp.im_filt_reproj;
        clear temp;
        % Calculate elevation for each pixel (this should work for all images now, since they should all be square and the same size)
        % For non-square images, BUGFIX THIS!
        mid    = [1+(size(ims{1}, 1)-1)/2; 1+(size(ims{1}, 2)-1)/2];        % centre of image
        r_full = 8 * size(ims{1}, 2) / 24;                                 % theoretical value for 24mm high chip that is fully covered by fisheye circular image
        [x, y] = meshgrid(1:size(ims{1}, 1), 1:size(ims{1}, 2));
        r      = sqrt((x-mid(1)).^2 + (y-mid(2)).^2);
        ele    = asind(r / 2 / r_full) * 2;%r/r_full * 90;
        ele_s  = ele; ele_s(y<mid(2)) = -ele_s(y<mid(2)); 
            % 1. extract images
        rotit = rotation(j);
        for ch = 1:4
            for sc = 1:4
                im = ims{sc};
                im(isnan(im)) = 0;
                if rotit
                    im = rot90(permute(im, [2 1 3]), 2);        % rotate image if necessary to place Milky Way horizontally %SHOULD BE POSSIBLE TO KEEP MW VERTICAL
                else
                    %im = fliplr(im);%BUGFIX I THINK THIS IS MAKING A MESS!
                    im = im;%Do nothing
                end
                if ch < 4
                    im(:, :, [1:ch-1 ch+1:3]) = 0; 
                    sumim{sc, ch}  = sum(im, 3);
                else
                    sumim{sc, ch}  = mean(im, 3);
                end   % if colour channel, set other channels to 0 (works for plotting and contrast calculation)
                plotim{sc, ch} = im;                            % save back into ims structure, as this will later be plotted
            end %sc = 1:4
        end %ch = 1:4
        wh_im = plotim{2, 4};% select slot 2, which contains the 4° filtered image
           figure(); 
       imshow((wh_im)./(nmax));% I don't think we want this flip% imshow(flipud(wh_im)./(nmax));
        %    imshow(flipud(wh_im)./(3*10^10)); [brightMW, dimMW] = ginput(2);  %USER SELECTS BRIGHT AND DIM ENDS OF MILKY WAY       
        %    px = %PIXELS AWAY FROM HORIZON
        % elevat = asind((Px - 500.5) * 3 / 2000) * 2; %THE ELEVATION
        xscale = sind((-60:30:60 +90) / 2) * 2000 / 3 +500.5;
        svpath = [ datasets{i} '_' 'scene' sprintf('%03d', j) '_norm_' sprintf('%g', nmax/10000)]; %why? m^2 = 10000 cm^2!
        pdfsave(gcf, fullfile(pdffolder,[ svpath '_specrad.pdf']));
        export_fig(fullfile(pdffolder,[ svpath '_specrad.png']), '-native');
    end
  end %i  = 1:length(datasets)
%%  Plot with a colourmap for log10 photons across the whole dynamic range
fhcmp = colourmapsky(para,datasets);% outputs the figure handle for errors
%% Make a histogram of measurements in each dataset
nmoptions = [3*10^10; 3*10^12; 1*10^11; 1*10^12; 3*10^14];
for i = 1:length(datasets)%[3 6]%
    hst = figure();
    hold on
    allfiles    = dir(fullfile(para.paths.root,datasets{i}, para.paths.filtfolder, '*filt.mat'));
    pdffolder   = fullfile(para.paths.root,datasets{i});%'E:\James data\04Jochen_pdfoutput';
    para.paths.datapath = fullfile(para.paths.root,datasets{i});
    rotation = zeros( 1,length(allfiles));%no rotation right now
     for j = 1:length(allfiles)
    %now plot them
        [~, f] = fileparts(allfiles(j).name);
        % I had been using the reprojected images, but it probably makes
        % more sense to use the spherical coordinates ones
        temp    = load(fullfile(allfiles(j).folder,f), 'im_filt_HDR');
        ims     = temp.im_filt_HDR;
        clear temp;

        % Calculate elevation for each pixel (this should work for all images now, since they should all be square and the same size)
    wh_im = ims{2};%1st slot is filter level: [2,4,8,16]
    sky_1 =  abs((wh_im(:,:,1)+ wh_im(:,:,2)+wh_im(:,:,3)));
    skypix = sky_1(sky_1 ~= 0);
    skyrad = log10(skypix/10000);
    % histogram(skypix/10000, 'EdgeColor','none');
    histogram(skyrad, 'EdgeColor','none');
    % set(gca,'xscale','log')
     end

    %set axes to ELF4LP default
    loax = round(min(log10(nmoptions)-4-2));
    upax = round(max(log10(nmoptions)-4)+3);
    goodax = linspace(loax,upax, abs(upax-loax)+1);
    ylim([0,1.5e4])
    xlim([min(goodax), max(goodax)])
    lbnm = strcat(repmat({'10^{'},1,length(goodax))',num2str(goodax'), repmat({'}'},1,length(goodax))');
    set(gca,'XTick',goodax, 'XTickLabels', lbnm);

    legend(allfiles.name)
    title(datasets{i})
    xlabel('Irradiance (photons/s/cm2/nm)')
    pdfsave(hst, fullfile(pdffolder,[ strcat(datasets{i}) '_hist.pdf']));
    hold off
end %i  = 1:length(datasets)
end