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
%---    --- --- --- 	---%%%
%% add all the functions needed
filePath = matlab.desktop.editor.getActiveFilename; %path to this file
foldPath = fileparts(filePath); %find the folder this file is in
%from that folder
addpath(genpath(fullfile(foldPath,'jochensmolka-elf-461d23388273')));%Specific ELF version used previously
addpath(genpath(fullfile(foldPath,'ELF4LP')));%ELF functions for studying light pollution
addpath(genpath(fullfile(foldPath,'kakearney-cptcmap-pkg-845bf83')));%colour map functions from File Exchange
addpath(genpath(fullfile(foldPath,'altmany-export_fig-3175417')));%PDF exporting functions from File Exchange
elf_paths;%run the function for finding ELF paths
%% select folder to process
%EVERY folder should contain a file called "brackets.info" with the number of
%the first and the last image in the bracket (I think in chronological order)
% so if the first bracket contained 5 images, and the second one was three
% brackets.info would have the contents
% 1 5
% 6 9
% the numbers in each row must be separated by tabs, not spaces
if ~exist('inputfolder','var')
    %if no inputfolder has been specified, get the user to select one
    para = elf_para('%');
else %if the user has specificed an inputfolder, make that the root directory
    para = elf_para(inputfolder);
end
%don't be surprised if you get asked for the top folder
%find datasets in the folder
[~, ~, datasets] = elf_checkdata(para);
nohor = NaN(1, 4);%which of these datasets should have the horizon (image edge) removed? (set to NaN)
%In the past I used 70 to include only 90--20° (i.e. exclude 0--20°)

%% convert to HDR scenes
warning('off','imageio:tiffmexutils:libtiffWarning')%I don't think this is relevant
for i = 1:length(datasets)% for each dataset
    elf_main1_HdrAndInt(datasets{i}, '*dng',false, 00)%Convert brackets to HDR. TODO check for rotation here
    night_filter(datasets{i}); close all %filter to 4 different spatial resolutions
end
warning('on','imageio:tiffmexutils:libtiffWarning')%turn it back on for other steps
%% Calculate irradiance and display (can be skipped)
for i = 1:length(datasets)% for each dataset
    allfiles    = dir(fullfile(para.paths.root,datasets{i}, para.paths.filtfolder, '*filt.mat')); %find all filtered files for this dataset
%     pdffolder   = fullfile(para.paths.root,datasets{i});%write out PDFs to the correct folder for this dataset
%     para.paths.datapath = fullfile(para.paths.root,datasets{i});%save the path to the folder for this dataset
    
    % Calculate total irradiance
    for j = 1:length(allfiles)% for each output file, calculate the cosine-corrected irradiance across the image
        irrad = night_irradiance(strcat(allfiles(j).folder, '/',allfiles(j).name)); % in photons / s / m2 / nm
        disp('Irradiance (photons/s/cm2/nm)');
        disp(irrad/10000)%display in photons / cm2 (m2 / cm2)
    end

    % Means of Johnsen et al., 2006 data for starlight:
    % 400-500 (34:66)   :    2.6e6
    % 500-600 (67:100)  :    4.9e6
    % 600-700 (101:133) :    5.1e6
    % 400-700: 4.21

end
%% Reproject (only needs to be done once)
for i = 1:length(datasets)%% for each dataset
    allfiles    = dir(fullfile(para.paths.root,datasets{i}, para.paths.filtfolder, '*filt.mat'));%find all filtered files for this dataset
    pdffolder   = fullfile(para.paths.root,datasets{i});%write out PDFs to the correct folder for this dataset
    para.paths.datapath = fullfile(para.paths.root,datasets{i});%save the path to the folder for this dataset
    rotation = zeros( 1,length(allfiles));%no rotation right now
    for j = 1:length(allfiles)% for each output file, reproject to fisheye shape
                        elf_support_logmsg('         Reprojecting Scene %03d\n', j);%inform the user
        ims             = elf_readwrite(para, 'loadfilt_mat', ['scene' sprintf('%03d', j)]); %load the file
        im_filt_reproj  = cell(size(ims)); %convert to cell
        roti = repmat(rotation(j), 1, length(ims)); %rotate if required
        for ii = 1:length(ims) %for all images
            im_filt_reproj{ii} = elf_project_reproject2fisheye_simple(ims{ii}, [], [], [1000 1000], roti(ii)); % im, azi, ele, imsize
        end
        [~, f] = fileparts(allfiles(j).name); %find the folder to save in 
        save(fullfile(allfiles(j).folder, [f '_reproj.mat']), 'im_filt_reproj'); %save the reprojected image in the same folder
    end
end

%%  Plot and save relative to an appropriate maximum
for i = 1:length(datasets)%% for each dataset
    allfiles    = dir(fullfile(para.paths.root,datasets{i}, para.paths.filtfolder, '*filt.mat'));%find all filtered files for this dataset
    pdffolder   = fullfile(para.paths.root,datasets{i});%write out PDFs to the correct folder for this dataset
    para.paths.datapath = fullfile(para.paths.root,datasets{i});%save the path to the folder for this dataset
    rotation = zeros( 1,length(allfiles));%no rotation right now
    % nmax is the maximum value in photons / s / m2 / nm
    % as used by Foster et al., 2018a: https://www.doi.org/10.1098/rspb.2017.2322
    % 3*10^10 is typical of clear moonless nights
    % 1*10^11 is for moderate light pollution or Aurora Borealis (norrsken)
    % 1*10^12 is for strong light pollution (e.g. Johannesburg)
    % 3*10^12 is good for full moon
    % 3*10^14 is good for security lighting (e.g. streetlights)
    nmoptions = [3*10^10; 3*10^12; 1*10^11; 1*10^12; 3*10^14];

    for j = 1:length(allfiles)% for each output file, plot as a linear image
        irrad = night_irradiance(strcat(allfiles(j).folder, '/',allfiles(j).name)); % in photons / s / m2 / nm
        [mini, indi] = min(abs(log10(mean(irrad(4,:))) - log10(nmoptions)));%find the closest maximum value
        if mini >1 %if out by more than a log unit, inform the user
            disp 'uh oh, no good display value found';
        end
        nmax = nmoptions(indi);%use the minimum
    %now plot them
        [~, f] = fileparts(allfiles(j).name); %find the folder to save in 
        temp    = load(fullfile(allfiles(j).folder,[f,'_reproj']), 'im_filt_reproj');%load the file
        ims     = temp.im_filt_reproj;%extract the reprojected image
        clear temp;%remove th rest of the file
        % Calculate elevation for each pixel (this should work for all images now, since they should all be square and the same size)
        % For non-square images, BUGFIX THIS!
        mid    = [1+(size(ims{1}, 1)-1)/2; 1+(size(ims{1}, 2)-1)/2];        % centre of image
        r_full = 8 * size(ims{1}, 2) / 24;                                 % theoretical value for 24mm high chip that is fully covered by fisheye circular image
        [x, y] = meshgrid(1:size(ims{1}, 1), 1:size(ims{1}, 2));
        r      = sqrt((x-mid(1)).^2 + (y-mid(2)).^2);
        ele    = asind(r / 2 / r_full) * 2;%r/r_full * 90;
        ele_s  = ele; ele_s(y<mid(2)) = -ele_s(y<mid(2)); %make negative values one side of the zenit
            % 1. extract images
        rotit = rotation(j);%find the rotation value
        for ch = 1:4 %for each channel (R, G, B, W)
            for sc = 1:4 %for each filter level (?)
                im = ims{sc}; %find the image
                im(isnan(im)) = 0;%replace any NaNs with zeroes
                if rotit
                    im = rot90(permute(im, [2 1 3]), 2);        % rotate image if necessary 
                else
                    %im = fliplr(im);%BUGFIX I THINK THIS IS MAKING A MESS!
                    im = im;%Do nothing
                end
                if ch < 4
                    im(:, :, [1:ch-1 ch+1:3]) = 0; 
                    sumim{sc, ch}  = sum(im, 3); %add all channels to make a total radiance channel
                else
                    sumim{sc, ch}  = mean(im, 3);
                end   % if colour channel, set other channels to 0 (works for plotting and contrast calculation)
                plotim{sc, ch} = im;                            % save back into ims structure, as this will later be plotted
            end %sc = 1:4
        end %ch = 1:4
        wh_im = plotim{2, 4};% select slot 2, which contains the 4° filtered image
           figure(); %open the figure
       imshow((wh_im)./(nmax));% I don't think we want this flip% imshow(flipud(wh_im)./(nmax));
        xscale = sind((-60:30:60 +90) / 2) * 2000 / 3 +500.5; % x scale degrees across the image (not used?)
        svpath = [ datasets{i} '_' 'scene' sprintf('%03d', j) '_norm_' sprintf('%g', nmax/10000)]; %make the pdf file name; m^2 = 10000 cm^2!
        pdfsave(gcf, fullfile(pdffolder,[ svpath '_specrad.pdf'])); %save the PDF (requires export_fig)
        export_fig(fullfile(pdffolder,[ svpath '_specrad.png']), '-native'); %save to native resolution if possible
    end
  end %i  = 1:length(datasets)
%%  Plot with a colourmap for log10 photons across the whole dynamic range
fhcmp = colourmapsky(para,datasets);% Convert radiance to a colourmap; outputs the figure handle for errors
%% Make a histogram of measurements in each dataset
nmoptions = [3*10^10; 3*10^12; 1*10^11; 1*10^12; 3*10^14]; %using the same maximisation values as before
for i = 1:length(datasets)%% for each dataset
    hst = figure(); %open the figure
    hold on %keep it open
    allfiles    = dir(fullfile(para.paths.root,datasets{i}, para.paths.filtfolder, '*filt.mat'));%find all filtered files for this dataset
    pdffolder   = fullfile(para.paths.root,datasets{i});%write out PDFs to the correct folder for this dataset
    para.paths.datapath = fullfile(para.paths.root,datasets{i});%save the path to the folder for this dataset
    %rotation = zeros( 1,length(allfiles));%no rotation right now
     for j = 1:length(allfiles)% for each output file, reproject to fisheye shape
    %now plot them
        [~, f] = fileparts(allfiles(j).name);%find the file name
        % more sense to use the spherical coordinates ones
        temp    = load(fullfile(allfiles(j).folder,f), 'im_filt_HDR');%load the filtered file
        ims     = temp.im_filt_HDR;%extract just the image
        clear temp;%remove the rest of the file

        % Calculate elevation for each pixel (this should work for all images now, since they should all be square and the same size)
        wh_im = ims{2};%1st slot is filter level: [2,4,8,16]
        sky_1 =  abs((wh_im(:,:,1)+ wh_im(:,:,2)+wh_im(:,:,3)));%make a white image (sum across channels)
        skypix = sky_1(sky_1 ~= 0); %only pixels with nonzero values are counted in the image histogram
        skyrad = log10(skypix/10000);%convert from values per m2 to per cm2
        histogram(skyrad, 'EdgeColor','none');
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