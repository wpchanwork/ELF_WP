function night_filter(dataset, horlimit, verbose)
% NIGHT_FILTER filters all HDR scenes for an environment for a 1
% degree and 10 degree resolution. The filtering algorithm (almost) fully corrects
% for the distortion introduced by the equirectangular projection.
%
% Uses: elf_support_logmsg, elf_paths, elf_para, elf_para_update, 
%       elf_readwrite, elf_support_formatA4, elf_filter 
%
% Loads files: HDR scenes as mat in scene folder
% Saves files: filtered images as mat in filt folder
% 
% Typical timing PER SCENE (on ELFPC):
%     192s for filtering both resolutions
%     2.5s to save jpg
%     1.5s rest
%
%   + 192s PER ENVIRONMENT to calculate initial correction images

%% check inputs
if nargin < 3, verbose = true; end
if nargin < 2 || isempty(horlimit), horlimit = NaN; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- nightELF Step 4: Filtering scenes -----\n');

%% Set up paths and file names; read info, infosum and para
elf_paths;
para        = elf_para('reset', dataset, '*.dng');
para        = elf_para_update(para);                                       % Combine old parameter file with potentially changed information in current elf_para
allfiles    = elf_io_dir(fullfile(para.paths.datapath, para.paths.scenefolder, '*.mat'));
fnames_im   = {allfiles.name};                                              % collect image names
infosum     = elf_readwrite(para, 'loadinfosum');                          % loads the old infosum file (which contains projection information)

                    elf_support_logmsg('      Processing %d scenes in environment %s\n', length(fnames_im), dataset);
                 
%% Set Filtering scales (specific for Milky Way)
para.ana.scales_deg = [2 4 8 3];
keyscale = 4; % the key scale that will be visualized and saved
                    
%% Process one scene at a time
it = 1;
for setnr = 1:length(fnames_im)
    %% Load HDR image
    im_HDR  = elf_readwrite(para, 'loadHDR_mat', fnames_im{setnr});
       
    %% Filter HDR image
    % panel 1: original image
    fh      = elf_support_formatA4(4);
    p1      = uipanel('Parent', fh, 'Position', [0 1/3 1 2/3]);
              elf_plot_image(im_HDR, infosum, p1, 'equirectangular', 1);
    ah1     = axes('Parent', p1, 'units', 'normalized', 'position', [0 0 .05 1]);
              text(0.5, 0.5, 'HDR image (delinearised for display)', 'units', 'normalized', 'position', [0 0.5], 'rotation', 90, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontweight', 'bold');
              axis(ah1, [0 1 0 1], 'off');
              
    % panel 2: 2 deg filtered image
    p2      = uipanel('Parent', fh, 'Position', [0 0 .5 1/3]);
    ah(1)   = axes('Parent', p2, 'Position', [0 0 1 1]); % the first scale in the para.ana.scales_deg
    
    % panel 3: 3 deg filtered image
    p3      = uipanel('Parent', fh, 'Position', [0.5 0 .5 1/3]);
    ah(2)   = axes('Parent', p3, 'Position', [0 0 1 1]);  % the second scale in the para.ana.scales_deg
            axis(ah(1:2), 'off');
            set(fh, 'Name', sprintf('Scene #%d of %d', setnr, length(allfiles)));
            drawnow;
    ah(3:4) = ah(1:2);
    
    %% Filter images
    if it==1
        % on first iteration, calculate correction image for borders, and save projection information (%TODO)
        [im_filt_HDR, para] = elf_filter(para, im_HDR, 'scene', verbose, ah, infosum, 1);
        % elf_readwrite(para, 'saveinfosum', [], infosum); % saves infosum AND para
    else
        im_filt_HDR = elf_filter(para, im_HDR, 'scene', verbose, ah, infosum);
    end
    it      = it+1;
    
    %% save projected and filtered images
    elf_readwrite(para, 'savefilt_mat', sprintf('scene%03d', setnr), im_filt_HDR);

    [~,f]       = fileparts(fnames_im{setnr});

    % save visualization
    % Get the displayed images
    panelID = 2 -mod(keyscale, 2);
    hImg10 = findobj(ah(panelID), 'Type', 'image');
    blur10 = hImg10(end).CData;     % the second scale in the para.ana.scales_deg filtered pixels
    
    % Save with no white border (PNG). Normalize if needed.
    imwrite(normalize8(blur10), fullfile(para.paths.datapath, para.paths.filtfolder, sprintf('%s_blur%g.png', f, para.ana.scales_deg(keyscale))));

    %Calculate summary
    im_filt_HDR = elf_readwrite(para, 'loadfilt_mat', sprintf('scene%03d', setnr));
    data = elf_wrap_rgb(im_filt_HDR{keyscale});                 % This directly reflect which scale in para.ana.scales_deg is picked
    intMean = elf_analysis_datasetmean(data, [], 1, 'logmean');
        
    outstat = fullfile(para.paths.datapath, para.paths.filtfolder, sprintf('%s_blur%g.csv', f, para.ana.scales_deg(keyscale)));
    elf_analysis_writestats(intMean, outstat);

    meanIm  = elf_io_readwrite(para, 'loadHDR_tif', fnames_im{setnr});
    h = elf_plot_intSummary_elev(intMean, meanIm, infosum);
    
    para.paths.fname_meanivep_jpg = fullfile(para.paths.outputfolder_pub, sprintf('%s_blur%g_summary.jpg', f, para.ana.scales_deg(keyscale)));
    elf_io_readwrite(para, 'savemeanivep_jpg', '', h.fh);
end
    
                    elf_support_logmsg('      Summary: All HDR scenes for environment %s have been filtered and saved to mat.\n\n', para.paths.dataset);




