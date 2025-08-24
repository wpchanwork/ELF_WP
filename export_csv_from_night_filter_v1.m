dataset = 'run';


para        = elf_para('reset', dataset, '*.dng');
para        = elf_para_update(para);                                       % Combine old parameter file with potentially changed information in current elf_para
allfiles    = elf_io_dir(fullfile(para.paths.datapath, para.paths.scenefolder, '*.mat'));
fnames_im   = {allfiles.name};                                              % collect image names
infosum     = elf_readwrite(para, 'loadinfosum');                          % loads the old infosum file (which contains projection information)

                    elf_support_logmsg('      Processing %d scenes in environment %s\n', length(fnames_im), dataset);

%% Process one scene at a time
for setnr = 1:length(fnames_im)
    im_filt_HDR = elf_readwrite(para, 'loadfilt_mat', sprintf('scene%03d', setnr));
    data = elf_wrap_rgb(im_filt_HDR{2});                 % This directly reflect which scale in para.ana.scales_deg is picked
    intMean = elf_analysis_datasetmean(data, [], 1, 'logmean');
    
    [~,f]       = fileparts(fnames_im{setnr}); 
    outstat = fullfile(para.paths.datapath, para.paths.filtfolder, [f '.csv']);
    elf_analysis_writestats(intMean, outstat);

    meanIm  = elf_io_readwrite(para, 'loadHDR_tif', fnames_im{setnr});
    h = elf_plot_intSummary_elev(intMean, meanIm, infosum);
    
    para.paths.fname_meanivep_jpg = fullfile(para.paths.outputfolder_pub, [f '_summary.jpg']);
    elf_io_readwrite(para, 'savemeanivep_jpg', '', h.fh);
end
