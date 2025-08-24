function S = elf_wrap_rgb(img)
% Convert an HxWx3 RGB image into one "data" struct expected by ELF.
% It fills:
%   S.int.(means/std/min/max/median/perc25/perc75/percmin/percmax) -> 4 x H
%      (channels in order: R,G,B,W), each row-wise stat along y (per row)
%   S.totalint.(mean/std/min/max/median/perc25/perc75/percmin/percmax) -> 1 x 4
%   S.totalint.hist -> 1 x B (dummy histogram, not used by writestats)
%   S.totalint.region_meanele -> Hx1 (row “elevation” index)

    img = double(img);
    assert(ndims(img)==3 && size(img,3)==3, 'img must be HxWx3');

    [H, W, ~] = size(img);
    R = img(:,:,1); G = img(:,:,2); B = img(:,:,3);
    Wht = (R + G + B) / 3;                 % “white” channel

    % Helper: row-wise stats (across columns, i.e., along y)
    rowstat = @(M, f) f(M, [], 2);         % min/max use [],2
    rmean   = @(M) mean(M, 2, 'omitnan');
    rstd    = @(M) std(M, 0, 2, 'omitnan');
    rmedian = @(M) median(M, 2, 'omitnan');
    rprc    = @(M,p) prctile(M, p, 2);     % Hx1

    % Build 4 x H arrays (R,G,B,W), then transpose into 4 x H
    cat4 = @(vR,vG,vB,vW) [vR vG vB vW]';

    S.int.means    = cat4(rmean(R),   rmean(G),   rmean(B),   rmean(Wht));
    S.int.std      = cat4(rstd(R),    rstd(G),    rstd(B),    rstd(Wht));
    S.int.min      = cat4(rowstat(R,@min), rowstat(G,@min), rowstat(B,@min), rowstat(Wht,@min));
    S.int.max      = cat4(rowstat(R,@max), rowstat(G,@max), rowstat(B,@max), rowstat(Wht,@max));
    S.int.median   = cat4(rmedian(R), rmedian(G), rmedian(B), rmedian(Wht));
    S.int.perc25   = cat4(rprc(R,25), rprc(G,25), rprc(B,25), rprc(Wht,25));
    S.int.perc75   = cat4(rprc(R,75), rprc(G,75), rprc(B,75), rprc(Wht,75));
    % Set “min percentile”/“max percentile” to 2.5/97.5 to match elf_analysis_writestats labels
    S.int.percmin  = cat4(rprc(R,2.5), rprc(G,2.5), rprc(B,2.5), rprc(Wht,2.5));
    S.int.percmax  = cat4(rprc(R,97.5),rprc(G,97.5),rprc(B,97.5),rprc(Wht,97.5));

    % Whole-image stats per channel (1x4)
    whole = @(M,f) f(M(:), 'omitnan');
    S.totalint.mean    = [whole(R,@mean),   whole(G,@mean),   whole(B,@mean),   whole(Wht,@mean)];
    S.totalint.std     = [std(R(:),0,'omitnan'), std(G(:),0,'omitnan'), std(B(:),0,'omitnan'), std(Wht(:),0,'omitnan')];
    S.totalint.min     = [min(R(:)), min(G(:)), min(B(:)), min(Wht(:))];
    S.totalint.max     = [max(R(:)), max(G(:)), max(B(:)), max(Wht(:))];
    S.totalint.median  = [median(R(:),'omitnan'), median(G(:),'omitnan'), median(B(:),'omitnan'), median(Wht(:),'omitnan')];
    S.totalint.perc25  = [prctile(R(:),25), prctile(G(:),25), prctile(B(:),25), prctile(Wht(:),25)];
    S.totalint.perc75  = [prctile(R(:),75), prctile(G(:),75), prctile(B(:),75), prctile(Wht(:),75)];
    S.totalint.percmin = [prctile(R(:),2.5), prctile(G(:),2.5), prctile(B(:),2.5), prctile(Wht(:),2.5)];
    S.totalint.percmax = [prctile(R(:),97.5),prctile(G(:),97.5),prctile(B(:),97.5),prctile(Wht(:),97.5)];

    % Region “elevation” (first column in writestats)
    S.totalint.region_meanele = (1:H)';     % Hx1

    % Provide a dummy histogram (same length across images). Not used by writestats but
    % elf_analysis_datasetmean expects totalint.hist to exist.
    % Here we use 256-bin histogram of the white channel:
    counts = histcounts(Wht(:), 0:256);     % 1x256, edges 0..255
    S.totalint.hist = counts;               % row vector
end