function [corr_2d, max_idx] = grid2dPrnPlane(curr_corr, numOfGrid, normFlag, disp_axis)
    corr_multid = reshape(curr_corr, numOfGrid, numOfGrid, numOfGrid);
    if strcmp(disp_axis, 'xy') || strcmp(disp_axis, 'yx')
        corr_multid = corr_multid(floor(numOfGrid/2), :, :);
    elseif strcmp(disp_axis, 'xz') || strcmp(disp_axis, 'zx')
        corr_multid = corr_multid(:, floor(numOfGrid/2), :);
    elseif strcmp(disp_axis, 'yz') || strcmp(disp_axis, 'zy')
        corr_multid = corr_multid(:, :, floor(numOfGrid/2));
    end
    corr_2d = squeeze(corr_multid);
    corr_2d = flipud(corr_2d);
    [max_corr, max_idx] = max(corr_2d,[],'all','linear'); 
    if(normFlag)
        corr_2d = corr_2d ./ max_corr; % πÈ“ªªØ
    end
end