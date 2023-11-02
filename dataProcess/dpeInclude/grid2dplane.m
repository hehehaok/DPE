function [corr_2d, max_idx] = grid2dplane(curr_corr, numOfGrid, normFlag, disp_axis)
    if strcmp(disp_axis, 'xy') || strcmp(disp_axis, 'yx')
        remove_axis = 1;
    elseif strcmp(disp_axis, 'xz') || strcmp(disp_axis, 'zx')
        remove_axis = 2;
    elseif strcmp(disp_axis, 'yz') || strcmp(disp_axis, 'zy')
        remove_axis = 3;
    end
    corr_multid = reshape(curr_corr, numOfGrid, numOfGrid, numOfGrid);
    corr_2d = max(corr_multid, [], remove_axis); 
    corr_2d = squeeze(corr_2d);
    corr_2d = flipud(corr_2d);
    [max_corr, max_idx] = max(corr_2d,[],'all','linear'); 
    if(normFlag)
        corr_2d = corr_2d ./ max_corr; % πÈ“ªªØ
    end
end