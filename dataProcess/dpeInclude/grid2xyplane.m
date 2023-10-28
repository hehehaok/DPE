function [curr_corr_xy_2d, max_idx] = grid2xyplane(curr_corr, numOfGrid, normFlag)
    corr_xy_zdt = reshape(curr_corr, [], numOfGrid^2); % 625*625
    corr_xy = max(corr_xy_zdt); % 1*625
    curr_corr_xy_2d = flipud(reshape(corr_xy, [numOfGrid, numOfGrid])); % 25*25
    [max_corr, max_idx] = max(curr_corr_xy_2d,[],'all','linear'); 
    if(normFlag)
        curr_corr_xy_2d = curr_corr_xy_2d ./ max_corr; % πÈ“ªªØ
    end
end