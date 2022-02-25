function [ ] = show_state(params, res_file, fr, cm, range)
%show_state(params, res_file, fr, cm, range)
%   Display linear drive in 3D
%   
%   Input:  params - model parameters (with cells' coordinates)
%           res_file - path to a file containing the simulation results
%           fr - simulation step to display
%           cm - colormap to use
%           range - linear drive dynamic range

ex = matfile(res_file);

for i=1:length(fr)
    Itot = ex.Itot(:, fr(i))';

    Itot = Itot - range(1);
    Itot = Itot./(range(2)-range(1));
    Itot(Itot<0) = 0;
    Itot(Itot>1) = 1;

    act = Itot.*256;
    act(act>256) = 256;
    act(act<1) = 1;
    clrs = cm(round(act), :);
    scatter3(params.wc(:, 1), params.wc(:, 2), params.wc(:, 3), 10, clrs, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
    daspect([1, 1, 1]);
    set(gca, 'ZDir', 'reverse')
    set(gca, 'YDir', 'reverse')
    
    axis off

    
end




