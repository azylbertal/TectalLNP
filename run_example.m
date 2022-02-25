%% run simulation and save results
load('optimized_params.mat');
glmprs = setup_sim(params);
[~,sps,Itot, Istm,] = run_sim(glmprs); % Simulate GLM response
save('sim_results.mat', 'sps', 'Itot', 'Istm', '-v7.3');

%% Animate results and save as an .avi file
load('optimized_params.mat');
figure
v = VideoWriter('burst_anim.avi');
v.FrameRate = 80;
open(v);
range = [-30, 30];
hot_map = hot(256);
basemap = zeros(256, 3);
basemap(1:127, 2) = ((1:127) - 1) / 127;
basemap(128:end, :) = hot_map((128:256)-10, :);

            for fr = 18000:20000
                show_state(params, 'sim_results.mat', fr, basemap, range)

                k = gca;

                if fr==51000
                    pos = k.CameraPosition(3);
                    tar = k.CameraTarget;
                end
                k.CameraViewAngleMode = 'manual';
                axis equal
                view(mod(fr./5, 360), 45);
                frm = getframe(gcf);
                writeVideo(v,frm);
                drawnow

            end
            close(v);
