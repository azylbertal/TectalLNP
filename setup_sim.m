function [gg] = setup_sim(params)
% [gg] = setup_sim(params)
%   Set up coupling matrix and filters according to model parameters

    
    wc = params.wc;
    frame_dur = params.frame_dur;
    steps_per_frame = params.steps_per_frame;
    ex_sigma = params.ex_sigma;
    in_sigma = params.in_sigma;
    ex_tau_sec = params.ex_tau_sec;
    in_tau_sec = params.in_tau_sec;
    ex_gain = params.ex_gain;
    in_gain = params.in_gain;
    run_steps = params.run_steps;
    hemisphere_correction = params.hemisphere_correction;
    dc = params.dc;
    random_seed = params.random_seed;
    transient_steps = params.transient_steps;
    
    dtSp = frame_dur/steps_per_frame;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
    rng(random_seed);
    
    gg.dtSp = dtSp;

    ex_tau_steps = ex_tau_sec./dtSp;
    in_tau_steps = in_tau_sec./dtSp;


   
    gg.inhibition_downsample = 100;
    gg.transient_steps = transient_steps;
    iht_ex=0:ex_tau_steps*5;
    iht_in=0:gg.inhibition_downsample:in_tau_steps*5;
    
    ihcpl_ex = ex_gain.*exp(-iht_ex/ex_tau_steps)';
    ihcpl_in = in_gain.*exp(-iht_in/in_tau_steps)';
    gg.ihcpl_ex = ihcpl_ex;
    gg.ihcpl_in = ihcpl_in;

    left_inds = wc(:, 2) > 309;
    right_inds = wc(:, 2) <= 309;

    nneur = size(wc, 1);
    slen = run_steps+transient_steps;
    gg.Istm = zeros(slen, nneur, 'single');
    initial_inhibition = (10.*(-exp(-(1:transient_steps)./(transient_steps/5))))';
    gg.Istm(1:transient_steps, :) = repmat(initial_inhibition, 1, nneur);

    gg.dc = dc;

    gg.iht_ex = iht_ex;
    gg.iht_in = iht_in;
    
    disp('Calculating cross connectivity');

    gg.w_ex = zeros(nneur, nneur);
    gg.w_in = zeros(nneur, nneur);
    
    for jj = 1:nneur
        cl1 = wc(jj, :);
        dists = sqrt((cl1(1)-wc(:, 1)).^2 + (cl1(2)-wc(:, 2)).^2 + (cl1(3)-wc(:, 3)).^2)';
        if left_inds(jj)
            same_hemisphere = left_inds';
        else
            same_hemisphere = right_inds';
        end
        hc(same_hemisphere) = 1;
        hc(~same_hemisphere) = hemisphere_correction;
        
        gg.w_ex(jj, :) = exp(-dists.^2./(2*ex_sigma^2)) .* hc;
        gg.w_in(jj, :) = exp(-dists.^2./(2*in_sigma^2)) .* hc;
    end
    gg.w_ex = single(gg.w_ex);
    gg.w_in = single(gg.w_in);
        
    gg.steps_per_frame = steps_per_frame;

end

