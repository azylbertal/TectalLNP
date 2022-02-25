% ========================================================================
function [tsp,sps,Itot,Istm] = run_sim(glmprs,verbose)
% [tsp,sps,Itot,Istm] = run_sim(glmprs,verbose)
% 
% Compute ongoing activity (multi-cell) coupled-glm.
%
% Uses time rescaling to sample conditionally Poisson process
% For speed, maintain downsampled version of the slow inhibitory input,
% which is upsampled to match spiking timebase as needed.
%
% Heavily based on code from Pillow et al, Nature 2008
% https://github.com/pillowlab/GLMspiketools
%
%   Input:  glmprs - struct with model paramteres
%           verbose - whether to display bin-specific information
%
%   Output: tsp - imaging time base
%           sps - spikes in imaging frame time base (without transient
%           response)
%           Itot - linear drive to each neuron (simulation time base)
%           Istm - external stimulus (simulation time base)
%   


if ~exist('verbose', 'var')
    verbose = true;
end

glmprs.nlfun = @exp;

ncells = size(glmprs.Istm,2);
rlen = size(glmprs.Istm, 1);

nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs.dtSp;

hlen = length(glmprs.ihcpl_ex);
hlen_in = length(glmprs.ihcpl_in);
% -------------  Compute filtered resp to signal ----------------
Istm = glmprs.Istm;
Istm = bsxfun(@plus,Istm,glmprs.dc); % add DC term
I_ex = Istm'; % total filter output
rlen_in = size(I_ex, 2)/glmprs.inhibition_downsample;
I_in = zeros(size(I_ex, 1), rlen_in+2, 'single');


% --------------- Set up simulation dynamics variables ---------------
Itot = I_ex;
tsp(1,1:ncells) = {zeros(round(rlen/25),1)};  % allocate space for spike times
nsp = zeros(1,ncells);
jbin = 1;  % counter ?
tspnext = exprnd(1,1,ncells); % time of next spike (in rescaled time) 
rprev = zeros(1,ncells); % Integrated rescaled time up to current point
ex_filt = single(glmprs.ihcpl_ex');
in_filt = single(glmprs.ihcpl_in');

spcells=[];
spnums = zeros(1, rlen);
% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    tic;
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen); % Bins to update in this iteration
    nii = length(iinxt);  % Number of bins
    next_weight = jbin/glmprs.inhibition_downsample - floor(jbin/glmprs.inhibition_downsample);
    interp_inhibition = (1-next_weight) .* I_in(:, floor(jbin/glmprs.inhibition_downsample) + 1) ...
    +(next_weight) .* I_in(:, floor(jbin/glmprs.inhibition_downsample) + 2);
    total_input = I_ex(:, iinxt) - interp_inhibition;
    Itot(:, iinxt) = total_input;
    rrnxt = (glmprs.nlfun(total_input)*dt)'; % Cond Intensity
    rrcum = cumsum(rrnxt+[rprev;zeros(nii-1,ncells)],1);  % Cumulative intensity
    if all(tspnext >= rrcum(end,:)) % No spike in this window
        spnums(jbin) = 0;
        jbin = iinxt(end)+1;
        rprev = rrcum(end,:);
    else   % Spike!
        [ispks,jspks] =  find(rrcum>=repmat(tspnext,nii,1));
        spcells = unique(jspks(ispks == min(ispks))); % cell number(s)
        ispk = iinxt(min(ispks)); % time bin of spike(s)
        rprev = rrcum(min(ispks),:); % grab accumulated history to here

        % Record this spike
        mxi = min(rlen, ispk+hlen); % determine bins for adding h current
        iiPostSpk = ispk+1:mxi;
        
        ispk_in = floor(ispk / glmprs.inhibition_downsample);
        mxi_in = min(rlen_in, ispk_in+hlen_in);
        iiPostSpk_in = ispk_in+1:mxi_in;
        
        ex_weights = glmprs.w_ex(:, spcells);
        in_weights = glmprs.w_in(:, spcells);
        a = repmat(ex_filt, length(spcells), 1);
        b = repmat(in_filt, length(spcells), 1);
        outgoing_total_ex = ex_weights*a;
        outgoing_total_in = in_weights*b;
        next_weight = ispk/glmprs.inhibition_downsample - floor(ispk/glmprs.inhibition_downsample);

        if ~isempty(iiPostSpk)
            I_ex(:, iiPostSpk) = I_ex(:, iiPostSpk)+outgoing_total_ex(:, 1:mxi-ispk);
            I_in(:, iiPostSpk_in) = I_in(:, iiPostSpk_in)+(1-next_weight).*outgoing_total_in(:, 1:mxi_in-ispk_in);
            I_in(:, iiPostSpk_in+1) = I_in(:, iiPostSpk_in+1)+(next_weight).*outgoing_total_in(:, 1:mxi_in-ispk_in);
        end


        for ic = 1:length(spcells)
            icell = spcells(ic);
            nsp(icell) = nsp(icell)+1;
            tsp{icell}(nsp(icell),1) = ispk*dt;
            rprev(icell) = 0;  % reset this cell's integral
            tspnext(icell) = exprnd(1); % draw RV for next spike in this cell
        end
        spnums(jbin) = length(spcells);
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/(sum(nsp));
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
    if verbose
        disp(['Bin: ', num2str(jbin), ' Spiking cells: ', num2str(length(spcells)), ' Time:', num2str(toc)]);
    end

end
    % Remove any extra bins from tsp and compute binned spike train 'sps'
sps = sparse(zeros(rlen/glmprs.steps_per_frame,ncells));
for jj = 1:ncells
    tsp{jj} = tsp{jj}(1:nsp(jj));
    if ~isempty(tsp{jj})
        sps(:, jj) = histcounts(tsp{jj}, 0:dt*glmprs.steps_per_frame:(rlen*dt));
    end
end

sps(1:round(glmprs.transient_steps/glmprs.steps_per_frame), :) = [];
end

