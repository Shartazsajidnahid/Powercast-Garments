% Run greedy segment fitting algorithm on the 4 time series.
% FitSegmentFn is
function [b_all, Irhat2, Iihat2, segs] = tensor_LinFitGreedy(Vr, Vi, Ir, Ii, FitSegmentFn, bic_param)

n_d = size(Vr,1);

N = size(Vr,2);
nseg = floor(N/2);
for i=1:nseg
    segs{i} = 2*i-1 : 2*i;
end
segs{nseg} = 2*nseg-1 : N;

segment_costs = zeros(1,nseg);
for i = 1: nseg
    for j = 1 : n_d
        segment_costs(i) = segment_costs(i) + get_segment_cost(transpose(Vr(j,:)), transpose(Vi(j,:)), transpose(Ir(j,:)), transpose(Ii(j,:)), segs{i}, FitSegmentFn, N, bic_param);
    end
end

for i = 1 : nseg-1
    merged_cost = tensor_cost_after_merge(Vr, Vi, Ir, Ii, i, i+1, segs, FitSegmentFn, N, bic_param);
    BIC_delta(i) = merged_cost - segment_costs(i) - segment_costs(i+1);
end

while 1
    nseg_now = length(segs);
    [best_delta, best_idx] = min(BIC_delta);
    if best_delta < 0 % merge
        
        segs{best_idx} = [segs{best_idx} segs{best_idx+1}];
        segs(best_idx+1) = [];
        
        segment_costs(best_idx) = 0;
        for j = 1 : n_d
            segment_costs(best_idx) = segment_costs(best_idx) + get_segment_cost(transpose(Vr(j,:)), transpose(Vi(j,:)), transpose(Ir(j,:)), transpose(Ii(j,:)), segs{best_idx}, FitSegmentFn, N, bic_param);
        end
        
        segment_costs(best_idx+1) = [];
        nseg_now = nseg_now - 1;
        
        % Recompute neighboring deltas, then delete the current delta
        if best_idx > 1
            cost_merge_left = tensor_cost_after_merge(Vr, Vi, Ir, Ii, best_idx-1, best_idx, segs, FitSegmentFn, N, bic_param);
            BIC_delta(best_idx-1) = cost_merge_left - segment_costs(best_idx-1) - segment_costs(best_idx);
        end
        if best_idx < nseg_now
            cost_merge_right = tensor_cost_after_merge(Vr, Vi, Ir, Ii, best_idx, best_idx+1, segs, FitSegmentFn, N, bic_param);
            BIC_delta(best_idx+1) = cost_merge_right - segment_costs(best_idx) - segment_costs(best_idx+1);
        end
        BIC_delta(best_idx) = [];
    else
        break
    end
end

Irhat = [];
Iihat = [];

b_all = cell(n_d, length(segs));
for j = 1: n_d
    for i=1:length(segs)
        seg = segs{i};
        [b_all{j,i}, Irhat_seg, Iihat_seg, ~] = FitSegmentFn(transpose(Vr(j,seg)), transpose(Vi(j,seg)), transpose(Ir(j,seg)), transpose(Ii(j,seg)), N, bic_param);
        Irhat = [Irhat; Irhat_seg];
        Iihat = [Iihat; Iihat_seg];
    end
    Irhat2{j,1} = Irhat';
    Iihat2{j,1} = Iihat';
end


% assert(length(Irhat) == N);
