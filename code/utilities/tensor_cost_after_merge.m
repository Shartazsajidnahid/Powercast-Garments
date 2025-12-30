function cost = tensor_cost_after_merge(Vr, Vi, Ir, Ii, ind1, ind2, segs, FitSegmentFn, N, bic_param)

cost = 0;
for j = 1: size(Vr,1)
    cost = cost + get_segment_cost(transpose(Vr(j,:)), transpose(Vi(j,:)), transpose(Ir(j,:)), transpose(Ii(j,:)),[segs{ind1} segs{ind2}], FitSegmentFn, N, bic_param);
end