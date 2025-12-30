function V_real_in = windowing_fn(V_real_in, tt, t_d, window_size, sigma)

half_window = floor(window_size/2);

if tt < ceil(window_size/2)
    V_real_in = V_real_in(1:tt+half_window);
    weight = normpdf([-half_window:1:half_window], 0, sigma);
    weight = weight(end-length(V_real_in)+1:end);
elseif tt > t_d - half_window
    V_real_in = V_real_in(tt-half_window:end);
    weight = normpdf([-half_window:1:half_window], 0, sigma);
    weight = weight(1:length(V_real_in));
else
    V_real_in = V_real_in(tt-half_window:tt+half_window);
    weight = normpdf([-half_window:1:half_window], 0, sigma);
end


% weight = ones(size(weight));


V_real_in = weight' .* V_real_in;