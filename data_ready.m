
% data = readmatrix('garments_raw.xlsx');
% save('garments_raw.mat', 'data');

a = load('garments_raw_train.mat');

data = a.data;

clean_data = [];

window_size = 30;

for i = 1 : floor(window_size): size(data,1) - window_size*2
    % i
    % F = extract_features(data(i:i+window_size*2,2),data(i:i+window_size*2,3),30,60);
    F = extract_features(data(i:i+window_size*2,5),data(i:i+window_size*2,6),25,50);

    clean_data = [clean_data;F];

end

data = clean_data;
save('clean_data_train.mat', 'data')

