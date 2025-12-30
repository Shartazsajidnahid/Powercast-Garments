function plot_retensor_res(myTensor, re_myTensor)

figure;
subplot(1,2,1); imagesc(myTensor(:,:,1));colorbar; title('true G');
subplot(1,2,2); imagesc(re_myTensor(:,:,1));colorbar;title('fitted G');

figure;
subplot(1,2,1); imagesc(myTensor(:,:,2));colorbar;title('true alpha_r');
subplot(1,2,2); imagesc(re_myTensor(:,:,2));colorbar;title('fitted alpha_r');

figure;
subplot(1,2,1); imagesc(myTensor(:,:,3));colorbar; title('true B');
subplot(1,2,2); imagesc(re_myTensor(:,:,3));colorbar;title('fitted B');

figure;
subplot(1,2,1); imagesc(myTensor(:,:,4));colorbar;title('true alpha_i');
subplot(1,2,2); imagesc(re_myTensor(:,:,4));colorbar;title('fitted alpha_i');
