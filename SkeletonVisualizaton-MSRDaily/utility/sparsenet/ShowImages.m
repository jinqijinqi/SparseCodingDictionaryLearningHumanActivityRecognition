%% Display the images
figure(3); colormap(gray);
for(n=1:num_images)
    this_image=reshape(IMAGES(:,n),image_size,image_size)';
    imagesc(this_image);
    pause;
end