function DisplayImage(Image, delay) %% Display Image
clims =  [min(Image, [], 'all') 2*std(reshape(Image, 1, []))];
figure
for i = 1:length(Image(1,1,:))
    imagesc(Image(:,:,i), clims)
    colorbar
    pause(delay)
end