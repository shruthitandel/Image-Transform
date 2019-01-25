%Project 2 Q :1.2

clc;
close all;

 src_img1=imread('paint-orig.png');
  %src_im = rgb2gray(src_img1);
  %src_img1 = im2double(src_img); 
  target_img=imread('image-rotate.png');
%  target_img2 = rgb2gray(target_img);
%  target_img1 = im2double(target_img);
% %src_img = checkerboard(4,2);
%src_img = checkerboard;
subplot(1,3,1);
imshow(src_img1);
[srcx, srcy] = getpts;
title('Original image');

subplot(1,3,2);
%figure;
imshow(target_img);
[trgx,trgy] = getpts;
title('Target image');

 X=[srcx(1,1),srcy(1,1),1 ,0,0,0; 
   0,0,0,srcx(1,1),srcy(1,1),1;
   srcx(2,1),srcx(2,1),1,0,0,0;
   0,0,0,srcx(2,1),srcy(2,1),1;
   srcx(3,1),srcy(3,1),1,0,0,0;
   0,0,0,srcx(3,1),srcy(3,1),1];


X2=[trgx(1,1) trgy(1,1) trgx(2,1) trgy(2,1) trgx(3,1) trgy(3,1)]';
affine_v= pinv(X)*X2;

aff_mat=[affine_v(1,1) affine_v(2,1) affine_v(3,1); affine_v(4,1) affine_v(5,1) affine_v(6,1); 0 0 1];
result = zeros([size(src_img1,1),size(src_img1,2)]);

for i=1:size(result,1)
    for j=1:size(result,2)
         current_v=[i;j;1];
         new_v=mtimes(inv(aff_mat),current_v);
            if(new_v(1,1)>=1 && new_v(2,1)>=1 && new_v(1,1)<=size(src_img1,1) && new_v(2)<=size(src_img1,2))
                 intensity = src_img1(round(new_v(1,1)),round(new_v(2,1)));
                 result(i,j) = intensity;
            end         
    end               
end
 
subplot(1,3,3);
imshow(result);
title('Resulting Image');
% 
% 
