%Ethan Nagelvoort
%821234668
%HW1

%1
I = imread('Flooded_house.jpg','jpg');
figure(1);
imshow(I);
title('RGB image');

%2
Red = I;
Green = I;
Blue = I;
for columns=1:704
    for rows=1:536
        Green(rows,columns,1)=0; 
        Green(rows,columns,3)=0;
        Red(rows,columns,2)=0; 
        Red(rows,columns,3)=0;
        Blue(rows,columns,1)=0; 
        Blue(rows,columns,2)=0;
    end
end
figure(2);
imshow(Red);
title('Red Band image');
figure(3);
imshow(Green);
title('Green Band image');
figure(4);
imshow(Blue);
title('Blue Band image');

%3.1
I = rgb2ycbcr(I);
figure(5);
imshow(I);
title('YCbCr color space');

%3.2
I = ycbcr2rgb(I);

%4
ycbcr = rgb2ycbcr(I);
figure(6);
imshow(ycbcr(:,:,1));
title('Y Band image');
figure(7);
imshow(ycbcr(:,:,3));
title('Cr Band image');
figure(8);
imshow(ycbcr(:,:,2));
title('Cb Band image');

%5
cb = ycbcr(:,:,2);
cr = ycbcr(:,:,3);
for row = 2:2:536
    for col= 2:2:704
        cb(row,col) = 0;
        cr(row,col) = 0;
    end
end
cbsub = cb(1:2:536,1:2:704);
crsub = cr(1:2:536,1:2:704);
figure(9);
subplot(1,2,1); subimage(cbsub); title('Cb band subsampled');
subplot(1,2,2); subimage(crsub); title('Cr band subsampled');

%6.1
linecb = cb;
linecr = cr;     
for row = 1:535
    for col = 1:703
        if (linecb(row, col) == 0)
            if(mod(row,2)==0)
                linecb(row, col) = 0.5*(linecb(row-1, col)) + 0.5*(linecb(row+1, col)); 
                linecr(row, col) = 0.5*(linecr(row-1, col)) + 0.5*(linecr(row+1, col)); 
            end
        end   
    end
end

for row = 1:536
    for col = 1:703
        if (linecb(row, col) == 0) 
            if(mod(col,2)==0)   
                linecb(row, col) = 0.5*(linecb(row, col-1)) + 0.5*(linecb(row, col+1)); 
                linecr(row, col) = 0.5*(linecr(row, col-1)) + 0.5*(linecr(row, col+1)); 
            end
        end
        
    end
end

figure(10);
subplot(1,2,1); subimage(linecb); title('Cb band linear interpolation');
subplot(1,2,2); subimage(linecr); title('Cr band linear interpolation');

%6.2
cbrc = cb;
crrc = cr;
for row = 1:536
    for col = 1:704
        if (cbrc(row, col) == 0)
            if(mod(row,2) == 0) 
                cbrc(row, col) = cbrc(row-1, col);
                crrc(row, col) = crrc(row-1, col);
            end
        end
    end
end

for row = 1:536
    for col = 1:704
        if (cbrc(row, col) == 0) 
            if(mod(row,2) ~= 0) 
                cbrc(row, col) = cbrc(row, col-1);
                crrc(row, col) = crrc(row, col-1);
            end
        end
    end
end

figure(11);
subplot(1,2,1); subimage(cbrc); title('Cb band row/col replication');
subplot(1,2,2); subimage(crrc); title('Cr band row/col replication');

%7
lineRegRBG = ycbcr2rgb(cat(3, ycbcr(:,:,1), linecb, linecr));
rowColRepRBG = ycbcr2rgb(cat(3, ycbcr(:,:,1), cbrc, crrc));

%8
figure(12);
subplot(1,2,1); subimage(lineRegRBG); title('RGB reconstructed from linear interpolation');
subplot(1,2,2); subimage(I); title('Original RGB');
figure(13);
subplot(1,2,1); subimage(rowColRepRBG); title('RGB reconstructed from row/col replication');
subplot(1,2,2); subimage(I); title('Original RGB');

%9
%{
The image from linear interpolation is of better quality than the one from
row/column replication. This is because in row/column replication you are
basically extending pixels from the previous row/column. This gives the
image more of a rougher look. In linear interpolation, the averages of
pixels around a missing pixel is used to fill in that missing pixel.
Because of that, the missing pixel is filled with a more smoother pixel and
so giving the image a more smoother look.
%}

%10
MSE = (sum(sum((I-lineRegRBG).^2)))/(536*704);
fprintf('Red MSE: %f\n', MSE(:,:,1));
fprintf('Green MSE: %f\n', MSE(:,:,2));
fprintf('Blue MSE: %f\n', MSE(:,:,3));
%{
Red MSE: 3.456507
Green MSE: 1.499754
Blue MSE: 6.014912
Blue has the maximum amount of distortion.
Green has the least amout of distortion.
Red has a modest amount of distortion, being in between the other two
colors.
%}

%11
fprintf('Original Image: %f\n', (size(I,1)*size(I,2)*3));
fprintf('Subsampled Image: %f\n', (size(ycbcr(:,:,1),1)*size(ycbcr(:,:,1),2) + size(cbsub,1)*size(cbsub,2) + size(crsub,1)*size(crsub,2)));
fprintf('Final Compression: %f\n', (size(I,1)*size(I,2)*3)/(size(ycbcr(:,:,1),1)*size(ycbcr(:,:,1),2) + size(cbsub,1)*size(cbsub,2) + size(crsub,1)*size(crsub,2)));
%{
Final Compression was 2.000000, so the original image was reduced to half
its size when it was being subsampled.
%}

