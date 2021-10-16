%Ethan Nagelvoort
%821234668
%HW2

%Encoder
I = imread('Flooded_house.jpg','jpg');
og = I;
I = rgb2ycbcr(I);
y = I(:,:,1);
cb = I(:,:,2);
cr = I(:,:,3);
for row = 2:2:536
    for col= 2:2:704
        cb(row,col) = 0;
        cr(row,col) = 0;
    end
end
cbsub = cb(1:2:536,1:2:704);
crsub = cr(1:2:536,1:2:704);
r = size(cbsub,1);
c = size(cbsub,2);

%a
y8 = y;
cb8 = cbsub;
cr8 = crsub;
pDCT = @dct2 ; 
y_DCT = blkproc (y8,[8 8],pDCT);
cb_DCT = blkproc (cb8,[8 8],pDCT);
cr_DCT = blkproc (cr8,[8 8],pDCT);
block1= y_DCT(41:48,1:8); %8*6 = 48, hence 6th row is 41:48
block2= y_DCT(41:48,9:16);
figure(1);
subplot(1,2,1); subimage(block1); title('First Block');
subplot(1,2,2); subimage(block2); title('Second Block');
fprintf("The DCT coefficient matrix of block 1:\n");
disp(block1);
fprintf("The DCT coefficient matrix of block 2:\n");
disp(block2);
my_DCT = fix(y_DCT);
mycb_DCT = fix(cb_DCT);
mycr_DCT = fix(cr_DCT);
figure(2);
subplot(1,3,1); subimage(my_DCT); title('DCT of Y');
subplot(1,3,2); subimage(cb_DCT); title('DCT of cb');
subplot(1,3,3); subimage(cr_DCT); title('DCT of cr');

%b
lum_quan = [16 11 10 16 24 40 51 61;
12 12 14 19 26 58 60 55;
14 13 16 24 40 57 69 56;
14 17 22 29 51 87 89 62;
18 22 37 56 68 109 103 77;
24 35 55 64 81 104 113 92;
49 64 78 87 108 121 120 101;
72 92 95 98 112 100 103 99];
chrom_quan = [17 18 24 47 99 99 99 99;
18 21 26 66 99 99 99 99;
24 26 56 99 99 99 99 99;
47 66 99 99 99 99 99 99;
99 99 99 99 99 99 99 99;
99 99 99 99 99 99 99 99;
99 99 99 99 99 99 99 99;
99 99 99 99 99 99 99 99];
lumin = @(y_DCT)round(y_DCT./lum_quan);
y_lumin = blkproc(y_DCT,[8 8],lumin);
block1  = y_lumin(41:48,1:8);
block2  = y_lumin(41:48,9:16);
figure(3);
subplot(2,2,1), subimage(block1), title("Quantized block 1.");
subplot(2,2,2), subimage(block2), title("Quantized block 2.");
cbq = @(cb_DCT)round(cb_DCT./chrom_quan);
crq = @(cr_DCT)round(cr_DCT./chrom_quan);
cb_quan = blkproc(cb_DCT,[8 8],cbq);
cr_quan = blkproc(cr_DCT,[8 8],crq);
blockcb1  = cb_quan(41:48,1:8);
blockcb2  = cb_quan(41:48,9:16);
blockcr1  = cr_quan(41:48,1:8);
blockcr2  = cr_quan(41:48,9:16);
%DC DCT coefficients
fprintf("The DC DCT coefficient of quantized block 1:\n");
fprintf("%d\n", block1(1,1));
fprintf("The DC DCT coefficient of quantized block 2:\n");
fprintf("%d\n", block2(1,1));
%ZIGZAG AC DCT coefficients
%zigzag through block1
currentPosition=1;
AC_DCT1 = [];
for col=1:15
    if col>8
        if mod(col,2)~=0
            k=mod(col,8);
            i=8;
            for currentCol=k+1:8
                AC_DCT1(currentPosition)=block1(i,currentCol);
                i=i-1;
                currentPosition=currentPosition+1;
            end
        else
           k=mod(col,8);
            currentCol=8;
            for i=k+1:8
                AC_DCT1(currentPosition)=block1(i,currentCol);
                currentPosition=currentPosition+1;
                currentCol=currentCol-1;
            end
        end
    else
        if mod(col,2)~=0
             i=col;
            for currentCol=1:col
                AC_DCT1(currentPosition)= block1(i,currentCol);
                i= i-1;
                currentPosition=currentPosition+1;
            end
        else
           currentCol=col;
            for i=1:col
                AC_DCT1(currentPosition)= block1(i,currentCol);
                currentPosition=currentPosition+1;
                currentCol=currentCol-1;
            end
        end
    end
end
%zigzag through block2
currentPosition=1;
AC_DCT2 = [];
for col=1:15
    if col>8
         if mod(col,2)~=0
           k=mod(col,8);i=8;
            for currentCol=k+1:8
                AC_DCT2(currentPosition)=block2(i,currentCol);
                i=i-1;
                currentPosition=currentPosition+1;
            end
        else
            currentCol=8;
            k=mod(col,8); 
            currentCol=8;
            for i=k+1:8
                AC_DCT2(currentPosition)=block2(i,currentCol);
                currentPosition=currentPosition+1;
                currentCol=currentCol-1;
            end
         end
    else
         if mod(col,2)~=0
            i=col;
            for currentCol=1:col
                AC_DCT2(currentPosition)= block2(i,currentCol);
                i=i-1;
                currentPosition=currentPosition+1;
            end
        else
           currentCol=col;
            for i=1:col
                AC_DCT2(currentPosition)= block2(i,currentCol);
                currentPosition=currentPosition+1;
                currentCol=currentCol-1;
            end
        end
    end
end
fprintf("The AC DCT coefficients of block1:\n");
display(AC_DCT1(2:end));
fprintf("The AC DCT coefficients of block2:\n");
display(AC_DCT2(2:end));

%Decoder

%c
inverse_Y = @(y_lumin)round(y_lumin.*lum_quan);
inverse_cb = @(cb_quan)round(cb_quan.*chrom_quan);
inverse_cr = @(cr_quan)round(cr_quan.*chrom_quan);
block1 = blkproc(y_lumin,[8 8],inverse_Y);
block2 = blkproc(cb_quan,[8 8],inverse_cb);
block3 = blkproc(cr_quan,[8 8],inverse_cr);
figure(4);
subplot(1,3,1); subimage(block1); title('Inverse Quantized of Y');
subplot(1,3,2); subimage(block2); title('Inverse Quantized of Cb');
subplot(1,3,3); subimage(block3); title('Inverse Quantized of Cr');

%d
pDCT = @idct2;
y_inv = blkproc(block1, [8 8], pDCT);
cb_inv = blkproc(block2, [8 8], pDCT);
cr_inv = blkproc(block3, [8 8], pDCT);
y_inv = fix(y_inv);
cb_inv = fix(cb_inv);
cr_inv = fix(cr_inv);
y_inv = uint8(y_inv);
rec_cb = uint8(zeros( 536, 704));
rec_cr = uint8(zeros( 536, 704));
rec_cb(1:2:end, 1:2:end) = cb_inv(1:end,1:end);
rec_cr(1:2:end, 1:2:end) = cr_inv(1:end,1:end);
for row = 2:2: 535 
    for col = 1:2:704
            rec_cb(row, col) = rec_cb(row-1, col)/2 + rec_cb(row+1, col)/2;
            rec_cr(row, col) = rec_cr(row-1, col)/2 + rec_cr(row+1, col)/2;
    end
end
for row = 1: 536
    for col = 2:2:703
            rec_cb(row, col) = rec_cb(row, col-1)/2 + rec_cb(row, col+1)/2;
            rec_cr(row, col) = rec_cr(row, col-1)/2 + rec_cr(row, col+1)/2;
    end
end
for row = 1:536
     for col = 704
         rec_cb(row, col) = (rec_cb(row, col-1)); 
         rec_cr(row, col) = (rec_cr(row, col-1)); 
     end 
 end
 for row = 536
     for col = 1:704
         rec_cb(row, col) = (rec_cb(row-1, col)); 
         rec_cr(row, col) = (rec_cr(row-1, col)); 
     end 
 end
reconst_ycbcr = cat(3,y_inv,rec_cb,rec_cr);
reconst_rgb = ycbcr2rgb(reconst_ycbcr);
figure(5)
subplot(1,2,1); subimage(og); title('Original RBG image');
subplot(1,2,2); subimage(reconst_rgb); title('Reconstructed image');

%error
error = y - reconst_ycbcr(:,:,1);
figure(6);
subplot(1,1,1); subimage(error); title('Error image');

%psnr
MSE = 1/(536*704)*(sum(sum((og-reconst_rgb).^2)));
PSNR = 10*log10(255.^2/MSE);
fprintf("The PSNR value of the luminance component is:\n");
display(PSNR(:,:,1));
figure(7)
imshow(reconst_rgb)





