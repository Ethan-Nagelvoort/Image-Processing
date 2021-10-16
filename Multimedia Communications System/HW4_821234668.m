%Ethan Nagelvoort
%821234668
%HW4

walk = VideoReader('walk_qcif.avi');
o = read(walk);
for i = 1:37
    f(:,:,:,i) = rgb2ycbcr(o(:,:,:, i));
end
vf=f(:,:,:,6:10);
y = f(:,:,1,6:10);
cb = f(1:2:144,1:2:176,2,6:10);
cr = f(1:2:144,1:2:176,3,6:10);
yac = zeros(5,(144/8)*(176/8)*(8*8-2));
ydc = zeros(5,(144/8)*(176/8));
rac = 1;
f1 = vf(:,:,:,1);
[ydct,cbdct,crdct,yq,cbq,crq] = dctquan(y(:,:,1),cb(:,:,1),cr(:,:,1));
counter = 1;
[yiq,cbiq,criq,yidct,cbidct,cridct] = inversedct(yq,cbq,crq);
[upsample] = row_col(yidct,cbidct,cridct);
decoded(:,:,:,counter) = ycbcr2rgb(upsample);
counter = counter+1;
[ydc,yac,dc,ac] = acdcblock(yq,ydc,yac,rac);
figure();
subplot(2,1,1), subimage(o(:,:,:,6)),title(['Orginal RGB Frame for I frame 6 ']);
subplot(2,1,2), subimage(decoded(:,:,:,1)),title(['reconstruction of frame 6 (I-frame)']);
rac = rac+1;
[res, recon] = motion_vect(vf);
for i = 2:5
    [ydct,cbdct,crdct,yq,cbq,crq] = dctquan(uint8(res(:,:,i)),0,0);
    [ydc,yac,dc,ac] = acdcblock(yq,ydc,yac,rac);
    [yiq,cbiq,criq,yidct,cbidct,cridct] = inversedct(yq,0,0);
    [ydct,cbdct,crdct,yq,cbq,crq] = dctquan(uint8(recon(:,:,i)),cb(:,:,i),cr(:,:,i));
    [yiq,cbiq,criq,yidct,cbidct,cridct] = inversedct(yq,cbq,crq);
    [upsample] = row_col(yidct,cbidct,cridct);
    decoded(:,:,:,counter) = ycbcr2rgb(upsample);
    counter = counter + 1;
    figure();
    subplot(3,1,1), subimage(o(:,:,:,i+5)),title(['Orginal RGB Frame for P frame ',num2str(i+5)]);
    subplot(3,1,2), subimage(uint8(res(:,:,i))),title(['Difference for P frame ',num2str(i+5)]);
    subplot(3,1,3), subimage(decoded(:,:,:,i)),title(['Reconstruction for P frame ',num2str(i+5)]);
    rac = rac+1;
end

function [ydct,cbdct,crdct,yq,cbq,crq] = dctquan(y,cbsub,crsub)
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
pdct = @dct2;
ydct = blkproc(y,[8 8],pdct);
lumin = @(ydct)round(ydct./lum_quan);
yq = blkproc(ydct,[8 8],lumin);
cbdct = blkproc(cbsub,[8 8],pdct);
crdct = blkproc(crsub,[8 8],pdct);
cbq = blkproc(cbdct,[8 8],@(cbdct)round(cbdct./chrom_quan));
crq = blkproc(crdct,[8 8],@(crdct)round(crdct./chrom_quan));
end

function [yiq,cbiq,criq,yidct,cbidct,cridct] = inversedct(yq,cbq,crq)
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
yiq = blkproc(yq,[8 8],@(yq)round(yq.*lum_quan));
cbiq = blkproc(cbq,[8 8],@(cbq)round(cbq.*chrom_quan));
criq = blkproc(crq,[8 8],@(crq)round(crq.*chrom_quan));
idct = @(block_struct)idct2(block_struct.data);
yidct  = uint8(blockproc(yiq,[8 8],idct));
cbidct = uint8(blockproc(cbiq,[8 8],idct));
cridct = uint8(blockproc(criq,[8 8],idct));
end

function [upsample] = row_col(y,ucb,ucr)

upsample = zeros(144,176,3);
upsample(1:2:144,1:2:176,2) = ucb(:,:);
upsample(1:2:144,1:2:176,3) = ucr(:,:);
rec_cb = uint8(zeros( 144, 176));
rec_cr = uint8(zeros( 144, 176));
rec_cb(1:2:end, 1:2:end) = ucb(1:end,1:end);
rec_cr(1:2:end, 1:2:end) = ucr(1:end,1:end);
for row = 2:2: 142 
        for col = 1:2:176
                rec_cb(row, col) = rec_cb(row-1, col)/2 + rec_cb(row+1, col)/2;
                rec_cr(row, col) = rec_cr(row-1, col)/2 + rec_cr(row+1, col)/2;
        end
end
    for row = 1: 143
        for col = 2:2:175
                rec_cb(row, col) = rec_cb(row, col-1)/2 + rec_cb(row, col+1)/2;
                rec_cr(row, col) = rec_cr(row, col-1)/2 + rec_cr(row, col+1)/2;
        end
    end
    for row = 1:144
         for col = 176
             rec_cb(row, col) = (rec_cb(row, col-1)); 
             rec_cr(row, col) = (rec_cr(row, col-1)); 
         end 
     end
     for row = 144
         for col = 1:176
             rec_cb(row, col) = (rec_cb(row-1, col)); 
             rec_cr(row, col) = (rec_cr(row-1, col)); 
         end 
     end
upsample = cat(3,y,rec_cb,rec_cr);
end


function [res, recon]=motion_vect(vf)
dx = zeros(9,11,2,5);
dy = zeros(9,11,2,5);
SAD = zeros(5,1);
res = int16(zeros(144,176,5));
recon= int16(zeros(144,176,5));
recon(:,:,1)  = vf(:,:,1,1);
for fnum=2:5
    b = vf(:,:,:,fnum);
    bb = vf(:,:,:,fnum-1);
    rf = bb(:,:,1);
    cf = b(:,:,1);
    rc = 1;
    cc = 1;
    numOfComp = 0;
    numOfAdd_Sub = 0;
    for rBLK=1:16:144
        for cBLK=1:16:176
            if ((rBLK - 8) > 0)
                tSide = rBLK - 8;
            else
                tSide = 1;
            end
            if ((cBLK - 8) > 0)
                lSide = cBLK - 8;
            else
                lSide = 1;
            end
            if ((rBLK + 24) < 145)
                bSide = rBLK + 24;
            else
                bSide = 144;
            end
            if ((cBLK + 24) < 177)
                rSide = cBLK + 24;
            else
                rSide = 176;
            end
            dx(rc,cc,1,fnum) = cBLK;
            dy(rc,cc,1,fnum) = rBLK;
            crntBLK = cf(rBLK:rBLK+15, cBLK:cBLK+15);
            ref = rf(tSide:tSide+15, lSide:lSide+15);
            min = 100000000000000000000000000000000000000000000;
            numOfAdd_Sub = numOfAdd_Sub + 2;
            cv = lSide;
            rv = tSide;
            for i = tSide:bSide-15    
                for ii = lSide:rSide-15
                    ref = rf(i:15+i, ii:ii+15);
                    sad = sum(abs(int16(crntBLK(:))- int16(ref(:))));
                    numOfComp = numOfComp + 1;
                    numOfAdd_Sub = numOfAdd_Sub + 2;
                    if(sad < min)
                        min = sad;
                        cv = ii;
                        rv = i;
                    end
                end
            end 
            SAD(fnum) = SAD(fnum) + min;
            dx(rc,cc,2,fnum) = cv;
            dy(rc,cc,2,fnum) = rv;
            res(rBLK:rBLK + 15,cBLK:cBLK + 15,fnum) = int16(cf(rBLK:rBLK+15, cBLK:cBLK+15)) - int16(rf(rv:rv + 15, cv:cv + 15));
            recon(rBLK:rBLK+15, cBLK:cBLK + 15, fnum) = recon(rv:rv + 15, cv:cv+15, fnum-1)+ res(rBLK:rBLK + 15, cBLK:cBLK + 15, fnum);
            cc = cc+1;
        end 
        rc = rc+1;
        cc = 1;
    end 
    figure();
    quiver(dx(:,:,1,fnum),dy(:,:,1,fnum),(dx(:,:,2,fnum)) - (dx(:,:,1,fnum)),(dy(:,:,2,fnum) - (dy(:,:,1,fnum))));
    axis([0 176 0 144]);
    title(['Motion Vector for Frames ',num2str(fnum+4),' and ',num2str(fnum+5)]);
end
end

function [ydc,yac,dc,ac] = acdcblock(f,ydc,yac,r)
c = zeros(24948);
coeff = size(c,1);
ac = 1;
dc = 1;
if (r <= 5)
    block = zeros(8,8);
    for row = 1:8:144
        for col = 1:8:176
            if ((row+7) <= 144) 
                if ((col+7) <= 176)
                    block = f(row:row+7, col:col+7);
                    currentPosition=1;
                    acdct = [];
                    for col=1:15
                        if col>8
                             if mod(col,2)~=0
                               k=mod(col,8);i=8;
                                for currentCol=k+1:8
                                    acdct(currentPosition)=block(i,currentCol);
                                    i=i-1;
                                    currentPosition=currentPosition+1;
                                end
                            else
                                currentCol=8;
                                k=mod(col,8); 
                                currentCol=8;
                                for i=k+1:8
                                    acdct(currentPosition)=block(i,currentCol);
                                    currentPosition=currentPosition+1;
                                    currentCol=currentCol-1;
                                end
                             end
                        else
                             if mod(col,2)~=0
                                i=col;
                                for currentCol=1:col
                                    acdct(currentPosition)= block(i,currentCol);
                                    i=i-1;
                                    currentPosition=currentPosition+1;
                                end
                            else
                               currentCol=col;
                                for i=1:col
                                    acdct(currentPosition)= block(i,currentCol);
                                    currentPosition=currentPosition+1;
                                    currentCol=currentCol-1;
                                end
                            end
                        end
                    end
                    yac(r,ac:ac+(8*8-2)) = acdct(2:end);
                    ydc(r,dc:dc+(8-1)) = acdct(1);
                    ac = ac + (8*8-1);
                    dc = dc + 1;
                end
            end
        end
    end
end
end