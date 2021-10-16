walk = VideoReader('walk_qcif.avi');
f = read(walk);
for i = 1:37
    f(:,:,:,i) = rgb2ycbcr(f(:,:,:, i));
end
vf=f(:,:,:,6:10);
dx = zeros(9,11,2,5);
dy = zeros(9,11,2,5);
SAD = zeros(5,1);
res = int16(zeros(144,176,5));
recon= int16(zeros(144,176,5));
recon(:,:,1)  = vf(:,:,1,1);
for fnum=2:5
    b = vf(:,:,:,fnum);
    bb = vf(:,:,:,fnum-1);
    cb = b(:,:,2);
    cr = b(:,:,3);
    for row = 2:2:144
        for col= 2:2:176
            cb(row,col) = 0;
            cr(row,col) = 0;
        end
    end
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
    for row = 2:2: 142 
        for col = 1:2:176
                cb(row, col) = cb(row-1, col)/2 + cb(row+1, col)/2;
                cr(row, col) = cr(row-1, col)/2 + cr(row+1, col)/2;
        end
    end
    for row = 1: 143
        for col = 2:2:175
                cb(row, col) = cb(row, col-1)/2 + cb(row, col+1)/2;
                cr(row, col) = cr(row, col-1)/2 + cr(row, col+1)/2;
        end
    end
    for row = 1:144
         for col = 176
             cb(row, col) = (cb(row, col-1)); 
             cr(row, col) = (cr(row, col-1)); 
         end 
     end
     for row = 144
         for col = 1:176
             cb(row, col) = (cb(row-1, col)); 
             cr(row, col) = (cr(row-1, col)); 
         end 
     end
    figure();
    quiver(dx(:,:,1,fnum),dy(:,:,1,fnum),(dx(:,:,2,fnum)) - (dx(:,:,1,fnum)),(dy(:,:,2,fnum) - (dy(:,:,1,fnum))));
    axis([0 176 0 144]);
    title(['Motion Vector for Frames ',num2str(fnum+4),' and ',num2str(fnum+5)]);
    figure();
    subplot(2,2,1),subimage(cf),title(['Original Y-component for P Frame ',num2str(fnum+5)]);
    subplot(2,2,2),subimage(uint8(recon(:,:,fnum))),title(['Reconstructed Y-component for P Frame ',num2str(fnum+5)])
    subplot(2,2,3),subimage(uint8(abs(res(:,:,fnum)))),title(['Error Image for P Frame ', num2str(fnum+5)])
    figure();
    subplot(1,2,1),subimage(ycbcr2rgb(vf(:,:,:,fnum)));title(['Original RGB Image for P Frame ',num2str(fnum+5)]);
    subplot(1,2,2),subimage(ycbcr2rgb(cat(3,uint8(recon(:,:,fnum)),uint8(cb),uint8(cr)))),title(['Reconstructed RGB Image for P Frame ',num2str(fnum+5)]);
    fprintf("Total number of addition/subtractions for frame %d: %d\nNumber of addition/subtractions per block in frame %d: %d\n", fnum+5,numOfAdd_Sub, fnum+5,round(numOfAdd_Sub/99));
    fprintf("Total number of comparasions for frame %d: %d\nNumber of comparasions per block for frame %d: %d\n",fnum+5,numOfComp,fnum+5, round(numOfComp/99));
end 