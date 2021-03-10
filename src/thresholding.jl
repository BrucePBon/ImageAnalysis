#global








#adaptative
function niblack( image::Array{T,2}, width::Integer; k = 0.1 ) where {T<:Real}

    intA  = integralArray( image, type=Float64 );
    intA2 = integralArray( image, type=Float64, fun=(x)->(x^2));

    imgc  = zeros( Float32, size(image) );

    h, w = size(image);
    npix = width*width;

    for col in 1:w
        for row in 1:h

            TL = max( 1, row - width ), max( 1, col - width )
            BR = min( h, row + width ), min( w, col + width )

            npix2 = (length(TL[1]:BR[1])+1)*(length(TL[2]:BR[2])+1);

            sumI  = integralArea(  intA, TL, BR );
            sumI2 = integralArea( intA2, TL, BR );

            mean = sumI/npix2;
            m2   =  mean*mean;

            std  = sqrt( (sumI2 + npix2*m2 - 2*sumI*mean)/npix2 );

            th = mean + k*std;

            imgc[row,col] = ( image[row,col] <= th ) ? 1 : 0;
        end
    end

    return imgc;
end

function sauvola( image::Array{T,2}, width::Integer; R=1, k=0.3 ) where {T<:Real}

    intA  = integralArray( image, type=Float64 );
    intA2 = integralArray( image, type=Float64, fun=(x)->(x^2));

    imgc  = zeros( Float32, size(image) );

    h, w = size(image);
    npix = width*width;

    for col in 1:w
        for row in 1:h

            TL = max( 1, row - width ), max( 1, col - width )
            BR = min( h, row + width ), min( w, col + width )

            npix2 = (length(TL[1]:BR[1])+1)*(length(TL[2]:BR[2])+1);

            sumI  = integralArea(  intA, TL, BR );
            sumI2 = integralArea( intA2, TL, BR );

            mean = sumI/npix2;
            m2   =  mean*mean;

            std  = sqrt( (sumI2 + npix2*m2 - 2*sumI*mean)/npix2 );

            th = mean + k*( std/R - 1 )*(mean-1);

            imgc[row,col] = ( image[row,col] <= th ) ? 1 : 0;
            #imgc[row,col] = th
        end
    end

    return imgc;
end

heavensideStep( num ) = ( num <= 0 ) ? 0 : 1;

function adaptOtsu( image::Array{T,2}, width::Integer, step::Integer, bindata::AbstractRange; R=0.1 ) where {T<:Real}

    imgc = zeros( Float32, size(image) )
    img  = copy(image)
    h, w = size(image)

    len  = 2*width + 1;
    #step = 2*width + 1;
    centers_r = [ row - width for row in len:step:h ]
    centers_c = [ col - width for col in len:step:w ]

    thresholds = zeros( Float32, length(len:step:h), length(len:step:w) );

    nhImg,_,_ = normhistogram( imgc, bindata );
    imgOt = otsu( nhImg );

    # computing threshold for each grid element
    contc = 1;
    for col in len:step:w
        contr = 1;
        for row in len:step:h

            TL = max( 1, row - width ), max( 1, col - width )
            BR = min( h, row + width ), min( w, col + width )

            TL = ( row - len + 1, col - len + 1 )
            BR = ( row, col )

            subimg = image[ TL[1]:BR[1], TL[2]:BR[2] ]
            nh,_,_ = normhistogram( subimg, bindata )
            subOt  = otsu( nh )

            th = heavensideStep( R/abs(imgOt - subOt) - 1 )*150 + subOt

            thresholds[contr,contc] = th;
            contr += 1;
        end
        contc += 1;
    end

    ch = length( centers_r );
    cw = length( centers_c );

    # Interpolation
    dist = step + 1;

    # centers
    for idxc in 1:length(centers_c)
        col = centers_c[idxc]
        for idxr in 1:length(centers_r)
            row = centers_r[idxr]
            #imgc[row,col] = thresholds[idxr,idxc]
            imgc[row,col] = img[row,col] > thresholds[idxr,idxc]
        end
    end

    # vertically in line
    for idxc in 1:length(centers_c)
        col = centers_c[idxc];
        for idxr in 1:length(centers_r)-1
            U   = centers_r[ idxr ];
            B   = centers_r[idxr+1];
            thU = thresholds[ idxr , idxc ]
            thB = thresholds[idxr+1, idxc ]
            for row in U+1:B-1
                thI = ( B - row )/dist * thU + ( row - U )/dist * thB
                #imgc[ row, col ] = ( B - row )/dist * thU + ( row - U )/dist * thB
                imgc[ row, col ] = img[ row, col ] > thI
            end
        end
    end

    # laterally in line
    for idxc in 1:length(centers_c)-1
        L = centers_c[idxc]
        R = centers_c[idxc+1];
        for col in L+1:R-1
            fL = ( R - col )/dist;
            fR = ( col - L )/dist;
            for idxr in 1:length(centers_r)
                row = centers_r[idxr]
                thL = thresholds[idxr,  idxc ]
                thR = thresholds[idxr, idxc+1]
                #imgc[ row, col ] = fL * thL + fR * thR;
                imgc[ row, col ] = img[ row, col ] > fL * thL + fR * thR
            end
        end
    end

    # middle regions
    for idxc in 1:length(centers_c)-1
        L = centers_c[idxc]
        R = centers_c[idxc+1]
        for idxr in 1:length(centers_r)-1
            U = centers_r[idxr]
            B = centers_r[idxr+1]
            thUL = thresholds[  idxr ,  idxc  ]
            thUR = thresholds[  idxr , idxc+1 ]
            thBL = thresholds[ idxr+1,  idxc  ]
            thBR = thresholds[ idxr+1, idxc+1 ]

            for col in L+1:R-1
                fR  = ( col - L )/dist;
                fL  = ( R - col )/dist;

                for row in U+1:B-1
                    fU = ( B - row )/dist
                    fB = ( row - U )/dist

                    itpL  = fB * thBL + fU * thUL
                    itpR  = fB * thBR + fU * thUR
                    itpTH = fR * itpR + fL * itpL
                    # imgc[row,col] = itpTH
                    imgc[ row, col ] = img[ row, col ] > itpTH
                end
            end
        end
    end

    # first row
    fr = 1:centers_r[1]-1
    for col in 1:centers_c[1]
        #imgc[fr,col] .= thresholds[1,1]
        imgc[fr,col] .= img[fr,col] .> thresholds[1,1]
    end
    for col in centers_c[end]:w
        #imgc[fr,col] .= thresholds[1,end]
        imgc[fr,col] .= img[fr,col] .> thresholds[1,end]
    end
    for cidx in 1:length(centers_c)-1
        L   = centers_c[cidx]
        R   = centers_c[cidx+1]
        thL = thresholds[1,cidx]
        thR = thresholds[1,cidx+1]

        #imgc[fr,L] .= thL
        imgc[fr,L] .= img[fr,L] .> thL

        for col in L+1:R-1
            fR = ( col - L )/dist;
            fL = ( R - col )/dist;
            #imgc[fr,col] .= fR * thR + fL * thL
            imgc[fr,col] .= img[fr,col] .> ( fR * thR + fL * thL )
        end
    end

    # first col
    fc = 1:centers_c[1]-1
    for row in 1:centers_r[1]
        #imgc[row,fc] .= thresholds[1,1]
        imgc[row,fc] .= img[row,fc] .> thresholds[1,1]
    end
    for row in centers_r[end]:h
        #imgc[row,fc] .= thresholds[end,1]
        imgc[row,fc] .= img[row,fc] .> thresholds[end,1]
    end
    for ridx in 1:length(centers_r)-1
        U = centers_r[ridx]
        B = centers_r[ridx+1]

        thU = thresholds[ ridx , 1]
        thB = thresholds[ridx+1, 1]

        #imgc[U,fc] .= thU
        imgc[U,fc] .= img[U,fc] .> thU

        for row in U+1:B-1
            fB = ( row - U )/dist;
            fU = ( B - row )/dist;
            #imgc[row,fc] .= fU * thU + fB * thB
            imgc[row,fc] .= img[row,fc] .> fU * thU + fB * thB
        end
    end

    # last row
    lr = centers_c[end]:w
    for col in 1:centers_c[1]
        #imgc[lr,col] .= thresholds[end,1]
        imgc[lr,col] .= img[lr,col] .> thresholds[end,1]
    end
    for col in centers_c[end]:w
        #imgc[lr,col] .= thresholds[end,end]
        imgc[lr,col] .= img[lr,col] .> thresholds[end,end]
    end
    for cidx in 1:length(centers_c)-1
        L = centers_c[cidx]
        R = centers_c[cidx+1]

        thL = thresholds[end,cidx]
        thR = thresholds[end,cidx+1]

        #imgc[lr,L] .= thL
        imgc[lr,L] .= img[lr,L] .> thL

        for col in L+1:R-1
            fR = ( col - L )/dist;
            fL = ( R - col )/dist;
            #imgc[ lr,col] .= fR * thR + fL * thL
            imgc[ lr, col ] .= img[ lr, col ] .> ( fR * thR + fL * thL )
        end
    end

    # last col
    lc = centers_r[end]:h

    for row in 1:centers_r[1]
        #imgc[row,lc] .= thresholds[1,end]
        imgc[row,lc] .= img[row,lc] .> thresholds[1,end]
    end
    for row in centers_r[end]:h
        #imgc[row,lc] .= thresholds[end,end]
        imgc[row,lc] .= img[row,lc] .> thresholds[end,end]
    end
    for ridx in 1:length(centers_r)-1
        U = centers_r[ridx]
        B = centers_r[ridx+1]

        thU = thresholds[ ridx , end]
        thB = thresholds[ridx+1, end]

        #imgc[U,lc] .= thU
        imgc[U,lc] .= img[U,lc] .> thU

        for row in U+1:B-1
            fB = ( row - U )/dist;
            fU = ( B - row )/dist;
            #imgc[row,lc] .= fU * thU + fB * thB
            imgc[row,lc] .= img[row,lc] .> ( fU * thU + fB * thB )
        end
    end

    return imgc, thresholds
end
