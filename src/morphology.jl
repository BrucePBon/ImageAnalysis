function dilate( mask, size; kern=nothing )

    if kern == nothing
        kern = ones( UInt8, size, size )
    end

    mask_ = convolution( UInt8.(mask), kern );

    for e in 1:length(mask_)
        mask_[e] = mask_[e] > 1 - 0.1
    end

    return mask_
end

function erode( mask, size; kern=nothing)

    if kern == nothing
        kern = ones( UInt8, size, size )
    end

    th = sum( kern );

    mask_ = convolution( UInt8.( mask ), kern );

    for e in 1:length(mask_)
        mask_[e] = ( mask_[e] > th - 0.1 )
    end

    return mask_
end

function distanceTransform( mask )

    distances = zeros( UInt16, size(mask) )

    kern_size = 3;
    keepgoing = true
    while keepgoing

        kern  = ones( UInt8, kern_size, kern_size );
        th    = sum( kern );

        mask_ = convolution( UInt8.( mask ), kern );
        h, w  = size( mask_ );
        off   = div( kern_size, 2 );

        trues = 0;
        for col in off:w-off
            for row in off:h-off
                if mask_[row,col] >= length(kern) - 0.1
                    distances[ row-off+1, col-off+1 ] += 1;
                    trues += 1;
                end
            end
        end

        if trues == 0
            keepgoing = false;
        end

        kern_size += 2;
    end

    println( "maximum kern_size: ", kern_size )

    return distances
end
