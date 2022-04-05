function dilate!( mask::Array{T,N}, dil::Array{Bool,N}, rad; kern=nothing, th=nothing ) where {T,N}

    if kern == nothing
        size = ones( Int64, N ) .* rad; 
        kern = ones( UInt8, size... )
    end

    conv = FFTConvolution_crop( mask, kern, typ=Float32 );

    th = ( th == nothing ) ? 1 - 0.1 : th;

    @inbounds @simd for e in 1:length(conv)
        dil[e] = conv[e] > th;
    end

    return dil
end

function dilate( mask::Array{T,N}, rad; kern=nothing, th=nothing ) where {T,N}

    return dilate!( mask, zeros(Bool,size(mask)), rad, kern=kern, th=th ); 
end

function erode!( mask::Array{T,N}, ero::Array{Bool,N}, rad; kern=nothing, th=nothing ) where {T,N}

    if kern == nothing
        kern = ones( UInt8, (ones( Int64, N ) .*  rad)... ) # ones( rad, rad ) or ones( rad, rad, rad )
    end

    th = ( th == nothing ) ? sum( kern ) - 0.1 : th;

    conv = FFTConvolution_crop( mask, kern, typ=Float32 );

    for e in 1:length(conv)
        ero[e] = ( conv[e] > th )
    end

    return ero

end

function erode( mask::Array{T,N}, rad; kern=nothing, th=nothing ) where {T,N}

    return erode!( mask, zeros(Bool,size(mask)), rad, kern=kern, th=th ); 
end

function closing( mask, size; kern=nothing )
    s1 = dilate( mask, size, kern=kern  )
    s2 = erode( s1, size, kern=kern )
    return s2
end

function opening( mask, size; kern=nothing )
    s1 = erode( mask, size, kern=kern )
    s2 = dilate( s1, size, kern=kern )
    return s2
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
