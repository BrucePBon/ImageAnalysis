function meanstd( img ) 

	T  = ( eltype(img) <: Integer ) ? Int64 : Float32; 

	I  = T(0)
	I2 = T(0)
	N  = length(img); 
	@inbounds @simd for idx in 1:N
		I  += img[idx]
		I2 += img[idx]*img[idx]
	end
	# mean = sumI / N 
	# std  : sum( ( I - mean )^2 ) = sum( I^2 + mean^2 - 2*I*mean ) = sumI2 + Nmean^2 - 2mean*sumI
	
	mean  = I/N
	mean2 = N*mean*mean
	std   = sqrt( ( I2 + mean2 - I*2*mean )/(N-1) ) 

	return mean*N, std	
end

function linearBinSearch( array::Array{T1,1}, query::T2 ) where {T1<:Number,T2<:Number}

    idx = 0;
    @inbounds for e in 1:length(array)
        query <= array[ e ] && break
        idx = e;
    end

    return idx
end

function binaryBinSearch( array::Array{T1,1}, query::T2 ) where {T1<:Number,T2<:Number}

    hi, lo, probe = length( array ) + 1, 0, 0;

    @inbounds while 1 < hi - lo
        probe = ( hi + lo ) >>> 1
        if array[probe] < query
            lo = probe
        else
            hi = probe
        end
    end

    return lo
end

function downsample( image, step )
    h, w = size(image)
    rowindices = 1:step:h
    colindices = 1:step:w

    subsampled_img = zeros( eltype(image), length(rowindices), length(colindices) )

    ccont = 1
    for col in colindices
        rcont = 1
        for row in rowindices

            subsampled_img[rcont,ccont] = image[ row, col ]
            rcont += 1
        end
        ccont += 1
    end

    return subsampled_img
end


function meanDownsample( image::Array{}, scale )

	h, w   = size( image )
	nh, nw = length( scale:scale:h ), length( scale:scale:w )
	N      = scale*scale
	imgd   = zeros( eltype(image), nh, nw )

	ccont = 1;
	for col in scale:scale:w
		rcont = 1;
		for row in scale:scale:h
			mean = sum( image[ row-scale+1:row, col-scale+1:col ] )/N
			imgd[ rcont, ccont ] = mean;

		rcont += 1;
		end
	ccont += 1;
	end

	return imgd
end

function reverse!( image::AbstractArray{T,2} ) where {T<:Real }

	h, w = size( image );
	halfh, halfw = div.( size( image ), 2 );

	for col in 1:halfw
		t = image[ :, col ]
		image[ :, col ] .= image[ :, w-col+1 ]
		image[ :, w-col+1 ] .= t
	end

	for row in 1:halfh
		t = image[ row, : ]
		image[ row, : ] .= image[ h-row+1, : ]
		image[  h-row+1, : ] .= t
	end
end
