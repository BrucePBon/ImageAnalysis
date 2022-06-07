""" Computing the mean and standard deviation in one loop"""
function meanstd( img; typ=Float32 ) 

	# Mean = sumI / N 
	# STD numerator = sum( ( I - mean )^2 ) 
	#               = sum( I^2 + mean^2 - 2*I*mean ) 
	#               = sum(I^2) + N*mean^2 - 2mean*sumI
	#               = sum(I^2) + (sumI)^2/N - 2(sumI)^2/N
	#               = sum(I^2) - (sumI)^2/N

	# Both sumI and sumI2 can be computed in one loop

	sumI  = typ(0); 
	sumI2 = typ(0);
	N = length(img); 
	@inbounds @simd for idx in 1:N
		I = convert( typ, img[idx] );
		sumI  += I;
		sumI2 += I*I;
	end
	
	mean = sumI/N
	std  = sqrt( ( sumI2 - sumI^2/N )/(N-1) ) 

	return mean, std	
end

function maxval( a )

	len = length(a)

	if len < 4
		return maximum( a ) 
	end

	simdstart = len%4 + 1

    m1 = a[1]
    m2 = a[2]
    m3 = a[3]
    m4 = a[4]
	@inbounds @simd for i in simdstart:4:len
        m1 = ( a[ i ] > m1 ) ? a[ i ] : m1
        m2 = ( a[i+1] > m2 ) ? a[i+1] : m2
        m3 = ( a[i+2] > m3 ) ? a[i+2] : m3
        m4 = ( a[i+3] > m4 ) ? a[i+3] : m4
	end

    return maximum( [ m1, m2, m3, m4, maximum(a[1:simdstart-1]) ] )
end


function linearBinSearch( array::Array{T1,1}, query::T2 ) where {T1<:Number,T2<:Number}

	if ( query > array[end] || query < array[1] )
		return 0
	end

    idx = 1;
    @inbounds for e in 1:length(array)
        query <= array[ e ] && break
        idx = e;
    end

    return idx
end

function binaryBinSearch( array::Array{T1,1}, query::T2 ) where {T1<:Number,T2<:Number}

	if ( query > array[end] || query < array[1] )
		return 0
	end

    hi, lo, probe = length( array ) + 1, 1, 0;

    @inbounds while 1 < ( hi - lo )
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
