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

# padding
function integralArray( image::AbstractArray{T,2}; type=Float32, fun=(x)->(x) ) where {T<:Real}
    intArr = zeros( type, size(image) .+ 1 )

	h, w = size( image )
	for col in 2:w+1
		for row in 2:h+1
			intArr[row,col] = fun(image[row-1,col-1]) + intArr[row-1,col] + intArr[row,col-1] - intArr[row-1,col-1]
		end
	end
	return intArr;
end

function integralArrayRev( image::AbstractArray{T,2}; type=Float32, fun=(x)->(x) ) where {T<:Real}

	intArr = zeros( type, size(image) .+ 1 )

	h, w = size( image )
	for col in 2:w+1
		for row in 2:h+1
			intArr[row,col] = fun(image[h-row+2,w-col+2]) + intArr[row-1,col] + intArr[row,col-1] - intArr[row-1,col-1]
		end
	end
	return intArr;


	#=
    intArr = integralArray( image, type=type, fun=fun );

	h, w = size( intArr ) .- 1;
	halfh, halfw = div.( ( h, w ), 2 );

	for col in 2:1+halfw
		t = intArr[ 2:end, col ]
		intArr[ 2:end, col ] .= intArr[ 2:end, w-col+2 ]
		intArr[ 2:end, w-col+2 ] .= t
	end

	for row in 2:1+halfh
		t = intArr[ row, 2:end ]
		intArr[ row, 2:end ] .= intArr[ h-row+2, 2:end ]
		intArr[  h-row+2, 2:end ] .= t
	end
	=#

	return intArr;
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

function integralArea( intArr, TL, BR )
	TL = TL .+ 1;
	BR = BR .+ 1;
	area = intArr[ BR[1], BR[2] ] - intArr[ BR[1], TL[2] ] - intArr[ TL[1], BR[2] ] + intArr[ TL[1], TL[2] ]
end

function integralAreaZNCC( sumS::Array{T,2}, sumS2::Array{T,2},
                           row, col, ss, si ) where {T<:Real}
    pad = 1;
    row = row + pad;
    col = col + pad;
    r0, r1 = max( 1, row - si[1] ), min( ss[1] + pad, row );
    c0, c1 = max( 1, col - si[2] ), min( ss[2] + pad, col );

    opS  = 0.0
    opS += sumS[ r1, c1 ];
    opS -= sumS[ r1, c0 ];
    opS -= sumS[ r0, c1 ];
    opS += sumS[ r0, c0 ];

    opS2  = 0.0
    opS2 += sumS2[ r1, c1 ];
    opS2 -= sumS2[ r1, c0 ];
    opS2 -= sumS2[ r0, c1 ];
    opS2 += sumS2[ r0, c0 ];

    return opS, opS2
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
