# Computing bins.
# A range is suitable when the min, max and step values are known ( min:step:max ).
# When the step is not know, the number of divisions is used to compute a step.
# In this case, the inputs should a tuple with ( min, max, divisions )

function bins( inputs::Tuple{T1,T2,T3} ) where {T1<:Real,T2<:Real,T3<:Real}
    min, max, divs = convert.( Float32, inputs )
    return bins( min, max, (max-min)/divs )
end

bins( min, max, step ) = bins( ( min, max, step ) )
bins( min::Float32, max::Float32, step::Float32 ) = collect( min:step:max )
bins( range::AbstractRange{T} ) where { T<:Real } = collect( convert(Float32, range.start):convert(Float32, range.step):convert( Float32, range.stop ) )

# Histogram implementation

function histogram( image, bin_data::Union{Tuple{T1,T2,T3}, AbstractRange{T4}} ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

    bin_intervals = bins( bin_data );
    histogram = zeros( Int32, length(bin_intervals) - 1 );
    outs = 0; # pixels that don't fall into any interval are added to "outs"

    @inbounds for pixel in image

        # check if pixel falls into any of the intervals
        idx = binaryBinSearch( bin_intervals, pixel )

        if idx == 0
            outs += 1;
        else
            histogram[idx] += 1;
        end
    end

    return histogram, outs, bin_intervals
end

function normhistogram( image, bin_data::Union{Tuple{T1,T2,T3}, AbstractRange{T4}} ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

    bin_intervals = bins( bin_data );
    histogram = zeros( Float32, length(bin_intervals) - 1 );
    onef = Float32(1.0)
    outs = 0;

    @inbounds for pixel in image
        idx = binaryBinSearch( bin_intervals, pixel )
        if idx == 0
            outs += 1;
        else
            histogram[idx] += onef;
        end
    end
	N = length( image );
	@simd for e in 1:length(histogram)
		@inbounds histogram[e] /= N
	end

    return histogram, outs, bin_intervals
end

# Cumulative density function from histogram

function cdf( image, bin_data::Union{Tuple{T1,T2,T3}, AbstractRange{T4}} ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

    bin_intervals = bins( bin_data );
    cdf = zeros( Float32, length(bin_intervals) - 1 );
    onef = Float32(1.0);
    outs = 0;

    @inbounds for pixel in image
        idx = binaryBinSearch( bin_intervals, pixel )
        if idx == 0
            outs += 1;
        else
            cdf[idx] += onef;
        end
    end

	N = length( image )
	cdf[1] = cdf[1]/N
	@inbounds for e in 2:length(cdf)
		cdf[e] = cdf[e]/N + cdf[e-1]
	end

    return cdf, outs, bin_intervals
end

# Histogram image equalization

function eqim( img, bin_data::Union{Tuple{T1,T2,T3}, AbstractRange{T4}} ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

	bs   = bins( bin_data );
	_cdf, _, _ = cdf( img, bin_data );
    eqim = copy( img );

    for e in 1:length(img)
        pixel = img[e]
        for idx in 1:length(bs)-1
            if pixel > bs[idx] && pixel <= bs[idx+1]
                eqim[e] = floor( eltype(img), length(bs)*_cdf[idx] );
            end
        end

    end
    return eqim
end

# Otsu thresholding
# Does it only work on integers...?

function o_varB( k, sumP, sumIP, meanT )
	return ( meanT*sumP[k] - sumIP[k] )^2/( sumP[k]*( 1 - sumP[k] ) )
end

function otsu( p::Array{Float32,1} )

	L = length( p )
	sumP  = copy( p )
	sumIP = copy( p )
	@simd for i in 2:L
		 @inbounds sumP[i] += sumP[ i - 1 ]
	end
	@simd for i in 2:L
		 @inbounds sumIP[i] = i*p[i] + sumIP[ i - 1 ]
	end
	meanT = sumIP[end]

	maxK = 2
	maxvarB = o_varB( maxK, sumP, sumIP, meanT )
	for k in 3:L-1
		varB = o_varB( k, sumP, sumIP, meanT )
		if ( varB > maxvarB )
			maxK = k
			maxvarB = varB
		end
	end

	return maxK
end

function otsu( img, bin_data::Union{Tuple{T1,T2,T3}, AbstractRange{T4}} ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

	nh, _, _ = normhistogram( img, bin_data )
	ot = otsu( nh );

	println( "Otsu threshold = ", ot );
	mask = zeros( Bool, size(img) );

	for e in 1:length(img)
		mask[e] = img[e] > ot
	end

	return mask
end

function otsuArr( p::Array{Float32,1} )

	arr   = copy( p )

	L = length( p )
	sumP  = copy( p )
	sumIP = copy( p )
	@simd for i in 2:L
		 @inbounds sumP[i] += sumP[ i - 1 ]
	end
	@simd for i in 2:L
		 @inbounds sumIP[i] = i*p[i] + sumIP[ i - 1 ]
	end
	meanT = sumIP[end]

	for k in 1:L
		arr[k] = o_varB( k, sumP, sumIP, meanT )
	end

	return arr
end
