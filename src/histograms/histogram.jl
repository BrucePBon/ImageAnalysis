""" 
	BINS
"""
# range
get_bins( bin_range::AbstractRange ) = Float32.( collect( bin_range ) )

# ( min, max, nbins )
function get_bins( bin_data )
	min, max, nbins = Float32.( bin_data );
	return get_bins( min:(max-min)/nbins:max  );
end

"""
	HISTOGRAMS
"""
function histogram( image, numbins::Integer ) 
	min, max = extrema( image );
	return histogram( image, (min,max,numbins) )
end

function histogram( image, bin_data::Union{Tuple{<:Real,<:Real,<:Real}, AbstractRange }; typ=Int64 )

    bins = get_bins( bin_data );
    hist = zeros( typ, length(bins) - 1 );
    outs = 0; # keeps count of pixels that don't fall into any interval

	oneT = typ(1);
    for pixel in image
        idx = binaryBinSearch( bins, pixel )

        if idx == 0
            outs += 1;
        else
            hist[idx] += oneT;
        end
    end

    return hist, outs, bins
end

function normhistogram( image, numbins::Integer ) 
	min, max = extrema( image );
	return normhistogram( image, (min,max,numbins) )
end

function normhistogram( image, bin_data::Union{Tuple{<:Real,<:Real,<:Real}, AbstractRange } )

	h, o, b = histogram( image, bin_data, typ=Float32 )

	N = length( image );
	@inbounds @simd for e in 1:length(h)
		h[e] /= N
	end

    return h, o, b
end

# Cumulative density function from histogram

function cdf( image, bin_data::Union{Tuple{T1,T2,T3}, AbstractRange{T4}} ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

    cdf, outs, bins  = histogram( image, bin_data, typ=Float32 )

	N = length( image )
	cdf[1] = cdf[1]/N
	@inbounds for e in 2:length(cdf)
		cdf[e] = cdf[e]/N + cdf[e-1]
	end

    return cdf, outs, bins
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

function eqim( img::Array{<:AbstractFloat,2}, bin_data::Union{Tuple{T1,T2,T3}, AbstractRange{T4}} ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

	bs   = bins( bin_data );
	_cdf, _, _ = cdf( img, bin_data );
    eqim = copy( img );

    for e in 1:length(img)
        pixel = img[e]
        for idx in 1:length(bs)-1
            if pixel > bs[idx] && pixel <= bs[idx+1]
                eqim[e] = floor( length(bs)*_cdf[idx] );
            end
        end

    end
    return eqim
end