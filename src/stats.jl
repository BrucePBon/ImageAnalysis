function mean( img::AbstractArray{T,2} ) where {T<:Real}

	ret = 0.0
	N = length( img );

	@simd for li in 1:N
		@inbounds ret += img[li]
	end

	return ret/N
end

function variance( img::AbstractArray{T,2} ) where {T<:Real}

	mean = mean( img )

	ret = 0.0
	N = length( img )

	@simd for li in 1:N
		@inbounds ret += ( img[li] - mean )*( img[li] - mean )
	end

	return ret/N
end

function variance( img::AbstractArray{T,2}, mean ) where {T<:Real}

	ret = 0.0
	N = length( img )

	@simd for li in 1:N
		@inbounds ret += ( img[li] - mean )*( img[li] - mean )
	end

	return ret/N
end

function std( img::AbstractArray{T,2} ) where {T<:Real}

	return sqrt( variance( img ) )
end

function std( img::AbstractArray{T,2}, mean ) where {T<:Real}

	return sqrt( variance( img, mean ) )
end

function cov( img1::Array{T1,2}, img2::Array{T2,2}) where {T1<:Real,T2<:Real}

	length( img1 ) !== length( img2 ) && ( error("Arrays have different lengths") )

	m1 = mean( img1 )
	s1 = std( img1, m1 )

	m2 = mean( img2 )
	s2 = std( img2, m2 )

	ret = 0.0
	N = length( img1 )

	@simd for li in 1:N
		@inbounds ret += ( img1[li] - m1 )*( img2[li] - m2 )
	end

	return ret/N
end

function correlation( img1::Array{T1,2}, img2::Array{T2,2} ) where {T1<:Real,T2<:Real}

	length( img1 ) !== length( img2 ) && ( error("Arrays have different lengths") )

	m1 = mean( img1 )
	s1 = std( img1, m1 )

	m2 = mean( img2 )
	s2 = std( img2, m2 )

	ret = 0.0
	N = length( img1 )

	@simd for li in 1:N
		@inbounds ret += ( img1[li] - m1 )*( img2[li] - m2 )
	end

	return ret/(N*s1*s2)
end

function crossCorrelationZNCC( f::Array{T,2}, g::Array{T,2} ) where {T<:Real}

		lN = div( length( f ), 16 );

		# statstics for ZNCC
		#=
		meanf  = mean( f );
		meanf2 = meanf * meanf;

		meang  = mean( g );
		meang2 = meang * meang;
		=#

		# integral Arrays to compute stdf at each possible translation
		sumf  = integralArray( f );
		sumf2 = integralArray( f, fun=(x)->(x*x) );
		sumg  = integralArray( g );
		sumg2 = integralArray( g, fun=(x)->(x*x) );

	    # correlation size = N + M - 1
	    csize = size(f) .+ size(g) .- 1;
	    h, w  = csize;
		num   = zeros( Float32, csize ); # we compute the numerator with FFT for each cross-correlation element
		corr  = zeros( Float32, csize );

		#=
		padf  = zeros( Complex{T}, csize );
		padf[ 1:size(f,1), 1:size(f,2) ] .= complex.( f .- meanf, 0.0 );

		padg  = zeros( Complex{T}, csize );
		padg[ 1:size(g,1), 1:size(g,2) ] .= complex.( g .- meang, 0.0 );

		plan  =  plan_fft!( padf );
		iplan = plan_ifft!( padf );

		# cross-correlation ( f - meanf ) and ( g - meang ) with FFT
		plan * padf;
		plan * padg;
		@inbounds @simd for e in 1:length(padf)
			padf[e] = conj( padf[e] ) * padg[e];
		end

		iplan * padf;
		@inbounds @simd for e in 1:length(padf)
			corr[e] = real( padf[e] )
		end

		shifts = div.( csize, 2 ).- 1
	    Base.circshift!( num, corr, shifts )
		=#

		# computing the denominator
		# Computing denominator through summed-area tables.
	    @inbounds for col in 1:csize[2]

			minCol_g = 1 + max( 0, col - size(f,2) )
			maxCol_g = min( col, size(g,2) )

			#
			minCol_f = max( 1, 1 + size(f,2) - col )
	        maxCol_f = size(f,2) - max( 0, col - size(g,2) )

			lenCol_g = length( minCol_g:maxCol_g );

	        for row in 1:csize[1]

				minRow_g = 1 + max( 0, row - size(f,1) );
            	maxRow_g = min( row, size(g,1) );

				#
				minRow_f = max( 1, 1 + size(f,1) - row )
	            maxRow_f = size(f,1) - max( 0, row - size(g,1) )

				lenRow_g = length( minRow_g:maxRow_g )
				N = lenRow_g*lenCol_g

				#
				F, F2 = integralAreaZNCC( sumf, sumf2, row, col, size(f), size(g) )
				G, G2 = integralAreaZNCC( sumg, sumg2, row, col, size(g), size(g) )

				#
				meanf  = F/N;
				meanf2 = meanf * meanf;
				stdf   = sqrt( F2 + N*meanf2 - 2*F*meanf );
				meang  = G/N;
				meang2 = meang * meang;
				stdg   = sqrt( G2 + N*meang2 - 2*G*meang );

				#
				subg = view( g, minRow_g:maxRow_g, minCol_g:maxCol_g )
				subf = view( f, minRow_f:maxRow_f, minCol_f:maxCol_f )

				#
				zncc = 0.0
				for e in 1:length(subg)
					zncc += ( subf[e] - meanf )*( subg[e] - meang )
				end

				#
				corr[ row, col ] = zncc/(stdg*stdf)

				#=
				if N < lN
					corr[ row, col ] = 0.0
				else
		            F, F2 = integralAreaZNCC( sumf, sumf2, row, col, size(f), size(g) )
					G, G2 = integralAreaZNCC( sumg, sumg2, row, col, size(g), size(g) )

		            corr[ row, col ] = num[ row, col ]/( sqrt( F2 + N*meanf2  - 2*F*meanf )*sqrt( G2 + N*meang2  - 2*G*meang ) );
				end
				=#
	        end
	    end

		return corr;
end
