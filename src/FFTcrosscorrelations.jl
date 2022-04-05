
function FFTCrossCorrelation( image::Array{<:Real,N}, kernel::Array{<:Real,N}; typ=Float32, shift=true ) where {N}
	return FFTCrossCorrelation!( image, kernel, zeros( typ, size(image) .+ size(kernel) .- 1 ), shift=shift ); 
end

function FFTCrossCorrelation!( image::Array{<:Real,N}, kernel::Array{<:Real,N}, corr::Array{T,N}; shift=true ) where {T,N}

	@assert all( size(corr) .== size(image) .+ size(kernel) .- 1 ); 

	paddedkernel = zeros( Complex{T}, size(corr) ); 
	paddedkernel[ 1:size(kernel,1), 1:size(kernel,2), 1:size(kernel,3) ] .= convert.( Complex{T}, kernel )

	paddedimage  = zeros( Complex{T}, size(corr) ); 
	paddedimage[  1:size( image,1), 1:size( image,2), 1:size( image,3) ] .= convert.( Complex{T},  image )

	FFTW.fft!( paddedimage  );
	FFTW.fft!( paddedkernel );

	@inbounds @simd for idx in 1:length(paddedimage)
		paddedimage[idx] = paddedimage[idx] * conj( paddedkernel[idx] )
	end

	FFTW.ifft!( paddedimage ); 

	@inbounds @simd for idx in 1:length(paddedimage)
		corr[idx] = real( paddedimage[idx] )
	end

	if shift
		println( "ImageAnalysis.FFTCrossCorrelation!.shifting" ) 
        shifted = copy( corr ); 
        shifts  = div.( size(kernel) .+ 1, 2 ); 
        Base.circshift!( corr, shifted, shifts );
    end

	return corr
end

function FFTZNCC( image::Array{<:Real,N}, kernel::Array{<:Real,N}; typ=Float32 ) where {N}

	corr = Array{typ,N}(undef,size(image) .+ size(kernel) .- 1 );

	# writing kernel - meank into paddedkernel

	# sum ( k - meank )^2 = sum( k^2 ) + sum( meank^2 ) - 2*meank*sum( k ) = sumk2 + N*meank*meank - 2*meank*sumk = 
	#                                                                      = sumk2 + meank*sumk - 2*meank*sumk 
	#                                                                      = sumk2 - meank*sumk
	sumk  = typ( 0.0 );
	sumk2 = typ( 0.0 ); 

	@inbounds @simd for idx in 1:length(kernel)
		kf = typ( kernel[idx] ); 
		sumk  += kf
		sumk2 += kf*kf
	end

	meank = sumk/length(kernel); 
	denK  = sqrt( ( sumk2 - meank*sumk ) ) # sqrt( sum( k - meank )^2 )

	paddedkernel = zeros( Complex{typ}, size(corr) ); 
	for z in 1:size(kernel,3), c in 1:size(kernel,2)
		@simd for r in 1:size(kernel,1)
			paddedkernel[r,c,z] = Complex{typ}(typ(kernel[r,c,z]) - meank,0);
		end
	end
	
	# writing image - meani into paddedimage
	meani = typ(0.0)
	@inbounds @simd for idx in 1:length(image)
		meani += typ(image[idx]); 
	end
	meani /= length(image); 
	paddedimage  = zeros( Complex{typ}, size(corr) ); 
	for z in 1:size(image,3), c in 1:size(image,2)
		@simd for r in 1:size(image,1)
			paddedimage[r,c,z] = Complex{typ}(typ(image[r,c,z]) - meank,0);
		end
	end

	sumI2 = zeros( Float32, size(image) .+ 1 )
	integralArrayZNCC!( paddedimage, sumI2 ); 

	FFTW.fft!( paddedimage  );
	FFTW.fft!( paddedkernel );

	@inbounds @simd for idx in 1:length(paddedimage)
		paddedimage[idx] = paddedimage[idx] * conj( paddedkernel[idx] )
	end

	FFTW.ifft!( paddedimage ); 

	@inbounds @simd for idx in 1:length(paddedimage)
		corr[idx] = real( paddedimage[idx] )
	end

	shifted = copy( corr ); 
	shifts  = size(kernel) .- 1 #div.( size(kernel), 2 ); 
	Base.circshift!( corr, shifted, shifts );

	for zet in 1:size( corr,3 ), col in 1:size( corr, 2 ), row in 1:size( corr, 1 )
		z0, z1 = max( zet - size(kernel,3), 1 ) - 1, min( zet, size(image,3) )
		c0, c1 = max( col - size(kernel,2), 1 ) - 1, min( col, size(image,2) )
		r0, r1 = max( row - size(kernel,1), 1 ) - 1, min( row, size(image,1) )

		denI = integralArea( sumI2, (r0,c0,z0), (r1,c1,z1 ) ); 
		denI = ( denI <= 0.0 ) ? Inf : denI; 
		corr[ row, col, zet ] = corr[ row,col, zet ]/( sqrt( denI )*denK ); 
	end

	return corr
end

function FFTZNCC_crop( image::Array{<:Real,N}, kernel::Array{<:Real,N}; typ=Float32 ) where {N}

	corr = FFTZNCC( image, kernel, typ=typ ); 
	ranges = UnitRange.( size(kernel), size(image) .+ size(kernel ) .- 1 ); 
	
	return corr[ 1:size(image,1), 1:size(image,2), 1:size(image,3) ]; 
end


function FFTCrossCorrelation_crop( image::Array{<:Real,N}, kernel::Array{<:Real,N}; typ=Float32, shift=true ) where {N}
	return FFTCrossCorrelation_crop!( image, kernel, zeros( typ, size(image) ), shift=shift ); 
end

function FFTCrossCorrelation_crop!( image::Array{<:Real,N}, kernel::Array{<:Real,N}, corr::Array{T,N}; shift=true ) where {T,N}

	@assert all( size(corr) .== size(image) ); 

	uncrop = FFTCrossCorrelation( image, kernel, shift=shift ); 
	for z in 1:size(corr,3), x in 1:size(corr,2)
		@simd for y in 1:size(corr,1)
			@inbounds corr[y,x,z] = uncrop[y,x,z]
		end
	end

	return corr
end