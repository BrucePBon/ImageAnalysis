
# Convolution and cross-correlation
# -------------------------------------------

# Convolutions

# 1-. Standard FFT convolution. 
# Returns a convolution with size equal to: size(image) .+ size(kernel) .- 1. 

""" General implementation """
function FFTConvolution( image::Array{<:Real,N}, kernel::Array{<:Real,N}; typ=Float32 ) where {N}

	return FFTConvolution!( image, kernel, zeros( typ, size(image).+size(kernel).-1) ); 
end

""" In-place implementation, allowing to reuse "conv" """
function FFTConvolution!( image::Array{<:Real,N}, kernel::Array{<:Real,N}, conv::Array{T,N} ) where {T,N}

	@assert all( size(conv) .== size(image) .+ size(kernel) .- 1 ); 

	padi = zeros(Complex{T},size(conv)); 
	padk = zeros(Complex{T},size(conv));
	padi[ 1:size(image,1), 1:size(image,2), 1:size(image,3) ] .= convert.( Complex{T}, image )
	padk[ 1:size(kernel,1), 1:size(kernel,2), 1:size(kernel,3) ] .= convert.( Complex{T}, kernel )

	FFTW.fft!( padi );
	FFTW.fft!( padk );

	@inbounds @simd for idx in 1:length(padi)
		padi[idx] = padi[idx] * padk[idx]
	end

	FFTW.ifft!( padi ); 

	@inbounds @simd for idx in 1:length(padi)
		conv[idx] = real( padi[idx] )
	end

	return conv
end

# 2-. Cropped FFT convolution: 
# Returns a convolution with the same size as the input image

""" General implementation """
function FFTConvolution_crop( image::Array{<:Real,N}, kernel::Array{<:Real,N}; typ=Float32 ) where {N}
	return FFTConvolution_crop!( image, kernel, zeros( typ, size(image) ) ); 
end

""" In-place implementation, allowing to reuse "conv" """
function FFTConvolution_crop!( image::Array{<:Real,N}, kernel::Array{<:Real,N}, conv::Array{T,N} ) where {T,N}

	@assert all( size(conv) .== size(image) ); 
	
	padk = zeros( Complex{T}, size(image) .+ size(kernel) .- 1 );
	padk[ 1:size(kernel,1), 1:size(kernel,2), 1:size(kernel,3) ] .= convert.( Complex{T}, kernel )

	padi  = zeros( Complex{T}, size(image) .+ size(kernel) .- 1 );
	padi[ 1:size(image,1), 1:size(image,2), 1:size(image,3) ] .= convert.( Complex{T}, image )

	FFTW.fft!( padi );
	FFTW.fft!( padk );

	@inbounds @simd for idx in 1:length(padi)
		padi[idx] = padi[idx] * padk[idx]
	end

	FFTW.ifft!( padi ); 

	offs = div.( size( kernel ), 2 ); 
	zoff = ( N < 3 ) ? 0 : offs[3]; 

	for zet in 1+zoff:size(image,3)+zoff
		for col in 1+offs[2]:size(image,2)+offs[2]
		@simd for row in 1+offs[1]:size(image,1)+offs[1]
			conv[row-offs[1],col-offs[2],zet-zoff] = real( padi[row,col,zet] )
	end end end

	return conv
end

# Spatial convolution.
# This is generally slower than the FFT convolution, unless for kernel of small size 

function spatialConvolution( image::Array{T1,2}, kernel::Array{T2,2}; typ=Float32 ) where {T1<:Real,T2<:Real}

	N1, N2 = size( image  );
	M1, M2 = size( kernel );

	paddedImage = similar(image, N1+2*M1-2, N2+2*M2-2);
	paddedImage[ M1:(M1+N1-1), M2:(M2+N2-1) ] .= image;

	corrmatrix = zeros( typ, N1 + M1 - 1, N2 + M2 - 1 )

	for col in 1:( N2 + M2 - 1)
		for row in 1:( N1 + M1 - 1 )
			corrmatrix[ row, col ] = sum( paddedImage[ row:row+M1-1, col:col+M2-1 ] .* kernel )
		end
	end

	return corrmatrix
end