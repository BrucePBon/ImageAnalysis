function distanceArray( array::Array{T,N} ) where {T,N}

	distances = Array{Float32,N}(undef,size(array))
	center    = (size(array,1),size(array,2),size(array,3)) ./ 2 .+ 0.5

	@inbounds begin 
	for dep in 1:size(array,3)
		dd2 = ( dep - center[3] )^2
		for col in 1:size(array,2)
			dc2 = ( col - center[2] )^2
			dcd = dc2 + dd2
			for row in 1:size(array,1)
				dr = ( row - center[1] )

				distances[row,col,dep] = sqrt( muladd( dr, dr, dcd ) )
			end
		end
	end end

	return distances	
end

function distoid( input::Array{<:Real,N} ) where {N}

	distances = distanceArray( input )

	# distoid = sum( distances .* ( input ./ sumI ) ) = sum( distances .* input ) / sumI
	distoid = 0.0
	sumI    = 0.0
	@inbounds @simd for idx in 1:length(input)
		sumI    += input[idx]
		distoid += distances[idx]*input[idx]
	end

	return distoid / sumI 
end

function distoid( input::Array{<:Real,N}, distances ) where {N}

	# distoid = sum( distances .* ( input ./ sumI ) ) = sum( distances .* input ) / sumI
	distoid = 0.0
	sumI    = 0.0
	@inbounds @simd for idx in 1:length(input)
		sumI    += input[idx]
		distoid += distances[idx]*input[idx]
	end

	return distoid / sumI 
end
