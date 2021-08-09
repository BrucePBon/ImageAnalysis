function entropy( img )
	sum = 0.0
	@inbounds @simd for idx in 1:length(img)
		sum += img[idx]
	end

	entrpy = 0.0
	@inbounds @simd for idx in 1:length(img)
		entrpy += (img[idx] == 0) ? 0.0 : log( 2, img[idx]/sum )*img[idx]/sum
	end

	return -entrpy
end
