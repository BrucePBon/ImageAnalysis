function localMinima( image::Array{T,2} ) where {T<:Real}

    minima_r = [];
    minima_c = [];

    prevGH = zeros( Bool, size(image,1) )
    for col in 2:size(image,2)

        column = view( image, :, col )

        gradH = image[ :, col ] .- image[ :, col - 1 ];

        prevGV = false;
        for row in 2:size(image,1)

            gradV = column[row] - column[row-1]

            if ( prevGV && gradV > 0 && prevGH[row] && gradH[row] > 0 )
                push!( minima_r, row-1 )
                push!( minima_c, col-1 )
            end

            prevGV = gradV < 0;
        end

        prevGH .= gradH .< 0
    end

    return( minima_r, minima_c )
end
