using LinearAlgebra


########################### 2D



function mean_structure_tensor( img::Array{<:Real,2}; square_size=(5,5), grid_size=(5,5), overlap=(0,0), threshold=0, mindist=5 )

    planarity = zeros( Float64, size(img) ); # planarity score for each (or a subset) of pixel in the input image

    sq_rad    = div.( square_size .- 1, 2 );       #  sq_rad  + 1 +  sq_rad  = ( 2 *  sq_rad  + 1 ) = square_size
    grid_rad  = div.(  grid_size  .- 1, 2 );       # grid_rad + 1 + grid_rad = ( 2 * grid_rad + 1 ) =  grid_size 
    pad_rad   = grid_rad .* square_size .+ sq_rad; # ( grid_rad * sq_size ) + sq_rad + 1 + sq_rad + ( grid_rad * sq_size ) = ( 2 * ( grid_rad * sq_size + sq_rad ) )
    pad_rad   = pad_rad .- overlap .* grid_rad;    # krita explained

    padsize   = size( img ) .+ 2 .* pad_rad ;
    intA      = zeros(Float64, padsize .+ 1);
    integralArray_pad!( img, intA, pad_rad );

    return mean_structure_tensor!( img, intA, planarity, sq_rad, grid_rad, overlap, pad_rad, threshold, mindist );
end

function mean_structure_tensor!( img::Array{<:Real,2}, intA, planarity, square_rad, grid_rad, overlap, pad_rad, threshold=0, mindist=5 )
    # convenient quantities
    lows  = (  1,  1  ) .+ ( pad_rad );
    highs = size( img ) .+ ( pad_rad );
    vert  = ( 2*square_rad[1] + 1 - overlap[1], 0  ); 
    horz  = ( 0, 2*square_rad[2] + 1  - overlap[2] );
    npix  = prod(( 2 .* square_rad .+ 1 )); 

    normv = [ (sqrt(y*y + x*x) <  mindist) ? (0.0, 0.0) : (y,x)./sqrt(y*y + x*x) for y in -grid_rad[1]:grid_rad[1], x in -grid_rad[2]:grid_rad[2] ];
    normv[grid_rad[1]+1,grid_rad[2]+1] = ( 0.0, 0.0 )
    println("testing", size(normv))


    # these will accomodate the means and signs of the 9x9 mean neighbourhood
    structure_tensor = zeros( Float32, 2, 2 ); 

    debug = false

    atans = zeros( Float32, size(planarity) ) ./ 0;
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl = ( y,x ) .- square_rad .- 1; 
        br = ( y,x ) .+ square_rad; 
        mean_center = ImageAnalysis.integralArea( intA, tl, br )/npix; 

        ( mean_center <= threshold ) && ( continue; )

        debug && println( (y, x), "\t", tl, "\t", br  )

        idx = 1;
        total_w = 0; 
        for xoff in -grid_rad[2]:grid_rad[2], yoff in -grid_rad[1]:grid_rad[1]

            off = vert .* yoff .+ horz .* xoff; 
            mean_each = ImageAnalysis.integralArea( intA, tl .+ off, br .+ off )/npix; 

            # local direction don't change, so they are initialized outside the loop
            direction = normv[idx];
            # weighting each direction by the mean intensity similarity
            dist     = mean_each #abs( mean_center - mean_each )*(mean_center*0.8 >= mean_each); 
            w        = dist #1/(1+dist); 
            total_w += w; 

            # constructing the structure tensor; 
            structure_tensor[1,1] += direction[1]*direction[1]*w
            structure_tensor[2,2] += direction[2]*direction[2]*w
            structure_tensor[2,1] += direction[1]*direction[2]*w
            idx += 1;   

            debug && println( "\t", off, ", ", tl .+ off, ", ", br .+ off, ", ", 
                              round.(direction,digits=3), ", ", round(dist,digits=3), ", ", round(w,digits=3), ", ", round.(structure_tensor,digits=3) )
        end            

        structure_tensor[1,2] = structure_tensor[2,1]; 
        structure_tensor ./= total_w; 

        e1, e2 = abs.( eigvals( structure_tensor ) ); 
        evecs  = eigvecs( structure_tensor ); 

        false &&  println( "\t\t", round.(structure_tensor,digits=3), "\t", (e1,e2) )

        e1e2score = 2/( e1/e2 + e2/e1 ); # [1,0] -> [:),:(]
        esmall = ( e1 < e2 ) ? e1 : e2; 
        vsmall = ( e1 < e2 ) ?  1 :  2; 

        ebig   = ( e1 < e2 ) ? e2 : e1; 
        escore = (1-esmall*2/(esmall+ebig))

        planarity[y-pad_rad[1],x-pad_rad[2]] = escore # ( e1/e2 + e2/e1 ) #1/(1+log10( e1/e2 + e2/e1 ))
        atans[y-pad_rad[1],x-pad_rad[2]] = atand( evecs[1,vsmall], evecs[2,vsmall] ) # ( e1/e2 + e2/e1 ) #1/(1+log10( e1/e2 + e2/e1 ))


        # reseting strucutre tensor 
        @inbounds @simd for idx in 1:length(structure_tensor)
            structure_tensor[idx] = 0.0
        end
    end
    return planarity, atans
end



##################### 3D



function mean_structure_tensor( vol::Array{<:Real,3}; cube_size=(5,5,5), grid_size=(5,5,5), overlap=(0,0,0), threshold=0, mindist=5 )

    planarity = zeros( Float64, size(vol) );       # planarity score for each (or a subset) of pixel in the input image

    cb_rad    = div.( cube_size .- 1, 2 );       #  cb_rad  + 1 +  cb_rad  = ( 2 *  cb_rad  + 1 ) = cube_size
    grid_rad  = div.( grid_size .- 1, 2 );       # grid_rad + 1 + grid_rad = ( 2 * grid_rad + 1 ) = grid_size 
    pad_rad   = grid_rad .* cube_size .+ cb_rad; # ( grid_rad * sq_size ) + cb_rad + 1 + cb_rad + ( grid_rad * sq_size ) = ( 2 * ( grid_rad * sq_size + cb_rad ) )
    pad_rad   =  pad_rad .- overlap .* grid_rad;    # krita explained

    padsize   = size( vol ) .+ 2 .* pad_rad ;
    intA      = zeros(Float64, padsize .+ 1);
    integralArray_pad!( vol, intA, pad_rad );

    return mean_structure_tensor!( vol, intA, planarity, cb_rad, grid_rad, overlap, pad_rad, threshold, mindist );
end

function mean_structure_tensor!( vol::Array{<:Real,3}, intA, planarity, cube_rad, grid_rad, overlap, pad_rad, threshold=0, mindist=5 )

    # convenient quantities
    npix  = prod((2 .* cube_rad .+ 1)); 
    lows  = ( 1, 1, 1 ) .+ ( pad_rad );
    highs = size( vol ) .+ ( pad_rad );
    vert  = ( 2*cube_rad[1] + 1 - overlap[1], 0,  0 ); 
    horz  = ( 0, 2*cube_rad[2] + 1  - overlap[2], 0 );
    deep  = ( 0, 0, 2*cube_rad[3] + 1  - overlap[3] );


    # computing normal vectors pointing to each cube in the grid
    grid_offs = UnitRange.( -1 .* grid_rad, grid_rad );     
    norm_vecs = Array{Tuple{Float32,Float32,Float32},3}(undef,length.(grid_offs)...); 
    tidx      = 1; 
    for z in grid_offs[3], x in grid_offs[2], y in grid_offs[1]
        norm_vecs[tidx] = (y,x,z)./sqrt(y*y + x*x + z*z)
        tidx += 1
    end
    norm_vecs[grid_rad[1]+1,grid_rad[2]+1, grid_rad[3]+1] = ( 0.0, 0.0, 0.0 )

    # these will accomodate the means and signs of the 9x9 mean neighbourhood
    structure_tensor = zeros( Float32, 3, 3 ); 

    debug = false
    @inbounds for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]

        if vol[y-pad_rad[1],x-pad_rad[2],z-pad_rad[3]] == 0
            continue
        end

        
        tl = ( y,x,z ) .- cube_rad .- 1; 
        br = ( y,x,z ) .+ cube_rad; 
        mean_center = ImageAnalysis.integralArea( intA, tl, br )/npix; 
        
        ( mean_center <= threshold ) && ( continue; )

        debug && println( (y, x, z), "\t", tl, "\t", br  )

        idx = 1;
        total_w = 0; 
        for zoff in grid_offs[3], xoff in grid_offs[2], yoff in grid_offs[1]

            off = vert .* yoff .+ horz .* xoff .+ deep .* zoff; 
            mean_each = ImageAnalysis.integralArea( intA, tl .+ off, br .+ off )/npix; 

            direction = norm_vecs[idx];

            # weighting each direction by the mean intensity similarity
            dist     = mean_each #abs( mean_center - mean_each )*(mean_center*0.8 >= mean_each); 
            w        = dist      #1/(1+dist); 
            total_w += w; 

            # constructing the structure tensor; 
            structure_tensor[1,1] += direction[1]*direction[1]*w
            structure_tensor[2,2] += direction[2]*direction[2]*w
            structure_tensor[3,3] += direction[3]*direction[3]*w
            structure_tensor[1,2] += direction[1]*direction[2]*w
            structure_tensor[1,3] += direction[1]*direction[3]*w
            structure_tensor[2,3] += direction[2]*direction[3]*w

            idx += 1;   

            debug && println( "\t", off, ", ", tl .+ off, ", ", br .+ off, ", ", 
                              round.(direction,digits=3), ", ", round(dist,digits=3), ", ", round(w,digits=3), ", ", round.(structure_tensor,digits=3) )
        end            

        structure_tensor[2,1] = structure_tensor[1,2]; 
        structure_tensor[3,1] = structure_tensor[1,3]; 
        structure_tensor[3,2] = structure_tensor[2,3]; 

        structure_tensor ./= total_w; 

        evals = abs.( eigvals( structure_tensor ) ); 
        evecs = eigvecs( structure_tensor ); 
        eperm = Base.sortperm( evals ); 
        e3, e2, e1 = evals[ eperm ]; 
        v1, v2, v3 = evecs[:,eperm[3]], evecs[:,eperm[2]], evecs[:,eperm[1]]; 
                                            #      plane              not plane
        smallratio   = (1-e3*3/sum(evals)); # ==1 (e[0]<<<e[:]), ==0 (e[3]=1/e[:]), 
        similarratio = 1/(e1/e2 + e2/e1);   # ==1 (e[1]==e[2] ), ==0 (e[1]!=e[2] )
        escore = similarratio*smallratio;   # the smaller the better

        false &&  println( "\t\t", round.(structure_tensor,digits=3), "\t", (e1,e2) )

        planarity[y-pad_rad[1],x-pad_rad[2],z-pad_rad[3]] = escore 

        # reseting structure tensor 
        @inbounds @simd for idx in 1:length(structure_tensor)
            structure_tensor[idx] = 0.0
        end
    end
    return planarity
end
