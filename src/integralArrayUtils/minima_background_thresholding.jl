# the idea is to apply this for bead removal, to provide an "automatic" threshold to remove background axima
function minima_background_th( vol::Array{<:Real,3}, rad, th; radmin=(2,2,2), fmin=1, radmax=(3,3,3), fmax=1, ovpmin=(0,0,0), ovpmax=(0,0,0) ); 

    volth = copy(vol);

    mn = mean_minima_pad( vol, radmin, 0, fmin=fmin, ovp=ovpmin );
    mx = mean_maxima_pad( vol, radmax, 0, fmax=fmax, ovp=ovpmax );
    
    # computing integral array of mn ( number of minima ) and 
    # of vol ( intensity value ), to be able to compute the
    # mean intensity of the minima in any square subregion of 
    # the volume

    mn_ia  = integralArray( Float32.( mn ) ); 
    vol_ia = integralArray_masked( vol, mn ); 
    
    # connected component labelling of the maxima

    out   = ImageComponentAnalysis.label_components( mx, trues(3,3,3) );
    Nlbls = maximum( out ); 

    centroids = [ [0.0,0.0,0.0] for x in 1:Nlbls ]; 
    numVoxels = [ 0 for x in 1:Nlbls ];
    meanValue = [ 0.0 for x in 1:Nlbls ];
    # O(N), computing centroid and mean intensity for each connected component
    @inbounds for dep in 1:size(mx,3), col in 1:size(mx,2), row in 1:size(mx,1)
        if ( out[row,col,dep] == 0 )
            continue
        end
        lbl_idx = out[row,col,dep];
        centroids[lbl_idx] .+= [row,col,dep]
        meanValue[lbl_idx]  += vol[row,col,dep]; 
        numVoxels[lbl_idx]  += 1; 
    end
    @simd for idx in 1:Nlbls
        @inbounds centroids[idx] ./= numVoxels[idx]
        @inbounds meanValue[idx]  /= numVoxels[idx]
    end

    # O(N), Storing the coordinates of each connected component
    coords = [ zeros( Int64, numVoxels[idx] ) for idx in 1:Nlbls] ; 
    counts = [ 0 for idx in 1:Nlbls ]; 
    for idx in 1:length(out)
        if out[idx] > 0
            counts[out[idx]] += 1
            coords[out[idx]][counts[out[idx]]] = idx; 
        end
    end

    # computing the mean intensity of minima around each centroid,
    # and discarting ccls whose intensity is not significantly higher
    # than the background (minima) around it

    for idx in 1:length(centroids)
        r, c, z = round.( Int64, centroids[idx] )

        needsfilt = false
        if any( (r,c,z) .- rad .< 1 ) || any( (r,c,z) .+ rad .> size(vol) )
            needsfilt = true
        else
            tl  = (r,c,z) .- rad .- 1
            br  = (r,c,z) .+ rad
            sum = integralArea( vol_ia, tl, br );
            N   = integralArea(  mn_ia, tl, br ); 

            if meanValue[idx] < sum/N*th

                needsfilt = true
            end
        end

        if needsfilt
            for coord in coords[idx]
                mx[coord...] = false
            end
        end
    end

    return mx
end