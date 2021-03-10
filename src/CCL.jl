function OAATccl( mask; connectivity=( (-1,0), (0,1), (1,0), (0,-1) ) )

    h, w  = size(mask);
    mask_ = copy(mask);

    components = [];
    comp_idx = 0;
    stack = [];

    for col in 1:w
        for row in 1:h

            # starting the process
            if mask_[ row, col ] > 0.5
                comp_idx += 1;
                push!( components, [ ] )
                push!( stack, (row,col) );
            end

            # Iterating
            while length(stack) > 0
                idx = pop!( stack )
                mask_[ idx... ] = 0
                push!( components[comp_idx], idx );

                for off in connectivity
                    r, c = min.( (h,w), max.( 1, idx .+ off ) )
                    if ( mask_[r,c] > 0.5 )
                        push!( stack, ( r, c ) );
                    end
                end
            end

        end
    end
    return components;
end
