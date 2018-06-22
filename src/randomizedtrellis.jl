function randomizedtrellis(cover,stegos,message,hhat,h,embpath,randomnodes::Float64)
    x=readlsb(cover,embpath);
    newstegos=deepcopy(stegos)
    stego=deepcopy(stegos[1])

    path=zeros(UInt8,2^h,length(x));
    numofblocks=length(message);

    wght=zeros(2^h);        #holds the costs of current "state"
    newwght=zeros(2^h);     #contained for holding costs of currently generated states
    wght[2:end]=Inf64;      
    indx=1  #this indexes the cover
    indm=1  #this indexes the message bits
    y=zeros(UInt8,length(x))
    calculated_costs = 0
    for i in 1:numofblocks #iterate over the number of subblocks, each encodes one message bit
        hh=min(h,length(message)-indm+1)  #this is effective height of the submatrix, which is needed
        for j in 1:length(hhat) #iterate over columns of the H matrix
            # println("pidx $(embpath[indx])")
            for k in 0:2^hh-1 # iterate over each state of the matrix
                # println()
                # there are two possibilities to reach the state k. Either we code zero, 
                # in which case the price (weight) is equal to weight of the same state in previous iteration + possible change,
                # or we code one and the state is reached from other state determined by the column of the hhat matrix.
                ind0=k+1 #index of the case, when the current LSB should be zero, i.e we are not adding current hhat[j]
                ind1=(k ⊻ hhat[j])&(2^hh-1)+1 #index to the state to which we come in case we are adding current hhat[j]

                #the determination of the cost is not entirely easy, since we arrive from two
                # different sources. The option here is either to use the current version of the image
                # and have the LSB cleared, which should correspond to w00, or come from the alternative and have
                # and have LSB set, which should correspond to w11

                #if we want to decide the way randomly
                w00 = isinf(wght[ind0])?wght[ind0]:rand()
                w11 = isinf(wght[ind1])?wght[ind1]:rand()

                if rand()>randomnodes
                    if Int(mod(stegos[ind0][embpath[indx]],2))==0
                        w00=wght[ind0]
                        w11=(isinf(wght[ind1]))?wght[ind1]:LSBcostfun(stegos[ind1],embpath[indx])[2] ;
                    else
                        w11=wght[ind1]
                        w00=(isinf(wght[ind0]))?wght[ind0]:LSBcostfun(stegos[ind0],embpath[indx])[1];
                    end
                    calculated_costs +=1
                end

                # choose the less costly path to the state, update weights
                # and update the current feature vector and sums
                if w00<w11
                    copy!(newstegos[ind0],stegos[ind0])
                    path[ind0,indx] = 0x0 
                    newwght[ind0] = (isinf(w00))?w00:matchlsb!(newstegos[ind0],0x0,embpath[indx])
                else
                    copy!(newstegos[ind0],stegos[ind1])
                    path[ind0,indx] = 0x1
                    newwght[ind0] = (isinf(w11))?w11:matchlsb!(newstegos[ind0],0x1,embpath[indx])
                end
            end
            indx+=1
            (stegos,newstegos)=swap(stegos,newstegos)
            (wght,newwght)=swap(wght,newwght)
        end

        # once we have traversed the trellis, we need to determine which path through trellises codes the message and
        # remove the others
        # println("reducing the number of states")
        for j = 0:2^(hh-1)-1
            midx=2*j + message[indm]+1
            wght[j+1] = wght[midx]
            copy!(stegos[j+1],stegos[midx])
         end
        wght[2^(hh-1)+1:end]=Inf64;
        newwght[2^(hh-1)+1:end]=Inf64;
        indm+=1
    end

    #there is no need to walkback the trellis, since we have already embedded image in stegos[1]
    embedding_cost = wght[1]
    info(@sprintf("embedding costs by stc = %g calculated_costs = %.3f",embedding_cost,Float64(calculated_costs)/(2^h*length(cover))))

    # check that the extracted message match the embedded message
    if !isinf(embedding_cost)
        extm=extract(readlsb(stegos[1],embpath),hhat,h)[1:length(message)];
        (sum(abs.(extm.!=message))==0)?info("message has been correctly extracted"):warn("message has not been correctly extracted");
    end
    return(deepcopy(stegos[1]),wght[1])
end

randomizedtrellis(cover,stegos,message,hhat,h,randomnodes::Float64) = randomizedtrellis(cover,stegos,message,hhat,h,randperm(length(cover)),randomnodes)
