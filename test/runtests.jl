using STC.SUniward
using Base.Test
using Images

function testreflectedsetindex()
  a=zeros(20,30)
  b= padarray(a,Pad(:symmetric,[16,16],[16,16]))
  @testset begin
    for j in 1:size(a,2)
        for i in 1:size(a,1)
            a[i,j]=1
            SUniward.reflectedsetindex!(b,i,j,1)
            c = padarray(a,Pad(:symmetric,[16,16],[16,16]))
            @test sum(abs.(c.parent - b.parent)) == 0
        end
    end
  end
end

function testsuniwardcosts()
  cover=load("suniward.jld","cover")
  rho=load("suniward.jld","rho")
  @testset begin
      @test sum(abs.(rho-SUniward.suniwardcosts(cover)))<1e-6
  end
end


function testoneitemconv2_1(flt)
    # a=rand(1:30,30,30)
    a=rand(30,30)
    # i,j,v=3,3,1
    (i,j,v)=rand(1:size(a,1)),rand(1:size(a,2)),rand()
    aa = padarray(SUniward.sameconv2(a,flt),Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))
    a[i,j]+=v
    sum(abs.(SUniward.sameconv2(a,flt) - SUniward.oneitemconv2!(aa,flt,i,j,v).parent[size(flt,1)+1:end-size(flt,1),size(flt,2)+1:end-size(flt,2)]))
end

function testoneitemconv2_2(flt)
    a=rand(30,30)
    aa = padarray(zeros(size(a)),Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))
    for i in 1:size(a,1)
      for j in 1:size(a,2)
        SUniward.oneitemconv2!(aa,flt,i,j,a[i,j])
      end
    end
    sum(abs.(SUniward.sameconv2(a,flt) - aa.parent[size(flt,1)+1:end-size(flt,1),size(flt,2)+1:end-size(flt,2)]))
end

function testweightedconv(flt)
    sigma=1
    cover=rand(30,30)
    stego=rand(30,30)
    diff=cover-stego
    
    paddedcover= padarray(cover,Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))
    paddeddiff= padarray(diff,Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))
    invr = padarray(zeros(size(diff)),Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))

    rc = SUniward.sameconv2(paddedcover.parent, flt)
    invr.parent.=1./(abs.(rc)+sigma)

    aa = padarray(zeros(size(diff)),Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))
    for i in 1:size(diff,1)
      for j in 1:size(diff,2)
        SUniward.oneitemconv2!(aa,invr,flt,i,j,diff[i,j])
      end
    end

    costs=SUniward.sameconv2(diff,flt).*invr.parent[size(flt,1)+1:end-size(flt,1),size(flt,2)+1:end-size(flt,2)]
    sum(abs.(costs - aa.parent[size(flt,1)+1:end-size(flt,1),size(flt,2)+1:end-size(flt,2)]))
end

function testincrementalcosts(flt)
    sigma=1
    cover=rand(30,30)
    stego=rand(30,30)
    diff=cover-stego
    
    paddedcover= padarray(cover,Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))
    paddeddiff= padarray(diff,Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))
    invr = padarray(zeros(size(diff)),Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))

    rc = SUniward.sameconv2(paddedcover.parent, flt)
    invr.parent.=1./(abs(rc)+sigma)

    aa = padarray(zeros(size(diff)),Pad(:symmetric,[size(flt,1),size(flt,2)],[size(flt,1),size(flt,2)]))
    distortion=zero(eltype(cover))
    (ul,vl)=(cld(size(flt,1),2),cld(size(flt,2),2))
    (up,vp)=(fld(size(flt,1),2),fld(size(flt,2),2))
    for i in 1:size(diff,1)
      for j in 1:size(diff,2)
        idxs=(max(1,i-up):min(size(diff,1),i+ul-1),max(1,j-vp):min(size(diff,2),j+vl-1))
        # println((i,j), " ",idxs)
        distortion-=sumabs(aa[idxs...])
        SUniward.oneitemconv2!(aa,invr,flt,i,j,diff[i,j])
        distortion+=sumabs(aa[idxs...])
      end
    end

    costs=SUniward.sameconv2(diff,flt).*invr.parent[size(flt,1)+1:end-size(flt,1),size(flt,2)+1:end-size(flt,2)]
    sum(abs.(sum(abs.(costs)-distortion)))
end

function testincrementalsuniward()
    cover=UInt8.(rand(0:255,30,30))
    stego=UInt8.(rand(0:255,30,30))
    image=SUniward.IncrementalSUniward(cover)
    for i in 1:size(cover,1)
      for j in 1:size(cover,2)
        image[i,j]=stego[i,j]
      end
    end
    return(abs(SUniward.suniwardcosts(Float64.(cover),Float64.(stego))-image.distortion))
end

function testincrementalsuniward2()
    cover=UInt8.(rand(0:255,30,30))
    stego=UInt8.(rand(0:255,30,30))
    image=SUniward.IncrementalSUniward(cover)
    for i in 1:size(cover,1)
      for j in 1:size(cover,2)
        image[sub2ind(size(cover),i,j)]=stego[i,j]
      end
    end
    return(abs(SUniward.suniwardcosts(Float64.(cover),Float64.(stego))-image.distortion))
end


function testLSBcostfun()
    cover=UInt8.(rand(0:255,30,30))
    pidx=rand(1:length(cover))
    cover[pidx]=8
    stego=deepcopy(cover)
    stego[pidx]=9

    image=SUniward.IncrementalSUniward(cover);
    w0,w1=SUniward.LSBcostfun(image,pidx)
    
    r=sum(abs.(w1-SUniward.suniwardcosts(Float64.(cover),Float64.(stego))))

    cover=UInt8.(rand(0:255,30,30))
    pidx=rand(1:length(cover))
    cover[pidx]=9
    stego=deepcopy(cover)
    stego[pidx]=8

    image=SUniward.IncrementalSUniward(cover);
    w0,w1=SUniward.LSBcostfun(image,pidx)

    r+=sum(abs.(w0-SUniward.suniwardcosts(Float64.(cover),Float64.(stego))))
    return(r)
end

function testtrypmone()
    cover=UInt8.(rand(0:255,30,30))
    image=SUniward.IncrementalSUniward(cover);

    d1=0.0
    d2=0.0
    for pidx in 1:length(cover)
        d1 += min(SUniward.tryvalue(image,cover[pidx]+1,pidx),SUniward.tryvalue(image,cover[pidx]-1,pidx))
        d2 += min(SUniward.trypmone(image,pidx)...)

        s = rand([-1,1]);
        s = (image[pidx] == 0)? 1 : s;
        s = (image[pidx] == 255)? -1 : s;
        image[pidx]=cover[pidx] + s
    end

    return(abs(d1-d2))
end

function testbinaryembedding(s=8,h=3)
    println("binary embedding with s = $s h = $h")
    payload = 0.4;
    cover=rawview(channelview(load("1.pgm")));
    l=size(cover);
    s=div(s,2);
    cover = cover[div(l[1],2)-s+1:div(l[1],2)+s,div(l[2],2)-s+1:div(l[2],2)+s];

    numofcols=Int(round(1/payload));
    hhat=rand(1:2^h-1,numofcols);
    numofblocks=Int(floor(length(cover)/numofcols));
    message=rand(0:1,numofblocks);
    embpath = randperm(length(cover));

    stegos = SUniward.IncrementalSUniward(cover,2^h;sigma=1.0,T=Float16);
    @time (stego,distortion)=STC.variablecoding(cover,stegos,message,hhat,h,embpath);
    δ = Int.(cover) - Int.(stego.image)
    println(@sprintf("stc cost = %g true cost = %g pixel increased / decreased %d / %d",distortion,SUniward.suniwardcosts(cover,stego.image),sum(δ .> 0),sum(δ .< 0)));

    stegos = SUniward.IncrementalSUniward(cover,2^h;sigma=1.0,T=Float16);
    @time (stego,distortion)=STC.variablecoding(cover,stegos,message,hhat,h,embpath);
    δ = Int.(cover) - Int.(stego.image)
    println(@sprintf("stc cost = %g true cost = %g pixel increased / decreased %d / %d",distortion,SUniward.suniwardcosts(cover,stego.image),sum(δ .> 0),sum(δ .< 0)));

    stegos = [SUniward.SUniwardAdd(cover) for i in 1:2^h];
    @time (stego,distortion)=STC.variablecoding(cover,stegos,message,hhat,h,embpath);
    δ = Int.(cover) - Int.(stego.image)
    println(@sprintf("stc cost = %g true cost = %g pixel increased / decreased %d / %d",distortion,SUniward.suniwardcosts(cover,stego.image),sum(δ .> 0),sum(δ .< 0)));
    println()
end

function testternaryembedding(s=8,h=3,)
    println("ternary embedding with s = $s h = $h")
    payload = 0.66*0.1;
    cover=rawview(channelview(load("1.pgm")));
    l=size(cover);
    s=div(s,2);
    cover = cover[div(l[1],2)-s+1:div(l[1],2)+s,div(l[2],2)-s+1:div(l[2],2)+s];

    numofcols=Int(round(1/payload));
    hhat=rand(1:3^h-1,numofcols);
    numofblocks=Int(floor(length(cover)/numofcols));
    message=rand(0:2,numofblocks);
    embpath = randperm(length(cover));

    stegos = SUniward.IncrementalSUniward(cover,3^h;sigma=1.0,T=Float16);
    @time (stego,distortion)=STC.variablecoding(cover,stegos,message,hhat,h,embpath,3);
    δ = Int.(cover) - Int.(stego.image)
    println(@sprintf("stc cost = %g true cost = %g pixel increased / decreased %d / %d",distortion,SUniward.suniwardcosts(cover,stego.image),sum(δ .> 0),sum(δ .< 0)));

    stegos = [SUniward.SUniwardAdd(cover) for i in 1:3^h];
    @time (stego,distortion)=STC.variablecoding(cover,stegos,message,hhat,h,embpath,3);
    δ = Int.(cover) - Int.(stego.image);
    println(@sprintf("stc cost = %g true cost = %g pixel increased / decreased %d / %d",distortion,SUniward.suniwardcosts(cover,stego.image),sum(δ .> 0),sum(δ .< 0)));
    println()
end

function testsaturation()
    println("testing saturation")
    for k in [1,2,254,255]
        h = 3;
        cover = UInt8.(k*ones(8,8))
        payload = 0.66*0.4;
        numofcols=Int(round(1/payload));
        hhat=rand(1:3^h-1,numofcols);
        numofblocks=Int(floor(length(cover)/numofcols));
        message=rand(0:2,numofblocks);
        embpath = randperm(length(cover));
        stegos = [SUniward.SUniwardAdd(cover) for i in 1:3^h];
        @time (stego,distortion)=STC.variablecoding(cover,stegos,message,hhat,h,embpath,3);
        δ = Int.(cover) - Int.(stego.image)
        println(@sprintf("stc cost = %g true cost = %g pixel increased / decreased %d / %d",distortion,SUniward.suniwardcosts(cover,stego.image),sum(δ .> 0),sum(δ .< 0)));
    end 
    println()
end

function testadditiveembedding(h=3,s=8,p=0.3)
    payload = p;
    cover=rawview(channelview(load("1.pgm")));
    l=size(cover);
    s=div(s,2);
    cover = cover[div(l[1],2)-s+1:div(l[1],2)+s,div(l[2],2)-s+1:div(l[2],2)+s];

    numofcols=Int(round(1/payload));
    numofblocks=Int(floor(length(cover)/numofcols));
    message=rand(0:1,numofblocks);
    embpath = randperm(length(cover));
    hhat=rand(1:2^h-1,numofcols);
    ρ = SUniward.suniwardcosts(Float64.(cover))
    costfun(img,j) = ((mod(img[j],2) == 0)? (0.0,ρ[j]) : (ρ[j],0.0))
    stego,distortion = STC.additivecoding(cover,message,hhat,h,costfun,embpath);
    δ = Int.(cover) - Int.(stego)
    println(@sprintf("stc cost = %g true cost = %g pixel increased / decreased %d / %d",distortion,SUniward.suniwardcosts(cover,stego),sum(δ .> 0),sum(δ .< 0)));

    stegos = [SUniward.SUniwardAdd(cover) for i in 1:2^h]
    @time (stego,distortion)=STC.variablecoding(cover,stegos,message,hhat,h,embpath,2);
    δ = Int.(cover) - Int.(stego.image)
    println(@sprintf("stc cost = %g true cost = %g pixel increased / decreased %d / %d",distortion,SUniward.suniwardcosts(cover,stego.image),sum(δ .> 0),sum(δ .< 0)));
end


testsuniwardcosts()
testreflectedsetindex()

@testset begin
    @test testoneitemconv2_1([1 2 3; 3 4 5; 6 7 8])<1e-6
    @test testoneitemconv2_1([1 2; 3 4; 6 7])<1e-6
    @test testoneitemconv2_1([1 2 3; 3 4 5])<1e-6
    @test testoneitemconv2_1(randn(16,16))<1e-6
end

@testset begin
    @test testoneitemconv2_2([1 2 3; 3 4 5; 6 7 8])<1e-6
    @test testoneitemconv2_2([1 2; 3 4; 6 7])<1e-6
    @test testoneitemconv2_2([1 2 3; 3 4 5])<1e-6
    @test testoneitemconv2_2(randn(16,16))<1e-6
end

@testset begin
    @test testweightedconv([1 2 3; 3 4 5; 6 7 8])<1e-6
    @test testweightedconv([1 2; 3 4; 6 7])<1e-6
    @test testweightedconv([1 2 3; 3 4 5])<1e-6
    @test testweightedconv(randn(16,16))<1e-6
end

@testset begin
    @test testincrementalsuniward()<1e-6
    @test testLSBcostfun()<1e-6
    @test testtrypmone()<1e-6
end

testbinaryembedding(32,3)
testternaryembedding(32,3)
testsaturation()