
using Images
using FileIO
using STC
using STC.SUniward
using JLD
using ArgParse
using DataFrames

function parsearguments(args)
    s = ArgParseSettings("experiment_1.jl")

    @add_arg_table s begin
        "-p"
            default = 0.3
            arg_type=Float64
            help = "payload to embed (default 0.3bpp"
        "--input"
            default = "test/1.pgm"
            arg_type=String
            help = "path to the input image (in pgm format)"
        "--output"
            default = "/tmp/1.pgm"
            arg_type=String
            help = "path to the output image (in pgm format)"
        "--height"
            default = 3
            arg_type=Int
            help = "height of the embedding matrix"
        "--method"
            default = "variable"
            help = "list of methods to use during embedding"
    end

    parse_args(args, s,as_symbols=true)
end

function saveandshow(cover,stego,d1,ofname)
    δ = Int.(cover) - Int.(stego.image)
    d2 = SUniward.suniwardcosts(cover,stego.image);
    cp = sum(δ .> 0);
    cm = sum(δ .< 0);
    println(@sprintf("stc cost = %g true cost = %g pixel increased / decreased %d / %d",d1,d2,cp,cm));
    !isdir(dirname(ofname)) && mkpath(dirname(ofname))
    save(ofname,Gray.(stego.image/255))
    open(ofname[1:end-3]*"meta","w") do fid
        println(fid,@sprintf("%g %g %d %d",d1,d2,cp,cm))
    end
    true
end

settings=parsearguments(ARGS)
h = settings[:height]
cover=rawview(channelview(load(settings[:input])));
numofcols=Int(round(1/(0.66*settings[:p])));
hhat=rand(1:3^h-1,numofcols);
numofblocks=Int(floor(length(cover)/numofcols));
message=rand(0:2,numofblocks);
embpath = randperm(length(cover));

if settings[:method] == "variable"
    stegos = SUniward.IncrementalSUniward(cover,3^h;sigma=1.0,T=Float16);
    (stego,d1)=STC.variablecoding(cover,stegos,message,hhat,h,embpath,3);
    saveandshow(cover,stego,d1,settings[:output])
elseif settings[:method] == "additive"
    ρ = SUniward.suniwardcosts(Float64.(cover));
    stegos = [SUniward.SUniwardAdd(cover,ρ) for i in 1:3^h];
    (stego,d1)=STC.variablecoding(cover,stegos,message,hhat,h,embpath,3);
    saveandshow(cover,stego,d1,settings[:output])
else 
    error("unknown embedding method ",settings[:method])
end