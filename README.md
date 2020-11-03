# STC

This code has been used to create experiments for the paper "Exploring Non-Additive Distortion in Steganography" by Tomas Pevny and Andrew Ker, presented at 6th ACM Workshop on Information Hiding and Multimedia Security in Innsbruck, Austria on June 20--22, 2018.

At the moment, the code has been tested with ~~Julia 0.6.2.~~ Julia 1.5.2

After installation, go to the test directory and execute `julia runtests.jl`, which verifies that everything is correct. Be sure that `julia` knows the path to the package, otherwise do `push!(LOAD_PATH,"/path/to/directory/with/STC`.

# Examples

To embed an image using S-Uniward, run

`julia examples/suniward.jl --input path/to/imagefile.pgm -- output path/to/outputfile.pgm --height 3 -p 0.3`

  - `method`    additive / variable use additive / non-additive S-Uniward
  - `--height`  controls the height of the sub-matrix in the parity-check matrix
  - `-p`        controls the payload
  - `-h`        shows the help


# Using your own distortion function

To use your own distortion function, you need should create a type extending  `AbstractImageDistortion` abstract type and define following methods: `matchlsb!`, `readlsb`, `readlsbq`, `LSBcostfun`, `trypmone`

For a minimal viable example see `src/SUniward/minchanges.jl` demonstrating minimization of a number of embedding changes. 
