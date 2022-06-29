using CSV, DataFrames

sample_files = readdir(joinpath(@__DIR__, "data"))
sample_files = sample_files[occursin.("short", sample_files)]

sample_file_list = map(x -> CSV.File(joinpath(@__DIR__, "data", x)) |> DataFrame, sample_files)
snvs_df = vcat(sample_file_list[1], sample_file_list[2])
snvs_df = sample_file_list[1]

for i in 2:length(sample_file_list)

    snvs_df = vcat(snvs_df, sample_file_list[i])

end

CSV.write(joinpath(@__DIR__, "data/snvs.csv"), snvs_df)
