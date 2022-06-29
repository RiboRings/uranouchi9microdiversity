using CSV, DataFrames

sample_files = readdir(joinpath(@__DIR__, "data"))
sample_files = sample_files[occursin.("snvs", sample_files)][2:end]

genome_list = CSV.File("top_genome_list.txt") |> DataFrame

for cur_sample in sample_files

    cur_df = CSV.File(joinpath("data", cur_sample)) |> DataFrame

    keep_rows = map(x -> x ∈ genome_list.selgen, cur_df.genome)
    out_df = cur_df[keep_rows, :]
    
    CSV.write(string(@__DIR__, "/data/", cur_sample, "_short.csv"), out_df)

end
