
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Data Collection

################################################################################

# Replace missings in a DataFrame.

function rep_missings(df::DataFrame)

    for col in names(df)
        df[ismissing.(df[:,col]), col] .= 0.0
    end

    return df

end


# Full join an array of dataframes by colname and return, removing missings.

function join_dfs(dfs::Array{DataFrame}, colname::String)

    length(dfs) >= 2 || error("Provide >1 DataFrame in the dfs array.")

    full_df = dfs[1]
    for i in 2:length(dfs)
        full_df = outerjoin(full_df, dfs[i], on = Symbol(colname))
    end
    full_df = rep_missings(full_df)

    return full_df

end


# Functions to pull out cell parameters from experiment outputs and convert
# into plotting-friendly data frames.

# Return the barcode counts of a Grow_Kill_Rec_Out structure.

function bc_counts(output::Grow_Kill_Rec_Out)

    # Create empty dataframe if no cells in output.
    if length(output.cells) == 0
        df = DataFrame(barcode = Float64[], count = Int64[])
    else
        df = DataFrame(barcode = map(x -> x.barcode, output.cells))
        # depreciated by 1.4.1
        # df = by(df, :barcode, nrow)
        df = combine(groupby(df, :barcode), nrow)
        colnames = [:barcode, :count]
        rename!(df, colnames)
    end

    return df

end


# Return the barcode counts for all replicates in an Experiment_Output structure.

function all_bc_counts(exp_output::Experiment_Output)

    all_dfs = Array{DataFrame}(undef, 0)
    for i in 1:length(exp_output.CO_outputs)
        temp_df = bc_counts(exp_output.CO_outputs[i])
        colnames = [:barcode, Symbol(exp_output.CO_outputs[i].rep_name)]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end
    for i in 1:length(exp_output.DT_outputs)
        temp_df = bc_counts(exp_output.DT_outputs[i])
        colnames = [:barcode, Symbol(exp_output.DT_outputs[i].rep_name)]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end

    full_df = join_dfs(all_dfs, "barcode")

    return full_df

end


# Return the mean R score of each barcode in an Experiment_Output structure.

function mean_bc_Rs(output::Grow_Kill_Rec_Out)

    # Create empty dataframe if no cells in output.
    if length(output.cells) == 0
        df = DataFrame(barcode = Float64[], mean_R = Float64[])
    else
        # First need function to return the mean R of a dfs R col.
        # function mean_R(df) ; mean(df[!,:R]) ; end
        df = DataFrame(barcode = map(x -> x.barcode, output.cells),
                             R = map(x -> x.R, output.cells))
        # depreciated by 1.4.1
        #df = by(df, :barcode, mean_R)
        df = combine(groupby(df, :barcode), :R => mean)
        colnames = [:barcode, :mean_R]
        rename!(df, colnames)
    end

    return df

end


# Return the sum of all R scores of each barcode in an Experiment_Output structure.

function sum_bc_Rs(output::Grow_Kill_Rec_Out)

    # Create empty dataframe if no cells in output.
    if length(output.cells) == 0
        df = DataFrame(barcode = Float64[], sum_R = Float64[])
    else
        # First need function to return the mean R of a dfs R col.
        # function mean_R(df) ; mean(df[!,:R]) ; end
        df = DataFrame(barcode = map(x -> x.barcode, output.cells),
                             R = map(x -> x.R, output.cells))
        # depreciated by 1.4.1
        #df = by(df, :barcode, mean_R)
        df = combine(groupby(df, :barcode), :R => sum)
        colnames = [:barcode, :sum_R]
        rename!(df, colnames)
    end

    return df

end


# R score by time

function R_by_t(output::Grow_Kill_Rec_Out)

    length(output.pulse_ts) == length(output.pulse_Rs) || error("pulse_ts should be the same length as pulse_Rs.")
    tvec = reduce(vcat, fill.(output.pulse_ts, length.(output.pulse_Rs)))
    Rvec = reduce(vcat, output.pulse_Rs)
    df = DataFrame(t = tvec, R = Rvec)

    return df

end


# N by time

function N_by_t(output::Grow_Kill_Rec_Out)

    length(output.Nvec) == length(output.tvec) || error("Nvec should be the same length as tvec.")
    df = DataFrame(t = output.tvec, N = output.Nvec)

    return df

end

# N, R and E by time

function NRE_by_t(output::Grow_Kill_Rec_Out)

    length(output.Nvec) == length(output.tvec) == length(output.Rvec) == length(output.Evec) || error("Nvec, Rvec and Evec should be the same length as tvec.")
    df = DataFrame(t = output.tvec, N = output.Nvec, R = output.Rvec, E = output.Evec)

    return df

end


# Return a dataframe of the N (Nvec) at each t_pulse (tvec) for an
# Experiment_Output structure.

function all_N_by_ts(exp_output::Experiment_Output)

    all_dfs = Array{DataFrame}(undef, 0)
    for i in 1:length(exp_output.CO_outputs)
        temp_df = N_by_t(exp_output.CO_outputs[i])
        temp_df[!, Symbol("Rep")] .= exp_output.CO_outputs[i].rep_name
        push!(all_dfs, temp_df)
    end
    for i in 1:length(exp_output.DT_outputs)
        temp_df = N_by_t(exp_output.DT_outputs[i])
        temp_df[!, Symbol("Rep")] .= exp_output.DT_outputs[i].rep_name
        push!(all_dfs, temp_df)
    end

    full_df = vcat(all_dfs...)

    return full_df

end

# Return a dataframe of the N (Nvec), R (Rvec) and E (Evec) at each
# t_pulse (tvec) for an Experiment_Output structure.

function all_NRE_by_ts(exp_output::Experiment_Output)

    all_dfs = Array{DataFrame}(undef, 0)
    for i in 1:length(exp_output.CO_outputs)
        temp_df = NRE_by_t(exp_output.CO_outputs[i])
        temp_df[!, Symbol("Rep")] .= exp_output.CO_outputs[i].rep_name
        push!(all_dfs, temp_df)
    end
    for i in 1:length(exp_output.DT_outputs)
        temp_df = NRE_by_t(exp_output.DT_outputs[i])
        temp_df[!, Symbol("Rep")] .= exp_output.DT_outputs[i].rep_name
        push!(all_dfs, temp_df)
    end

    full_df = vcat(all_dfs...)

    return full_df

end



# Return a DataFrame with all the cell parameters of all cells from an
# experiment output.

function all_cell_params(exp_output::Experiment_Output)

    all_dfs = Array{DataFrame}(undef, 0)
    for i in 1:length(exp_output.CO_outputs)
        temp_df = DataFrame(
        Cell_ID = map(x -> x.cell_ID, exp_output.CO_outputs[i].cells),
        barcode = map(x -> x.barcode, exp_output.CO_outputs[i].cells),
              b = map(x -> x.b, exp_output.CO_outputs[i].cells),
              d = map(x -> x.d, exp_output.CO_outputs[i].cells),
              R = map(x -> x.R, exp_output.CO_outputs[i].cells),
              )
        temp_df[!, Symbol("Rep")] .= exp_output.CO_outputs[i].rep_name
        push!(all_dfs, temp_df)
    end
    for i in 1:length(exp_output.DT_outputs)
        temp_df = DataFrame(
        Cell_ID = map(x -> x.cell_ID, exp_output.DT_outputs[i].cells),
        barcode = map(x -> x.barcode, exp_output.DT_outputs[i].cells),
              b = map(x -> x.b, exp_output.DT_outputs[i].cells),
              d = map(x -> x.d, exp_output.DT_outputs[i].cells),
              R = map(x -> x.R, exp_output.DT_outputs[i].cells),
              )
        temp_df[!, Symbol("Rep")] .= exp_output.DT_outputs[i].rep_name
        push!(all_dfs, temp_df)
    end

    full_df = vcat(all_dfs...)

    return full_df

end


# Return the barcode counts and corresponding number of resistant cells (R>0.0)
# in each barcode lineage across all replicates.

function all_bc_counts_tot_R(output)

    # Check output is one of the correct structure types.
    typeof(output) == Expanded_Split_Cells || typeof(output) == Experiment_Output || error("The output provided is not one of the supported data structures.")
    # Collect CO and DT outputs accordingly.
    CO_outputs = Array{CancerCell}[]
    DT_outputs = Array{CancerCell}[]

    if typeof(output) == Expanded_Split_Cells
        for i in 1:size(output.CO_flasks)[2]
            push!(CO_outputs, output.CO_flasks[:,i])
        end
        for i in 1:size(output.DT_flasks)[2]
            push!(DT_outputs, output.DT_flasks[:,i])
        end
    end
    if typeof(output) == Experiment_Output
        for i in 1:size(output.CO_outputs)[1]
            push!(CO_outputs, output.CO_outputs[i].cells)
        end
        for i in 1:size(output.DT_outputs)[1]
            push!(DT_outputs, output.DT_outputs[i].cells)
        end
    end

    all_dfs = Array{DataFrame}(undef, 0)

    for i in 1:length(CO_outputs)

        temp_bc_df = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]))
        temp_bc_df = combine(groupby(temp_bc_df, :barcode), nrow)
        colnames = [:barcode, :count]
        rename!(temp_bc_df, colnames)

        # Pull out the binary resistant phenotypic state of cells:
        R_bin_pheno = map(x -> x.R, CO_outputs[i])
        # Convert to 1 if >0 and 0 if =0.
        R_bin_pheno[R_bin_pheno .> 0] .= 1.0
        temp_R_df = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]),
                             R = R_bin_pheno)
        temp_R_df = combine(groupby(temp_R_df, :barcode), :R => sum)
        colnames = [:barcode, :n_R_cells]
        rename!(temp_R_df, colnames)

        rep_name = string("CO", i)
        temp_df = outerjoin(temp_bc_df, temp_R_df, on = :barcode)
        colnames = [:barcode, Symbol(rep_name), Symbol(string(rep_name, "_R"))]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end

    for i in 1:length(DT_outputs)

        temp_bc_df = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]))
        temp_bc_df = combine(groupby(temp_bc_df, :barcode), nrow)
        colnames = [:barcode, :count]
        rename!(temp_bc_df, colnames)

        # Pull out the binary resistant phenotypic state of cells:
        R_bin_pheno = map(x -> x.R, DT_outputs[i])
        # Convert to 1 if >0 and 0 if =0.
        R_bin_pheno[R_bin_pheno .> 0] .= 1.0
        temp_R_df = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]),
                             R = R_bin_pheno)
        temp_R_df = combine(groupby(temp_R_df, :barcode), :R => sum)
        colnames = [:barcode, :n_R_cells]
        rename!(temp_R_df, colnames)

        rep_name = string("DT", i)
        temp_df = outerjoin(temp_bc_df, temp_R_df, on = :barcode)
        colnames = [:barcode, Symbol(rep_name), Symbol(string(rep_name, "_R"))]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end

    full_df = join_dfs(all_dfs, "barcode")

    return full_df

end


# Return the barcode counts and corresponding mean resistant phenotype (R)
# in each barcode lineage across all replicates.

function all_bc_counts_mean_R(output)

    # Check output is one of the correct structure types.
    typeof(output) == Expanded_Split_Cells || typeof(output) == Experiment_Output || error("The output provided is not one of the supported data structures.")
    # Collect CO and DT outputs accordingly.
    CO_outputs = Array{CancerCell}[]
    DT_outputs = Array{CancerCell}[]

    if typeof(output) == Expanded_Split_Cells
        for i in 1:size(output.CO_flasks)[2]
            push!(CO_outputs, output.CO_flasks[:,i])
        end
        for i in 1:size(output.DT_flasks)[2]
            push!(DT_outputs, output.DT_flasks[:,i])
        end
    end
    if typeof(output) == Experiment_Output
        for i in 1:size(output.CO_outputs)[1]
            push!(CO_outputs, output.CO_outputs[i].cells)
        end
        for i in 1:size(output.DT_outputs)[1]
            push!(DT_outputs, output.DT_outputs[i].cells)
        end
    end

    all_dfs = Array{DataFrame}(undef, 0)

    for i in 1:length(CO_outputs)

        temp_bc_df = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]))
        temp_bc_df = combine(groupby(temp_bc_df, :barcode), nrow)
        colnames = [:barcode, :count]
        rename!(temp_bc_df, colnames)


        temp_R_df = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]),
                             R = map(x -> x.R, CO_outputs[i]))
        temp_R_df = combine(groupby(temp_R_df, :barcode), :R => mean)
        colnames = [:barcode, :mean_R]
        rename!(temp_R_df, colnames)

        rep_name = string("CO", i)
        temp_df = outerjoin(temp_bc_df, temp_R_df, on = :barcode)
        colnames = [:barcode, Symbol(rep_name), Symbol(string(rep_name, "_mR"))]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end

    for i in 1:length(DT_outputs)

        temp_bc_df = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]))
        temp_bc_df = combine(groupby(temp_bc_df, :barcode), nrow)
        colnames = [:barcode, :count]
        rename!(temp_bc_df, colnames)

        temp_R_df = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]),
                             R = map(x -> x.R, DT_outputs[i]))
        temp_R_df = combine(groupby(temp_R_df, :barcode), :R => mean)
        colnames = [:barcode, :mean_R]
        rename!(temp_R_df, colnames)

        rep_name = string("DT", i)
        temp_df = outerjoin(temp_bc_df, temp_R_df, on = :barcode)
        colnames = [:barcode, Symbol(rep_name), Symbol(string(rep_name, "_mR"))]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end

    full_df = join_dfs(all_dfs, "barcode")

    return full_df

end


# Return the barcode counts and corresponding:
#  i) total number of resistant (R>0.0) cells, and
# ii) mean resistant phenotype (R)
# in each barcode lineage across all replicates.

function all_bc_counts_tot_and_mean_R(output)

    # Check output is one of the correct structure types.
    typeof(output) == Expanded_Split_Cells || typeof(output) == Experiment_Output || error("The output provided is not one of the supported data structures.")
    # Collect CO and DT outputs accordingly.
    CO_outputs = Array{CancerCell}[]
    DT_outputs = Array{CancerCell}[]

    if typeof(output) == Expanded_Split_Cells
        for i in 1:size(output.CO_flasks)[2]
            push!(CO_outputs, output.CO_flasks[:,i])
        end
        for i in 1:size(output.DT_flasks)[2]
            push!(DT_outputs, output.DT_flasks[:,i])
        end
    end
    if typeof(output) == Experiment_Output
        for i in 1:size(output.CO_outputs)[1]
            push!(CO_outputs, output.CO_outputs[i].cells)
        end
        for i in 1:size(output.DT_outputs)[1]
            push!(DT_outputs, output.DT_outputs[i].cells)
        end
    end


    all_dfs = Array{DataFrame}(undef, 0)

    for i in 1:length(CO_outputs)

        temp_bc_df = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]))
        temp_bc_df = combine(groupby(temp_bc_df, :barcode), nrow)
        colnames = [:barcode, :count]
        rename!(temp_bc_df, colnames)

        # Pull out the binary resistant phenotypic state of cells:
        R_bin_pheno = map(x -> x.R, CO_outputs[i])
        # Convert to 1 if >0 and 0 if =0.
        R_bin_pheno[R_bin_pheno .> 0] .= 1.0

        temp_R_df1 = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]),
                                     R = R_bin_pheno)
        temp_R_df1 = combine(groupby(temp_R_df1, :barcode), :R => sum)
        temp_R_df2 = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]),
                                     R = map(x -> x.R, CO_outputs[i]))
        temp_R_df2 = combine(groupby(temp_R_df2, :barcode), :R => mean)
        colnames1 = [:barcode, :n_R_cells]
        colnames2 = [:barcode, :mean_R]
        rename!(temp_R_df1, colnames1)
        rename!(temp_R_df2, colnames2)
        # Join the total number of Resistant cells with mean Resistance
        # phenotype dataframes.
        temp_R_df = outerjoin(temp_R_df1, temp_R_df2, on = :barcode)
        # Now add to counts.
        rep_name = string("CO", i)
        temp_df = outerjoin(temp_bc_df, temp_R_df, on = :barcode)
        colnames = [:barcode, Symbol(rep_name), Symbol(string(rep_name, "_R")), Symbol(string(rep_name, "_mR"))]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end

    for i in 1:length(DT_outputs)

        temp_bc_df = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]))
        temp_bc_df = combine(groupby(temp_bc_df, :barcode), nrow)
        colnames = [:barcode, :count]
        rename!(temp_bc_df, colnames)

        # Pull out the binary resistant phenotypic state of cells:
        R_bin_pheno = map(x -> x.R, DT_outputs[i])
        # Convert to 1 if >0 and 0 if =0.
        R_bin_pheno[R_bin_pheno .> 0] .= 1.0

        temp_R_df1 = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]),
                                     R = R_bin_pheno)
        temp_R_df1 = combine(groupby(temp_R_df1, :barcode), :R => sum)
        temp_R_df2 = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]),
                                     R = map(x -> x.R, DT_outputs[i]))
        temp_R_df2 = combine(groupby(temp_R_df2, :barcode), :R => mean)
        colnames1 = [:barcode, :n_R_cells]
        colnames2 = [:barcode, :mean_R]
        rename!(temp_R_df1, colnames1)
        rename!(temp_R_df2, colnames2)
        # Join the total number of Resistant cells with mean Resistance
        # phenotype dataframes.
        temp_R_df = outerjoin(temp_R_df1, temp_R_df2, on = :barcode)
        # Now add to counts.
        rep_name = string("DT", i)
        temp_df = outerjoin(temp_bc_df, temp_R_df, on = :barcode)
        colnames = [:barcode, Symbol(rep_name), Symbol(string(rep_name, "_R")), Symbol(string(rep_name, "_mR"))]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end

    full_df = join_dfs(all_dfs, "barcode")

    return full_df

end


# Return the barcode counts and corresponding:
#  i) total number of resistant (R>0.0) cells, and
# ii) total number of escape mutations (E>0.0)
# in each barcode lineage across all replicates.

function all_bc_counts_tot_R_tot_E(output)

    # Check output is one of the correct structure types.
    typeof(output) == Expanded_Split_Cells || typeof(output) == Experiment_Output || error("The output provided is not one of the supported data structures.")
    # Collect CO and DT outputs accordingly.
    CO_outputs = Array{CancerCell}[]
    DT_outputs = Array{CancerCell}[]

    if typeof(output) == Expanded_Split_Cells
        for i in 1:size(output.CO_flasks)[2]
            push!(CO_outputs, output.CO_flasks[:,i])
        end
        for i in 1:size(output.DT_flasks)[2]
            push!(DT_outputs, output.DT_flasks[:,i])
        end
    end
    if typeof(output) == Experiment_Output
        for i in 1:size(output.CO_outputs)[1]
            push!(CO_outputs, output.CO_outputs[i].cells)
        end
        for i in 1:size(output.DT_outputs)[1]
            push!(DT_outputs, output.DT_outputs[i].cells)
        end
    end


    all_dfs = Array{DataFrame}(undef, 0)

    for i in 1:length(CO_outputs)

        temp_bc_df = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]))
        temp_bc_df = combine(groupby(temp_bc_df, :barcode), nrow)
        colnames = [:barcode, :count]
        rename!(temp_bc_df, colnames)

        # Pull out the binary resistant phenotypic state of cells:
        R_bin_pheno = map(x -> x.R, CO_outputs[i])
        # Convert to 1 if >0 and 0 if =0.
        R_bin_pheno[R_bin_pheno .> 0] .= 1.0

        temp_R_df = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]),
                                     R = R_bin_pheno)
        temp_R_df = combine(groupby(temp_R_df1, :barcode), :R => sum)
        temp_E_df = DataFrame(barcode = map(x -> x.barcode, CO_outputs[i]),
                                     E = map(x -> x.E, CO_outputs[i]))
        temp_E_df = combine(groupby(temp_E_df, :barcode), :E => sum)
        colnames1 = [:barcode, :n_R_cells]
        colnames2 = [:barcode, :n_E_cells]
        rename!(temp_R_df, colnames1)
        rename!(temp_e_df, colnames2)
        # Join the total number of Resistant cells with mean Resistance
        # phenotype dataframes.
        temp_RE_df = outerjoin(temp_R_df, temp_E_df, on = :barcode)
        # Now add to counts.
        rep_name = string("CO", i)
        temp_df = outerjoin(temp_bc_df, temp_RE_df, on = :barcode)
        colnames = [:barcode, Symbol(rep_name), Symbol(string(rep_name, "_R")), Symbol(string(rep_name, "_E"))]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end

    for i in 1:length(DT_outputs)

        temp_bc_df = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]))
        temp_bc_df = combine(groupby(temp_bc_df, :barcode), nrow)
        colnames = [:barcode, :count]
        rename!(temp_bc_df, colnames)

        # Pull out the binary resistant phenotypic state of cells:
        R_bin_pheno = map(x -> x.R, DT_outputs[i])
        # Convert to 1 if >0 and 0 if =0.
        R_bin_pheno[R_bin_pheno .> 0] .= 1.0

        temp_R_df = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]),
                                     R = R_bin_pheno)
        temp_R_df = combine(groupby(temp_R_df, :barcode), :R => sum)
        temp_E_df = DataFrame(barcode = map(x -> x.barcode, DT_outputs[i]),
                                     R = map(x -> x.E, DT_outputs[i]))
        temp_E_df = combine(groupby(temp_E_df, :barcode), :E => mean)
        colnames1 = [:barcode, :n_R_cells]
        colnames2 = [:barcode, :n_E_cells]
        rename!(temp_R_df, colnames1)
        rename!(temp_E_df, colnames2)
        # Join the total number of Resistant cells with mean Resistance
        # phenotype dataframes.
        temp_RE_df = outerjoin(temp_R_df, temp_E_df, on = :barcode)
        # Now add to counts.
        rep_name = string("DT", i)
        temp_df = outerjoin(temp_bc_df, temp_RE_df, on = :barcode)
        colnames = [:barcode, Symbol(rep_name), Symbol(string(rep_name, "_R")), Symbol(string(rep_name, "_E"))]
        rename!(temp_df, colnames)
        push!(all_dfs, temp_df)
    end

    full_df = join_dfs(all_dfs, "barcode")

    return full_df

end


# Simply return the barcode counts and sum of the number of resistant mutations
# in a CancerCell Array.

function bc_counts_tot_R_mut(cells::Array{CancerCell})

    bc_df = DataFrame(barcode = map(x -> x.barcode, cells))
    bc_df = combine(groupby(bc_df, :barcode), nrow)
    colnames = [:barcode, :count]
    rename!(bc_df, colnames)

    # Pull out the binary resistant phenotypic state of cells:
    R_bin_pheno = map(x -> x.R, cells)
    # Convert to 1 if >0 and 0 if =0.
    R_bin_pheno[R_bin_pheno .> 0] .= 1.0
    R_df = DataFrame(barcode = map(x -> x.barcode, cells),
                     R = R_bin_pheno)
    R_df = combine(groupby(R_df, :barcode), :R => sum)
    colnames = [:barcode, :n_R_cells]
    rename!(R_df, colnames)

    bc_R_df = outerjoin(bc_df, R_df, on = :barcode)
    return bc_R_df
end


# Just barcode counts given a rep name from cells.

function get_counts(cells::Array{CancerCell}, rep_name::String)
    df = DataFrame(barcode = map(x -> x.barcode, cells))
    # depreciated by 1.4.1
    # df = by(df, :barcode, nrow)
    df = combine(groupby(df, :barcode), nrow)
    colnames = [:barcode, Symbol(rep_name)]
    rename!(df, colnames)
    return df
end

# Same again but for barcode ids.

function get_cid_counts(cells::Array{CancerCell}, rep_name::String)
    df = DataFrame(cell_ID = map(x -> x.cell_ID, cells))
    # depreciated by 1.4.1
    # df = by(df, :barcode, nrow)
    df = combine(groupby(df, :cell_ID), nrow)
    colnames = [:cell_ID, Symbol(rep_name)]
    rename!(df, colnames)
    return df
end

# Only extract barcodes from a list of chosen lineages and return with a
# given rep name as the colname.

function get_chosen_counts(cells::Array{CancerCell}, chosen_bcs::Array{Float64},
    rep_name::String)

    df = DataFrame(barcode = map(x -> x.barcode, cells))
    # Only keep barcodes that are in the chosen_bcs array.
    df = filter(row -> row.barcode ∈ chosen_bcs, df)
    df = combine(groupby(df, :barcode), nrow)
    colnames = [:barcode, Symbol(rep_name)]
    rename!(df, colnames)
    return df

end



###############################################################################

# More detailed data collation functions that take an 'Expanded_Split_Cells'
# structure and return barcode resistance distributions of interest.

# Return the number of resistant lineages - those with at least 1 resistance
# mutation - found in i of I flasks following a growth expansions period.

function n_res_found_in_i(expand_split_out::Expanded_Split_Cells)

    df = all_bc_counts_tot_R(expand_split_out)

    @rput df

    R"""
    library(tidyr)
    library(plyr)
    library(dplyr)
    # Only look at DT replicates.
    df <- df[, c(grep("barcode|DT", colnames(df)))]

    # Only look at barcodes with a positive count in at least one DT flask.
    df <- df[rowSums(df[, c(grep("DT\\d{1}\\>", colnames(df)))] == 0) < length(grep("DT\\d{1}\\>", colnames(df))), ]

    # Want the number of barcodes found in i flasks (DTi > 0),
    # and  the number of resistant barcodes found in i flasks (DTi_R > 0).
    df["all_bc_found_in"] <- rowSums(df[, c(grep("DT\\d{1}\\>", colnames(df)))] > 0)
    df["res_bc_found_in"] <- rowSums(df[, c(grep("DT\\d{1}_R", colnames(df)))] > 0)

    # Convert to counts of each.
    all_bc_fi_counts <- data.frame(table(df$all_bc_found_in))
    res_bc_fi_counts <- data.frame(table(df$res_bc_found_in))
    colnames(all_bc_fi_counts) <- c("found_in", "all_count")
    colnames(res_bc_fi_counts) <- c("found_in", "res_count")
    fi_counts <- join(all_bc_fi_counts, res_bc_fi_counts)
    fi_counts[is.na(fi_counts)] <- 0

    """
    @rget fi_counts
    return fi_counts

end


# For each drug-treatment replicate flask from a simulation, return:
# i)  the total number of resistant lineages and resistant cells (R >= 1.0).
# ii) the number of resistant lineages and resistant cells unique to that flask.

function tot_and_unique_res_per_flask(expand_split_out::Expanded_Split_Cells)

    df = all_bc_counts_tot_R(expand_split_out)

    @rput df

    R"""
    cat("0.")
    library(tidyr)
    library(plyr)
    library(dplyr)
    # Only look at DT replicates.
    df <- df[, c(grep("barcode|DT", colnames(df)))]

    # Only look at barcodes with a positive count in at least one DT flask.
    df <- df[rowSums(df[, c(grep("DT\\d{1}\\>", colnames(df)))] == 0) < length(grep("DT\\d{1}\\>", colnames(df))), ]

    # Here we are only interested in resistant barcode lineages.
    df <- df[rowSums(df[, c(grep("DT\\d{1}_R", colnames(df)))]) > 0 ,]
    # Make empty dataframe if no count.
    if(nrow(df) == 0){
        df[1 ,] <- 0
        }
    cat(paste0("There are ", nrow(df), " resistant barcodes.\n"))
    # Create a second dataframe that only includes those unique to a flask.
    df_u <- df[rowSums(df[, c(grep("DT\\d{1}_R", colnames(df)))] > 0) == 1 ,]
    cat(paste0("There are ", nrow(df_u), " unique resistant barcodes.\n"))

    cat("1.")

    # Make empty dataframe if no count.
    if(nrow(df_u) == 0){
        df_u[1 ,] <- 0
        }

    cat("2.")

    # Split each into a df for cell counts and resistance mutation counts, then
    # turn each from wide to long, saving n_cels and n_res_muts per flask.

    dfc <- df[, c("barcode", grep("DT\\d{1}\\>", colnames(df), value = T))]
    dfR <- df[, c("barcode", grep("DT\\d{1}_R", colnames(df), value = T))]
    colnames(dfR) <- sub("_R", "", colnames(dfR))

    dfc <- gather(dfc, key = "flask", value = "n_cells", grep("barcode", colnames(dfc), invert = T))
    dfR <- gather(dfR, key = "flask", value = "n_res_muts", grep("barcode", colnames(dfR), invert = T))
    df <- join(dfc, dfR)

    cat("3.")

    df_uc <- df_u[, c("barcode", grep("DT\\d{1}\\>", colnames(df_u), value = T))]
    df_uR <- df_u[, c("barcode", grep("DT\\d{1}_R", colnames(df_u), value = T))]
    colnames(df_uR) <- sub("_R", "", colnames(df_uR))

    df_uc <- gather(df_uc, key = "flask", value = "n_cells", grep("barcode", colnames(df_uc), invert = T))
    df_uR <- gather(df_uR, key = "flask", value = "n_res_muts", grep("barcode", colnames(df_uR), invert = T))
    df_u <- join(df_uc, df_uR)

    cat("4.")

    # Want the total number of resistant barcode lineages and cells, per flask
    #  and the total number of resistant barcode lineages and cells unique to a
    #  flask.

    all_n_res_bc_lins <- plyr::ddply(df, .(flask), summarise, n_res_bcs=sum(n_res_muts > 0))
    all_n_res_cells <- plyr::ddply(df[df$n_res_muts > 0 ,], .(flask), summarise, n_res_cells=sum(n_cells))
    # Again, fill in n cells if no counts to sum.
    if(nrow(all_n_res_cells) == 0){
        all_n_res_cells <- data.frame(flask = all_n_res_bc_lins$flask, n_res_cells = 0)
    }

    cat("5.")

    uniq_n_res_bc_lins <- plyr::ddply(df_u, .(flask), summarise, uniq_n_res_bcs=sum(n_res_muts > 0))
    uniq_n_res_cells <- plyr::ddply(df_u[df_u$n_res_muts > 0 ,], .(flask), summarise, uniq_res_cells=sum(n_cells))
    # Again, fill in n cells if no counts to sum.
    if(nrow(uniq_n_res_cells) == 0){
        uniq_n_res_cells <- data.frame(flask = uniq_n_res_bc_lins$flask, uniq_res_cells = 0)
    }

    cat("6.")

    n_res_bc_lins_cells <- list(all_n_res_bc_lins, all_n_res_cells, uniq_n_res_bc_lins, uniq_n_res_cells)
    n_res_bc_lins_cells <-  Reduce(function(x, y, ...) join(x, y, ...), n_res_bc_lins_cells)


    # Add columns that are the proportion of a flask's resistant barcodes and
    # cells that are unique to that flask.
    n_res_bc_lins_cells["prop_uniq_res_bcs"] <- n_res_bc_lins_cells$uniq_n_res_bcs/n_res_bc_lins_cells$n_res_bcs
    n_res_bc_lins_cells["prop_uniq_res_cells"] <- n_res_bc_lins_cells$uniq_res_cells/n_res_bc_lins_cells$n_res_cells
    # Remove NAs
    n_res_bc_lins_cells[is.na(n_res_bc_lins_cells)] <- 0
    """

    @rget n_res_bc_lins_cells
    return n_res_bc_lins_cells

end


###############################################################################

# Collect the summary and correlation statistics directly from a simulation
# output.

# Summary Statistics:
# Save the following given a count dataframe where columns =
# samples and rows = barcode lineage, values = counts...
# NB here have to manually provide Passage information.

#   i) Gini coefficient.
#  ii) q-Div value (where q = Hill number).
# iii) Number of barcodes found in i flasks.
#  iv) The q-Div values but re-calculated using the combined
#      relative frequencies of lineages found in i of I flasks.
#   v) The qD_dissimilarity calculated using qD_diss function.

function collect_summ_stats(count_df::DataFrame, Passage::Int64)

    @rput count_df
    @rput Passage

    R"""
    #source("C:/Users/Freddie/Google Drive/Barcode_Simulations/Cancer_Barcode_Sim/CBC_Analysis_Custom_Functions.R")
    source("/data/BCI-EvoCa2/freddie/Barcode_Simulations/Cancer_Barcode_Sim/CBC_Analysis_Custom_Functions.R")
    colnames(count_df) <- paste0(colnames(count_df), "_P", Passage)
    summ_stats_df <- collect_summ_stats(count_df)
    """
    @rget summ_stats_df
    return summ_stats_df

end



# Correlation Statistics:
# Save the following given a count dataframe where columns =
# samples and rows = barcode lineage, values = counts...
# NB here have to manually provide Passage information.

#   i) Pearson's correlation co-efficient (r).
#  ii) Spearman's rank correlation (rho).
# iii) Pearson's, but only on shared lineages.
#  iv) Spearman's, but only on shared lineages.
#   v) Pearson's, but only on barcodes that make up
#      the cumulative frequency to X (i.e. top X of sample).
#  vi) Spearman's, but only on barcodes that make up
#      the cumulative frequency to X (i.e. top X of sample).

function collect_corr_stats(count_df::DataFrame, Passage::Int64)

    @rput count_df
    @rput Passage

    R"""
    #source("C:/Users/Freddie/Google Drive/Barcode_Simulations/Cancer_Barcode_Sim/CBC_Analysis_Custom_Functions.R")
    source("/data/BCI-EvoCa2/freddie/Barcode_Simulations/Cancer_Barcode_Sim/CBC_Analysis_Custom_Functions.R")
    colnames(count_df) <- paste0(colnames(count_df), "_P", Passage)
    corr_stats_df <- collect_corr_stats(count_df, 0.25)
    """
    @rget corr_stats_df
    return corr_stats_df

end


# Now, given the output of an 'expand_split_cells', return the summary
# statistics for:
#  i) all CO cells in the I CO flasks (to give a 'null' expectation), and
# ii) DT cells with a resistant phenotype (R > 0.0) in the I DT flasks.

function split_res_summ_stats(expand_split_out::Expanded_Split_Cells)

    # Need to split into CO and DT groups, and account for 0 resistant
    # mutations.
    df_out = all_bc_counts_tot_and_mean_R(expand_split_out)
    @rput df_out
    R"""
    df_CO <- df_out[, c(grep("CO\\d{1}\\>", colnames(df_out)))]
    # Pull out all the lineages with >1 Resistance 'mutation'.
    df_DT <- df_out[c(rowSums(df_out[, c(grep("DT\\d{1}_R\\>", colnames(df_out)))]) > 0) ,]
    df_DT <- df_DT[, c(grep("DT\\d{1}\\>|DT\\d{1}_R\\>", colnames(df_DT)))]
    # Need to replace all counts that don't have a resistance mutation (_R = 0)
    # with a 0 - only interested in resistant counts.
    # First replace any number of resistance mutations with a 1:
    df_DT[, c(grep("_R", colnames(df_DT)))][df_DT[, c(grep("_R", colnames(df_DT)))] > 0] <- 1
    # Now * count by _R - will convert all non-resistant lineages to 0.
    for(i in seq(1, ncol(df_DT)-1, 2)){
        df_DT[, i] <- df_DT[, i]*df_DT[, i+1]
    }
    # Remove _R columns
    df_DT <- df_DT[, c(grep("DT\\d\\>", colnames(df_DT)))]
    # Make empty dataframe if no counts.
    if(nrow(df_DT) == 0){
        df_DT[1 ,] <- 0
    }
    """
    @rget df_CO
    @rget df_DT

    summ_df_CO = collect_summ_stats(df_CO, 0)
    summ_df_DT = collect_summ_stats(df_DT, 0)

    summ_df = vcat(summ_df_CO, summ_df_DT)
    return summ_df

end

# Now, given the output of an 'expand_split_cells', return the summary
# statistics for:
#  i) all CO cells in the I CO flasks (to give a 'null' expectation), and
# ii) DT cells with a resistant phenotype (R > 0.0) in the I DT flasks.

function split_res_corr_stats(expand_split_out::Expanded_Split_Cells)

    # Need to split into CO and DT groups, and account for 0 resistant
    # mutations.
    df_out = all_bc_counts_tot_and_mean_R(expand_split_out)
    @rput df_out
    R"""
    df_CO <- df_out[, c(grep("CO\\d{1}\\>", colnames(df_out)))]
    # Pull out all the lineages with >1 Resistance 'mutation'.
    df_DT <- df_out[c(rowSums(df_out[, c(grep("DT\\d{1}_R\\>", colnames(df_out)))]) > 0) ,]
    df_DT <- df_DT[, c(grep("DT\\d{1}\\>|DT\\d{1}_R\\>", colnames(df_DT)))]
    # Need to replace all counts that don't have a resistance mutation (_R = 0)
    # with a 0 - only interested in resistant counts.
    # First replace any number of resistance mutations with a 1:
    df_DT[, c(grep("_R", colnames(df_DT)))][df_DT[, c(grep("_R", colnames(df_DT)))] > 0] <- 1
    # Now * count by _R - will convert all non-resistant lineages to 0.
    for(i in seq(1, ncol(df_DT)-1, 2)){
        df_DT[, i] <- df_DT[, i]*df_DT[, i+1]
    }
    # Remove _R columns
    df_DT <- df_DT[, c(grep("DT\\d\\>", colnames(df_DT)))]
    # Make empty dataframe if no counts.
    if(nrow(df_DT) == 0){
        df_DT[1 ,] <- 0
    }
    """
    @rget df_CO
    @rget df_DT

    corr_df_CO = collect_corr_stats(df_CO, 0)
    corr_df_DT = collect_corr_stats(df_DT, 0)

    corr_df = vcat(corr_df_CO, corr_df_DT)
    return corr_df

end
