
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Experiment Functions

# A collection of functions that run pre-made experimental set-ups.

# To distinguish from other functions, all those that run an experiment should
# begin with 'Run_Exp...'.

################################################################################

# Import Threads functions for multi-threading.
import Base.Threads.@threads
import Base.Threads.@spawn


# Experimental set-up.

# 4 Replicates per treatment (Control vs Drug-Treatment).
# If Passage number > 1, can provide the output of a previous Passage as
# an 'Experiment_Input' object to continue with a N_seed randomly sampled
# cells from the previous flask.

function Run_Exp(N::Int64, b::Float64, d::Float64, p::Float64, mu::Float64,
    sig::Float64, del::Float64, R_real::String,
    n_pulse::Int64, Nmax::Int64, N_seed::Int64,
    t_CO::Float64, t_DT::Float64, Passage::Int64,
    insta_kill::Bool, lim_probs::Bool;
    al=0.0::Float64, psi=0.0::Float64,
    Exp_Input::Experiment_Input=Experiment_Input())

    if Passage == 1

        # First expand and split the cells for t = 6.0.
        #exp_split_cells = expand_split_cells(N, b, d, 6.0, N_seed, 4, p, mu, sig, del, psi=psi, al=al, R_real=R_real, use_lim_probs=lim_probs)
        # First expand and split the cells for t = 7.0: for QR_HCTbc.
        # exp_split_cells = expand_split_cells(N, b, d, 7.0, N_seed, 4, p, mu, sig, del, psi=psi, al=al, R_real=R_real, use_lim_probs=lim_probs)
        # First expand and split the cells for t = 9.0: for QR_SW6bc.
        exp_split_cells = expand_split_cells(N, b, d, 9.0, N_seed, 4, p, mu, sig, del, psi=psi, al=al, R_real=R_real, use_lim_probs=lim_probs)


        Exp_Input.CO_inputs = exp_split_cells.CO_flasks
        Exp_Input.DT_inputs = exp_split_cells.DT_flasks
    end

    CO_Outputs = Array{Grow_Kill_Rec_Out}(undef, 0)
    DT_Outputs = Array{Grow_Kill_Rec_Out}(undef, 0)

    c1 = @spawn push!(CO_Outputs, grow_kill_rec_cells(Exp_Input.CO_inputs[:,1], t_CO, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, rep_name = string("CO", 1, "P", Passage)))
    c2 = @spawn push!(CO_Outputs, grow_kill_rec_cells(Exp_Input.CO_inputs[:,2], t_CO, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, rep_name = string("CO", 2, "P", Passage)))
    c3 = @spawn push!(CO_Outputs, grow_kill_rec_cells(Exp_Input.CO_inputs[:,3], t_CO, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, rep_name = string("CO", 3, "P", Passage)))
    c4 = @spawn push!(CO_Outputs, grow_kill_rec_cells(Exp_Input.CO_inputs[:,4], t_CO, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, rep_name = string("CO", 4, "P", Passage)))
    d1 = @spawn push!(DT_Outputs, grow_kill_rec_cells(Exp_Input.DT_inputs[:,1], t_DT, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, drug_kill=true, insta_kill=insta_kill, rep_name = string("DT", 1, "P", Passage)))
    d2 = @spawn push!(DT_Outputs, grow_kill_rec_cells(Exp_Input.DT_inputs[:,2], t_DT, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, drug_kill=true, insta_kill=insta_kill, rep_name = string("DT", 2, "P", Passage)))
    d3 = @spawn push!(DT_Outputs, grow_kill_rec_cells(Exp_Input.DT_inputs[:,3], t_DT, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, drug_kill=true, insta_kill=insta_kill, rep_name = string("DT", 3, "P", Passage)))
    d4 = @spawn push!(DT_Outputs, grow_kill_rec_cells(Exp_Input.DT_inputs[:,4], t_DT, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, drug_kill=true, insta_kill=insta_kill, rep_name = string("DT", 4, "P", Passage)))

    wait(c1); wait(c2); wait(c3); wait(c4); wait(d1); wait(d2); wait(d3); wait(d4);

    return Experiment_Output(CO_Outputs, DT_Outputs)

end


# Run the experiment and save the output.

function Run_Exp_save_output(N::Int64, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64,
    R_real::String, n_pulse::Int64, Nmax::Int64,
    N_seed::Int64, t_CO::Float64, t_DT::Float64, Nsim::Int64, Passage::Int64,
    insta_kill::Bool, lim_probs::Bool;
    psi=0.0::Float64, al=0.0::Float64)

    cd("Outputs")

    # Save insta_kill and lim_probs as either 0/1.

    if insta_kill == true
        insta_kill_state = "1"
    elseif insta_kill == false
        insta_kill_state = "0"
    end
    if lim_probs == true
        lim_probs_state = "1"
    elseif lim_probs == false
        lim_probs_state = "0"
    end


    # Check for dir, and create if not
    dir_name = string("CBC_Exp_out_N-", N, "_b-", round(b, digits=3),
    "_d-", round(d, digits = 3), "_p-", round(p, digits = 10),
    "_mu-", round(mu, digits = 10), "_sig-", round(sig, digits = 10),
    "_del-", round(del, digits = 10),
    "_psi-", round(psi, digits = 10), "_al-", round(al, digits = 10),
    "_R_real-", R_real,
    "_n_pulse-", n_pulse, "_Nmax-", Nmax,
    "_N_seed-", N_seed, "_t_CO-", t_CO, "_t_DT-", t_DT,
    "_ik-", insta_kill_state, "_lp-", lim_probs_state)

    if isdir(dir_name) == true
        cd(dir_name)
    else
        mkdir(dir_name)
        cd(dir_name)
    end

    NameSim = string("Sim_", Nsim)
    # Make a sub-dir for the different sim iterations
    if isdir(NameSim) == true
        cd(NameSim)
    else
        mkdir(NameSim)
        cd(NameSim)
    end

    # Call outside the for loop first for scoping.
    exp_in = Experiment_Input()

    for P in 1:Passage
        if P == 1
            exp_out = Run_Exp(N, b, d, p, mu, sig, del, R_real, n_pulse, Nmax,
            N_seed, t_CO, t_DT, P, insta_kill, lim_probs,
            psi=psi, al=al)
        else

            exp_out = Run_Exp(N, b, d, p, mu, sig, del, R_real, n_pulse, Nmax,
            N_seed, t_CO, t_DT, P, insta_kill, lim_probs,
            psi=psi, al=al,
            Exp_Input=exp_in)
        end

        # If the current Passage has a replicate that has no cells
        # remaining, kill this sim.
        if sum(map(x -> length(x.cells), exp_out.CO_outputs) .== 0) > 0
            println("A control replicate has 0 cells remaining. Ending sim.")
            break
        end
        if sum(map(x -> length(x.cells), exp_out.DT_outputs) .== 0) > 0
            println("A drug-treatment replicate has 0 cells remaining. Ending sim.")
            break
        end
        # Save both the total number of resistant cells (R > 0.0) and the mean
        # resistant phenotype (mean(R)) as well as the counts per barcode
        # lineage.
        bc_R_df = all_bc_counts_tot_and_mean_R(exp_out)
        #NRE_df = all_NRE_by_ts(exp_out)
        N_df = all_N_by_ts(exp_out)

        #@rput bc_R_df ; @rput NRE_df ; @rput P
        @rput bc_R_df; @rput N_df; @rput P
        R"""
        write.csv(bc_R_df, file=paste0("bc_counts_tot_mean_R_P", P, ".csv"), row.names = F)
        #write.csv(NRE_df, file=paste0("NRE_by_t_P", P, ".csv"), row.names = F)
        write.csv(N_df, file=paste0("N_by_t_P", P, ".csv"), row.names = F)
        """

        # Might have had enough to save the barcode distributions and N vector,
        # but not enough cells to draw the next Passage's flasks from. Check
        # for that now.
        if sum(map(x -> length(x.cells), exp_out.CO_outputs) .< N_seed) > 0
            println("Not enough cells in the control replicates to seed the next Passages flask. Ending sim.")
            break
        end
        if sum(map(x -> length(x.cells), exp_out.DT_outputs) .< N_seed) > 0
            println("Not enough cells in the drug-treatment replicates to seed the next Passages flask. Ending sim.")
            break
        end


        exp_in = Experiment_Input()
        exp_in.CO_inputs=Array{CancerCell,2}(undef, N_seed, 4)
        exp_in.DT_inputs=Array{CancerCell,2}(undef, N_seed, 4)
        for i in 1:4
            exp_in.CO_inputs[:,i] = sample(exp_out.CO_outputs[i].cells, N_seed, replace=false)
            exp_in.DT_inputs[:,i] = sample(exp_out.DT_outputs[i].cells, N_seed, replace=false)
        end

    end

    cd("../../../")

end


################################################################################

# We want to just run the growth expansion of barcoded cells and subsequent
# splitting into i flasks to look at which combinations of parameters influence:
# i) the number of resistant barcode lineages (has at least one resistance
#    mutation) that are found in i of I flasks.
# ii) the number and proportion of resistant cells/barcode lineages that are
#    only found in 1 of the I total flasks.
# Repeat for a total of Nsim times, saving the barcode distributions of
# interest.

function Run_Exp_Split_Res_Dist(N::Int64, b::Float64, d::Float64,
    t_exp::Float64, N_seed::Int64, N_reps::Int64, p::Float64,
    mu::Float64, sig::Float64, del::Float64, lim_probs::Bool, Nsim::Int64;
    psi=0.0, al=0.0)

    cd("Outputs_Split_Res")

    if lim_probs == true
        lim_probs_state = "1"
    elseif lim_probs == false
        lim_probs_state = "0"
    end

    found_in_dfs = Array{DataFrame}(undef, 0)
    n_uniq_res_bc_lins_dfs = Array{DataFrame}(undef, 0)
    summ_dfs = Array{DataFrame}(undef, 0)
    #corr_dfs = Array{DataFrame}(undef, 0)

    for sim_rep in 1:Nsim

        param_df = DataFrame(N = N, b = round(b, digits = 3),
        d = round(d, digits = 3), t_exp = t_exp, N_seed = N_seed,
        N_reps = N_reps, p = round(p, digits = 10), mu = round(mu, digits = 10),
        sig = round(sig, digits = 10), del = round(del, digits = 10),
        psi = round(psi, digits = 10), al = round(al, digits = 2),
        lim_probs = lim_probs_state, Nsim = sim_rep)

        out = expand_split_cells(N, b, d, t_exp, N_seed, N_reps, p, mu,
        sig, del, psi=psi, al=al, use_lim_probs=lim_probs)

        found_in_df = n_res_found_in_i(out)
        n_uniq_res_bc_lins_df = tot_and_unique_res_per_flask(out)
        # Collect summary statistics.
        summ_df = split_res_summ_stats(out)
        #corr_df = split_res_corr_stats(out)
        # Add run parameter information
        found_in_df = hcat(found_in_df, repeat(param_df, nrow(found_in_df)))
        n_uniq_res_bc_lins_df = hcat(n_uniq_res_bc_lins_df, repeat(param_df, nrow(n_uniq_res_bc_lins_df)))
        summ_df = hcat(summ_df, repeat(param_df, nrow(summ_df)))
        #corr_df = hcat(corr_df, repeat(param_df, nrow(corr_df)))
        # Save to array of dataframes
        push!(found_in_dfs, found_in_df)
        push!(n_uniq_res_bc_lins_dfs, n_uniq_res_bc_lins_df)
        push!(summ_dfs, summ_df)
        #push!(corr_dfs, corr_df)

    end

    found_in_df = vcat(found_in_dfs...)
    n_uniq_res_bc_lins_df = vcat(n_uniq_res_bc_lins_dfs...)
    summ_df = vcat(summ_dfs...)
    #corr_df = vcat(corr_dfs...)

    found_in_name = string("CBC_Exp_Split_Res_Dist_found-in_N-", N, "_b-", round(b, digits=3),
    "_d-", d, "_t_exp-", t_exp, "_N_seed-", N_seed, "_N_reps-", N_reps,
    "_p-", round(p, digits = 10), "_mu-", round(mu, digits = 10),
    "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
    "_psi-", round(psi, digits = 10), "_al-", round(al, digits = 2),
    "_lp-", lim_probs_state, "_Nsim-", Nsim, "_output.csv")

    n_uniq_res_bc_lins_name = string("CBC_Exp_Split_Res_Dist_unique-per-flask_N-", N, "_b-", round(b, digits=3),
    "_d-", d, "_t_exp-", t_exp, "_N_seed-", N_seed, "_N_reps-", N_reps,
    "_p-", round(p, digits = 10), "_mu-", round(mu, digits = 10),
    "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
    "_psi-", round(psi, digits = 10), "_al-", round(al, digits = 2),
    "_lp-", lim_probs_state, "_Nsim-", Nsim, "_output.csv")

    summ_df_name = string("CBC_Exp_Split_Res_Dist_summ_stats-", N, "_b-", round(b, digits=3),
    "_d-", d, "_t_exp-", t_exp, "_N_seed-", N_seed, "_N_reps-", N_reps,
    "_p-", round(p, digits = 10), "_mu-", round(mu, digits = 10),
    "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
    "_psi-", round(psi, digits = 10), "_al-", round(al, digits = 2),
    "_lp-", lim_probs_state, "_Nsim-", Nsim, "_output.csv")

    corr_df_name = string("CBC_Exp_Split_Res_Dist_corr_stats-", N, "_b-", round(b, digits=3),
    "_d-", d, "_t_exp-", t_exp, "_N_seed-", N_seed, "_N_reps-", N_reps,
    "_p-", round(p, digits = 10), "_mu-", round(mu, digits = 10),
    "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
    "_psi-", round(psi, digits = 10), "_al-", round(al, digits = 2),
    "_lp-", lim_probs_state, "_Nsim-", Nsim, "_output.csv")

    CSV.write(found_in_name, found_in_df)
    CSV.write(n_uniq_res_bc_lins_name, n_uniq_res_bc_lins_df)
    CSV.write(summ_df_name, summ_df)
    #CSV.write(corr_df_name, corr_df)

    cd("../")

end


################################################################################

# To investigate what happens to the proportion of resistant cells over time,
# given a set of simulation parameters, we wan't to answer the following
# questions:
# i)  What happens to the proportion of resistant cells in a population if we
#     start with N0 sensitive (R = 0) and 0 resistant cells and continuously
#     bottleneck and grow them for total time t?
# ii) What happens to the proportion of resistant cells in a population if we
#     start with a N0 sensitive cells (R = 0) and then grow these cells for
#     a total time t?

function Run_Exp_Res_Prop_Over_Time(N0::Int64, b::Float64, d::Float64,
    mu::Float64, sig::Float64, del::Float64,
    Nmax::Int64, tmax::Float64, cont_growth::Bool;
    psi=0.0::Float64, al=0.0::Float64,
    nbottleneck=1::Int64, dt_rec=0.0::Float64)

    # Use neither limiting probabilities.
    init_cells = seed_cells(N0, b, d, 0.0, mu, sig, del,
    use_lim_probs=false, psi=psi, al=al)

    if cont_growth == true

        if nbottleneck > 1
            error("Can not use an nbottleneck > 1 and use the continuous growth version.")
        end

        prop_Rs = []
        tvec = []
        Nvec = []
        t = 0.0
        cells = init_cells
        push!(tvec, t)
        push!(Nvec, length(cells))
        n_Rs = sum(map(x -> x.R > 0, cells))
        prop_R = n_Rs/length(cells)
        push!(prop_Rs, prop_R)
        while t < tmax
            cells = grow_cells(cells, dt_rec,
                Nmax, mu, sig, del, psi=psi, al=al,
                R_real="b", drug_presence=0)
            # Update time and save recording vectors.
            t += dt_rec
            push!(tvec, t)
            push!(Nvec, length(cells))
            n_Rs = sum(map(x -> x.R > 0, cells))
            prop_R = n_Rs/length(cells)
            push!(prop_Rs, prop_R)
            if length(cells) > Nmax
                break
            end
        end

        return DataFrame(t = tvec, N = Nvec, prop_R = prop_Rs)

    else

        if dt_rec > 0.0
            error("dt_rec can only be used with the continuous growth version.")
        end

        prop_Rs = []
        tvec = []
        Nvec = []
        cells = init_cells
        push!(Nvec, length(cells))
        n_Rs = sum(map(x -> x.R > 0, cells))
        prop_R = n_Rs/length(cells)
        push!(prop_Rs, prop_R)

        for i in 1:nbottleneck
            cells = grow_cells(cells, tmax,
                Nmax, mu, sig, del, al=al, psi=psi,
                R_real="b", drug_presence=0).cells
            # Update time and save recording vectors.
            push!(Nvec, length(cells))
            n_Rs = sum(map(x -> x.R > 0, cells))
            prop_R = n_Rs/length(cells)
            push!(prop_Rs, prop_R)
            cells = sample(cells, N0, replace=false)
        end

        return DataFrame(nbottleneck = 0:nbottleneck, N = Nvec, prop_R = prop_Rs)

    end

end


# Run the experiment and save the output. Set Nmax to 10^10 - just use tmax.

function Run_Exp_Res_Prop_Over_Time_save_output(N0::Int64, b::Float64, d::Float64,
    mu::Float64, sig::Float64, del::Float64,
    Nmax::Int64, tmax::Float64, cont_growth::Bool, nsim::Int64;
    al=0.0::Float64, psi=0.0::Float64,
    nbottleneck=1::Int64, dt_rec=0.0::Float64)

    df = Run_Exp_Res_Prop_Over_Time(N0, b, d, mu, sig, del,
    Nmax, tmax, cont_growth, psi=psi, al=al,
    nbottleneck=nbottleneck, dt_rec=dt_rec)

    if cont_growth == true
        cg_state = "1.0"
    elseif cont_growth == false
        cg_state = "0.0"
    end

    cd("Outputs_Res_Prop")

    df_name = string("CBC_Exp_Res_Prop_Over_Time_N0-", N0,
    "_b-", round(b, digits=3),"_d-", d, "_mu-", round(mu, digits = 10),
    "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
    "_psi-", round(psi, digits = 10), "_al-", round(al, digits = 2),
    "_tmax-", tmax, "_cont-growth-", cg_state, "_nbottleneck-", nbottleneck,
    "_dt-rec-", dt_rec, "_nsim-", nsim, "_output.csv")

    CSV.write(df_name, df)

    cd("../")

end

################################################################################

# This experiment investigates whether we can identify the cost incurred by
# resistant cells by comparing successful lineages in the drug-treatment
# flasks with those found in the control flasks.
# Hypothesis: because cells with a resistant phenotype (R > 0.0) incur a
# growth cost according to δ (del), lineages that perform well in the
# drug-treatment flasks should perform less well in the control flask than the
# average lineage.
# With this difference, it might be possible to identify differences between
# the simulations where resistance evolves due to non-genetic sources of
# phenotypic variability (sig > 0.0) and those where there is a cost of
# resistance (del > 0.0).

# N.b. only ever run until the 2nd Passage and 4 replicates.

# Run the experiment and save the output.

function Run_Exp_CO_DT_diff_save_output(N::Int64, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64,
    R_real::String, n_pulse::Int64, Nmax::Int64,
    N_seed::Int64, t_CO::Float64, t_DT::Float64, Nsim::Int64,
    insta_kill::Bool, lim_probs::Bool;
    psi=0.0::Float64, al=0.0::Float64)

    #cd("C:\\Barcode_Output_Copies\\Outputs_CO-DT_Diff")
    cd("Outputs_CO-DT_Diff")

    # Save insta_kill and lim_probs as either 0/1.

    if insta_kill == true
        insta_kill_state = "1"
    elseif insta_kill == false
        insta_kill_state = "0"
    end
    if lim_probs == true
        lim_probs_state = "1"
    elseif lim_probs == false
        lim_probs_state = "0"
    end

    # Save df name according to params.
    df_name = string("CBC_Exp_CO_DT_diff_out_N-", N, "_b-", round(b, digits=3),
    "_d-", round(d, digits = 3), "_p-", round(p, digits = 10),
    "_mu-", round(mu, digits = 10), "_sig-", round(sig, digits = 10),
    "_del-", round(del, digits = 10),
    "_psi-", round(psi, digits = 10), "_al-", round(al, digits = 2),
    "_R_real-", R_real, "_n_pulse-", n_pulse, "_Nmax-", Nmax,
    "_N_seed-", N_seed, "_t_CO-", t_CO, "_t_DT-", t_DT,
    "_ik-", insta_kill_state, "_lp-", lim_probs_state, "_output.csv")

    # List of output dataframes.
    lin_dfs = Array{DataFrame}(undef, 0)

    for sim_rep in 1:Nsim

        # Call outside the for loop first for scoping.
        exp_in = Experiment_Input()
        bc_df = DataFrame()

        # Only ever run for 2 Passages.
        for P in 1:2
            if P == 1
                exp_out = Run_Exp(N, b, d, p, mu, sig, del, R_real, n_pulse, Nmax,
                N_seed, t_CO, t_DT, P, insta_kill, lim_probs,
                psi=psi, al=al)
            else
                exp_out = Run_Exp(N, b, d, p, mu, sig, del, R_real, n_pulse, Nmax,
                N_seed, t_CO, t_DT, P, insta_kill, lim_probs,
                psi=psi, al=al,
                Exp_Input=exp_in)
            end

            # If the current Passage has a replicate that has no cells
            # remaining, kill this sim.
            if sum(map(x -> length(x.cells), exp_out.CO_outputs) .== 0) > 0
                println("A control replicate has 0 cells remaining. Ending sim.")
                break
            end
            if sum(map(x -> length(x.cells), exp_out.DT_outputs) .== 0) > 0
                println("A drug-treatment replicate has 0 cells remaining. Ending sim.")
                break
            end

            exp_in = Experiment_Input()
            exp_in.CO_inputs=Array{CancerCell,2}(undef, N_seed, 4)
            exp_in.DT_inputs=Array{CancerCell,2}(undef, N_seed, 4)
            for i in 1:4
                exp_in.CO_inputs[:,i] = sample(exp_out.CO_outputs[i].cells, N_seed, replace=false)
                exp_in.DT_inputs[:,i] = sample(exp_out.DT_outputs[i].cells, N_seed, replace=false)
            end

            # Save counts and phenotype identities.
            bc_df = all_bc_counts_tot_and_mean_R(exp_out)
        end

        # Now need to save the chosen DT and CO lineages...
        @rput bc_df

        R"""
        # Just keep barcode ID and count columns.
        bc_df <- bc_df[, c(1, grep("\\D{2}\\d{1}\\>", colnames(bc_df)))]

        # Pull out the top 10 barcodes in DT2-4.
        bc_dfD <- bc_df[, c(1, grep("DT2|DT3|DT4", colnames(bc_df)))]
        bc_dfD <- data.frame(barcode = bc_dfD$barcode, DT = rowSums(bc_dfD[, 2:ncol(bc_dfD)]))
        bc_dfD <- bc_dfD[order(bc_dfD$DT, decreasing = T) ,]
        # Collect the top 10 barcode IDs (n.b. might include 0's).
        dt_bcs <- bc_dfD$barcode[1:10]

        # Pull out the top 10 barcodes in CO2-4.
        bc_dfC <- bc_df[, c(1, grep("CO2|CO3|CO4", colnames(bc_df)))]
        bc_dfC <- data.frame(barcode = bc_dfC$barcode, CO = rowSums(bc_dfC[, 2:ncol(bc_dfC)]))
        bc_dfC <- bc_dfC[order(bc_dfC$CO, decreasing = T) ,]
        # Collect the top 10 barcode IDs (n.b. might include 0's).
        co_bcs <- bc_dfC$barcode[1:10]

        # Collect control counts from CO1 for each group.
        dt_df <- bc_df[which(bc_df$barcode %in% dt_bcs) ,]
        co_df <- bc_df[which(bc_df$barcode %in% co_bcs) ,]
        dt_counts <- dt_df$CO1
        co_counts <- co_df$CO1
        """

        @rget dt_counts
        @rget co_counts
        count_df = DataFrame(dt_count = dt_counts, co_count = co_counts, nsim = sim_rep)
        push!(lin_dfs, count_df)

    end

    lin_df = vcat(lin_dfs...)

    # Join and merge output dataframes with sim param information.

    param_df = DataFrame(N = N, b = round(b, digits = 3),
    d = round(d, digits = 3), p = round(p, digits = 10),
    mu = round(mu, digits = 10), sig = round(sig, digits = 10),
    del = round(del, digits = 10),
    psi = round(psi, digits = 10), al = round(al, digits = 2),
    R_real = R_real, n_pulse = n_pulse, Nmax = Nmax, N_seed = N_seed,
    t_CO = t_CO, t_DT = t_DT, ik = insta_kill_state,
    lim_probs = lim_probs_state, Nsim = Nsim)

    param_df = repeat(param_df, nrow(lin_df))
    out_df = hcat(lin_df, param_df)

    CSV.write(df_name, out_df)

    #cd("C:\\Users\\Freddie\\Google Drive\\Barcode_Simulations\\Cancer_Barcode_Sim")
    cd("../")

end


################################################################################

# Assessing the distribution of resistant lineages and the relationship between
# resistance and lineage identity after sampling a POT sample, expanding for
# a given time (which i'm trying to optimise), and then sampling K cells W
# times (without replacement) - to simulate sampling into a 96-well plate.
# Then want to consider different values of time for the further expansion step
# (delt_2) and values of K. Ideally, want a few wells to be enriched for
# the cells that are resistant in the POT sample, and want these lineages to
# still be (the majority) resistant.

# Start with N0 cells, expand for t_exp (this will be cell-line specific), then
# take POT_N cells, expand these for delt_2 (to be optimised), sample K cells
# from this expanded pool (Nt_2) W times without replacement. Resistance
# evolution parameters are assigned as previously.

# function Run_Exp_Exp_Well_Dist(N0::Int64, b::Float64, d::Float64, p::Float64,
#     mu::Float64, sig::Float64, del::Float64, R_real::String, lim_probs::Bool,
#     t_exp::Float64, POT_N::Int64, delt_2::Float64, K::Int64, W::Int64;
#     al=0.0::Float64, psi=0.0::Float64)
#
#     # Seed uniquely barcoded cells.
#     init_cells = seed_cells(N0, b, d, p, mu, sig, del, use_lim_probs=lim_probs,
#                             psi=psi, al=al)
#
#     # Expand cells for initial expansion time in absence of drug.
#     exp_cells = grow_cells(init_cells, t_exp, 10^10, mu, sig, del, psi=psi,
#                            al=al, R_real=R_real, drug_presence=0, t_frac=0.025)
#
#     # Save the barcode counts and resistant phenotypes at this point.
#     fin_dfs = Array{DataFrame}(undef, 0)
#     exp_df = bc_counts_tot_R_mut(exp_cells.cells)
#     rename!(exp_df, [:barcode, :count_exp, :n_R_cells_exp])
#     push!(fin_dfs, exp_df)
#
#     # Sample POT_N of these expanded cells.
#     POT_cells = sample(exp_cells.cells, POT_N, replace = false)
#
#     # Save the barcode counts and resistant phenotypes at this point.
#     POT_df = bc_counts_tot_R_mut(POT_cells)
#     rename!(POT_df, [:barcode, :count_POT, :n_R_cells_POT])
#     push!(fin_dfs, POT_df)
#
#     # Expand the POT cells for an additional delt_2 time (this will be 'EP1').
#     EP1_cells = grow_cells(POT_cells, delt_2, 10^10, mu, sig, del, psi=psi,
#                            al=al, R_real=R_real, drug_presence=0, t_frac=0.025)
#
#     # Save the barcode counts and resistant phenotypes at this point.
#     EP1_df = bc_counts_tot_R_mut(EP1_cells.cells)
#     rename!(EP1_df, [:barcode, :count_EP1, :n_R_cells_EP1])
#     push!(fin_dfs, EP1_df)
#
#     # Sample K cells W times without replacement.
#     K * W < length(EP1_cells.cells) || error("There are not enough cells to split into the chosen number of wells. Need at least K * W.")
#
#     # Sample all the cells at once, and then reshape into replicate wells.
#     well_cells = sample(EP1_cells.cells, (K * W), replace = false)
#     well_cells = reshape(well_cells, (K, W))
#
#     # Go through and save the barcode counts and total number of resistant cells
#     # per replicate well.
#     well_dfs = Array{DataFrame}(undef, 0)
#     for i in 1:size(well_cells)[2]
#         temp_df = bc_counts_tot_R_mut(well_cells[:,i])
#         colnames = [:barcode, Symbol(string("count_", i)), Symbol(string("n_R_cells_", i))]
#         rename!(temp_df, colnames)
#         push!(well_dfs, temp_df)
#     end
#
#     well_df = join_dfs(well_dfs, "barcode")
#     push!(fin_dfs, well_df)
#
#     # Now join the dataframes from each point.
#     fin_df = join_dfs(fin_dfs, "barcode")
#
#     # Save df name according to params.
#     df_name = string("CBC_Exp_Well_Dist_N0-", N0,
#     "_b-", round(b, digits=3), "_d-", round(d, digits = 3),
#     "_p-", round(p, digits = 10), "_mu-", round(mu, digits = 10),
#     "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
#     "_R_real-", R_real, "_lp-", lim_probs, "_t_exp-", t_exp,
#     "_POT_N-", POT_N, "_delt_2-", delt_2, "_K-", K, "_W-", W, "_output.csv")
#
#     # Actually only interested in the resistant lineages, so pull these out
#     # in R before saving. Can infer the number of non-resistant lineages from
#     # the simulation parameters.
#
#     @rput fin_df; @rput df_name
#     R"""
#     Rbc_df <- fin_df[rowSums(fin_df[, grep("n_R_cells.*", colnames(fin_df))]) > 0 ,]
#     write.csv(Rbc_df, df_name)
#     """
#
# end


# # HCTbc parameters:
#
# delt_2s = [1.0, 2.0, 3.0]
# Ks = [1000, 2000, 5000, 10000]
#
# for delt_2 in delt_2s
#     for K in Ks
#
#         Run_Exp_Exp_Well_Dist(10^6, 0.693, 0.07, 0.0, 10^-6, 0.1, 0.0, "l", true,
#         7.0, 10^6, delt_2, K, 96)
#
#     end
# end
#
# # SW6bc parameters:
#
# delt_2s = [1.0, 2.0, 3.0]
# Ks = [1000, 2000, 5000, 10000]
#
# for delt_2 in delt_2s
#     for K in Ks
#
#         Run_Exp_Exp_Well_Dist(10^6, 0.565, 0.026, 0.0, 10^-6, 0.01, 0.0, "l", true,
#         9.0, 10^6, delt_2, K, 96)
#
#     end
# end


################################################################################

# Same as above, but don't perform subsequent expansion and sampling into
# replicate 96-wells... instead, do the initial expansion, sample n x replicates
# of K cells, and save the resistant lineages and cell counts for the POT, and
# these sub-samples. I can then use these to illustrate how different parameter
# values that lead to the same equilibrium frequencies of resistance can
# still lead to different lineage distributions...
# (I is the number of replicates to sample).

function Run_Exp_Samp_I_Reps(N0::Int64, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64, R_real::String, lim_probs::Bool,
    t_exp::Float64, K::Int64, I::Int64;
    al=0.0::Float64, psi=0.0::Float64)

    # Seed uniquely barcoded cells.
    init_cells = seed_cells(N0, b, d, p, mu, sig, del, use_lim_probs=lim_probs,
                            psi=psi, al=al)

    # Expand cells for initial expansion time in absence of drug.
    exp_cells = grow_cells(init_cells, t_exp, 10^10, mu, sig, del, psi=psi,
                           al=al, R_real=R_real, drug_presence=0, t_frac=0.025)

    # Save the barcode counts and resistant phenotypes at this point.
    fin_dfs = Array{DataFrame}(undef, 0)
    exp_df = bc_counts_tot_R_mut(exp_cells.cells)
    rename!(exp_df, [:barcode, :count_exp, :n_R_cells_exp])
    push!(fin_dfs, exp_df)

    # Sample K cells for I replicates:
    K * I < length(exp_cells.cells) || error("There are not enough cells to split into the chosen number of replicates. Need at least K * I.")

    # Sample all the cells at once, and then reshape into replicate wells.
    rep_cells = sample(exp_cells.cells, (K * I), replace = false)
    rep_cells = reshape(rep_cells, (K, I))

    # Go through and save the barcode counts and total number of resistant cells
    # per replicate well.
    rep_dfs = Array{DataFrame}(undef, 0)
    for i in 1:size(rep_cells)[2]
        temp_df = bc_counts_tot_R_mut(rep_cells[:,i])
        colnames = [:barcode, Symbol(string("count_", i)), Symbol(string("n_R_cells_", i))]
        rename!(temp_df, colnames)
        push!(rep_dfs, temp_df)
    end

    rep_df = join_dfs(rep_dfs, "barcode")
    push!(fin_dfs, rep_df)

    # Now join the dataframes from each point.
    fin_df = join_dfs(fin_dfs, "barcode")

    # Save df name according to params.
    df_name = string("CBC_Exp_Well_Dist_N0-", N0,
    "_b-", round(b, digits=3), "_d-", round(d, digits = 3),
    "_p-", round(p, digits = 10), "_mu-", round(mu, digits = 10),
    "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
    "_R_real-", R_real, "_lp-", lim_probs, "_t_exp-", t_exp,
    "_K-", K, "_I-", I, "_output.csv")

    # Actually only interested in the resistant lineages, so pull these out
    # in R before saving. Can infer the number of non-resistant lineages from
    # the simulation parameters.

    @rput fin_df; @rput df_name
    R"""
    Rbc_df <- fin_df[rowSums(fin_df[, grep("n_R_cells.*", colnames(fin_df))]) > 0 ,]
    write.csv(Rbc_df, df_name)
    """

end



# Run_Exp_Samp_I_Reps(2*10^5, 0.893, 0.200, 0.0, 10^-7, 0.001, 0.0, "l", true,
#          6.0, 2*10^5, 2)
#
# Run_Exp_Samp_I_Reps(2*10^5, 0.893, 0.200, 0.0, 10^-7, 0.0, 0.0025, "l", true,
#          6.0, 2*10^5, 2)
#
# Run_Exp_Samp_I_Reps(2*10^5, 0.893, 0.200, 0.0, 10^-5, 0.1, 0.0, "l", true,
#          6.0, 2*10^5, 2)
#
# Run_Exp_Samp_I_Reps(2*10^5, 0.893, 0.200, 0.0, 10^-5, 0.0, 0.25, "l", true,
#          6.0, 2*10^5, 2)



################################################################################

# Want to recreate a potential in vitro experiment that aims to distinguish
# between the phenotypic switching (mu + sig) and cost of resistance (mu + del)
# evolutionary scenarios.
# Does this by expanding the cells for delt_1 (to recreate the beginning of all
# the barcode experiments), taking a sample of these cells (POT_N), adding
# CFSE to track the number of divisions (so set all cells' Ndiv variable to
# 0 at this point), expanding for a second amount of time (delt_2), sorting
# K cells into W wells (96, 364...), expanding these isolated populations
# for another period of time (delt_3) before finally comparing the distribution
# of resistance amongst all of the wells. Also aim to compare the CFSE signal
# at the time of splitting single cells into wells (which can quantify the
# number of divisions up to a maximum of 8) with the resistance identity of the
# wells post the final expansions (delt_3).
# Decided to filter the cells after the second expansion (delt_2) by two
# chosen number of divisions (div_1, and div_2). Then use these to sample K
# cells into W wells - but keeping the div_1 and div_2 populations seperate.
# (try and look at the Poisson distribution of expected
# births beforehand to inform choice of these two integers).


function Run_Exp_Exp_Well_Dist(N0::Int64, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64, R_real::String, lim_probs::Bool,
    delt_1::Float64, POT_N::Int64, delt_2::Float64, div_1::Int64, div_2::Int64,
    K::Int64, W::Int64, delt_3::Float64, nsims::Int64;
    al=0.0::Float64, psi=0.0::Float64)

    cd("Outputs_Well_Dist")

    output_dfs = Array{DataFrame}(undef, 0)

    # Run nsims times.
    for nsim in 1:nsims

        # Seed uniquely barcoded cells.
        init_cells = seed_cells(N0, b, d, p, mu, sig, del, use_lim_probs=lim_probs,
                                psi=psi, al=al)

        # Expand cells for initial expansion time in absence of drug.
        exp_cells = grow_cells(init_cells, delt_1, 10^10, mu, sig, del, psi=psi,
                               al=al, R_real=R_real, drug_presence=0, t_frac=0.025)

        # Sample POT_N of these expanded cells.
        POT_cells = sample(exp_cells.cells, POT_N, replace = false)

        # To simulate the CFSE staining step, set the number of divisions variable
        # for all of these sampled cells to 0.
        for i in 1:length(POT_cells)
            POT_cells[i].Ndiv = 0
        end

        # Expand the POT cells for an additional delt_2 time (this will be 'EP1').
        EP1_cells = grow_cells(POT_cells, delt_2, 10^10, mu, sig, del, psi=psi,
                               al=al, R_real=R_real, drug_presence=0, t_frac=0.025)

        # Decided to simulate filtering cells according to the number of divisions
        # they've undertaken using div_1 and div_2 (assume that i can filter cells
        # for this using the CFSE marker).
        EP1_dv1 = EP1_cells.cells[map(x -> x.Ndiv, EP1_cells.cells) .== div_1]
        EP1_dv2 = EP1_cells.cells[map(x -> x.Ndiv, EP1_cells.cells) .== div_2]

        # Sample K cells W times without replacement.
        K * W < length(EP1_dv1) || error("There are not enough EP1_dv1 cells to split into the chosen number of wells. Need at least K * W.")
        K * W < length(EP1_dv2) || error("There are not enough EP1_dv2 cells to split into the chosen number of wells. Need at least K * W.")

        # Sample all the cells at once, and then reshape into replicate wells.

        well_dv1_cells = sample(EP1_dv1, (K * W), replace = false)
        well_dv1_cells = reshape(well_dv1_cells, (K, W))

        well_dv2_cells = sample(EP1_dv2, (K * W), replace = false)
        well_dv2_cells = reshape(well_dv2_cells, (K, W))

        exp_dv1_cells = Array{Array{CancerCell}}(undef, 0)
        exp_dv2_cells = Array{Array{CancerCell}}(undef, 0)

        for i in 1:W
            temp_exp_dv1 = grow_cells(well_dv1_cells[:,i], delt_3, 10^10, mu, sig, del, psi=psi, al=al, R_real=R_real).cells
            temp_exp_dv2 = grow_cells(well_dv2_cells[:,i], delt_3, 10^10, mu, sig, del, psi=psi, al=al, R_real=R_real).cells
            push!(exp_dv1_cells, temp_exp_dv1)
            push!(exp_dv2_cells, temp_exp_dv2)
        end

        # Finally, i'm interested in:
        # i) the distribution of resistance amongst the wells following this final
        # expansion step, and
        # ii) if there is a difference between the group that was filtered for
        # lower proliferative capacity (dv1) vs those for higher (dv2).
        nR_dv1s = []
        nR_dv2s = []
        for i in 1:W
            # Number of resistant cells in first group
            nR_dv1 = sum(map(x -> x.R, exp_dv1_cells[i,:][1]))
            push!(nR_dv1s, nR_dv1)
            # and the second group
            nR_dv2 = sum(map(x -> x.R, exp_dv2_cells[i,:][1]))
            push!(nR_dv2s, nR_dv2)
        end

        # Save as dataframe.
        nR_df = DataFrame(nR_dv1 = nR_dv1s, nR_dv2 = nR_dv2s, nsim = nsim,
        N0 = N0, b = b, d = d, p = p, mu = mu, sig = sig, del = del,
        R_real = R_real, lim_probs = Int64(lim_probs), delt_1 = delt_1,
        POT_N = POT_N, delt_2 = delt_2, div_1 = div_1, div_2 = div_2,
        K = K, W = W, delt_3 = delt_3, psi = psi)

        push!(output_dfs, nR_df)

    end

    # Join all of the output dataframes.
    output_df = vcat(output_dfs...)

    # Save df name according to params.
    df_name = string("CBC_Exp_Well_Dist_N0-", N0,
    "_b-", round(b, digits=3), "_d-", round(d, digits = 3),
    "_p-", round(p, digits = 10), "_mu-", round(mu, digits = 10),
    "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
    "_R_real-", R_real, "_lp-", Int64(lim_probs), "_delt_1-", delt_1,
    "_POT_N-", POT_N, "_delt_2-", delt_2, "_div_1-", div_1, "_div_2-", div_2,
    "_K-", K, "_W-", W, "_delt_3-", delt_3, "_nsims-", nsims, "_psi-", psi,
    "_output.csv")


    @rput output_df; @rput df_name
    R"""
    output_df$nR_dv1 <- as.numeric(output_df$nR_dv1)
    output_df$nR_dv2 <- as.numeric(output_df$nR_dv2)
    write.csv(output_df, df_name, row.names = F)
    """

end



################################################################################

# Below are experimental functions for recreating potential experiments following
# the full QR in vitro experiment. Experiment numbers correspond to the numbers
# they were concieved in. Not all experiments are being simulated, so the
# numbers here are incomplete.


# Experiment 5:
# Taking CO and DT samples (so either R assigned using eq freq, or R = 1.0,
# respectively), expand for delt_1, and seeding into 3x groups, each consisting
# of W wells with K cells/well.
# Then run for 3x different recovery times: del_C1, del_C2, del_C3
# After each of these expansions, record
#  i) the total number of cells/well (if del > 0.0, potentially statistically
#    lower for the resistant wells).
# ii) the total number of resistant cells/well (if sig > 0.0, expect there to
#     be fewer resistant cells, as cells transition back to sensitive with
#     probability sig/division).

function Run_Exp_Experiment_5(N0::Int64, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64, R_real::String, delt_1::Float64,
    K::Int64, W::Int64, del_C1::Float64, del_C2::Float64, del_C3::Float64,
    nsims::Int64;
    al=0.0::Float64, psi=0.0::Float64, lim_probs=true::Bool)

    cd("Outputs_Experiment_5")

    output_dfs = Array{DataFrame}(undef, 0)

    for nsim in 1:nsims

        # Seed uniquely barcoded cells.
        # Assign the control cells using the equilibrium frequncies.
        CO_init_cells = seed_cells(N0, b, d, p, mu, sig, del, use_lim_probs=lim_probs,
                                psi=psi, al=al)
        DT_init_cells = seed_cells(N0, b, d, p, mu, sig, del, use_lim_probs=lim_probs,
                                psi=psi, al=al)

        # Assuming that all the DT cells have reached resistance.
        for i in 1:length(DT_init_cells)
            DT_init_cells[i].R = 1.0
        end

        # Expand cells for initial expansion time in absence of drug.
        CO_exp_cells = grow_cells(CO_init_cells, delt_1, 10^10, mu, sig, del, psi=psi,
                               al=al, R_real=R_real, drug_presence=0, t_frac=0.025)

        DT_exp_cells = grow_cells(DT_init_cells, delt_1, 10^10, mu, sig, del, psi=psi,
                              al=al, R_real=R_real, drug_presence=0, t_frac=0.025)

        # Sample K cells W times without replacement for each 3x del_Ci groups.
        K * W * 3 < length(CO_exp_cells.cells) || error("There are not enough CO_exp_cells cells to split into the chosen number of wells. Need at least K * W * 3.")
        K * W * 3 < length(DT_exp_cells.cells) || error("There are not enough DT_exp_cells cells to split into the chosen number of wells. Need at least K * W * 3.")

        # Sample all the cells at once, and then reshape into replicate wells.

        well_CO_cells = sample(CO_exp_cells.cells, (K * W * 3), replace = false)
        well_CO_cells = reshape(well_CO_cells, (K, W, 3))

        well_DT_cells = sample(DT_exp_cells.cells, (K * W * 3), replace = false)
        well_DT_cells = reshape(well_DT_cells, (K, W, 3))

        CO_exp_wells = Array{Array{CancerCell}}(undef, 0)
        DT_exp_wells = Array{Array{CancerCell}}(undef, 0)

        del_Cs = [del_C1, del_C2, del_C3]

        for j in 1:3
            for i in 1:W

                temp_CO = grow_cells(well_CO_cells[:,i,j], del_Cs[j], 10^10, mu, sig, del, psi=psi, al=al, R_real=R_real).cells
                temp_DT = grow_cells(well_DT_cells[:,i,j], del_Cs[j], 10^10, mu, sig, del, psi=psi, al=al, R_real=R_real).cells

                push!(CO_exp_wells, temp_CO)
                push!(DT_exp_wells, temp_DT)

            end
        end

        # I'm just interested in the total number of cells, and the total number of
        # resistant cells per well. In practice would reveal the resistant phenotype
        # distribution with a drug assay.

        del_Ci_N_CO = map(x -> length(x), CO_exp_wells)
        del_Ci_N_DT = map(x -> length(x), DT_exp_wells)

        del_Ci_NR_CO = map(x -> sum(map(y -> y.R, x)), CO_exp_wells)
        del_Ci_NR_DT = map(x -> sum(map(y -> y.R, x)), DT_exp_wells)

        # Save to a master_df
        temp_df = DataFrame(CO_N = del_Ci_N_CO, DT_N = del_Ci_N_DT,
                           CO_NR = del_Ci_NR_CO, DT_NR = del_Ci_NR_DT,
                           del_Ci = repeat(del_Cs, inner = W),
                           nsim = nsim, N0 = N0, b = b, d = d, p = p, mu = mu,
                           sig = sig, del = del, R_real = R_real,
                           delt_1 = delt_1, K = K, W = W, del_C1 = del_C1,
                           del_C2 = del_C2, del_C3 = del_C3, al = al,
                           psi = psi, lp = Int64(lim_probs))

        push!(output_dfs, temp_df)

    end

    # Merge output dfs

    output_df = vcat(output_dfs...)

    # Save df name according to params.
    df_name = string("CBC_Exp_experiment_5_N0-", N0,
    "_b-", round(b, digits=3), "_d-", round(d, digits = 3),
    "_p-", round(p, digits = 10), "_mu-", round(mu, digits = 10),
    "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
    "_R_real-", R_real, "_delt_1-", delt_1, "_K-", K, "_W-", W,
    "_del_C1-", del_C1, "_del_C2-", del_C2, "_del_C3-", del_C3, "_nsims-", nsims,
    "_output.csv")


    @rput output_df; @rput df_name
    R"""
    output_df$CO_N <- as.numeric(output_df$CO_N)
    output_df$DT_N <- as.numeric(output_df$DT_N)

    output_df$CO_NR <- as.numeric(output_df$CO_NR)
    output_df$DT_NR <- as.numeric(output_df$DT_NR)

    write.csv(output_df, df_name, row.names = F)
    """

end



# Experiment 6:
# Seed N0 cells. Expand for delt_1 - these produce the POT samples.
# Then sample POT_N cells, set all n_divisions to 0 and expand for delt_2
# - split into a 'high proliferation' (HP) and 'low proliferation' (LP) groups
# based on the number of divisions they've experienced.
# Then want to sample 4* groups of W cells, with K cells per well per each
# proliferation group (so 8* groups in total).
# Set an Wmax for the carrying capacity of each well (around 40,000 for 96-well)
# Grow a pair of HP-LP wells as the control - until Wmax.
# For the remaining three pairs, expose to pulse drug treatment:
# del_D1, del_D2 and del_D3 dictate the time to grow the three pairs for (unless
# Wmax is reached for).
# np_D1, np_D2 and np_D3 dictate the number of drug-pulses experienced by each
# of the remaining pair - need to ensure each del_Di/np_Di is equal: f
# for example would do del_D1 = 21, np_D1 = 7, del_D2 = 36, np_D2 = 12 and
# del_D3 = 51, np_D3 = 17.

# At this stage, want to save the following:
#   i) The number of cells per well
#  ii) The number of resistant cells per well
# iii) The number of escape mutations per well (when alpha > 0.0)
#  iv) The cell number over time until np_Di or Wmax (this be too many data
#      points?).

# Now select all of the wells that haven't gone extinct - then take K2 (if
# cell number > K2) and grow for delt_3 in the absence of drug (so no
# drug-killing). Make sure to keep track of which replicate is which (LP/HP and
# CO, or np_Di).
# Now, finally, save the following:
#   i) The number of cells per well
#  ii) The number of resistant cells per well
# iii) The number of escape mutations per well (when alpha > 0.0)

function Run_Exp_Experiment_6(N0::Int64, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64, al::Float64, R_real::String,
    lim_probs::Bool, psi::Float64, delt_1::Float64, POT_N::Int64,
    delt_2::Float64, K::Int64, W::Int64, Wmax::Int64,
    del_D1::Float64, del_D2::Float64, del_D3::Float64,
    np_D1::Int64, np_D2::Int64, np_D3::Int64, K2::Int64, delt_3::Float64,
    nsims::Int64;
    insta_kill=true::Bool)


    cd("Outputs_Experiment_6")

    output_dfs = Array{DataFrame}(undef, 0)

    # Run nsims times.
    for nsim in 1:nsims

        # Seed uniquely barcoded cells.
        init_cells = seed_cells(N0, b, d, p, mu, sig, del, use_lim_probs=lim_probs,
                                psi=psi, al=al)

        # Expand cells for initial expansion time in absence of drug.
        exp_cells = grow_cells(init_cells, delt_1, 10^10, mu, sig, del, psi=psi,
                               al=al, R_real=R_real, drug_presence=0, t_frac=0.025)

        # Sample POT_N of these expanded cells.
        POT_cells = sample(exp_cells.cells, POT_N, replace = false)

        # To simulate the CFSE staining step, set the number of divisions variable
        # for all of these sampled cells to 0.
        for i in 1:length(POT_cells)
            POT_cells[i].Ndiv = 0
        end

        # Expand the POT cells for an additional delt_2 time (this will be 'EP1').
        EP1_cells = grow_cells(POT_cells, delt_2, 10^10, mu, sig, del, psi=psi,
                               al=al, R_real=R_real, drug_presence=0, t_frac=0.025)

        # Number of births experienced by cells should ~Pois with rate (2bt),
        # therefore split LP and HP around this point.
        EP1_HP = EP1_cells.cells[map(x -> x.Ndiv, EP1_cells.cells) .>= (2 * b * delt_2)]
        EP1_LP = EP1_cells.cells[map(x -> x.Ndiv, EP1_cells.cells) .<  (2 * b * delt_2)]

        # Sample K cells W times without replacement for each 3x del_Di groups
        # and 1x CO group (4* groups in total).
        K * W * 4 < length(EP1_HP) || error("There are not enough EP1_HP cells to split into the chosen number of wells. Need at least K * W * 3.")
        K * W * 4 < length(EP1_LP) || error("There are not enough EP1_LP cells to split into the chosen number of wells. Need at least K * W * 3.")

        # Sample all the cells at once, and then reshape into replicate wells.

        well_EP1_HP = sample(EP1_HP, (K * W * 4), replace = false)
        well_EP1_HP = reshape(well_EP1_HP, (K, W, 4))

        well_EP1_LP = sample(EP1_LP, (K * W * 4), replace = false)
        well_EP1_LP = reshape(well_EP1_LP, (K, W, 4))

        HP_DT_out = Array{Array{CancerCell}}(undef, 0)
        LP_DT_out = Array{Array{CancerCell}}(undef, 0)

        HP_CO_out = Array{Array{CancerCell}}(undef, 0)
        LP_CO_out = Array{Array{CancerCell}}(undef, 0)

        prolif_groups = ["HP", "LP"]
        del_Dis = [del_D1, del_D2, del_D3]
        np_Dis = [np_D1, np_D2, np_D3]

        Di_HP_outs = Array{Grow_Kill_Rec_Out}(undef, 0)
        Di_LP_outs = Array{Grow_Kill_Rec_Out}(undef, 0)
        CO_HP_outs = Array{Grow_Kill_Rec_Out}(undef, 0)
        CO_LP_outs = Array{Grow_Kill_Rec_Out}(undef, 0)

        # First do the del_Di drug-treatments:
        for j in 1:3
            for i in 1:W

                temp_Di_HP_out = grow_kill_rec_cells(well_EP1_HP[:, i, j],
                del_Dis[j], mu, sig, del, np_Dis[j], Wmax, R_real=R_real, psi=psi,
                al=al, drug_kill=true, insta_kill=insta_kill,
                rep_name = string("D", j, "_HP_", i))

                push!(Di_HP_outs, temp_Di_HP_out)

                temp_Di_LP_out = grow_kill_rec_cells(well_EP1_LP[:, i, j],
                del_Dis[j], mu, sig, del, np_Dis[j], Wmax, R_real=R_real, psi=psi,
                al=al, drug_kill=true, insta_kill=insta_kill,
                rep_name = string("D", j, "_LP_", i))

                push!(Di_LP_outs, temp_Di_LP_out)

            end
        end
        # And now the CO treatment:
        for j in 4
            for i in 1:W

                temp_C_HP_out = grow_kill_rec_cells(well_EP1_HP[:, i, j],
                1000.0, mu, sig, del, 1000, Wmax, R_real=R_real, psi=psi,
                al=al,
                rep_name = string("C_HP_", i))

                push!(CO_HP_outs, temp_C_HP_out)

                temp_C_LP_out = grow_kill_rec_cells(well_EP1_LP[:, i, j],
                1000.0, mu, sig, del, 1000, Wmax, R_real=R_real, psi=psi,
                al=al,
                rep_name = string("C_LP_", i))

                push!(CO_LP_outs, temp_C_LP_out)

            end
        end

        # Want to save the number of cells, the number of resistant cells, and
        # the number of escape mutation cells per well, per condition.

        out_dfs_1 = Array{DataFrame}(undef, 0)

        for i in 1:(3*W)

            temp_df = DataFrame(N = length(Di_HP_outs[i].cells),
                                NR = sum(map(x -> x.R, Di_HP_outs[i].cells)),
                                NE = sum(map(x -> x.E, Di_HP_outs[i].cells)),
                                fin_t = Di_HP_outs[i].fin_t,
                                samp = Di_HP_outs[i].rep_name)

            push!(out_dfs_1, temp_df)

            temp_df = DataFrame(N = length(Di_LP_outs[i].cells),
                                NR = sum(map(x -> x.R, Di_LP_outs[i].cells)),
                                NE = sum(map(x -> x.E, Di_LP_outs[i].cells)),
                                fin_t = Di_LP_outs[i].fin_t,
                                samp = Di_LP_outs[i].rep_name)

            push!(out_dfs_1, temp_df)

        end

        # and for the control treatments....

        for i in 1:W

            temp_df = DataFrame(N = length(CO_HP_outs[i].cells),
                                NR = sum(map(x -> x.R, CO_HP_outs[i].cells)),
                                NE = sum(map(x -> x.E, CO_HP_outs[i].cells)),
                                fin_t = CO_HP_outs[i].fin_t,
                                samp = CO_HP_outs[i].rep_name)

            push!(out_dfs_1, temp_df)

            temp_df = DataFrame(N = length(CO_LP_outs[i].cells),
                                NR = sum(map(x -> x.R, CO_LP_outs[i].cells)),
                                NE = sum(map(x -> x.E, CO_LP_outs[i].cells)),
                                fin_t = CO_LP_outs[i].fin_t,
                                samp = CO_LP_outs[i].rep_name)

            push!(out_dfs_1, temp_df)

        end

        # And combine...

        out_df_1 = vcat(out_dfs_1...)

        # And now take K2 cells from each well and grow them all in the absence
        # of drug for delt_3 time (or until Wmax is reached).

        out_dfs_2 = Array{DataFrame}(undef, 0)

        for i in 1:(3*W)

            # First, need to account for if the cells drifted to extinction.
            if length(Di_HP_outs[i].cells) == 0

                temp_df = DataFrame(N = 0,
                                    NR = 0,
                                    NE = 0,
                                    fin_t = 0.0,
                                    samp = Di_HP_outs[i].rep_name)

            else

                if length(Di_HP_outs[i].cells) > K2
                    temp_cells = sample(Di_HP_outs[i].cells,
                    K2, replace = false)
                else
                    temp_cells = Di_HP_outs[i].cells
                end

                temp_cells = grow_cells(temp_cells, delt_3, Wmax, mu, sig, del,
                psi=psi, al=al, R_real=R_real)

                temp_df = DataFrame(N = length(temp_cells.cells),
                                    NR = sum(map(x -> x.R, temp_cells.cells)),
                                    NE = sum(map(x -> x.E, temp_cells.cells)),
                                    fin_t = temp_cells.fin_t,
                                    samp = Di_HP_outs[i].rep_name)

            end

            push!(out_dfs_2, temp_df)

            if length(Di_LP_outs[i].cells) == 0

                temp_df = DataFrame(N = 0,
                                    NR = 0,
                                    NE = 0,
                                    fin_t = 0.0,
                                    samp = Di_LP_outs[i].rep_name)

            else

                if length(Di_LP_outs[i].cells) > K2
                    temp_cells = sample(Di_LP_outs[i].cells,
                    K2, replace = false)
                else
                    temp_cells = Di_LP_outs[i].cells
                end

                temp_cells = grow_cells(temp_cells, delt_3, Wmax, mu, sig, del,
                psi=psi, al=al, R_real=R_real)

                temp_df = DataFrame(N = length(temp_cells.cells),
                                    NR = sum(map(x -> x.R, temp_cells.cells)),
                                    NE = sum(map(x -> x.E, temp_cells.cells)),
                                    fin_t = temp_cells.fin_t,
                                    samp = Di_LP_outs[i].rep_name)

            end

            push!(out_dfs_2, temp_df)

        end

        # Do the same for the control samples...

        for i in 1:W

            if length(CO_HP_outs[i].cells) > K2
                temp_cells = sample(CO_HP_outs[i].cells,
                K2, replace = false)
            else
                temp_cells = CO_HP_outs[i].cells
            end

            temp_cells = grow_cells(temp_cells, delt_3, Wmax, mu, sig, del,
            psi=psi, al=al, R_real=R_real)

            temp_df = DataFrame(N = length(temp_cells.cells),
                                NR = sum(map(x -> x.R, temp_cells.cells)),
                                NE = sum(map(x -> x.E, temp_cells.cells)),
                                fin_t = temp_cells.fin_t,
                                samp = CO_HP_outs[i].rep_name)

            push!(out_dfs_2, temp_df)

            if length(CO_LP_outs[i].cells) > K2
                temp_cells = sample(CO_LP_outs[i].cells,
                K2, replace = false)
            else
                temp_cells = CO_LP_outs[i].cells
            end

            temp_cells = grow_cells(temp_cells, delt_3, Wmax, mu, sig, del,
            psi=psi, al=al, R_real=R_real)

            temp_df = DataFrame(N = length(temp_cells.cells),
                                NR = sum(map(x -> x.R, temp_cells.cells)),
                                NE = sum(map(x -> x.E, temp_cells.cells)),
                                fin_t = temp_cells.fin_t,
                                samp = CO_LP_outs[i].rep_name)

            push!(out_dfs_2, temp_df)

        end

        out_df_2 = vcat(out_dfs_2...)

        insertcols!(out_df_1, ncol(out_df_1)+1, :out => 1)
        insertcols!(out_df_2, ncol(out_df_2)+1, :out => 2)

        out_df = vcat(out_df_1, out_df_2)
        insertcols!(out_df, ncol(out_df)+1, :nsim => nsim)
        push!(output_dfs, out_df)

    end

    # Merge output dfs

    output_df = vcat(output_dfs...)

    param_df = DataFrame(N0 = N0, b = b, d = d, p = p, mu = mu, sig = sig,
    del = del, al = al, R_real = R_real, lp = Int64(lim_probs), psi = psi,
    delt_1 = delt_1, POT_N = POT_N, delt_2 = delt_2, K = K, W = W, Wmax = Wmax,
    del_D1 = del_D1, del_D2 = del_D2, del_D3 = del_D3, np_D1 = np_D1,
    np_D2 = np_D2, np_D3 = np_D3, K2 = K2, delt_3 = delt_3, ik = insta_kill)

    output_df = hcat(output_df, repeat(param_df, nrow(output_df)))

    # Save df name according to params. n.b. that the file name is too long
    # if i try and save all of them in the file name - so just save the main
    # simulation variable parameters, and the rest can be saved in the output
    # dataframes.
    df_name = string("CBC_Exp_experiment_6_N0-", N0,
    "_b-", round(b, digits=3), "_d-", round(d, digits = 3),
    "_p-", round(p, digits = 10), "_mu-", round(mu, digits = 10),
    "_sig-", round(sig, digits = 10), "_del-", round(del, digits = 10),
    "_al-", round(al, digits = 10), "_R_real-", R_real,
    "_lp-", Int64(lim_probs), "_psi-", round(psi, digits = 10),
    "_nsims-", nsims,
    "_output.csv")

    @rput output_df; @rput df_name
    R"""
    output_df$N <- as.numeric(output_df$N)
    output_df$NR <- as.numeric(output_df$NR)
    output_df$NE <- as.numeric(output_df$NE)

    write.csv(output_df, df_name, row.names = F)
    """


end
