
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

function Run_Exp(N::Int64, b::Float64, d::Float64, t_exp::Float64,
    p::Float64, mu::Float64, sig::Float64, del::Float64,
    R_real::String, n_pulse::Int64, Nmax::Int64, N_seed::Int64,
    t_CO::Float64, t_DT::Float64, Passage::Int64,
    insta_kill::Bool, lim_probs::Bool;
    al=0.0::Float64, psi=0.0::Float64,
    Exp_Input::Experiment_Input=Experiment_Input())

    if Passage == 1

        # First expand and split the cells for t = t_exp
        exp_split_cells = expand_split_cells(N, b, d, t_exp, N_seed, 4, p, mu, sig, del, psi=psi, al=al, R_real=R_real, use_lim_probs=lim_probs)


        Exp_Input.CO_inputs = exp_split_cells.CO_flasks
        Exp_Input.DT_inputs = exp_split_cells.DT_flasks
    end

    CO_Outputs = Array{Grow_Kill_Rec_Out}(undef, 4)
    DT_Outputs = Array{Grow_Kill_Rec_Out}(undef, 4)

    # Run the bulk of the simulation on the threads. All scheduling done
    # automatically by @threads.
    Threads.@threads for i in 1:4
        CO_Outputs[i] = grow_kill_rec_cells(Exp_Input.CO_inputs[:,i], t_CO, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, rep_name = string("CO", i, "P", Passage), must_Nmax=true)
        DT_Outputs[i] = grow_kill_rec_cells(Exp_Input.DT_inputs[:,i], t_DT, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, drug_kill=true, insta_kill=insta_kill, rep_name = string("DT", i, "P", Passage), must_Nmax=true)
    end

    return Experiment_Output(CO_Outputs, DT_Outputs)

end


# Run the experiment and save the output.

function Run_Exp_save_output(N::Int64, b::Float64, d::Float64, t_exp::Float64,
    p::Float64, mu::Float64, sig::Float64, del::Float64,
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
    dir_name = string("CBC_Exp_out_N-", N,
    "_b-", round(b, digits=3), "_d-", round(d, digits = 3),
    "_t_exp-", round(t_exp, digits = 3), "_p-", round(p, digits = 10),
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

    # Create vectors to hold output dataframes until all Passages have run.
    bc_R_dfs = Array{DataFrame}(undef, 0)
    N_dfs = Array{DataFrame}(undef, 0)
    Ps = Array{Int64}(undef, 0)

    for P in 1:Passage
        if P == 1
            exp_out = Run_Exp(N, b, d, t_exp, p, mu, sig, del, R_real, n_pulse,
            Nmax, N_seed, t_CO, t_DT, P, insta_kill, lim_probs,
            psi=psi, al=al)
        else

            exp_out = Run_Exp(N, b, d, t_exp, p, mu, sig, del, R_real, n_pulse,
            Nmax, N_seed, t_CO, t_DT, P, insta_kill, lim_probs,
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
        # Saving the dataframes until all Passages have run.
        push!(bc_R_dfs, bc_R_df)
        push!(N_dfs, N_df)
        push!(Ps, P)

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

    # Now save the outputs:
    for i in 1:length(Ps)

        bc_R_df = bc_R_dfs[i]
        N_df = N_dfs[i]
        Pass = Ps[i]

        @rput bc_R_df; @rput N_df; @rput Pass
        R"""
        write.csv(bc_R_df, file=paste0("bc_counts_tot_mean_R_P", Pass, ".csv"), row.names = F)
        #write.csv(NRE_df, file=paste0("NRE_by_t_P", Pass, ".csv"), row.names = F)
        write.csv(N_df, file=paste0("N_by_t_P", Pass, ".csv"), row.names = F)
        """

    end

    cd("../../../")

end


################################################################################

# Track a subset of lineages through the grow_kill_rec function several times
# to illustrate how the cell trajectories differ over time in the simulation
# when:
# psi = 0.0 (drug-killing is deterministic)
# psi > 0.0 (drug-killing has a stochastic component)


function Run_Exp_Sub_Lins(N::Int, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64, use_lim_probs::Bool, n_lin::Int,
    tmax_1::Float64, Nmax_1::Int, tmax_2::Float64, Nmax_2::Int, n_recs::Int,
    n_pulse::Int, psi::Float64)

    init_cells = seed_cells(N, b, d, p, mu, sig, del, use_lim_probs=use_lim_probs)

    # Choose n_lin lineages that are resistant (R = 1.0) and n_lin that are
    # sensitive.

    Rlins = sample(map(x -> x.barcode, init_cells[map(x -> x.R == 1.0, init_cells)]),
        n_lin, replace=false)

    Slins = sample(map(x -> x.barcode, init_cells[map(x -> x.R == 0.0, init_cells)]),
        n_lin, replace=false)

    # Expand these cells:

    output = grow_cells(init_cells, tmax_1, Nmax_1, mu, sig, del)

    # Store the count of each resistant and sensitive lineages in seperate vectors.

    Rlin_vec = Array{Array{Int64}}(undef, 0)
    Slin_vec = Array{Array{Int64}}(undef, 0)

    push!(Rlin_vec, map(x -> single_bc_count(output.cells, x), Rlins))
    push!(Slin_vec, map(x -> single_bc_count(output.cells, x), Slins))

    trec = tmax_2/n_recs
    dt_t = tmax_2/n_pulse
    curr_dt_t = dt_t
    curr_t = 0.0
    tvec = [curr_t]

    while curr_t <= tmax_2

        output = grow_cells(output.cells, trec, Nmax_2, mu, sig, del)

        push!(Rlin_vec, map(x -> single_bc_count(output.cells, x), Rlins))
        push!(Slin_vec, map(x -> single_bc_count(output.cells, x), Slins))

        curr_t += trec
        push!(tvec, curr_t)

        if curr_t >= curr_dt_t && curr_dt_t < tmax_2
            curr_dt_t += dt_t

            # Kill section:
            kill_vec = Int64[]
            for i in 1:length(output.cells)
                # Random unif to determine kill.
                ran_k = rand(Uniform(0, 1))
                # if resistant...
                if output.cells[i].R > 0.0
                    if ran_k >= (1.0 - psi)
                        push!(kill_vec, i)
                    end
                # if sensitive...
                elseif output.cells[i].R == 0.0
                    if ran_k >= (0.0 + psi)
                        push!(kill_vec, i)
                    end
                end
            end

            # Now kill the cells according to kill_vec
            deleteat!(output.cells, kill_vec)
            # And save updated lineage counts.
            push!(Rlin_vec, map(x -> single_bc_count(output.cells, x), Rlins))
            push!(Slin_vec, map(x -> single_bc_count(output.cells, x), Slins))
            push!(tvec, curr_t)

        end

    end

    # Merge the resistant and sensitive lienage vectors

    Rlin_df = DataFrame(hcat(Rlin_vec...))
    Slin_df = DataFrame(hcat(Slin_vec...))

    @rput Rlin_df; @rput Slin_df; @rput tvec; @rput psi;

    R"""
    library(tidyr)
    library(ggplot2)
    library(grid)
    library(ggthemes)

    Rlin_df["lin"] <- 1:nrow(Rlin_df)
    Slin_df["lin"] <- 1:nrow(Slin_df)
    Rlin_df <- gather(Rlin_df, key = "time", value = "count",
                    grep("lin", colnames(Rlin_df), invert = T))
    Rlin_df$time <- rep(tvec, each = nrow(Rlin_df)/length(tvec))
    Rlin_df["Phenotype"] <- "Resistant"
    Slin_df <- gather(Slin_df, key = "time", value = "count",
                    grep("lin", colnames(Slin_df), invert = T))
    Slin_df$time <- rep(tvec, each = nrow(Slin_df)/length(tvec))
    Slin_df["Phenotype"] <- "Sensitive"
    lin_df <- rbind(Rlin_df, Slin_df)
    if(psi == 0.0){
        plot_title <- "\U03A8 = 0.0"
    } else {
        plot_title <- "\U03A8 = 0.1"
    }

    p <- ggplot(data = lin_df, aes(x = time, y = count, colour = Phenotype,
                group = interaction(lin, Phenotype))) +
        geom_line(size = 1.0, alpha = 0.6) +
        scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        annotation_logticks(base = 10, sides = "l") +
        xlab("\nTime") +
        ylab("Lineage Cell Count\n") +
        theme_minimal() +
        theme(panel.background = element_rect(colour = "grey80"),
            text = element_text(size = 18), element_line(size = 0.8),
            panel.border = element_rect(colour = "black", fill=NA, size=2)
        ) +
        scale_colour_manual(values = c("blue1", "orangered")) +
        ggtitle(plot_title)

    setwd("Plots")
    ggsave(paste0("Chosen_lin_example_psi-", psi, ".jpeg"), p,
    width = 10, height = 5)
    setwd("../")
    """


end
