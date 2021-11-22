
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Functions

################################################################################

###########################
# Growth Functions
###########################

# Description of this version of the Cancer Barcode Simulation (CBC).

# Seed N cells with given birth (b) and death (d) rates, and resistance (R)
# scores all initially as 0 (sensitive).

# Pre-existing Resistance.
# p - probability of cells being resistant (R = 1).
# bc_lib and use_lib are optional arguments - if not proivded the function
# will default to assigning barcodes that are unitue for 1:N.

# De-novo Resistance Mutation and Plasticity.
# μ (mu) - probability of resistance 'mutation' (R + 1) per cell-division.
# σ (sig) - the probability of a sensitive 'mutation' (R - 1) per cell-division.
# ... hence sig must be > 0. mu always refers to the probability of (R += 1).
# plastic - (Bool) if false, proceeds as simple mutation sim, where R increases
# by Poisson(mu) per division. If true, R is only ever 0 or 1, with the
# direction of movement conditional on its current state, and a mutation
# only ever leading to a change of +/-1.

# Implementing the Cost of Resistance:
# Cells are grown with their underlying birth and death rates.
# When cells have an R score > 0, these cells grow with:
# a fitness disadvantage in the absence of 'drug' (drug_presence = 0).
# a fitness advantage in the presence of 'drug' (drug_presence = 1).
# 'R_real' - dictates whether these differences due to R are realised as
# b - increase or decreases in birth rate, or
# d -    "     "      "     "  death rate.
# A cell's R score mutates with probability μ (mu), per division.
# The value of R keeps track of the number of resistance mutations.
# A second value, δ (del), quantifies the relative fitness cost/benefit of the mutation.
# The wild-type (no R mutations) has a fitness of 1.0
# In the absence of drug, the resistant cells (R > 0) have fitness of 1.0 - δ
# In the presence of drug, the resistant cells (R > 0) have fitness 1.0 + δ
# Whilst this framework can be implemented, currently, the cost is permanently
# invoked, and cells are killed instantaneously, with probability 0.0 if R > 0,
# or 1.0 if R = 0.
# Cells are grown either until t reaches tmax or N reaches Nmax.
# By default, the cost of resistance is realised in the birth rate (R_real = b).

# Implementing variability in the resistance phenotype:
# If al == be == 0.0, the simulation reverts to the binary resistance phenotype
# outlined above. If al != 0.0 | be != 0.0, the resistance phenotype is now
# drawn from a Beta distribution: R = Beta(al, be). The probability of death
# during the drug-kill step is still (1 - R).
# 'Mutation' means a sensitive cell gains a resistant phenotype, with
# probability μ (mu) per division - by drawing from the Beta(al, be).
# The cost of resistance - δ (del) - still functions on a binary basis:
# resistant cells (R > 0.0) incurr the cost, whereas sensitive (R = 0.0) do not.
# The 'plasticity' functionality - σ (sig) - still controls the rate at which
# cells transition from resistant (R > 0.0) to sensitive.


############################
# Custom Library Function.
############################

# Load a preset barcode library. Takes an array of floats - can be used to
# introduce bias in initial plasmid pool distribution - could be determined,
# for example, by sequencing to high depth.

function load_lib(bcs)
    bc_lib = BarcodeLibrary(bcs)
    return bc_lib
end

#############################################
# Resistance Equilibrium Proportion Function.
#############################################

# Given values of mu, sig and del, and the sensitive birth and death rates
# (br and dr), return the analytical solution to the stable proportion of
# resistance.

function stable_R_pheno(mu::Float64, sig::Float64, del::Float64,
    br::Float64, dr::Float64)
    @rput mu; @rput sig; @rput del; @rput br; @rput dr;
    R"""
    # Constructing Quadratic Formula
    result <- function(a,b,c){
        if(delta(a,b,c) > 0){ # first case D>0
            x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
            x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
            result = c(x_1,x_2)
        }
        else if(delta(a,b,c) == 0){ # second case D=0
            x = -b/(2*a)
        }
        else {"There are no real roots."} # third case D<0
    }

    # Constructing delta
    delta<-function(a,b,c){
        b^2-4*a*c
    }

    if(mu == 0.0 & sig == 0.0 & del == 0.0){
        p <- 0.0
    } else if(del == 0.0){
        p = mu/(mu + sig)
    } else {
        lam_r <- (br - dr)
        lam_R <- lam_r - (lam_r * del)
        # Assume splitting cost between birth and death rates.
        # i.e. R_real = "l"
        bR <- br * (1-del)

        # Convert into quadratic terms;
        a <- (-lam_R + lam_r)
        b <- (lam_R - (2*mu*br) - (2*sig*bR) - lam_r)
        c <- (2*mu*br)

        # Solve the quadratic;
        p <- result(a, b, c)[[2]]
    }
    p <- round(p, digits = 6)
    """
    @rget p; return p
end



############################
# Seed Cells Function.
############################


function seed_cells(N::Int64, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64;
    use_lim_probs=false::Bool,
    psi=0.0::Float64, al=0.0::Float64,
    bc_lib=BarcodeLibrary(1:N)::BarcodeLibrary, use_lib=false::Bool)

    cells = Array{CancerCell}(undef, N)
    Rs = zeros(N)
    Es = zeros(N)
    # Use lib to decide whether to use a BarcodeLibrary or default to 1:N.
    if use_lib
        bcs = sample(bc_lib.barcodes, N, replace = true)
    else
        bcs = bc_lib.barcodes
    end
    # Assign all the fields to the CancerCell accordingly.
    for i in 1:length(cells)
        cells[i] = CancerCell(i, bcs[i], b, d, Rs[i], Es[i], 0)
    end

    # Assign resistance to the initial cells for the 'pre-existing' scenario
    # using a Poisson r.v. with mean (n cells * p). If using binary resistance
    # phenotype, all = 1.0, otherwise use the Beta(al, be) distribution.
    if p > 0.0
        np = rand(Poisson(p * length(cells)))
        p_cells = sample(1:length(cells), np, replace = false)
        for i in 1:np
            #if res_phen == "bin"
                cells[p_cells[i]].R = 1.0
            #elseif res_phen == "var"
            #    R_phenos = round.(rand(Beta(al, be), np), digits = 4)
            #    cells[p_cells[i]].R = R_phenos[i]
            #end
        end
    end

    # Assign resistance to the initial cells according to the equilibrium
    # probabilities. If using binary resistance phenotype, all = 1.0,
    # otherwise use the Beta(al, be) distribution.
    if use_lim_probs == true
        # Expected proportion of resistance...
        exp_R = stable_R_pheno(mu, sig, del, b, d)
        # As a number of resistant cells...
        nR = Int64(round(exp_R*length(cells)))
        if nR > length(cells)
            nR = length(cells)
        end
        R_cells = sample(1:length(cells), nR, replace = false)
        for i in 1:nR
            #if res_phen == "bin"
                cells[R_cells[i]].R = 1.0
            #elseif res_phen == "var"
            #    R_phenos = round.(rand(Beta(al, be), nR), digits = 4)
            #    cells[R_cells[i]].R = R_phenos[i]
            #end
        end
    end

    # NB always assume that all cells start with 0 'escape mutations'.
    # Therefore all cells' 'E' values remain at 0.0.

    return cells

end


# Opposed to deepcopying the cell following a birth every time, it is much
# quicker to create a seperate copycell function, and call this within
# the birth-death function.

function copycell(orig_cell::CancerCell)
    new_cell::CancerCell = CancerCell(
    copy(orig_cell.cell_ID),
    copy(orig_cell.barcode),
    copy(orig_cell.b),
    copy(orig_cell.d),
    copy(orig_cell.R),
    copy(orig_cell.E),
    copy(orig_cell.Ndiv)
    )
end


############################
# Growth Function.
############################


function grow_cells(cells::Array{CancerCell}, tmax::Float64,
    Nmax::Int64, mu::Float64, sig::Float64, del::Float64;
    psi=0.0::Float64, al=0.0::Float64,
    R_real="b"::String, drug_presence=0::Int64, t_frac=0.050::Float64)

    out_cells = deepcopy.(cells)

    0 <= mu <= 1.0 || error("mu must be between 0 and 1.")
    0 <= sig <= 1.0 || error("sig must be between 0 and 1.")
    0 <= del <= 1.0 || error("del must be between 0 and 1.")
    0 <= al <= 1.0 || error("al must be between 0 and 1.")
    0 <= t_frac <= 1.0 || error("t_frac must be between 0 and 1.")

    length(unique(map(x -> x.b, out_cells))) == 1 || error("For this simulation, cells should have a uniform base birth and death rate.")
    length(unique(map(x -> x.d, out_cells))) == 1 || error("For this simulation, cells should have a uniform base birth and death rate.")
    R_real == "b" || R_real == "d" || R_real == "l" || error("R_real can only take 'b', 'd' or 'l' as a value.")
    drug_presence == 1 || drug_presence == 0 || error("drug_presence must only take either '0' (absence) or 1 (presence) as values.")


    # Maximum cell birth and death rates need to be changed for KMC steps if
    # using the cost of resistance version is being used.

    # if the drug is absent (drug_presence = 0):
    # the maximum death rate is d + (lam*δ) (if R_real = d)
    # if the drug is present (drug_presence = 1):
    # the maximum birth rate is b + (lam*δ) (if R_real = b).
    # whilst the maximum birth and death rates are:
    #                           b * (1+del),
    #                       and d * (1+del), (if R_real = l).


    # For now, regardless of the number of resistance mutations (R), a cell
    # will only incur the benefit/cost as δ.
    # R_real determines if the cost is incurred in the birth or death rates, or
    # shared proportionately amongst both rates; "b", "d" and "l", respectively.

    # Finally, implementing 'jackpot' mutations that occur with rate α (al) per
    # division when a cell is resistant where, if it occurs, a cell's 'E'
    # variable becomes 1.0, and that cell (and ancestors) no longer incur the
    # fitness cost according to δ (del).
    # Following a resistant cell dividing, switching back to sensitive
    # (controlled by σ (sig)) takes precedence over the jackpot mutations - i.e.
    # this probability is calculated first.
    # If a resistant cell with an escape mutation (R > 0.0 & E > 0.0) reverts
    # to sensitive (controlled by σ (sig)), both the resistant and escape
    # status are lost (R = 0.0, E = 0.0).

    bmax = maximum(map(x -> x.b, out_cells))
    dmax = maximum(map(x -> x.d, out_cells))
    # Save the population's pre-cost, net growth rate, λ (lam).
    lam = bmax - dmax
    # Save the normalised birth and death rates for R_real = l,
    #bnorm = bmax/(bmax + dmax)
    #dnorm = dmax/(bmax + dmax)

    # If using the cost of resistance simulation, need to update bmax and dmax
    # as the highest population b/d rates can now be higher if drug present.
    # The cost of resistance is only implemented if δ > 0.0.
    if del > 0.0
        if     drug_presence == 0 && R_real == "d"
            dmax = dmax + (lam*del)
            bdmax = bmax + dmax
        elseif drug_presence == 1 && R_real == "b"
            bmax = bmax + (lam*del)
            bdmax = bmax + dmax
        elseif drug_presence == 0 && R_real == "l"
            # Both are now smaller, so don't need to change either b/dmax.
            bdmax = bmax + dmax
        elseif drug_presence == 1 && R_real == "l"
            bmax = bmax * (1+del)
            dmax = dmax * (1+del)
            bdmax = bmax + dmax
        else
            bdmax = bmax + dmax
        end
    else
        bdmax = bmax + dmax
    end

    # Create vector to hold time, and initiate at 0.
    t = 0.0
    tvec = Float64[t]
    # Record the population size (Nt) every tmax*t_frac.
    t_rec_change = tmax*t_frac
    t_rec = t_rec_change

    # Have a seperate counter that keeps track of the number of live cells,
    # as we can't now use length(out_cells) - includes the dead cells until the
    # end.
    Nt = length(out_cells)
    # Have a vector to save the number of live cells every t_rec.
    Nvec = Int64[Nt]
    # Also have vectors to save the number of cells that have a resistance
    # mutation (Rvec) and the number of cells that have an escape mutation
    # (Evec). Keep track of the number only when an event happens with
    # count vectors to avoid repeatedely using map (expensive).
    Rcount = sum(map(x -> x.R, out_cells))
    Ecount = sum(map(x -> x.E, out_cells))
    Rvec  = Int64[Rcount]
    Evec = Int64[Ecount]

    # Opposed to sampling directly from out_cells, have a 'samp_vec'. Update
    # this following a birth or death accordingly:
    # if a birth, +1 to the length of the vector.
    # if a death, change this index in the sampling vec to 0, and add the index
    # to the 'kill_vec'.
    # Keep sampling from samp_vec until ran_cell_pos > 0
    # This way, don't have to use deleteat! (which is expensive) until the very
    # end.
    samp_vec = collect(1:length(out_cells))
    kill_vec = Array{Int64}(undef, 0)

    # Grow cells until t = tmax

    while t <= tmax

        ran_samp_pos = rand(1:length(samp_vec))
        ran_cell_pos = samp_vec[ran_samp_pos]
        # Keep sampling from samp_vec until hit a +ve number.
        if ran_cell_pos == 0
            while ran_cell_pos == 0
                ran_samp_pos = rand(1:length(samp_vec))
                ran_cell_pos = samp_vec[ran_samp_pos]
            end
        end

        ran = rand(Uniform(0, bdmax))

        # Update the time pre-event, to prevent including events that will
        # have occured post-tmax.

        dt = -1/(bdmax * Nt) * log(rand())
        t += dt

        # If the time-change has crossed t_rec, then update the tvec and
        # Nvec tracking vectors, and update t_rec by t_rec_change.
        if t >= t_rec
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Rvec, Rcount)
            push!(Evec, Ecount)
            t_rec += t_rec_change
        end

        # If the time-change has crossed tmax, then there wasn't time for
        # the birth or event to happen, t is capped at tmax, and the
        # simulation ends.
        if t > tmax
            t = tmax
            break
        end
        # If N is now >= Nmax, also end the simulation, whilst updating the
        # tracking vectors.
        if Nt >= Nmax
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Rvec, Rcount)
            push!(Evec, Ecount)
            break
        end
        # Otherwise, continue with the birth and death steps.

        # NB - the trade-off of resistance is realied as follows,
        # where λ = (b - d).

        # the drug is absent:
            # trade-off realised in birth rate.
                # b = b - λ*δ
                # d = d
            # trade-off realised in death rate.
                # b = b
                # d = d + λ*δ
            # trade-off realised in both rates.
                # b = b * (1-δ)
                # d = d * (1-δ)
        # the drug is present:
            # trade-off realised in birth rate.
                # b = b + λ*δ
                # d = d
            # trade-off realised in death rate.
                # b = b
                # d = d - λ*δ
            # trade-off realised in both rates.
                # b = b * (1+δ)
                # d = d * (1+δ)

        # NB that this means if using the 'l' cost version (both rates), then
        # both birth AND deat rates are lower. It is cell turnover that is
        # either reduced or increased depending on drug-presence.

        # Set cell_b and cell_d using the trade-off if performing the cost of
        # resistance simulation.

        # If the chosen cell has gained an 'escape' mutation (CancerCell
        # variable 'E') then the cell no longer experiences the cost of
        # resistance and the rates simply revert to the 'wild-type' birth
        # and death rates.

        if del > 0.0
            if out_cells[ran_cell_pos].R > 0

                if out_cells[ran_cell_pos].E > 0
                    # Escapes cost of resistance.
                    cell_b = out_cells[ran_cell_pos].b
                    cell_d = out_cells[ran_cell_pos].d
                else
                    # Otherwise, calculate the cost incurred according to 'R_real'.
                    if drug_presence == 0
                        if R_real == "b"
                            cell_b = out_cells[ran_cell_pos].b - (lam*del)
                            cell_d = out_cells[ran_cell_pos].d
                        end
                        if R_real == "d"
                            cell_b = out_cells[ran_cell_pos].b
                            cell_d = out_cells[ran_cell_pos].d + (lam*del)
                        end
                        if R_real == "l"
                            cell_b = out_cells[ran_cell_pos].b * (1-del)
                            cell_d = out_cells[ran_cell_pos].d * (1-del)
                        end
                    end

                    if drug_presence == 1
                        if R_real == "b"
                            cell_b = out_cells[ran_cell_pos].b + (lam*del)
                            cell_d = out_cells[ran_cell_pos].d
                        end
                        if R_real == "d"
                            cell_b = out_cells[ran_cell_pos].b
                            cell_d = out_cells[ran_cell_pos].d - (lam*del)
                        end
                        if R_real == "l"
                            cell_b = out_cells[ran_cell_pos].b * (1+del)
                            cell_d = out_cells[ran_cell_pos].d * (1+del)
                        end
                    end
                end
            else
                # NB If the cell has no resistance mutations, use the 'wild-type'
                # birth and death rates, regardless of drug_presence... doesn't
                # make biological sense as the drug should reduce fitness of the
                # wild-type... but okay for now because we're looking at relative
                # fitness... will have to reconsider if implement drug presence.
                cell_b = out_cells[ran_cell_pos].b
                cell_d = out_cells[ran_cell_pos].d
            end

            # Set to 0 if < 0
            if cell_b < 0
                cell_b = 0
            end
            if cell_d < 0
                cell_d = 0
            end
        # If not using cost of resistance simulation, simply assign rates from
        # the randomly selected cell.
        else
            cell_b = out_cells[ran_cell_pos].b
            cell_d = out_cells[ran_cell_pos].d
        end

        # Birth step
        if ran < cell_b

            # Update the cells number of divisions.
            out_cells[ran_cell_pos].Ndiv += 1
            ran_cell = copycell(out_cells[ran_cell_pos])

            # NB Here, BOTH DAUGHTER CELLS INHERIT THE MUTATION.
            # 'Grey -> Black + Black'

            # Calculate if a transition to or from the resistance phenotype
            # occurs - mu corresponds to the probability of a mutation to
            # resistance, whereas sig corresponds to a mutation from resistance,
            # per division (therefore <= 1.0).
            #if res_phen == "bin"
                ran_cell.R == 0 || ran_cell.R == 1 || error("A cell's R is not 0 or 1 in the binary resistance sim.")
            #end
            #if res_phen == "var"
            #    ran_cell.R >= 0 || error("A cell's R is negative in the variable resistance sim: al | be > 0.0.")
            #end
            ran_p = rand()
            # Draw a SECOND random number to do the 'escape' mutation. NB that
            # if you _don't_ do this, you can never get an escape mutation when
            # sig > al.
            ran_pE = rand()
            # Sensitive to Resistant with rate mu. Update the resistance
            # phenotype according to either the binary or variable model:
            if ran_cell.R == 0.0
                if mu > ran_p
                    #if res_phen == "bin"
                        ran_cell.R += 1
                        out_cells[ran_cell_pos].R += 1
                        # Also update Rcount with number of resistant cells.
                        Rcount += 2
                    #end
                    #if res_phen == "var"
                    #    new_R = round(rand(Beta(al, be)), digits = 4)
                    #    ran_cell.R = new_R
                    #    out_cells[ran_cell_pos].R = new_R
                    #end
                end
            # Resistant to Sensitive with rate sig. Sensitive is always R = 0.0,
            # regardless of the resistance phenotype model type.
            # This switching step takes precedent over the 'escape mutation'
            # rate (al), which occurs if the cell remains resistant with
            # rate al. Again, both daughter cells inherent the 'mutation'.
            elseif ran_cell.R > 0.0
                if sig > ran_p
                    ran_cell.R = 0.0
                    out_cells[ran_cell_pos].R = 0.0
                    # Update R count vector accordingly.
                    Rcount -= 1
                    if ran_cell.E == 1.0
                        # Cell also loses 'escape mutations'.
                        ran_cell.E = 0.0
                        out_cells[ran_cell_pos].E = 0.0
                        # Update E count vector accordingly.
                        Ecount -= 1
                    end
                # If the current cell does not have an escape mutation,
                # now let daughter cells acquire escape with probability al.
                elseif ran_cell.E == 0.0
                    if al > ran_pE
                        ran_cell.E = 1.0
                        out_cells[ran_cell_pos].E = 1.0
                        # Also update Ecount with number of resistant cells.
                        Ecount += 2
                    end
                else
                    # No mutations - update R and E count vectors accordingly.
                    Rcount += 1
                    Ecount += 1
                end
            end

            push!(out_cells, ran_cell)
            # Update samp_vec with the new length of the output cells.
            push!(samp_vec, (length(out_cells)))
            # Update number of cells to normalise dt.
            Nt += 1
        end

        # Death step
        if bmax <= ran < bmax + cell_d
            # Remove this chosen cell's index from the samp_vec.
            samp_vec[ran_samp_pos] = 0
            push!(kill_vec, ran_cell_pos)
            # Update number of cells to normalise dt.
            Nt -= 1
            # Update the R and E count vectors accordingly.
            if out_cells[ran_cell_pos].R == 1.0
                Rcount -= 1
                if out_cells[ran_cell_pos].E == 1.0
                    Ecount -= 1
                end
            end
        end

        # Break if no cells remain. Make use of the fact that:
        # samp_vec = N0 + nbirths
        # kill_vec = ndeaths
        # therefore if length(kill_vec) >= length(samp_vec), there are no cells
        # remaining. Also save this outcome to the record vectors.
        if length(kill_vec) >= length(samp_vec)
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Rvec, Rcount)
            push!(Evec, Ecount)
            break
        end

    end

    # Now perform all the kill steps using the indexes stored in kill_vec.
    deleteat!(out_cells, sort(kill_vec))

    fin_t = round(t, digits = 4)

    return Grow_Out(out_cells, Nvec, tvec, Rvec, Evec, fin_t)

end

###################################
# Grow, Kill Record Cells Function.
###################################


# Function that grows cells for a given time, tmax, using the kmc method.
# In this Binary version of Resistance Evolution, cells are either completely
# resistant (R >= 0), or sensitive (R = 0).
# The total time, tmax, is split into n_pulse equal periods - when time passes
# these points, cell parameters are recorded. Depnding on the optional drug_kill
# argument, cells are killed at these times with probability 0 (if sensitive)
# or 1 (if resistant).
# Have to decide an Nmax - if chosen, this will override tmax, and the
# simulation will end after a growth following t_pulse is > Nmax. The times at
# which t_pulse occurs will still be at tmax/n_pulse however.
# If store_pulse=true, stores cell_IDs, barcodes and R at each t_pulse.
# If insta_kill=true, perform the first kill-step at t=0.0. Continue with all
# future kill steps at each (tmax/n_pulse).

# Introducing the non-deterministic death of cells during the drug-kill
# step. Use the single parameter, psi, to control the magnitude of
# stochasticity.
# psi limited between 0.0 and 0.5.
# For some given psi, the probability of death during the drug-kill step:
# Resistant cells (R = 1.0): with probability (0.0 + psi).
# Sensitive cells (R = 0.0): with probability (1.0 - psi).

function grow_kill_rec_cells(cells::Array{CancerCell}, tmax::Float64,
    mu::Float64, sig::Float64, del::Float64, n_pulse::Int64, Nmax::Int64;
    psi=0.0::Float64, al=0.0::Float64,
    R_real="b"::String, drug_kill::Bool=false, insta_kill::Bool=false,
    rep_name::String="Unassigned",
    store_pulse::Bool=false, drug_presence=0::Int64)

    # Can only insta_kill if drug_kill = true
    if insta_kill == true && drug_kill == false
        error("Can only inst_kill if drug_kill=true")
    end

    # Initiate time.
    t = 0.0

    # Create vector of times for recording (and drug) pulses.
    t_pulse_change = tmax/n_pulse
    t_pulse = t_pulse_change

    # Create a vector of times, corresponding population size and mean
    # population R score to track the barcoded cell growth progress.
    # Also keep track of the number of resistant cells (R>0.0) and cells
    # with an escape mutation (E>0.0).
    Nvec = Int64[]
    tvec = Float64[]
    Rvec = Int64[]
    Evec = Int64[]

    # Also store the cell_IDs, barcodes and R scores: at the beginning, just
    # prior to each drug pulse and at the end of the growth period.
    pulse_ts = Float64[]
    pulse_cids = Array{Array{Int64, 1}, 1}(undef, 0)
    pulse_bcs = Array{Array{Float64, 1}, 1}(undef, 0)
    pulse_Rs = Array{Array{Float64, 1}, 1}(undef, 0)

    if store_pulse == true
        # Store each for t = 0
        push!(pulse_ts, t)
        push!(pulse_cids, map(x -> x.cell_ID, cells))
        push!(pulse_bcs, map(x -> x.barcode, cells))
        push!(pulse_Rs, map(x -> x.R, cells))
    end

    # Record time and population size and number of resistant and escape
    # mutations.
    push!(Nvec, length(cells))
    push!(tvec, t)
    push!(Rvec, sum(map(x -> x.R, cells)))
    push!(Evec, sum(map(x -> x.E, cells)))

    # If insta_kill=true and drug_kill=true, perform a drug-kill step before
    # continuing with the growth stage.

    if insta_kill == true
        # Go through the cells, and kill them with:
        # Resistant (R = 1.0): p(1 - R) + psi.
        # Sensitive (R = 0.0): p(1 - R) - psi.
        kill_vec = Int64[]
        for i in 1:length(cells)
            # Random unif to determine kill.
            ran_k = rand(Uniform(0, 1))
            # if resistant...
            if cells[i].R > 0.0
                if ran_k >= (1.0 - psi)
                    push!(kill_vec, i)
                end
            # if sensitive...
            elseif cells[i].R == 0.0
                if ran_k >= (0.0 + psi)
                    push!(kill_vec, i)
                end
            end
            # (deterministic version pre- using psi).
            # if ran_k >= cells[i].R
            #    push!(kill_vec, i)
            # end
        end
        # Now kill the cells according to kill_vec
        deleteat!(cells, kill_vec)
    end

    if store_pulse == true
        # Store each for t = 0
        push!(pulse_ts, t)
        push!(pulse_cids, map(x -> x.cell_ID, cells))
        push!(pulse_bcs, map(x -> x.barcode, cells))
        push!(pulse_Rs, map(x -> x.R, cells))
    end

    # Grow cells until t >=  tmax

    while round(t, digits = 2) < tmax

        # Perform drug-kill prior to birth-death stage.

        # If the time has reached the next pulse interval...
        if t >= t_pulse

            # Optional drug-kill step:
            if drug_kill == true
                # Go through the cells, and kill them with:
                # Resistant (R = 1.0): p(1 - R) + psi.
                # Sensitive (R = 0.0): p(1 - R) - psi.
                kill_vec = Int64[]
                for i in 1:length(cells)
                    # Random unif to determine kill.
                    ran_k = rand(Uniform(0, 1))
                    # if resistant...
                    if cells[i].R > 0.0
                        if ran_k >= (1.0 - psi)
                            push!(kill_vec, i)
                        end
                    # if sensitive...
                    elseif cells[i].R == 0.0
                        if ran_k >= (0.0 + psi)
                            push!(kill_vec, i)
                        end
                    end
                    # (deterministic version pre- using psi).
                    # if ran_k >= cells[i].R
                    #    push!(kill_vec, i)
                    # end
                end
                # Now kill the cells according to kill_vec
                deleteat!(cells, kill_vec)
            end

            if store_pulse == true
                # Store time, and R & barcode distributions post-drug-killing.
                push!(pulse_ts, t)
                push!(pulse_cids, map(x -> x.cell_ID, cells))
                push!(pulse_bcs, map(x -> x.barcode, cells))
                push!(pulse_Rs, map(x -> x.R, cells))
            end

            # And increase t_pulse by t_pulse_change for next record/kill.
            t_pulse += t_pulse_change
        end

        # Break if all cells have died, recording in the tracking vectors.
        if length(cells) == 0
            push!(Nvec, length(cells))
            push!(tvec, t)
            push!(Rvec, sum(map(x -> x.R, cells)))
            push!(Evec, sum(map(x -> x.E, cells)))
            break
        end

        # Record time and population size post-drugkill,
        # and number of resistant and escape-mutation cells.
        push!(Nvec, length(cells))
        push!(tvec, t)
        push!(Rvec, sum(map(x -> x.R, cells)))
        push!(Evec, sum(map(x -> x.E, cells)))

        # Grow cells
        grow_out = grow_cells(cells, t_pulse_change, Nmax, mu, sig,
        del, psi=psi, al=al, R_real=R_real, drug_presence=drug_presence)
        # Assign cells.
        cells = grow_out.cells
        # Add the grow_output's Nvec and tvec to the total
        # time and population size vectors for grow_kill_rec.
        # Also extract the output structures Rvec and Evec.
        append!(Nvec, grow_out.Nvec)
        append!(Rvec, grow_out.Rvec)
        append!(Evec, grow_out.Evec)
        # Need to account for the fact that the grow_out's t_vec has started
        # from t=0 again.
        grow_out.tvec .+= t
        append!(tvec, grow_out.tvec)

        # Break if N is now >= Nmax.
        if length(cells) >= Nmax
            break
        end

        # Break if no cells remain.
        if length(cells) == 0
            break
        end

        # Update the time. Can perform after the growth events, as grow_cells
        # breaks if cells > Nmax.
        t += t_pulse_change

    end

    # Store the final values.
    if store_pulse == true
        push!(pulse_ts, t)
        push!(pulse_cids, map(x -> x.cell_ID, cells))
        push!(pulse_bcs, map(x -> x.barcode, cells))
        push!(pulse_Rs, map(x -> x.R, cells))
    end

    fin_t = round(tvec[length(tvec)], digits = 4)

    return Grow_Kill_Rec_Out(cells, Nvec, tvec, Rvec, Evec,
    pulse_ts, pulse_cids, pulse_bcs, pulse_Rs,
    rep_name, fin_t)

end


##################################
# Expand, Split Cells Function.
##################################


# Expand N cells, split equally into 2 x N_reps replicates, each containing
# N_seed cells. Then store each of these expanded groups in a 'flask' object,
# which are assigned to either a control or drug-treatment group.
# Can then repeatedely pass this object to experimental functions to repeat
# evolution in parallel with chosen parameter values.

# p - pre-existing resistance. Cells at t=0 are resistant (R = 1) with prob(p).
# mu - probability of resistance 'mutation' (R + 1) per cell-division.
# sig - the probability of a sensitive 'mutation' (R - 1) per cell-division.
# ... hence sig must be > 0. mu always refers to the probability of (R += 1).
# plastic - (Bool) if false, proceeds as simple mutation sim, where R increases
# by Poisson(mu) per division (hence, double hits). If true, R is only ever 0
# or 1, with the direction of movement conditional on its current state, and
# a mutation only ever leading to a change of +/-1.

function expand_split_cells(N::Int64, b::Float64, d::Float64,
    t_exp::Float64, N_seed::Int64, N_reps::Int64, p::Float64,
    mu::Float64, sig::Float64, del::Float64;
    psi=0.0::Float64, al=0.0::Float64,
    R_real="b"::String, drug_presence=0::Int64, use_lim_probs=true::Bool)

    # Use limiting probabilities according to use_lim_probs, and plastic/res_cost
    if use_lim_probs == false
        init_cells = seed_cells(N, b, d, p, mu, sig, del, psi=psi, al=al,
        use_lim_probs=false)
    elseif use_lim_probs == true
        init_cells = seed_cells(N, b, d, p, mu, sig, del, psi=psi, al=al,
        use_lim_probs=true)
    end

    # Expand the cells for t_exp. Set the Nmax to >> than can simulate.
    exp_cells = grow_cells(init_cells, t_exp, 10^10, mu, sig, del,
    psi=psi, al=al, R_real=R_real, drug_presence=drug_presence).cells

    # Need to split into N_reps * 2 - N_reps x control (CO) 'flasks'.
    #                               - N_reps x drug-treatment (DT) 'flasks'.
    N_reps * 2 * N_seed < length(exp_cells) || error("There are not enough cells to split into the chosen number of replicates. Need at least (2 * N_reps * Nseed) cells post-expansion.")

    # Sample all the cells at once, and then reshape into replicates.
    rep_cells = sample(exp_cells, (N_reps * 2 * N_seed), replace = false)
    rep_cells = reshape(rep_cells, (N_seed, 2 * N_reps))

    # Group into control and drug-treatment groups:
    CO_flasks = rep_cells[:, 1:N_reps]
    DT_flasks = rep_cells[:, (N_reps + 1):(2*N_reps)]

    # Extract the original cell_IDs, barcodes and R scores from the expanded
    # pool of cells.
    orig_cids = map(x -> x.cell_ID, exp_cells)
    orig_bcs = map(x -> x.barcode, exp_cells)
    orig_Rs = map(x -> x.R, exp_cells)

    # Return all as an Expanded_Split_Cells data structure:
    return Expanded_Split_Cells(CO_flasks, DT_flasks, orig_cids, orig_bcs,
    orig_Rs)

end


################################################################################
