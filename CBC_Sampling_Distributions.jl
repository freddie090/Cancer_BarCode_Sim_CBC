
# Cancer Barcode Simulation
# Freddie Whiting - 2020

# Analytical Approach and Comparison to Expanding and Splitting Barcodes
# into N flasks. Also investigates the distribution of resistance phenotypes
# (R scores) in the subsequent flasks, given some R assignment with a Beta
# distribution.

using Distributions
using DataFrames
using RCall
using CSV
using SpecialFunctions
using Base.Threads

e = Base.MathConstants.e


# The birth-death process (geometric distribution)
# with
# Sampling without replacement (hypergeometric distribution)

# birth-death function:
# probability of seeing a lineage of size n at time t given birth/death rates.
function p_bdn(n::Int64, t::Float64, b::Float64, d::Float64)

    alph_n = d*(e^((b - d)*t) - 1)
    alph_d = (b*e^((b - d)*t) - d)
    alph = alph_n/alph_d

    beta_n = b*(e^((b - d)*t) - 1)
    beta_d = (b*e^((b - d)*t) - d)
    beta = beta_n/beta_d

    pn = (1 - alph)*(1 - beta)*(beta^(n - 1))
    if n > 0
        return pn
else
        return alph
    end
end

# birth-death function:
# probability of seeing a lineage of size n at time t given birth/death rates,
# conditioned on survival (n > 0).
function p_bdn_csurv(n::Int64, t::Float64, b::Float64, d::Float64)

    beta_n = b*(e^((b - d)*t) - 1)
    beta_d = (b*e^((b - d)*t) - d)
    beta = beta_n/beta_d

    pn = (1 - beta)*(beta^(n - 1))
    if n == 0
        error("n must be > 0.")
else
        return pn
    end
end

R"""
ggplot(data = bcs, aes(x = count, y = ..density.., fill = "simulation")) +
geom_histogram(bins = 100, alpha = 1.0) +
geom_line(data = dft, aes(x = x, y = y, colour = "probability"), size = 0.8) +
scale_colour_manual(values = "black") +
scale_fill_manual(values = "grey") +
xlab("Count") +
ylab("") +
theme_minimal() +
theme(panel.background = element_rect(colour = "grey80"),
      legend.title = element_blank())
"""

# sampling without replacement function:
# N = the total population size sampling from,
# K = the number of cells sampled from the entire population,
# n = the number of cells in the barcode lineage of interest,
# k = the number of successes (we are calculating the probability of this),
# i.e. Pr(X = k).

function p_swr(k::Int64, n::Int64, K::Int64, N::Int64; j=false)

    if j == false
        pk = pdf(Hypergeometric(n, (N - n), K), k)
    end
    if j == true
        # if using the p(j) form, instead we have the probability of j
        # individuals remaining following sampling k without replacement...
        # this is just equivalent to asking what is p(k) in the (N - K)
        # cells...
        K = (N - K)
        pk = pdf(Hypergeometric(n, (N - n), K), k)
    end
    return pk
end

# Analagous function for sampling with replacement

function p_bin(k::Int64, n::Int64, K::Int64, N::Int64)

    pk = pdf(Binomial(K, n/N), k)
    return pk

end

# Now need to marginalise over bd_n to get the probability of
# k | b, d, t, N, K
# Because bd_n gives p(n | b, d, t) but we need N, take the expected value of
# N | b, d, t.To get this, need to choose how many unique barcode lineages
# there are at t = 0, prior to expansion, = n_bc
# To speed up, can manually set the maximum N and K we expect to realistically
# see/sample following expansion.
# Change binom=true to get the corresponding dataframe for sampling with
# replacement.

function p_k(n_bc::Int64, K::Int64, b::Float64, d::Float64, t::Float64,
    nmax::Int64, Kmax::Int64; binom=true)

    # Get deterministic approximation of N
    N = n_bc*(e^((b - d)*t)) |> round |> Int64
    # Get the probability distribution of having a sub-pop of size n at time t
    ns = collect(0:nmax)
    # If using the default p(n) vector (empty), assign the p(n)s using p_bdn.
    p_ns = map(x -> p_bdn(x, t, b, d), ns)
    # Get the probability of sampling k of any given barcode, given they
    # are found n times in the total population - N - and also given we are
    # sampling K total individuals.
    # Repeat for k = 0:K
    ks = collect(0:Kmax)
    p_ks = Array{Float64}(undef, length(ks))
    for i in 1:length(ks)
        # p(k|n)
        if binom == false
            pkns = map(x -> p_swr(ks[i], x, K, N), ns)
        else
            pkns = map(x -> p_bin(ks[i], x, K, N), ns)
        end
        # p(n|b, d, t)
        pns = p_ns
        # Sum over all ((p(k|n)*p(n|b,d,t)))
        pk = sum(pkns .* pns)
        # Assign to vector
        p_ks[i] = pk
    end
    # Dataframe for plotting
    df = DataFrame(k = ks, p_k = p_ks)
    return df
end


# Plot for increasing orders of magnitude of n_bc:K...

function comp_p_ks(xlim::Float64)
    dfs = Array{DataFrame}(undef, 0)
    for i in 1:4
        # Number infected barcode cells/Number samples into flasks.
        Ninf = 100*(10^i)
        df_temp = p_k(Ninf, Ninf, log(2)+0.2, 0.2, 6.0, 10000, 1000)
        df_temp[:Ninf] = Ninf
        push!(dfs, df_temp)
    end
    df = reduce(vcat, dfs)
    @rput df
    @rput xlim
    R"""
    library(ggplot2)
    df$Ninf <- as.factor(df$Ninf)
    ggplot(data = df, aes(x = k, y = p_k, colour = Ninf)) +
    geom_line() +
    scale_x_continuous(limits = c(0, xlim))
    """
end


# Compare binomial and hypergeometric sampling... with the aim of showing
# that the binomial analogue is a sufficiently accurate replacement.

function comp_hypge_binom(n_bc::Int64, K::Int64, b::Float64, d::Float64, t::Float64,
    nmax::Int64, Kmax::Int64)

    hdf = p_k(n_bc, K, b, d, t, nmax, Kmax)
    bdf = p_k(n_bc, K, b, d, t, nmax, Kmax, binom=true)

    @rput hdf
    @rput bdf

    R"""
    library(ggplot2)
    ggplot() +
    geom_area(data = hdf, aes(x = k, y = p_k, fill = "hdf"), alpha = 0.6) +
    geom_area(data = bdf, aes(x = k, y = p_k, fill = "bdf"), alpha = 0.6)
    """
end

# Function to perform actual expansion and sampling to compare to analytical
# distributions

function grow_samp(n_bc, b, d ,t, K)

    bcs = collect(1:n_bc)
    # Get deterministic approximation of N
    N = n_bc*(e^((b - d)*t)) |> round |> Int64
    # Run gillespie on bcs
    bdmax = b + d
    while length(bcs) <= N
        ran_pos = rand(1:length(bcs))
        ran_cell = bcs[ran_pos]
        ran = rand(Uniform(0, bdmax))
        if ran < b
            push!(bcs, ran_cell)
        end
        if b <= ran < bdmax
            deleteat!(bcs, ran_pos)
        end
        length(bcs) > 0 || error("Sim drifted to extinction.")
    end
    K < length(bcs) || error("K < expanded cell count.")
    # Sample K from the expanded bcs
    sbcs = sample(bcs, K, replace = false)
    # Convert bcs to count df
    df = DataFrame(bc = sbcs)
    df = combine(groupby(df, :bc), nrow)
    rename!(df, [:bc, :count])
    # Calculate 0 counts of bcs not in df
    df0 = DataFrame(bc = setdiff(collect(1:n_bc), df[:bc]), count = 0)
    # Join dfs
    df = vcat(df0, df)
    return df
end


# Compare each method with physical simulation...

function comp_p_k_grow(n_bc, K, b, d, t, Nmax, Kmax)

    df1 = grow_samp(n_bc, b, d, t, K)
    df2 = p_k(n_bc, K, b, d, t, Nmax, Kmax)
    df3 = p_k(n_bc, K, b, d, t, Nmax, Kmax, binom=true)
    @rput df1
    @rput df2
    @rput df3
    @rput Kmax
    R"""
    library(ggplot2)
    ggplot() +
    geom_histogram(data = df1, aes(x = count, y = ..density..), bins = Kmax, alpha = 0.6) +
    geom_area(data = df2, aes(x = k, y = p_k, fill = "without-replacement"), alpha = 0.6) +
    geom_area(data = df3, aes(x = k, y = p_k, fill = "with-replacement"), alpha = 0.6) +
    scale_x_continuous(limits = c(0, Kmax))
    """
end

# As we're using the binomial estimate, for i flasks, each p.m.f. for p(k) is
# i.i.d. Therefore, the probability of seeing the same barcode in flasks 1 and
# 2, for example, k times = p(k | n)^i * p(n), summed over all ns.

# Given i flasks, function that returns the probability of seeing the same
# barcode:
# p(k1 = x1, k2 = x2, ... , ki = xi), if greater_than=false,
# provided xvec, where = [x1, x2, ... , xi]

function p_k_Xi(n_bc::Int64, K::Int64, b::Float64, d::Float64, t::Float64,
    i::Int64, xvec::Array{Int64}, nmax::Int64, Kmax::Int64)

    i > 0 || error("i must be > 0")
    length(xvec) == i || error("Length of xvec must = i, the number of flasks.")
    # Get deterministic approximation of N
    N = n_bc*(e^((b - d)*t)) |> round |> Int64
    # ... and all n's to calculate p(n) for.
    ns = collect(0:nmax)
    # p(n) doesn't change so can calculate outside of p_kis function.
    p_ns = map(x -> p_bdn(x, t, b, d), ns)

    # Pre-assign vector to hold p(k) for each i...
    pk_arr = zeros(length(ns), i)
    # For each specified ki = xi, calculate p(ki = xi | n) for each n
    for y in 1:i
        pk_arr[:,y] = map(x -> p_bin(xvec[y], x, K, N), ns)
    end
    # Take the product for the different ks for each i
    # p(k1 = x1 | n, k2 = x2 | n, ... , ki = xi | n)
    pks = mapslices(prod, pk_arr, dims = 2)
    # Then take the product with all the p_ns
    # p(k1 = x1 | n, k2 = x2 | n, ... , ki = xi | n) * p(n)
    pkns = pks .* p_ns
    # And take sum
    pk = sum(pkns)
    return pk

end

# We can also ask what is the probability that a barcode makes it into i flasks at least
# k times. For example, we might be interested in the probability that a barcode is found in
# each flask at least once.

# Make use of the compliment probabilities:
# p(x1 ≥ k, x2 ≥ k, ... , xi ≥ k) = 1 - p(x1 = 0, x2 = 0, ... , xi = 0) +
#                                       p(x1 = 1, x2 = 1, ... , xi = 1) + (...) +
#                                       p(x1 = (k - 1), x2 = (k - 1), ... , xi = (k - 1))


function p_alk_Xi(n_bc::Int64, K::Int64, b::Float64, d::Float64, t::Float64,
    i::Int64, k::Int64, nmax::Int64, Kmax::Int64)

    N = n_bc*(e^((b - d)*t)) |> round |> Int64
    # Get the probability distribution of having a sub-pop of size n at time t
    ns = collect(0:nmax)
    # Assign the p(n)s using p_bdn.
    p_ns = map(x -> p_bdn(x, t, b, d), ns)

    pk = 0.0
    # get p(k|n) for each < k...
    pkns = zeros(length(ns))
    for y in 0:(k - 1)
        pkns = pkns .+ map(x -> p_swr(y, x, K, N), ns)
    end

    # as we've calculated the complement set, p(xi > k) = 1 - pk
    pkns = 1 .- pkns
    # for i flasks, p(k|n)^i
    pkns = pkns .^ i
    # p(k) = p(k|n) * p(n)
    pks = pkns .* p_ns
    # summed for all ns
    pk += sum(pks)

    return pk
end


# Write a simulation to calculate the expansion and splitting into flasks to
# compare to the analytical distributions.

function grow_samp_split(n_bc, K, b, d, t, n_flasks; binom=false)

    bcs = collect(1:n_bc)
    # Get deterministic approximation of N
    N = n_bc*(e^((b - d)*t)) |> round |> Int64
    # Run gillespie on bcs
    bdmax = b + d
    while length(bcs) < N
        ran_pos = rand(1:length(bcs))
        ran_cell = bcs[ran_pos]
        ran = rand(Uniform(0, bdmax))
        if ran < b
            push!(bcs, ran_cell)
        end
        if b <= ran < bdmax
            deleteat!(bcs, ran_pos)
        end
        length(bcs) > 0 || error("Sim drifted to extinction.")
    end
    K < length(bcs) || error("K < expanded cell count.")

    # Sample K from the expanded bcs for each flask...
    flask_bcs = Array{DataFrame}(undef, 0)

    for i in 1:n_flasks
        if binom == false
            samp_ind = sample(1:length(bcs), K, replace = false)
            samp_cells = bcs[samp_ind]
            deleteat!(bcs, sort(samp_ind))
        else
            samp_cells = bcs[sample(1:length(bcs), K, replace = true)]
        end
        # Convert sampled bcs to count df
        df = DataFrame(bc = samp_cells)
        df = combine(groupby(df, :bc), nrow)
        rename!(df, [:bc, Symbol("Count_", i)])
        push!(flask_bcs, df)
    end

    # Merge all of the flasks...
    full_df = flask_bcs[1]
    for i in 2:length(flask_bcs)
        full_df = join(full_df, flask_bcs[i], on = Symbol("bc"), kind = :outer)
    end
    # Remove missings...
    for col in names(full_df)
        full_df[ismissing.(full_df[col]), col] = 0
    end

    # Calculate 0 counts of bcs not in df  - i.e. lost to sampling/drift.
    no_counts = setdiff(collect(1:n_bc), full_df[:bc])
    df0 = DataFrame(zeros(Int64, length(no_counts), n_flasks))
    colnames = map(x -> Symbol("Count_", x), 1:n_flasks)
    rename!(df0, colnames)
    df0 = hcat(DataFrame(bc = no_counts), df0)
    # Join dfs
    full_df = vcat(full_df, df0)

    return full_df
end


# Compare flasks 1-4, given sampling parameters and compare to analytical
# expectations of distributions.

function comp_flasks(n_bc::Int64, K::Int64, b::Float64, d::Float64, t::Float64,
    nmax::Int64)

    df1 = p_k_Xi(n_bc, K, b, d, t, 1, nmax)
    df2 = p_k_Xi(n_bc, K, b, d, t, 2, nmax)
    df3 = p_k_Xi(n_bc, K, b, d, t, 3, nmax)
    df4 = p_k_Xi(n_bc, K, b, d, t, 4, nmax)

    grow_df = grow_samp_split(n_bc, K, b, d, t, 4)

    @rput df1
    @rput df2
    @rput df3
    @rput df4

    @rput grow_df

    R"""
    library(ggplot2)
    library(tidyr)

    df1["flask"] <- 1
    df2["flask"] <- 2
    df3["flask"] <- 3
    df4["flask"] <- 4

    df <- rbind(df1, df2); df <- rbind(df, df3); df <- rbind(df, df4)
    df$flask <- as.factor(df$flask)

    gdf <- gather(grow_df, key = flask, value = count, 2:ncol(grow_df))
    gdf$flask <- sub("Count_([1-4])", replacement = "\\1", gdf$flask) %>% as.numeric %>% as.factor

    ggplot() +
        geom_histogram(data = gdf, aes(x = count, y = ..density.., fill = flask), bins = 100, colour = "black", position = "identity", alpha = 0.4) +
        #geom_density(data = gdf, aes(x = count, fill = flask), colour = "black", alpha = 0.4) +
        scale_y_sqrt() +
        geom_line(data = df, aes(x = k, y = p_k, colour = flask), size = 1.5) +
        facet_wrap(~flask) +
        scale_x_continuous(limits = c(0, 100))

    """
end


# The next step is to combine the probability distributions for barcode counts
# per flask with the resistance phenotpye distribution.

# We use the Beta distribution to assign resistance (R) scores at the
# beginning of the simulation.
# As we have the probability of any given barcode making it > x times into
# i flasks, we can also now ask -
# What is the probability of seeing a barcode with a resistance score between
# R1 and R2 > x times in i flasks...

# We can make use of the cdf of the Beta distribution.
# beta_inc() in SpecialFunctions which returns a tuple (j,k), where,
# given beta_inc(al, be, x), j = the probability of seeing Beta(al, be) <= x and
# k = the probability of seeing Beta(a, b) >= x

# We can combine this with our p.m.f for seeing a chosen barcode lineage
# >= k times in all i flasks.
# We are then interested in one of two distributions:
# i)  p(R > r , [x1, x2, ..., xi] ≥ k) - the probability that any one of the
# initial barcodes has a Resistance phenotype > r AND is seen in all i
# flasks at least k times.
# or
# ii) p(R > r | [x1, x2, ..., xi] ≥ k) - the conditional probability, given
# any one of the initial barcodes is seen at least k times in all i flasks,
# what is the probability it has a Resistance phenotype > r.
# As p(A|B) = p(AB)/p(B) and p(R,k) = p(R) * p(k), p(R|k) = p(R).

function p_alk_Xi_gtR(n_bc::Int64, K::Int64, b::Float64, d::Float64, t::Float64,
    i::Int64, k::Int64, al::Float64, be::Float64, R::Float64,
    nmax::Int64, Kmax::Int64; cond_prob=false)

    # First get the probability of seeing the barcode >= k times in i flasks
    pk = p_alk_Xi(n_bc, K, b, d, t, i, k, nmax, Kmax)
    # And now the probability of seeing a barcode lineage with an R score >=
    # R given the barcodes were assigned the resistance phenotype with a
    # Beta distribution Beta(al, be).
    pR = beta_inc(al, be, R)[2]

    if cond_prob == true
        return pR
    else
        return pk * pR
    end

end

# Need to now perform a simulation growth expansion split whilst also assigning
# R scores to the barcodes using the Beta(al, be) distribution.

function grow_samp_split_R(n_bc, K, b, d, t, n_flasks, al, be)

    bcs = collect(1:n_bc)
    Rs = rand(Beta(al, be), n_bc)
    # Save original bcs and Rs
    orig_bcs_Rs = DataFrame(bc = bcs, R = Rs)
    # Get deterministic approximation of N
    N = n_bc*(e^((b - d)*t)) |> round |> Int64
    # Run gillespie on bcs
    bdmax = b + d
    while length(bcs) < N
        ran_pos = rand(1:length(bcs))
        ran_cell = bcs[ran_pos]
        ran_R = Rs[ran_pos]
        ran = rand(Uniform(0, bdmax))
        if ran < b
            push!(bcs, ran_cell)
            push!(Rs, ran_R)
        end
        if b <= ran < bdmax
            deleteat!(bcs, ran_pos)
            deleteat!(Rs, ran_pos)
        end
        length(bcs) > 0 || error("Sim drifted to extinction.")
    end
    K < length(bcs) || error("K < expanded cell count.")

    # Sample K from the expanded bcs for each flask...
    flask_bcs = Array{DataFrame}(undef, 0)

    for i in 1:n_flasks
        samp_ind = sample(1:length(bcs), K, replace = false)
        samp_cells = bcs[samp_ind]
        deleteat!(bcs, sort(samp_ind))
        # Convert sampled bcs to count df
        df = DataFrame(bc = samp_cells)
        df = combine(groupby(df, :bc), nrow)
        rename!(df, [:bc, Symbol("Count_", i)])
        push!(flask_bcs, df)
    end

    # Merge all of the flasks...
    full_df = flask_bcs[1]
    for i in 2:length(flask_bcs)
        full_df = join(full_df, flask_bcs[i], on = [:bc], kind = :outer)
    end
    # Remove missings...
    for col in names(full_df)
        full_df[ismissing.(full_df[col]), col] = 0
    end
    # Calculate 0 counts of bcs not in df  - i.e. lost to sampling/drift.
    no_counts = setdiff(collect(1:n_bc), full_df[:bc])
    df0 = DataFrame(zeros(Int64, length(no_counts), n_flasks))
    colnames = map(x -> Symbol("Count_", x), 1:n_flasks)
    rename!(df0, colnames)
    df0 = hcat(DataFrame(bc = no_counts), df0)
    # Join dfs
    full_df = vcat(full_df, df0)
    # Add R column from original dataframe
    full_df = join(full_df, orig_bcs_Rs, on = [:bc], kind = :outer)

    return full_df
end


# Want to simulate the barcoding, expansion, splitting and killing of cells
# many times - then ask, of these simulations, how many barcode lineages are
# present in 1, 2, ..., i flasks? Then plot the distribution.
# Can then look at attempting to find an analytical solution to this problem.
# Because we're not mutating R, we can sample from the inverse cumulative
# distribution function to speed up the iterations.

function grow_samp_split_insta_kill(n_bc, K, b, d, t, n_flasks, al, be)

    # Given a probability between 0 and 1 (Cn), what is the corresponding size
    # of the barcode at time t?

    function bc_n(p::Float64, t::Float64, b::Float64, d::Float64)

        alph_n = d*(e^((b - d)*t) - 1)
        alph_d = (b*e^((b - d)*t) - d)
        alph = alph_n/alph_d

        beta_n = b*(e^((b - d)*t) - 1)
        beta_d = (b*e^((b - d)*t) - d)
        beta = beta_n/beta_d

        n = (log((1 - p)/(1 - alph)))/(log(beta))

        if n <= 0
            n = 0.0
        end
        return n
    end


    # Now we want to expand our barcodes and corresponding R scores
    # given b, d, t...
    bcs = collect(1:n_bc)
    Rs = rand(Beta(al, be), n_bc)
    # Save original bcs and Rs
    orig_bcs_Rs = DataFrame(bc = bcs, R = Rs)
    # Calculate a distribution of barcodes following growth for t
    freqs = map(x -> Int64(ceil(bc_n(x, t, b, d))), rand(length(bcs)))
    # Now repeat the barcodes and Rs by their updated frequency
    exp_bcs = vcat(fill.(bcs, freqs)...)
    exp_Rs = vcat(fill.(Rs, freqs)...)
    # Sample into the n_flasks...
    flask_samp_bcs = Array{DataFrame}(undef, 0)
    flask_kill_bcs = Array{DataFrame}(undef, 0)
    for i in 1:n_flasks
        length(exp_bcs) == length(exp_Rs) || error("Length of bcs and Rs should be the same.")
        # Sample K cells' barcodes and R scores without replacement.
        samp_ind = sample(1:length(exp_bcs), K, replace = false)
        samp_cells = exp_bcs[samp_ind]
        samp_Rs = exp_Rs[samp_ind]
        # Without replacement, so delete originals.
        deleteat!(exp_bcs, sort(samp_ind))
        deleteat!(exp_Rs, sort(samp_ind))
        # Convert sampled bcs to count df
        df = DataFrame(bc = samp_cells)
        df = combine(groupby(df, :bc), nrow)
        rename!(df, [:bc, Symbol("Count_", i)])
        push!(flask_samp_bcs, df)
        # Now perform the kill step based on R score.
        kill_ind = Array{Int64}(undef, 0)
        for r in 1:length(samp_Rs)
            randn = rand()
            if randn > samp_Rs[r]
                push!(kill_ind, r)
            end
        end
        # Delete the flask's cells who were killed in the kill-step.
        deleteat!(samp_cells, sort(kill_ind))
        deleteat!(samp_Rs, sort(kill_ind))
        df = DataFrame(bc = samp_cells)
        df = combine(groupby(df, :bc), nrow)
        rename!(df, [:bc, Symbol("Count_", i)])
        push!(flask_kill_bcs, df)
    end

    function convert_df(dfs::Array{DataFrame})
        # Merge all of the flasks...
        full_df = dfs[1]
        for i in 2:length(dfs)
            full_df = join(full_df, dfs[i], on = [:bc], kind = :outer)
        end
        # Remove missings...
        for col in names(full_df)
            full_df[ismissing.(full_df[col]), col] = 0
        end
        # Calculate 0 counts of bcs not in df  - i.e. lost to sampling/drift.
        no_counts = setdiff(collect(1:n_bc), full_df[:bc])
        df0 = DataFrame(zeros(Int64, length(no_counts), n_flasks))
        colnames = map(x -> Symbol("Count_", x), 1:n_flasks)
        rename!(df0, colnames)
        df0 = hcat(DataFrame(bc = no_counts), df0)
        # Join dfs
        full_df = vcat(full_df, df0)
        # Add R column from original dataframe
        full_df = join(full_df, orig_bcs_Rs, on = [:bc], kind = :outer)
        return full_df
    end

    full_samp_df = convert_df(flask_samp_bcs)
    full_kill_df = convert_df(flask_kill_bcs)

    return full_samp_df, full_kill_df

end


# Have a version that runs Nsim times, and returns a data frame of:
# the total number of barcodes found in 1, 2, ..., i flasks n times after
# and their corresponding R scores,
# running Nsim total times,
#  i) pre-kill: before doing the instant drug-kill step based on R scores.
# ii) post-kill: after doing the instant drug-kill step based on R scores.

function grow_samp_split_insta_kill_Nsim(n_bc, K, b, d, t, n_flasks, al, be, Nsim)

    # And pre- and post-data frames with counts and R scores.
    samp_dfs = Array{DataFrame}(undef, 0)
    kill_dfs = Array{DataFrame}(undef, 0)

    # Repeat following Nsim times...
    for n in 1:Nsim

        output = grow_samp_split_insta_kill(n_bc, K, b, d, t, n_flasks, al, be)
        full_samp_df = output[1]
        full_kill_df = output[2]

        samp_df = full_samp_df
        samp_df[:Nsim] = n
        samp_df[:found_in] = sum.(eachrow(samp_df[:,2:(n_flasks + 1)] .>= 1))
        samp_df[:kill] = "pre"
        push!(samp_dfs, samp_df)

        kill_df = full_kill_df
        kill_df[:Nsim] = n
        kill_df[:found_in] = sum.(eachrow(kill_df[:,2:(n_flasks + 1)] .>= 1))
        kill_df[:kill] = "post"
        push!(kill_dfs, kill_df)

    end

    # Merge barcode counts and R scores
    samp_df = reduce(vcat, samp_dfs)
    kill_df = reduce(vcat, kill_dfs)

    df = vcat(samp_df, kill_df)

    return df

end



# Plot the R densities from a grow_samp_split_insta_kill output

function plot_grow_samp_split_insta_kill(df)

    @rput df

    R"""
    # Load libraries
    library(ggplot2)
    library(tidyr)

    df_ext <- df[df$found_in > 0 ,]
    # Change to factors for plotting
    df_ext$found_in <- as.factor(df_ext$found_in)
    df_ext$kill <- as.factor(df_ext$kill)
    df_ext$kill <- factor(df_ext$kill, levels = c("pre", "post"))

    ggplot(data = df_ext, aes(x = R, y = ..count.., fill = found_in)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~kill, ncol = 1) +
    scale_y_sqrt()

    """

end


################################################################################

# Dec 2020 - Binary R Simulations

# Pre-Existing Resistance.

# Want to find the analytical expectation of seeing a resistant barcode in i
# flasks > 0 times given barcode lineages have a probability p of being r
# resistant at t=0.

# Beacuse each barcode lineage has probability p of being resistant in the
# pre-existing case, the probability of a resistant barcode lineage making it
# into i flasks is simply the product of the probabilitites. First we need to
# find the probability a barcode makes it into j of i flasks.

# Q: What is the probability that a barcode is found only in j out of i flasks
# at least once (p(k >= 1)) given it has a probability p(n) of being at size n
# in the expanded population of size N, and given we are sampling K individuals
# from N j times?
# Well we can think of this Q as a second binomial sampling step for a total
# of i (the number of flasks) draws. Seeing a barcode in j of the i flasks
# is equivalent to asking the probability of seeing i successes, where the
# probability of success is p(k >= 1 | n). We then need to marginalsie over
# all ns.
# The probability of seeing a barcode in j of i flasks is therefore:
# Sum(n=0 to N){(i choose j)(p(k >= 1 | n^j) * (1 - p(k >= 1 | n)^(i - j))*p(n))}

function p_alk_in_j_of_i(n_bc::Int64, K::Int64, b::Float64, d::Float64,
    t::Float64, i::Int64, j::Int64, nmax::Int64, Kmax::Int64)

    # Get deterministic approximation of N
    N = n_bc*(e^((b - d)*t)) |> round |> Int64
    # ... and all n's to calculate p(n) for.
    ns = collect(0:nmax)
    # calculate the p(n)s.
    p_ns = map(x -> p_bdn(x, t, b, d), ns)
    # For p(k >= 1), we only need to know p(k = 0). p(k >= 1) is then simply
    # 1 - p(k = 0). Use this as prob(success) in binomial.
    pk_0s = map(x -> p_bin(0, x, K, N), ns)
    pk_1ms = 1 .- pk_0s
    # Plug into binomial for j of i flasks.
    pk_bins = map(x -> pdf(Binomial(i, x), j), pk_1ms)
    # Now just need to marginalise over p(ns)....
    pjs = pk_bins .* p_ns
    pj = sum(pjs)
    return pj

end

# In the pre-existing case, the probability of seeing a resistant barcode
# lineage in j of i flasks is simply p * p(j).
# Return the expected number of barcodes to directly compare to the simulation.

function p_alk_in_j_of_i_pre_ex(n_bc::Int64, K::Int64, b::Float64, d::Float64,
    t::Float64, i::Int64, j::Int64, nmax::Int64, Kmax::Int64, p::Float64)

    pr_res = p_alk_in_j_of_i(n_bc, K, b, d, t, i, j, nmax, Kmax) * p * n_bc

    return pr_res

end

# For the de-novo mutation case,...

function prob_res_de_novo_mut()
    # Expected number of lineages with >= 1 resistant cell using
    #the Pr(Resistance) from Iwasa et al 2006.
    exp_pRs = 1 .- exp.((-Nmus .* (b/d)) .* log(b/(b - d)))
    # A second estimate of the number of lineages with >= 1 resistant mutation.
    # This is assuming the deterministic growth function, and integrating
    # betwen 0 and t.
    exp_pRs2 = (b.*mus) .* (((1/(b - d))*e^((b - d)*t)) - (1/(b - d)))
    # Can see from plotting that the Iwasa estimate is more accurate.

    # Can also calculate the expected number of extant resistant cells.
    #exp_nRs = log.(N)/((b/d - 1)*(log(b/b - d)))
    # This gives negative value?...

    # Instead, can get the expected number of resistant cells in a lineage -
    # this should just be the product of the number of births each cell has
    # experienced, the mutation rate per division, and the expected number of
    # cells at time t.
    exp_nRs = mus .* (N*2*b*t)
    # Can see from plotting that this matches with the TOTAL mean number of
    # resistant cells, across multiple simulation runs (which can also be
    # thought of as simply cells growing in parallel).

    # Can also calculate the expected number of extant cells - nR/pR
    exp_enRs = exp_nRs ./ exp_pRs

    # If we are conditioning on survival, we should find that we get the
    # same mean number of resistant lineages, independent of mu. We should
    # see more of them in higher mu values, but the independent lineages
    # should be the same size.
end


# Now have a function to simulate the cell growth/splitting steps.

# Write a simulation to calculate the expansion and splitting into flasks to
# compare to the analytical distributions. Use the functions from the Binary R
# code.

include("CBC_Binary_R_Simulations.jl")

function bin_R_grow_samp_split(n_bc, K, b, d, t, n_flasks, p, mu)

    if p > 0 && mu > 0
        error("Can not have p and mu > 0.")
    end

    out = expand_split_cells_cost_bin_R(n_bc, b, d, t, K, n_flasks, p, mu,
    0.0, 0.0, false, false, use_lim_probs=false)

    out_df = all_bc_counts_tot_R_mut(out)
    # Just pull out counts and number of R mutations
    bc_R_df = out_df[:, map(x->match(r"DT", x), names(out_df)) .!= nothing]

    return bc_R_df

end

# Turn a bc_R_df into
# i) the total number of barcodes found in i flasks, and
# ii) the number of resisatnt barcode lineages found in i flasks.

function bc_R_found_ins(bc_R_df::DataFrame)

    @rput bc_R_df
    R"""
    options(warn=-1)

    library(tidyr)
    library(plyr)
    library(dplyr)
    # Pull out resistant barcodes
    Rbc_R_df <- bc_R_df[rowSums(bc_R_df[, grep("DT\\d{1}_R", colnames(bc_R_df))]) > 0 ,]
    # If none, manually add 0 counts.
    if(nrow(Rbc_R_df) == 0){
        Rbc_R_df[1 ,] <- 0
    }
    # Now just pull out counts
    bc_R_df <- bc_R_df[, grep("DT\\d{1}\\>", colnames(bc_R_df))]
    Rbc_R_df <- Rbc_R_df[, grep("DT\\d{1}\\>", colnames(Rbc_R_df))]
    # Add found_in column to each
    bc_R_df["found_in"] <- rowSums(bc_R_df > 0)
    Rbc_R_df["found_in"] <- rowSums(Rbc_R_df > 0)
    bc_fi <- data.frame(table(bc_R_df$found_in))
    Rbc_fi <- data.frame(table(Rbc_R_df$found_in))
    colnames(bc_fi) <- c("found_in", "count")
    colnames(Rbc_fi) <- c("found_in", "Rcount")
    df_fi <- join(bc_fi, Rbc_fi)
    df_fi[is.na(df_fi)] <- 0
    """
    @rget df_fi
    return df_fi

end


# Run the simulation for the pre-existing case and compare to the theoretical
# expectations for different values of p.

function comp_th_sim_pre_ex_res(n_bc::Int64, K::Int64, b::Float64, d::Float64,
    t::Float64, ps::Array{Float64}, i::Int64, nmax::Int64, Kmax::Int64, nsim::Int64)

    sim_dfs = Array{DataFrame}(undef, 0)
    for x in 1:length(ps)
        for ns in 1:nsim
            temp_df = bc_R_found_ins(bin_R_grow_samp_split(n_bc, K, b, d, t, i, ps[x], 0.0))
            temp_df[:nsim] = ns
            temp_df[:p] = ps[x]
            push!(sim_dfs, temp_df)
        end
    end
    sim_df = vcat(sim_dfs...)
    # Theory expectations.
    th_dfs = Array{DataFrame}(undef, 0)
    js = collect(1:i)
    for p in ps
        ex_fis = []
        for j in js
            push!(ex_fis, p_alk_in_j_of_i_pre_ex(n_bc, K, b, d, t, i, j, nmax, Kmax, p))
        end
        push!(th_dfs, DataFrame(found_in=js, th_exp=ex_fis, p=p))
    end
    th_df = vcat(th_dfs...)

    return sim_df, th_df

end

function plot_th_sim_pre_ex_res(output)

    sim_df = output[1] ; th_df = output[2]
    @rput sim_df ; @rput th_df

    R"""
    library(ggplot2)
    sim_df <- sim_df[sim_df$found_in != 0 ,]
    sim_df$p <- as.factor(sim_df$p)
    th_df <- data.frame(sapply(th_df, as.numeric))
    th_df$p <- as.factor(th_df$p)
    named_rhos <- paste0("rho = ", unique(as.character(th_df$p)))
    names(named_rhos) <- unique(as.character(th_df$p))
    ggplot(data = sim_df, aes(x = found_in, y = Rcount, fill = p)) +
    geom_violin(alpha = 0.6) +
    geom_point(data = th_df, aes(x = found_in, y = th_exp, colour = "theory"), size = 3) +
    scale_colour_manual(values = "black") +
    xlab("Found in i flasks") +
    ylab("Number of Barcode Lineages") +
    theme_minimal() +
    facet_wrap(~p, ncol = 1, labeller = labeller(p = as_labeller(named_rhos)))
    """

end

# For the de-novo mutation case, first confirm that my estimates of:
# pR (the probability of resistance)
# nR (the total mean number of resistant cells within a lineage)
# eR (the mean number of resistance cells conditioned on resistance arising)
# match closely with simulations.

function comp_th_sim_de_nov_mut(b, d, t, nmax, mus::Array{Float64}, nsim, nrep)
    cell = seed_cells_cost_bin_R(1, b, d, 0.0, 0.0, 0.0, 0.0)
    mu_dfs = Array{DataFrame}(undef, 0)
    for mu in mus
        for j in 1:nrep
            ns = []
            muts = []
            for i in 1:nsim
                # grow cells.
                out = grow_cells_cost_kmc_bin_R(cell, t, nmax, mu, 0.0, 0.0)
                # collect the counts and number of reistant mutations.
                if length(out) == 0
                    push!(ns, 0)
                    push!(muts, 0)
                else
                    push!(ns, length(out))
                    push!(muts, sum(map(x -> x.R > 0, out)))
                end
            end
            temp_df = DataFrame(n = ns, n_mut = muts, nsim = 1:nsim, nrep = j, mu = mu)
            push!(mu_dfs, temp_df)
        end
    end
    mu_df = vcat(mu_dfs...)

    # Now get the theoretical expectations.
    N = 1*(e^((b - d)*t)) |> round |> Int64
    Nmus = N.*mus
    exp_pRs = 1 .- exp.((-Nmus .* (b/d)) .* log(b/(b - d)))
    exp_nRs = mus .* (N*2*b*t)
    exp_enRs = exp_nRs ./ exp_pRs

    th_df = DataFrame(mu = mus, exp_pR = exp_pRs, exp_nR = exp_nRs, exp_enR = exp_enRs)

    return mu_df, th_df

end

function plot_th_sim_de_nov_mut(output)

    mu_df = output[1]
    th_df = output[2]

    @rput mu_df
    @rput th_df

    R"""
    library(ggplot2)
    library(plyr)
    library(scales)

    mylog10_trans <- function (base = 10) {
      trans <- function(x) log(x + 1e-07, base)
      inv <- function(x) base^x
      trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base),
                domain = c(1e-100, Inf))
    }
    mybreaks <- c(0, (10^(-10:0) - 1e-07))

    mu_df <- data.frame(sapply(mu_df, as.numeric))
    mu_df$mu <- as.factor(mu_df$mu)
    th_df$mu <- as.factor(th_df$mu)
    sim_p_Rs <- plyr::ddply(mu_df, .(mu, nrep), function(x){nrow(x[x$n_mut > 0 ,])/nrow(x)})
    sim_n_Rs <- plyr::ddply(mu_df, .(mu, nrep), function(x){mean(x$n_mut)})
    sim_en_Rs <- plyr::ddply(mu_df, .(mu, nrep), function(x){mean(x[x$n_mut > 0 ,]$n_mut)})
    sim_en_Rs[is.na(sim_en_Rs)] <- 0
    sim_p_Rs$mu <- as.factor(sim_p_Rs$mu)

    p1 <- ggplot(data = sim_p_Rs, aes(x = mu, y = V1)) +
    geom_violin() +
    geom_point(data = th_df, aes(x = mu, y = exp_pR)) +
    scale_y_continuous(trans = "mylog10", breaks = mybreaks) +
    ylab("pR: Proportion Resistant")

    p2 <- ggplot(data = sim_n_Rs, aes(x = mu, y = V1)) +
    geom_violin() +
    geom_point(data = th_df, aes(x = mu, y = exp_nR)) +
    scale_y_continuous(trans = "mylog10", breaks = mybreaks) +
    ylab("nR: Mean number of Resistant Cells per Lineage")

    p3 <- ggplot(data = sim_en_Rs, aes(x = mu, y = V1)) +
    geom_violin() +
    geom_point(data = th_df, aes(x = mu, y = exp_enR)) +
    scale_y_continuous(trans = "pseudo_log", breaks = c(0, seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100))) +
    ylab("eR: Mean number of Extant Resistant Cells per Lineage")

    cowplot::plot_grid(p1, p2, p3, nrow = 1)

    """
end


################################################################################

# Probability of technical replicate sampling, assuming no amplification noise.
# i.e. assume underlying birth-death p.m.f (p(n)).

# j2 - dictates p(j1 = X | j2 = j2), for plotting the joint pmf.

function tech_rep_samp_dist_plot(n_bc, b, d, t, nmax, K, kmax, J, jmax, j2)

    # Get deterministic approximation of N
    N = n_bc*(e^((b - d)*t)) |> round |> Int64
    # ... and all n's to calculate p(n) for.
    ns = collect(0:nmax)
    # calculate the p(n)s.
    p_ns = map(x -> p_bdn(x, t, b, d), ns)

    # First, the extraction sampling step (excluding amplification noise for now).
    # calculate the p(ks).
    ks = collect(0:kmax)
    # p(k | n)... well prob = n/N, for K rounds of sampling.
    pkns = []
    for k in ks
        push!(pkns, map(x -> pdf(Binomial(K, x/N), k), ns))
    end
    # now to get p(k) = p(k | n) * p(n) for summed over all ns 0-> N.
    p_ks = sum.(map(x -> x.*p_ns, pkns))

    # Now for the sequencing sampling step.
    # calculate the p(js).
    js = collect(0:jmax)
    # p(j | k)... well prob = k/K, for J rounds of sampling.
    pjns = []
    for j in js
        push!(pjns, map(x -> pdf(Binomial(J, x/K), j), ks))
    end
    # now to get p(j) = p(j | k) * p(k) for summed over all ks 0-> K.
    p_js = sum.(map(x -> x.*p_ks, pjns))

    # Can now put all into dfs for plotting comparison.
    pn_df = DataFrame(n = ns, p_n = p_ns)
    pk_df = DataFrame(k = ks, p_k = p_ks)
    pj_df = DataFrame(j = js, p_j = p_js)
    @rput pn_df; @rput pk_df; @rput pj_df

    pjn_as = []
    pjn_bs = []
    for j in js
        push!(pjn_as, map(x -> pdf(Binomial(J, x/K), j), ks))
    end
    for j in j2
        push!(pjn_bs, map(x -> pdf(Binomial(J, x/K), j), ks))
    end
    # p(j1 = X | j2 = 10) =
    # p(j = X|k) * p(j = 10|k) * p(k), summed over all ks 0-> K.
    p_j2s = []
    for i in 1:length(pjn_as)
        push!(p_j2s, sum(pjn_as[i] .* pjn_bs[1] .* p_ks))
    end

    pj2_df = DataFrame(j = js, p_j = p_j2s)
    @rput pj2_df; @rput j2


    R"""
    library(ggplot2)
    library(plyr)
    library(scales)
    pn_df <- data.frame(sapply(pn_df, as.numeric))
    pk_df <- data.frame(sapply(pk_df, as.numeric))
    pj_df <- data.frame(sapply(pj_df, as.numeric))
    pj2_df <- data.frame(sapply(pj2_df, as.numeric))

    mylog10_trans <- function (base = 10) {
      trans <- function(x) log(x + 1e-07, base)
      inv <- function(x) base^x
      trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base),
                domain = c(1e-100, Inf))
    }
    mybreaks <- c(0, (10^(-10:0) - 1e-07))

    p1 <- ggplot() +
        geom_point(data = pn_df, aes(x = n, y = p_n, fill = "n"), colour = "black", shape = 21, size = 3, alpha = 0.4) +
        geom_point(data = pk_df, aes(x = k, y = p_k, fill = "k"), colour = "black", shape = 21, size = 3, alpha = 0.4) +
        geom_point(data = pj_df, aes(x = j, y = p_j, fill = "j"), colour = "black", shape = 21, size = 3, alpha = 0.4) +
        geom_point(data = pj2_df, aes(x = j, y = p_j, fill = "j2"), colour = "black", shape = 21, size = 3, alpha = 0.4) +
        scale_x_log10(limits = c(1, 1000)) +
        scale_y_continuous(limits = c(0, 0.25), trans = "mylog10", breaks = mybreaks)

    # Can compare to the simulated seq tech rep.
    bc_csv <- read.csv("C:/Users/whitin01/Google Drive/Barcode_Output_Copy/GC-FJHW-8163/Har_Sim_N-10^6_b_log2+0.2_d_0.2_texp_6.0_K-1x10^6_J-5x10^6.csv", stringsAsFactors=F)
    seq_tr1 <- data.frame(table(subset(bc_csv, sim_a_Har1 == j2)$sim_b_Har1))
    seq_tr1$Var1 <- as.numeric(seq_tr1$Var1)
    seq_tr2 <- data.frame(table(subset(bc_csv, sim_a_Har2 == j2)$sim_b_Har2))
    seq_tr2$Var1 <- as.numeric(seq_tr2$Var1)
    p2 <- ggplot() +
    geom_line(data = pj2_df, aes(x = j, y = p_j*1000000, colour = "j2"), alpha = 0.6, size = 2) +
    geom_point(data = seq_tr1, aes(x = Var1, y = Freq, colour = "tr1"), size = 3) +
    geom_point(data = seq_tr2, aes(x = Var1, y = Freq, colour = "tr2"), size = 3) +
    scale_x_log10(limits = c(1, 1000)) +
    scale_y_continuous(trans = "pseudo_log")

    cowplot::plot_grid(p1, p2, ncol = 1)

    """

end



###############################################################################

# Jan 2021 - CBC Noise Model.
# Exploring the Beta-Binomial distribution for modelling the 'unknown noise'
# during the sampling stems.

# Can reparametise the beta distribution so it is in terms of its mean (mu) and
# 'sample size' (v), which controls the 'spread' on the beta distribution...
# α = (mu * v)
# β = (1 - mu) * v

function muv_BetaBinomial(n::Int64, mu::Float64, v)

    α = (mu * v)
    β = (1 - mu) * v

    return BetaBinomial(n, α, β)

end

# So we can plot a beta binomial with a constant mean, but vary the spread to
# see how this influences the final distribution.

function plot_comb_muv_betabin(n::Int64, p::Float64, vs, nmax::Int64)

    xs = collect(0:nmax)
    bb_dfs = Array{DataFrame}(undef, 0)
    for i in 1:length(vs)
        tdf = DataFrame(x = xs, y = map(x -> pdf(muv_BetaBinomial(n, p, vs[i]), x), xs), v = vs[i])
        push!(bb_dfs, tdf)
    end
    bb_df = vcat(bb_dfs...)
    @rput bb_df
    R"""
    library(ggplot2)
    bb_df[is.na(bb_df)] <- 0
    bb_df$v <- as.factor(as.character(bb_df$v))
    ggplot(data = bb_df, aes(x = x, y = y, colour = v)) +
    geom_line()
    """
end
