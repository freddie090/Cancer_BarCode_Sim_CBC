
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Mutable Structures

################################################################################

# Barcode library - vector of unique barcodes to sample from when 'infecting'
# cells with barcodes.

mutable struct BarcodeLibrary
    barcodes::Array{Float64}
    #BarcodeLibrary() = new()
end

# Cancer Cell - cell_ID, barcode, birth and death rates and resistance score.

mutable struct CancerCell
    cell_ID::Int64 # ID if want to distinguish unique lineage from barcode ID.
    barcode::Float64 # Barcode identity.
    b::Float64 # birth rate.
    d::Float64 # death rate.
    R::Float64 # Resistant phenotype (binary).
    #E::Float64 # Escape mutation (binary) - if 1.0, R no longer incurs cost.
    Ndiv::Int64 # Number of divisions a given cell lineage has experienced.
    #CancerCell() = new()
end

# Output of grow_cells.

mutable struct Grow_Out
    cells::Array{CancerCell}
    Nvec::Array{Int64}
    tvec::Array{Float64}
    Rvec::Array{Int64}
    #Evec::Array{Int64}
    fin_t::Float64
end

# Output of grow_kill_rec_cells.

mutable struct Grow_Kill_Rec_Out
    cells::Array{CancerCell}
    Nvec::Array{Int64}
    tvec::Array{Float64}
    Rvec::Array{Int64}
    #Evec::Array{Int64}
    pulse_ts::Array{Float64}
    pulse_cids::Array{Array{Int64}}
    pulse_bcs::Array{Array{Float64}}
    pulse_Rs::Array{Array{Float64}}
    rep_name::String
    fin_t::Float64
end

# Expanded and split cells in replicate 'flasks' - either control or
# drug-treatment.
# Also holds the original distributions of cell_IDs, barcodes and R scores.

mutable struct Expanded_Split_Cells
    CO_flasks::Array{CancerCell,2}
    DT_flasks::Array{CancerCell,2}
    orig_cids::Array{Int64}
    orig_bcs::Array{Float64}
    orig_Rs::Array{Float64}
end

# Experimental output - an array of Grow_Kill_Rec_Out structures split
# by control (CO) and drug-treatment (DT) replicates. Also takes the original
# cell_IDs, barcodes and R scores from the expanded cells.

mutable struct Experiment_Output
    CO_outputs::Array{Grow_Kill_Rec_Out}
    DT_outputs::Array{Grow_Kill_Rec_Out}
    #orig_cids::Array{Int64}
    #orig_bcs::Array{Float64}
    #orig_Rs::Array{Float64}
end

# Experimental input -  An Array{CancerCell} divided into replicates which
# is either collected directly from expanded and split cells, or from a
# previous passage. Means previously grown cells can be 're-seeded' into
# a subsequent replicate flask.

mutable struct Experiment_Input
    CO_inputs::Array{CancerCell, 2}
    DT_inputs::Array{CancerCell, 2}
    Experiment_Input()=new()
end


################################################################################
