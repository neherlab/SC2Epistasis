# Load structure
load data/PDB/Spike/7krr.pdb

# Show as cartoon
hide everything, all
show cartoon, all

# Custom grey color
set_color pale_grey, [0.85, 0.85, 0.85]

# Color all in pale grey
color pale_grey, all

# -------------------
# S1 subunit domains
# -------------------
select NTD, resi 14-305
color cyan, NTD

select RBD, resi 319-541
color blue, RBD

select CTD1, resi 542-590
color orange, CTD1

select CTD2, resi 591-690
color yellow, CTD2

# -------------------
# Specific functional motifs in S2
# -------------------
select FP, resi 788-806
color red, FP

#select FPPR, resi 816-834
#color magenta, FPPR

select HR1, resi 910-984
color green, HR1

select CH, resi 985-1034
color lightgreen, CH

select CD, resi 1035-1067
color hotpink, CD

select HR2, resi 1163-1212
color violet, HR2

# Orient structure according to saved view
set_view (\
     0.760770679,    0.023988612,    0.648576975,\
    -0.648804545,    0.002344940,    0.760951936,\
     0.016733591,   -0.999708891,    0.017347621,\
     0.000000000,   -0.000000000, -654.684326172,\
   197.759201050,  198.350357056,  192.401489258,\
   516.158081055,  793.210571289,  -20.000000000 )

# Center and zoom for good view
zoom all

# Render high-resolution image
ray 1600,1200
png results/figures/spike_trimer_domains.png, dpi=600
