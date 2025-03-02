mol load pdb ../../step6.0_minimization.pdb
mol addfile ../../run1_fit-rot+trans_300-500ns.xtc first 1 last -1 step 1 waitfor all;

# set two groups to calculate salt bridges
set group1_start 33
set group1_end 44
set group2_start 90
set group2_end 422

# generate list of residues in two groups
set group1_res [list]
for {set i $group1_start} {$i <= $group1_end} {incr i} {lappend group1_res $i}
set group2_res [list]
for {set i $group2_start} {$i <= $group2_end} {incr i} {lappend group2_res $i}

# salt bridge cutoff
set distance_cutoff 4.0 ;# unit Å

# get charged residues
proc get_charged_residues {reslist} {
    set charged [list]
    foreach resid $reslist {
        set sel [atomselect top "resid $resid"]
        set resname [lindex [$sel get resname] 0]
        set charge ""
        set atoms ""
        switch $resname {
            ARG {set charge "+"; set atoms "NH1 NH2 NE"}
            LYS {set charge "+"; set atoms "NZ"} 
            ASP {set charge "-"; set atoms "OD1 OD2"}
            GLU {set charge "-"; set atoms "OE1 OE2"}
        }
        if {$charge ne ""} {lappend charged [list $resid $charge $atoms]}
        $sel delete
    }
    return $charged
}

# get resid of charged residues
set group1_charged [get_charged_residues $group1_res]
set group2_charged [get_charged_residues $group2_res]

# generate effective residue pairs
array set valid_pairs {}
foreach res1 $group1_charged {
    foreach res2 $group2_charged {
        if {[lindex $res1 1] ne [lindex $res2 1]} {
            set key "[lindex $res1 0]-[lindex $res2 0]"
            set valid_pairs($key) 1
        }
    }
}

#
array set saltbridge_counts {}
foreach pair [array names valid_pairs] {set saltbridge_counts($pair) 0}

# select atoms
set sel_group1_pos [atomselect top "(resid $group1_res) and (resname ARG LYS) and name NH1 NH2 NE NZ"]
set sel_group2_neg [atomselect top "(resid $group2_res) and (resname ASP GLU) and name OD1 OD2 OE1 OE2"]
set sel_group1_neg [atomselect top "(resid $group1_res) and (resname ASP GLU) and name OD1 OD2 OE1 OE2"]
set sel_group2_pos [atomselect top "(resid $group2_res) and (resname ARG LYS) and name NH1 NH2 NE NZ"]

# atoms to resid 
array set atom_to_resid {}
set all_atoms [atomselect top all]
foreach idx [$all_atoms get index] resid [$all_atoms get resid] {
    set atom_to_resid($idx) $resid
}
$all_atoms delete

# total frames
set total_frames [molinfo top get numframes]
set all [atomselect top all];

# calculation
for {set frame 0} {$frame < $total_frames} {incr frame} {
    animate goto $frame
    $all update;
    puts "Processing frame $frame/$total_frames"
    
    # positive-negative contacts
    set contacts [measure contacts $distance_cutoff $sel_group1_pos $sel_group2_neg]
    foreach a1 [lindex $contacts 0] a2 [lindex $contacts 1] {
        set pair "$atom_to_resid($a1)-$atom_to_resid($a2)"
        if {[info exists valid_pairs($pair)]} {incr saltbridge_counts($pair)}
    }
    
    # negative-positive contacts
    set contacts [measure contacts $distance_cutoff $sel_group1_neg $sel_group2_pos]
    foreach a1 [lindex $contacts 0] a2 [lindex $contacts 1] {
        set pair "$atom_to_resid($a1)-$atom_to_resid($a2)"
        if {[info exists valid_pairs($pair)]} {incr saltbridge_counts($pair)}
    }
}

# output
set out [open sb_probabilities.dat w]

puts $out "\nSalt Bridge Probabilities:"
foreach pair [lsort [array names saltbridge_counts]] {
    set prob [expr {double($saltbridge_counts($pair))/$total_frames*100}]
    puts $out "[format "Residue %-5s - %-5s: %6.2f%%" [lindex [split $pair -] 0] [lindex [split $pair -] 1] $prob]"
}

close $out

# output2
set out2 [open sb_probabilities_simple.dat w]

foreach pair [lsort [array names saltbridge_counts]] {
    set prob [expr {double($saltbridge_counts($pair))/$total_frames*100}]
    puts $out "[format "%-5s %-5s %6.2f%%" [lindex [split $pair -] 0] [lindex [split $pair -] 1] $prob]"
}

close $out2


$sel_group1_pos delete
$sel_group2_neg delete
$sel_group1_neg delete
$sel_group2_pos delete


exit

