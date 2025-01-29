# get average
proc mean  L {
    expr ([join $L +])/[llength $L].
}
# get standard deviation
proc mean2 list {
    set sum 0
    foreach i $list {set sum [expr {$sum+$i*$i}]}
    expr {double($sum)/[llength $list]}
}
proc stddev list {
    set m [mean $list]
    expr {sqrt([mean2 $list]-$m*$m)}

#trajectory type, CHANGE
set type "netcdf"; #dcd, netcdf

#define path, CHANGE
set path  "$/home/data/joch8192/GIL/dynamics/4kly_pD97" ## change based on your path   

# input and output, MM
set outfolder "${path}/bla/data"
set system    "4kly_pD97"
set dcdlisten "${path}/data.txt"

# create folders if not existing
file mkdir "${outfolder}"

# holder of averages and std. dev.
set averages {}
set stds     {}

#looping
for {set hallaq 0} {$hallaq < [llength $system]} {incr hallaq} {
    set input "[lindex $dcdlisten  $hallaq]"

    #variables for holding bond length alternation
    set bla {}

    # read psf/trajectory filenames
    set Eingabe [open "$input" r]
    set psffile "$[gets $Eingabe]"
    set dcdfiles {}
    while {true} {
        set dcdfiles "$dcdfiles [gets $Eingabe]"
        # check if eof reached
        if {[eof $Eingabe]} {
            break
        }
    }
    close $Eingabe

    #open trajectory psf
    set Mol [mol new "$psffile"]

    #get atom indeces for LYR296
    set sel5  [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C5"]
    set sel6  [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C6"]
    set sel7  [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C7"]
    set sel8  [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C8"]
    set sel9  [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C9"]
    set sel10 [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C10"]
    set sel11 [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C11"]
    set sel12 [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C12"]
    set sel13 [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C13"]
    set sel14 [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C14"]
    set sel15 [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name C15"]
    set sel16 [atomselect $Mol "segname DIM[string index [lindex $system $hallaq] end] and resid 248 and name N16"]
    set atom5  [$sel5  get index]
    set atom6  [$sel6  get index]
    set atom7  [$sel7  get index]
    set atom8  [$sel8  get index]
    set atom9  [$sel9  get index]
    set atom10 [$sel10 get index]
    set atom11 [$sel11 get index]
    set atom12 [$sel12 get index]
    set atom13 [$sel13 get index]
    set atom14 [$sel14 get index]
    set atom15 [$sel15 get index]
    set atom16 [$sel16 get index]
    $sel5  delete
    $sel6  delete
    $sel7  delete
    $sel8  delete
    $sel9  delete
    $sel10 delete
    $sel11 delete
    $sel12 delete
    $sel13 delete
    $sel14 delete
    $sel15 delete
    $sel16 delete

    foreach dcd $dcdfiles {
        set frame 0
        animate read $type $dcd beg $frame end $frame $Mol
        while {[molinfo $Mol get numframes] == 1} {
            #measure get average single/double
            set tmp1 "[measure bond "$atom6 $atom7"] [measure bond "$atom8 $atom9" ] [measure bond "$atom10 $atom11"] [measure bond "$atom12 $atom13"] [measure bond "$atom14 $atom15"]"
            set tmp2 "[measure bond "$atom5 $atom6"] [measure bond "$atom7 $atom8"] [measure bond "$atom9 $atom10"] [measure bond "$atom11 $atom12"] [measure bond "$atom13 $atom14"]"
            lappend bla  "[expr [mean $tmp1] - [mean $tmp2]]"

            animate delete all
            incr frame
            animate read $type $dcd beg $frame end $frame $Mol
        }
        animate delete all
        #puts "Done with ${dcd}."
    }
    mol delete $Mol

    #write distances
    set Ausgabe [open "${outfolder}/[lindex $system $hallaq].dat" w]
    puts $Ausgabe "BLA"
    set I 0
    foreach i $bla {
        puts $Ausgabe "$i"
        incr I
    }
    close $Ausgabe

    lappend averages "[lindex $system $hallaq] [mean $bla]"
    lappend stds     "[lindex $system $hallaq] [stddev $bla]"

}; #no more looping :(

#write averages
set Ausgabe [open "${outfolder}/000_averages.out" w]
puts $Ausgabe "System BLA"
foreach i $averages {
    puts $Ausgabe $i
}
close $Ausgabe

#write standard deviations
set Ausgabe [open "${outfolder}/000_stdevs.out" w]
puts $Ausgabe "System BLA"
foreach i $stds {
    puts $Ausgabe $i
}
close $Ausgabe

exit
