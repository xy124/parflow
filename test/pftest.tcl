set sig_digits 12

proc pftestFile {file message sig_digets} {
    set correct [pfload correct_output/$file]
    set new     [pfload                $file]
    set diff [pfmdiff $new $correct $sig_digets]
    if {[string length $diff] != 0 } {
	puts "FAILED : $message $diff"
	return 0
    } {
	return 1
    }
}



