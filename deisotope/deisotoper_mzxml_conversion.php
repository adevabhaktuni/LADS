
<?php
// Mod by krawson
// 2014.04.11
/*
	All extra file dependencies removed. All function calls are local and
	call flat files (.dta) rather than DB table. 
	Option added to analyze data from provided mzxml file. Need to specify 
	file type using -M options. Creates .dta files from extracted data
	from mzxml scan tags. 
	Accepts command line option to execute deisotpe method or simply unpack
	the peaks. 
	Command line option to define which levels you want to pay attention to
	between MS1, MS2, and MS3. defined by -l # to specify depth. 
*/

// Profiling call - delete before production!
//xdebug_start_trace('/home/krawson/Documents/Stanford/Projects/Deisotoper/PHP_debug/cachegrind.out.%t-%s');
//xdebug_start_trace();

//require_once("/var/www/html/gfy/lib/main_lib.php");
define('RATIO_LOG', getcwd().'/ratios.log');
//define('GREP_PATH', '/bin/grep');
//require_once(FILE_BASE_PATH."/www/modules/deisotope/lib/deisotope_lib.php");

/*
	MAIN
*/

// parse out the commandline options
$cli_options = getopt("hnf:o:dDl:MFu");

$noise_flag					= FALSE;
$debug_mode					= FALSE;
$major_debug				= FALSE;
$ms1_enable					= FALSE;
$ms2_enable					= TRUE;
$ms3_enable					= FALSE;
$mzxml_enable				= FALSE;
$unpack_only				= FALSE;

#Windows
$ppmWindow = 25; #ppm
$leftIsoWindow = 4.4; #amu
$rightIsoWindow = 2.2; #amu
$multWindow = $ppmWindow; #ppm

#Constants
$proton = 1.00727782;
$neutron = 13.00335483 - 12.00000000;

#Stop thresholds
$stopIntThres = 1; #multiple of median
$stopCountThres = 5; #how many more failures than successes to stop

#Print threshold
$min_count = 10;

#Ratio array
$open_ratio = array();

foreach($cli_options as $flag => $value) {
	switch($flag) {
		case 'h':
		print <<<EOF
		
deisotope.php -f <filename>
(c) 2010 by Konrad Karczewski
Adapted from Sean Beausoleil & The President and Fellows Of Harvard University

De-charges and de-isotopes dta file.

Command line options:
	h            Shows this help.
	n			 Add noise peaks back.
	f <filename> Filename of dta file
	o <filename> Output filename (defaults to input.deiso)
	d            Debug Mode
	p <ppm>      PPMWindow (default 25)
	D            Major Debug Mode

EOF;
		exit;
		break;
		case 'f':
			$filename = $value;
		break;
		case 'o':
			$out_file = $value;
		break;
		case 'n':
			$noise_flag = TRUE;
		break;
		case 'd':
			$debug_mode = TRUE;
		break;
		case 'D':
			$debug_mode = TRUE;
			$major_debug = TRUE;
		break;
		case 'l':
			if ($value == 1){
				$ms1_enable = TRUE;
				$ms2_enable = FALSE;
				$ms3_enable = FALSE;
			} elseif($value == 2){
				$ms2_enable = TRUE;
				$ms1_enable = FALSE;
				$ms3_enable = FALSE;
			} elseif($value == 3){
				$ms1_enable = TRUE;
				$ms2_enable = TRUE;
				$ms3_enable = FALSE;
			} elseif($value == 4){
				$ms3_enable = TRUE;
				$ms1_enable = FALSE;
				$ms2_enable = FALSE;
			} elseif($value == 5){
				$ms1_enable = TRUE;
				$ms2_enable = FALSE;
				$ms3_enable	= TRUE;
			} elseif($value == 6){
				$ms1_enable = FALSE;
				$ms2_enable = TRUE;
				$ms3_enable = TRUE;
			} elseif($value == 7){
				$ms1_enable = TRUE;
				$ms2_enable = TRUE;
				$ms3_enable = TRUE;
			} else {
				$ms1_enable = FALSE;
				$ms2_enable	= TRUE;
				$ms3_enable	= FALSE;
			}
		break;
        case 'p':
	         $ppmWindow = floatval($value);
		break;
		case 'M':
			$mzxml_enable = TRUE;
		break;
		case 'u':
			$unpack_only = TRUE;
		break;
		default:
		break;
	}
}


if($filename == NULL) {
		print <<<EOF
ERROR: Filename was not specified.

deisotope.php -f <filename>
(c) 2010 by Konrad Karczewski
Adapted from Sean Beausoleil & The President and Fellows Of Harvard University

EOF;
	exit;
}
if($out_file == NULL) {
	$out_file = $filename . ".deiso";
}

$path = dirname($out_file);
if(!is_dir($path)) {
	if (!(mkdir($path))) {
		client_error("Cannot create directory for storing unpacked dtas: " . $path . "\n", true, 21);
	}
}

# Load Ratio.log data
$open_ratio = loadRatios(RATIO_LOG);

# Determine file type and analysis type
$path = pathinfo($filename);

#Load mz/int data
if(is_readable($filename) && $path['extension'] == 'mzXML'){
	list($inMzInt, $precursor, $charge) = load_convert($filename);
} elseif(is_readable($filename) && $path['extension'] == 'dta') {
	list($inMzInt, $precursor, $charge) = loadData($filename);
	if ($major_debug){
		print "Started with: \n";
		print $precursor . "\t" . $charge . "\n";
		foreach ($inMzInt as $mz => $int){
			print $mz . "\t\t" . $int . "\n";
		}
	}
	$deIsoMzInt = deisotope($inMzInt, $precursor, $charge);
	if ($major_debug){
		print "Ended with: \n";
		print $precursor . "\t" . $charge . "\n";
		print "\n";
		foreach ($deIsoMzInt as $mz => $int){
			print $mz . "\t\t" . $int . "\n";
		}
	}
	writeData($deIsoMzInt, $precursor, $out_file, $min_count);
} elseif(is_dir($filename)){
	// print "DETECTED DTA DIRECTORY $filename\n";
	$files = array();
	foreach(glob($filename."/*.dta") as $file){
		$files[] = $file;
	}
	print count($files)." files detected\n";
	foreach($files as $curr){
		$path_dir = pathinfo($curr);
		// print $curr." ".$path_dir['extension']."\n";
		if(is_readable($curr) && $path_dir['extension'] == 'dta'){	
			list($inMzInt, $precursor, $charge) = loadData($curr);
			if ($major_debug){
				print "Started with: \n";
				print $precursor . "\t" . $charge . "\n";
				foreach ($inMzInt as $mz => $int){
					print $mz . "\t\t" . $int . "\n";
				}
			}
			$deIsoMzInt = deisotope($inMzInt, $precursor, $charge);
			if ($major_debug){
				print "Ended with: \n";
				print $precursor . "\t" . $charge . "\n";
				print "\n";
				foreach ($deIsoMzInt as $mz => $int){
					print $mz . "\t\t" . $int . "\n";
				}
			}
			// $out_file = $out_file."/".$curr;
			writeData($deIsoMzInt, $precursor, $curr, $min_count);
		
		} else {
			print "Unable to read $curr\n";
		}
	}
} else {
	print "ERROR: Cannot open provided file type.\n";
}

/*
	END OF MAIN
*/

////////////////////////////////////////////////////////////////////////////////
// ERROR METHODS
////////////////////////////////////////////////////////////////////////////////

/**
 *	register_error_ini_directives : registers directives for error reporting from PHP errors
 *
 *	@return	void	No return
 *
 */
function register_error_ini_directives() {
	// Turn off Display Errors for a production site
	ini_set("display_errors",TRUE);
	// Turn off Log Errors for maximum code security
	ini_set("log_errors",TRUE);
	// Log all of the errors in GFY to the GFY Logfile
	ini_set("error_log",LOGFILE);
}

/**
 *	client_error : raises an error to the user and halts execution if set
 *
 *	Shows an error message to the user and halts execution
 *  Typical call:  client_error("Error Message", TRUE);
 *
 *  @todo	client_error should be able to handle CLI and web-based errors
 *  @todo	client_error should make HTML-formatted errors too
 *
 *	@param	string	$error_string	Text of the error to raise
 *  @param	boolean	$flg_halt		Flag to halt execution if set. (FALSE if not set)
 *	@param	int		$errno			Error number to return to a calling program; not required (will default to FALSE)
 *	@return	void	Prints error out to stderr if halting, or stdout if not.
 *
 */
function client_error( $error_string, $flg_halt = FALSE, $errno = FALSE ) {
	log_message("GFY CLIENT ERROR: ".$error_string);
	if ($flg_halt) {
		log_message("GFY FATAL CLIENT ERROR: ".$error_string);
		print "FATAL CLIENT ERROR: ".$error_string."\n";
		if($errno && (is_numeric($errno) && $errno >= 1 && $errno <= 254)) {
			exit($errno); // exit with a nonzero status supplied, so that other programs can pick this up.
		} else {
			exit(127); // exit with a nonzero status, so that other programs can pick this up.
		}
	} else {
		log_message("GFY CLIENT WARNING: ".$error_string);
		$message = sprintf(<<<EOT

	<table width="90%%" align="center" cellpadding="0" cellspacing="0" class="warning-box">
			<tr>
				<td><img src="%s/images/icons/messagebox_warning.png"></td>
				<td>%s</td>
			</tr>
	</table>
			
EOT
, HTTP_BASE_PATH, $error_string );
		print $message;
	}
}

/**
 *	log_message : logs a message to the msd/gfy log file
 *
 *	Logs a message into the msd/gfy master log file
 *
 *  @todo	log_message should be configurable as to where it can send the error
 *  @todo	log_message needs to be written!
 *
 *	@param	string	$log_string	Text of the error to raise
 *	@return	boolean	True/False if logging was successful
 */
//if(!function_exists("log_message")) {
	function log_message( $log_string ) {
		// this needs to be better used!
		// probably use openlog and the ilk
		
		// uncomment the next line to display log strings in the current stdout
		//print $log_string ."\n" ;
		$server_name = trim(`/bin/hostname`);

		if(defined('LOGFILE') && is_writable(LOGFILE)){
			// grab some info for the logging
			$datetime 	= date("d-M-Y H:i:s");
			$user_id	= (isset($_SESSION['username']) ? $_SESSION['username'] : "X");
			$ip_addr 	= (isset($_SERVER['REMOTE_ADDR']) ? $_SERVER['REMOTE_ADDR'] : "localhost");
			$pid		= (function_exists("posix_getpid") ? posix_getpid() : "NoPID");

			error_log( "[".$datetime . "] @$server_name ".$ip_addr." - ".$pid." - ".$user_id." - ".$log_string."\n", 3, LOGFILE );
		}
		return TRUE;
	}
//}

////////////////////////////////////////////////////////////////////////////////
// Deisotoper methods
////////////////////////////////////////////////////////////////////////////////

/*
	DATA RETRIEVAL METHODS
*/

// Load and convert mzxml immediately
// Complete analysis on given levels, execute deisotope if option has been selected
function load_convert($filename){

	$path = pathinfo($filename);
	global $min_count, $ms1_enable, $ms2_enable, $ms3_enable, $mzxml_enable, $out_file, $unpack_only;
	print $filename."\n";
	$scanned = simplexml_load_file($filename);
	$preIntensity = -1;
	$preCharge = -1;
	$mzLoad = array();

	$count = 0;
	foreach($scanned->msRun->children() as $a){
		$preIntensity = -1;
		$preCharge = -1;
		$mH = 0;
		//print "Curr $count\n";
		if($a->getName() == 'scan'){
			$attrs = $a->attributes();

			if($attrs['msLevel'] == 2 && $ms2_enable){
				if($a->children()->getName() == 'precursorMz'){
					$b = $a->precursorMz;
					$battrs = $b->attributes();

					if($preIntensity == -1){
						$preIntensity = $battrs['precursorIntensity'];
					}	

					if($preCharge == -1){
						$preCharge = $battrs['precursorCharge'];
					}	

					$mH = (float)$b * $preCharge;
					$mH = round($mH, 9);			
				}
			} elseif($attrs['msLevel'] == 1 && $ms1_enable){
				$preIntensity = 0;
				$preCharge = 0;
			}

			$high = FALSE;
			if($ms2_enable || $ms3_enable){
				$high = TRUE;
			}

			$ms_id = (int)$attrs['num'];

			if($mH > 0 && $high){
				$label = $path['filename'].".".$ms_id.".".$count.".".$preCharge;
			} else {
				$label = $path['filename'].".".$ms_id."."."ms1";
			}

			if($a->peaks){
				$base = base64_decode($a->peaks);
				if($attrs['msLevel'] == 1 && $ms1_enable){
					//print "Decoding MS1 peaks\n";
					$mzLoad = mzxml_decode_peaks($a->peaks);
				} 
				if ($attrs['msLevel'] == 2 && $ms2_enable) {
					//print "Decoding MS2 peaks\n";
					$mzLoad = mzxml_decode_peaks($a->peaks);
				} 
				if ($attrs['msLevel'] == 3 && $ms3_enable){
					print "Given MS level ".$attrs['msLevel']."not currently supported";
				}

				if(!empty($mzLoad) || $mzLoad == NULL){
					if($ms1_enable && $mH == 0 && $attrs['msLevel'] == 1){
						//print "Unpacking MS1 scan\n";	
						$mzLoad = refine_ms1($mzLoad);
						if(!$unpack_only){
							$mzLoad = deisotope($mzLoad, $preIntensity, $preCharge);
						}
						writeData($mzLoad, "ms1", $out_file.$label.".dta", $min_count);
						$count++;
					} 
					if ($ms2_enable && $mH > 0 && $attrs['msLevel'] == 2){
						//print "Unpacking MS2 scan\n";
						$mzLoad = refine($mzLoad);
						if(!$unpack_only){
							$mzLoad = deisotope($mzLoad, $preIntensity, $preCharge);
						}
						writeData($mzLoad, $mH, $out_file.$label.".dta", $min_count);
						$count++;
					}
					if ($ms3_enable && $mH > 0 && $attrs['msLevel'] == 3){
						$mzLoad = refine($mzLoad);
						if(!$unpack_only){
							$mzLoad = deisotope($mzLoad, $preIntensity, $preCharge);
						}
						writeData($mzLoad, $mH, $out_file.$label.".ms3.dta", $min_count);
						$count++;
					}
				} else {
					print "ERROR: Peaks section could not be decoded.\n";
				}
			} else {
				print "ERROR: No peaks provided in current scan element.\n";
			}		
		}
	}
	print "Number of scans reported $count\n";
}

// Refine the current mzLoad array. Reduces decimal places on both key and values
function refine($mzLoad){
	$ret = array();
	while ($curr = current($mzLoad)) {
        $curr_key = number_format((float)key($mzLoad), 5, ".", "");
        $curr_val = number_format((float)$curr, 2, ".", "");
    	$ret["$curr_key"] = $curr_val;
    	next($mzLoad);
	}
	return $ret;
}

// Refine current mzLoad array similar to one above, specifically for ms1 data
function refine_ms1($mzLoad){
	$ret = array();
	while ($curr = current($mzLoad)) {
    
        $curr_key = number_format((float)key($mzLoad), 9, ".", "");
        $curr_val = number_format((float)$curr, 3, ".", "");
    
    	$ret["$curr_key"] = $curr_val;
    	next($mzLoad);
	}
	return $ret;
}

// read data from dta file, check to make sure that file is dta extension
function loadData($fn){
  $fh = fopen($fn, 'r');
  $mzLoad = array();
  $first = explode("\t", trim(fgets($fh)));
  $prec = $first[0];
  $prec_charge = $first[1];
  while(!feof($fh)){
    $line = fgets($fh);
    if (preg_match("/\d/", $line)){
      $split = explode("\t", trim($line));
      $mzLoad[$split[0]] = $split[1];
    }
  }
  fclose($fh);
  return array($mzLoad, $prec, $prec_charge);
}

// write de-isotoped, de-charged dta file
function writeData($out, $precursor, $fn, $min_count){
  if (count($out) >= $min_count){
    print $fn . "\n";
    $fh = fopen($fn, 'w') or die("Cannot write to: " . $fn . "\n");
    $firstline = $precursor . "\t1\n";
    fwrite($fh, $firstline);
    ksort($out);
    foreach ($out as $mz => $int){
      $str = $mz . "\t" . $int . "\n";
      fwrite($fh, $str);
    }
    //print "$fn\n";
    fclose($fh);
  }
}

// Open the ratio log and store in array for faster analysis in deisotope stage
function loadRatios($ratio_log){
	$opened_ratios = array();
	$file = fopen($ratio_log, 'r');
	while(!feof($file)){
		if(empty($file)){
			print "EMPTY LINE";
		}
		$line = fgets($file);
		$opened_ratios[] = $line;
	}
	return $opened_ratios;
}

/*
	DECODE METHODS
*/
/**
     * Decodes peak information input to function (ms2_peaks) from base64 to binary and saves them in the object
     *
     * @access public
     * @see mzxml_decode_peaks_to_binary
     */  
  	function decode_ms2_peaks( $base64_peak_information ) {
	  /*
	  print "decoding peaks: input = \n";var_dump($base64_peak_information);
	  */
  		$ms2 = mzxml_decode_peaks_to_binary($base64_peak_information);

		return $ms2;
	}

	/**
 *	mzxml_decode_peaks_to_binary : take in a base64 string from an mzXML file and convert to a binary representation
 *
 *	Load in base64 string and then convert this to a binary representation in network byte order
 *	(Essentially, just base64_decode the string) -- to decode this, simply pick up after the base64_decode
 *	This function is called to store data in the dd, and non-base64 decoded files are 33% smaller.
 *
 *	@param	STRING	$base64_data		Bas64 encoded data representing full scan information
 *	@return	INT (BINARY)				Binary peak data representation
 *
 */
function mzxml_decode_peaks_to_binary($base64_data) {

    return base64_decode($base64_data);
        
}

/**
 *	mzxml_decode_peaks : take in a base64 string from an mzXML file and convert to mz => int pairs
 *
 *	Load in base64 string and then convert this to mz => int pairs
 *	Encoding and decoding are precise to the second to last number
 *
 *	@param	STRING	$base64_data		Bas64 encoded data representing full scan information
 *	@return	ARRAY						array of peaks as {m/z => int}
 *
 */
function mzxml_decode_peaks ($base64_data) {

	$host_order_32_array = unpack("N*", base64_decode($base64_data));

	$result_array = array();

	while(TRUE) {

		// try to speed up
		reset($host_order_32_array);
		list($mz_key, $shifted_mz_val) = each($host_order_32_array);
		list($int_key, $shifted_int_val) = each($host_order_32_array);
		unset($host_order_32_array[$mz_key], $host_order_32_array[$int_key]);


		$mz_array = unpack("f", pack("I", $shifted_mz_val));

		$intensity_array = unpack("f", pack("I", $shifted_int_val));

		$result_array[(string)$mz_array[1]] = $intensity_array[1];

		if(count($host_order_32_array) == 0) {
			break;
		}

	}

	return $result_array;
}

/**
 *	mz_to_mh : take m/z value and convert it to an M+H value
 *
 *	@author Corey Bakalarski
 *	@param	FLOAT	$mz		            m/z value
 *	@param	INT 	$charge	            charge
 *	@return	FLOAT						M+H value
 *
 */

function mz_to_mh ($mz, $charge) {
    global $g_elemental_mass;
    
    //this function converts an m+H value to an m/z value
    return sprintf("%.05f",(($mz * $charge) + ($g_elemental_mass['proton']) - ($charge * $g_elemental_mass['proton'])));

}

/*
	DEISOTOPE ANALYSIS METHODS
*/

// Get Charge State
function get_charge_state($avg_peaks, $accurate_parent, $ms2_charge, $max_charge, $debug_mode = FALSE) {
	//spacing
	$c13_spacing = 13.00335483 - 12.00000000;

	$chi_thresh = 1; // if best chi square has this value or greater, make .dta files for both 2,3;

	$mono_arr = array();
	$charge_arr = array();

	if ($debug_mode) {
		print "starting with these averaged peaks surrounding $accurate_parent:\n";
		var_dump($avg_peaks);
	}

	$median_int = median(array_values($avg_peaks));
	//loop through to remove peaks < 20 % of median
	foreach ($avg_peaks as $mz => $int) {
		if ($int/$median_int < 0.01 ) {
			unset($avg_peaks[$mz]);
		}
	}

	if ($debug_mode) {
		print "took out extraneous peaks:\n";
		var_dump($avg_peaks);
	}

	//create theoretical isotopic envelope peaks from parent for each charge state considered
	$env_arr = array();
	/*
	 print "filling envelope array:\n";
	 */
	//only fill with known charge state if known, otherwise do 2-6
	if ($ms2_charge == NULL) {
		//print "charge wasn't set in xml file so trying from 2-6\n";
		$low_charge = 1;
		$high_charge = $max_charge;
	}
	else {
		//print "charge set in xml file so using $ms2_charge\n";
		$low_charge = $ms2_charge;
		$high_charge = $ms2_charge;
	}
	for ($z=$low_charge; $z <=$high_charge; $z++) {
		/*
		 print "working on charge = $z\n";
		 */
		$theo_iso_arr = array();
		for ($offset = -5; $offset <= 5; $offset++) {
			$theo_iso_arr[$offset] = $accurate_parent + (($c13_spacing*$offset)/$z);
		}
		/*
		 print "theoretical array:\n";
		 var_dump($theo_iso_arr);
		 */

		// JE:  want to knit theo_iso_arr with avg_peaks together such that there is an entry for all envelope offests.

		$offset_ind = 0;
		$avg_peak_ind = 0;
		$offset_arr = array_keys($theo_iso_arr);
		$avg_peak_arr = array_keys($avg_peaks);
		$add_zero_int_flg = 0; // flag denoting if gaps between first peak and subsequent peaks should be filled with zero intensity values

		/*
		 print "count offset array: ".count($offset_arr)."\n";
		 print "count avg peak array: ".count($avg_peak_arr)."\n";
		 */

		while ($offset_ind < count($offset_arr) && $avg_peak_ind < count($avg_peak_arr)) {

			$offset = $offset_arr[$offset_ind];
			$theo_mz = $theo_iso_arr[$offset];
			/*
			 print "offset index: ".$offset_ind."\n";
			 print "avg peak index: ".$avg_peak_ind."\n";
			 */

			$obs_mz = $avg_peak_arr[$avg_peak_ind];
			$obs_int = $avg_peaks[$obs_mz];

			/*
			 print "obs m/z: ".$obs_mz."\n";
			 print "obs int: ".$obs_int."\n";
			 print "theo mz: ".$theo_mz."\n";
			 */

			$min = $theo_mz - calc_diff($theo_mz,25);
			$max = 2*$theo_mz-$min;
			/*
			 print "matching $obs_mz to $min and $max around $theo_mz\n";
			 */
			if ($obs_mz < $min) {
				/*
				 print "$obs_mz < $min...skipping\n";
				 */
				$avg_peak_ind++;
			}else // observed m/z >= min
			if ($obs_mz <= $max) { // match
				/*
				 print "found match\n";
				 */

				if (
				(isset($env_arr[$z][$offset]["mz"]) && $obs_int > $env_arr[$z][$offset]["int"]) // mz already found and observed int > prior
				||
				(!(isset($env_arr[$z][$offset]["mz"]))) // mz not yet found
				) { //update
					$env_arr[$z][$offset]["mz"] = (float)$obs_mz;
					$env_arr[$z][$offset]["int"] = (int)$obs_int;

					if ($debug_mode) {
						print "Matched $obs_mz with $theo_mz for charge $z\n";
					}
				}
				$avg_peak_ind++;

			}else{ // observed > max
				/*
				 print "no match...";
				 */
				if (isset($env_arr[$z][$offset]["mz"])) { // found a value for this offset
					// do nothing; just go to the next offset
					/*
					 print "already found a match...skipping\n";
					 */
				}else{ // add zero intensity placeholder
					/*
					 print "no match found...adding zero intenzity for $theo_mz\n";
					 */
					$env_arr[$z][$offset]["mz"] = $theo_mz;
					$env_arr[$z][$offset]["int"] = 0;
				}

				$offset_ind++;
			}
		}
	}

	//trim off trailing 0 intensity values from env_arr//
	foreach ($env_arr as $z => $array) {
		/*
		 print "charge: $z\n";
		 print "array:\n";
		 var_dump($array);
		 */

		$offset_arr = array_keys($array);

		/*
		 print "offsets:\n";
		 var_dump($offset_arr);
		 */

		if (count($offset_arr)>1) {
			$offset = array_pop($offset_arr);
			$last_offset = array_shift($offset_arr);

			$stop_flag = 0;
			while ($stop_flag == 0 && $offset > $last_offset) {
				if ($env_arr[$z][$offset]["int"] == 0) {
					array_pop($env_arr[$z]);
					$offset--;
				}else{
					$stop_flag = 1;
				}
			}
		}
	}

	///////////////////////////////////////////////////////

	/*
	 if ($debug_mode) {
	 print "envelope array:\n";
	 var_dump($env_arr);
	 }
	 */
	$count = $max = $series = $max_series = 0;
	$monoscore_arr = array();
	//loop through env_arr to determine most likely charge state
	if ($ms2_charge == NULL) {
		foreach ($env_arr as $z => $pk_arr) {
			//sort so the peaks are ascending
			ksort($pk_arr);

			// count consecutive non-zero intensity peaks surrounding parent ion in series given charge state
			$count= 0;

			if ($debug_mode) {
				print "counting for charge $z:\n";
				/*
				 print "peaks array is\n";
				 var_dump($pk_arr);
				*/
			}

			$max_offset = 5; // look for envelope peaks up to 5 away from parent ion
			$min_range = $max_range = 0;
			foreach (array(-1,1) as $sign) {
				$done_flag = 0;
				/*
				 print "looking in $sign direction\n";
				 */
				$offset=1;
				while ($offset <= $max_offset && $done_flag == 0 ) {
					/*
					 print "offset is $offset\n";
					 */
					if (isset($pk_arr[$offset*$sign])) {
						$pks = $pk_arr[$offset*$sign];
						/*
						 print "pks: \n";
						 var_dump($pks);
						 */
					}
					if (isset($pks) && $pks["int"] >0 && isset($pk_arr[$offset*$sign])) {
						$offset++;
					}else{
						/*
						 print "found last item...";
						 */
						// found zero intensity; we're done
						if ($sign == 1) {
							$max_range = max(0,$offset-1);
							/*
							 print "max range is $max_range\n";
							 */
						}else{ // sign = -1
							$min_range = min(0,-$offset+1);
							/*
							 print "min range is $min_range\n";
							 */
						}
						$done_flag = 1;
					}
				}
				if ($done_flag == 0) { //reached end of array without finding empty element or zero intensity
					if ($sign == 1) {
						$max_range = max(0,$offset-1);
						/*
						 print "max range is $max_range\n";
						 */
					}else{ // sign = -1
						$min_range = min(0,-$offset+1);
						/*
						 print "min range is $min_range\n";
						 */
					}
				}
			}
			$count = $max_range - $min_range;

			//check tests -- count >= z-1, series of ions

			if ($debug_mode) {
				print "count for charge $z is $count\n";
			}

			// @todo: move magic numbers to constants file
			$mw = ($accurate_parent * $z) - (1.00728 * $z);
			$num_residues = $mw / 111;
			$num_carbons = (int)( $num_residues * 5.1 );
			$peaks_arr = $env_arr[$z];
			$ratio_arr = fetch_isotopic_envelope_ratios($num_carbons);

			//Relax stringency a little, but don't allow 1+ to win just bc there's 1+
			if ($count >= $z-2 && $count > 0) { // successful identification of potential charge state; find most likely mono given charge
				if ($debug_mode) {
					print "finding mono for charge = $z\n";
				}
				list ($mono,$score) = find_mono($peaks_arr, $ratio_arr, $debug_mode);
				/*
				 print "received mono $mono with score $score\n";
				 */
				$monoscore_arr[(string)$score]=array($z=>$mono);
				/*
				 print "current monoscore_arr:\n";
				 var_dump($monoscore_arr);
				 */
			}
		}
	}
	else {
		$z = $ms2_charge;
		$mw = ($accurate_parent * $z) - (1.00728 * $z);
		$num_residues = $mw / 111;
		$num_carbons = (int)( $num_residues * 5.1 );
		$peaks_arr = $env_arr[$z];
		$ratio_arr = fetch_isotopic_envelope_ratios($num_carbons);
		list ($mono,$score) = find_mono($peaks_arr, $ratio_arr, $debug_mode);
		$monoscore_arr[(string)$score]=array($z=>$mono);
	}
	/*
	 print "final monoscore array:\n";
	 var_dump($monoscore_arr);
	 print "found scores for ".count($monoscore_arr)." charge states\n";

	 if (count($monoscore_arr) > 0) {print "best score was ".min(array_keys($monoscore_arr))."\n";}else{print "no score yet\n";}
	 */
	if (count($monoscore_arr) == 0 || min(array_keys($monoscore_arr)) > $chi_thresh) { // not enough peaks to conclusively determine charge; make both 2 and 3;
		if ($debug_mode) {
			print "making dtas for default charge states\n";
		}
		$charge_arr[] = 0;
		$mono_arr[] = 0;
//		if ($ms2_charge == NULL) {
//			foreach (array(2,3) as $z) {
//				$mw = ($accurate_parent * $z) - (1.00728 * $z);
//				$num_residues = $mw / 111;
//				$num_carbons = (int)( $num_residues * 5.1 );
//				$peaks_arr = $env_arr[$z];
//				$ratio_arr = fetch_isotopic_envelope_ratios($num_carbons);
//
//				list($mono, $score) = find_mono($peaks_arr, $ratio_arr, $debug_mode);
//				$charge_arr[] = $z;
//				$mono_arr[] = $mono;
//			}
//		}
//		else {
//			//just make dta for the charge we know
//			$z = $ms2_charge;
//			$mw = ($accurate_parent * $z) - (1.00728 * $z);
//			$num_residues = $mw / 111;
//			$num_carbons = (int)( $num_residues * 5.1 );
//			$peaks_arr = $env_arr[$z];
//			$ratio_arr = fetch_isotopic_envelope_ratios($num_carbons);
//
//			list($mono, $score) = find_mono($peaks_arr, $ratio_arr, $debug_mode);
//			$charge_arr[] = $z;
//			$mono_arr[] = $mono;
//		}
	} else { // made at least one mono assignment

		krsort($monoscore_arr);
		$best = -1;
		foreach ($monoscore_arr as $score => $z_mono_arr) {
			if ($best == -1 && $score != -1) {
				$best = $score;
				/*
				 print "best score is $best\n";
				 */
			}
			if (abs($best - $score)/$best <= 0.2 || $score > 0.8) { // make dta for all charge states if the score is within 20% of the best, or score is above 0.8 (good match between observed and expected)
				foreach ($z_mono_arr as $z => $mono) {
					$charge_arr[]=$z;
					$mono_arr[]=$mono;
				}
			}
		}
	}
	if ($debug_mode) {

		print "making .dta files for the following mono mass(es):\n";
		var_dump($mono_arr);
		print "mono masses correspond to following charge state(s):\n";
		var_dump($charge_arr);
	}

	return array($mono_arr, $charge_arr);
}

// Median
function median ($list) {
	if (count($list)>2) {
		sort ($list);
		$mid = count($list)/2;
		if ($mid == floor($mid)) {// even; use average of two middle numbers
			return (($list[$mid]+$list[$mid+1])/2);
		}else{
			return ($list[$mid]);
		}
	}else if (count($list) == 2) {
		return (array_sum($list)/2);
	}else if (count($list) == 1) {
		return $list[0];
	}else{
		client_error("Cannot take median of nothing",FALSE);
		return NULL;
	}
}

// Fetch Isotopic Envelope Ratios
function fetch_isotopic_envelope_ratios($peptide_carbons) {


	// @todo: add more carbons to the ratios file
	global $open_ratio;
	//print "Open ratios array size ".count($open_ratio)."\n";

	$max_carbons = 400;
	$peptide_carbons = min($peptide_carbons,$max_carbons);

	//print "peptide_carbons $peptide_carbons\n";
	// load in the appropriate line from the file
	//$command = GREP_PATH." '^$peptide_carbons,' ".RATIO_LOG;
	//print "command $command\n";
	//$ratio_file_line = `$command`;
	$ratio_file_line = $open_ratio[$peptide_carbons-1];
	//print "ratio_file_line $ratio_file_line\n";
	//print "Retrieved ration from array ".$open_ratio[$peptide_carbons-1]."\n";
	$line_array = explode(",",$ratio_file_line);

	if (isset($line_array[5])) {
		$ratios_array = explode(" ",$line_array[5]);
	}else{
		client_error("grep command ".$command." failed to return the expected ratios array",TRUE);
	}
	return $ratios_array;
}

/*
 Function which compares the the expected isotopic ratios to the observered ratios
 Returns mono
 */

/*
 Function which compares the the expected isotopic ratios to the observered ratios
 Returns mono
 */

function find_mono($peaks_arr, $ratio_arr, $debug_mode = FALSE) {

	$mono_arr = array();
	$chi_arr = array();

	//make array of observed intensities
	$obs_arr = array();
	foreach ($peaks_arr as $offset=> $mz_int_arr) {
		$obs_arr[$offset] = $mz_int_arr["int"];
	}
	$exp_arr = $ratio_arr;

	$less_than = 5; //number of peaks to the left of selected ion to consider as mono
	$greater_than = 1; // number of peaks to the right of selected ion to consider as mono
	$starting_max_peaks = 5; // maximum numbers of peaks to consider in isotopic envelope at once


	// sum all peaks in frame with selected ion for normalization purposes (Willi's idea)

	/*
	 $sum_obs_int = 0;
	 foreach ($peaks_arr as $offset => $mz_int_arr) {
	 $sum_obs_int += $mz_int_arr["int"];
	 }
	 */

	for ($i = -$less_than; $i <= $greater_than; $i++) {

		/*
		 print "working on new peaks array item ".$i."\n";
		 */

		/*
		 print "offset is ".$i.", corresponding with m/z = ".$peaks_arr[$i]['mz']."\n";
		 */
		if (isset($peaks_arr[$i])) {
			if ($peaks_arr[$i]["int"] == 0) {
				/*
				 print "item ".$i." does not exist\n";
				 */

			}
			else { //generate slice of obs_array for calculating scores



				//with long distributions, sometimes extraneous peaks can creap into the frame set by the selected ion and charge state.  If there are more than three possible observed peaks, find mono given 3, 4 and 5 peaks, select best, as long as selected ion is within distribution
				$max_peaks = $starting_max_peaks;
				$proceed_flag = 0;
				while ($max_peaks >= 3 && $proceed_flag == 0 && $i + $max_peaks-1 >= 0) {
					/*
					 print "looking for up to $max_peaks peaks in distribution\n";
					 */

					$sub_obs_arr = array();
					$sub_exp_arr = array();
					$willi_chi_obs_arr = array();
					$chi_obs_arr = array();

					$offset = $i;
					$exp_ind = 0;
					$zero_flag = 0; // allow at most one zero at end of observed distribution

					while (isset($obs_arr[$offset]) && count($sub_obs_arr) < $max_peaks && $zero_flag == 0) { // keep observed array size to five, to avoid unreliable assignments


						if ($obs_arr[$offset] == 0) {
							$zero_flag++;
						}


						if (($zero_flag == 1 && count($sub_obs_arr)<=2) || $zero_flag == 0) { // only add zero if there's 2 or fewer peaks
							$sub_obs_arr[] = $obs_arr[$offset];

							$sub_exp_arr[] = $exp_arr[$exp_ind];

						}
						$offset++;
						$exp_ind++;

					}

					// pad observed array with one zero if there's only two or fewer peaks
					if (count($sub_obs_arr) <= 2 && $sub_obs_arr[count($sub_obs_arr)-1] != 0) {
						$sub_obs_arr[] = 0;
						$sub_exp_arr[] = $exp_arr[$exp_ind];


					}
					//note that we should proceed if there's 3 or fewer peaks prior to padding
					if (count($sub_obs_arr) <= 3 || (count($sub_obs_arr) == 4 && $sub_obs_arr[count($sub_obs_arr)-1] == 0)) {
						$proceed_flag=1;
					}else if (count($sub_obs_arr) < $starting_max_peaks && $max_peaks == $starting_max_peaks) {
						$max_peaks = count($sub_obs_arr); // i.e., 4; don't want to repeat score finding if we don't have to
					}


					// normalize observed array to all considered peaks in frame with selected ion (Willi's idea) (excluding trailing peaks (My idea) or by sum of sub-array (original idea); prepend observed array with observed peaks not assumed to be in distribution, given offset and mono; prepend expected distribution with zero
					$sum_sub_obs_int = array_sum($sub_obs_arr);

					foreach ($sub_obs_arr as $idx=>$int) {
						$willi_chi_obs_arr[]=$int;
						$chi_obs_arr[] = $int/$sum_sub_obs_int;
					}
					$idx = $i-1;
					$willi_sub_exp_arr = $sub_exp_arr;
					$stop_flg = 0;
					//tack back on the unexpected peaks in the front
					while ($idx >= -$less_than && $stop_flg == 0) {
						if (isset($peaks_arr[$idx])) {
							array_unshift($willi_chi_obs_arr,$peaks_arr[$idx]["int"]);
							array_unshift($willi_sub_exp_arr,0);
							$idx--;
						}else{
							$stop_flg=1;
						}
					}

					// normalize willi chi
					$sum_willi_sub_obs_int = array_sum($willi_chi_obs_arr);
					foreach ($willi_chi_obs_arr as $idx => $int) {
						$willi_chi_obs_arr[$idx] = $int/$sum_willi_sub_obs_int;
					}


					if ($debug_mode) {
						print "proceeding with \n";
						print "expected array:\n";
						var_dump($sub_exp_arr);
						print "observed array:\n";
						var_dump($sub_obs_arr);
						print "willi norm array\n";
						var_dump($willi_chi_obs_arr);
						print "willi expected array\n";
						var_dump($willi_sub_exp_arr);
						print "orig norm array\n";
						var_dump($chi_obs_arr);
					}

					$chi = calc_chi($chi_obs_arr, $sub_exp_arr, $debug_mode);
					// $correl = correl($sub_obs_arr, $sub_exp_arr);
					// $w_correl = correl($willi_chi_obs_arr,$willi_sub_exp_arr);
					// list($score, $m, $r2) = calc_slope($chi_obs_arr, $sub_exp_arr, $debug_mode);
					$willi_chi = calc_chi($willi_chi_obs_arr, $willi_sub_exp_arr, $debug_mode);
					list($w_score, $w_m, $w_r2) = calc_slope($willi_chi_obs_arr, $willi_sub_exp_arr, $debug_mode);

					$mz = $peaks_arr[$i]["mz"];
					$int = $peaks_arr[$i]["int"];
					// if ($debug_mode) {
					// 	print "\nchi for $mz was $chi, correl was $correl, score was $score with slope $m and r2 $r2\n";
					// 	print "willichi was $willi_chi, w_correl was $w_correl, w_score was $w_score, w_m was $w_m, w_r2 was $w_r2\n";
					// }

					$mono_arr[(string)$w_score][(string)$mz]["peaks"] = count($chi_obs_arr)-$zero_flag;
					// $mono_arr[(string)$w_score][(string)$mz]["correl"] = $correl;
					$mono_arr[(string)$w_score][(string)$mz]["chi"] = $chi;
					$mono_arr[(string)$w_score][(string)$mz]["w_chi"] = $willi_chi;
					// $mono_arr[(string)$w_score][(string)$mz]["score"] = $score;
					// $mono_arr[(string)$w_score][(string)$mz]["w_correl"]=$w_correl;

					// $chi_arr[(string)$w_score][(string)$mz]["peaks"] = count($chi_obs_arr)-$zero_flag;
					// $chi_arr[(string)$w_score][(string)$mz]["chi"] = $chi;
					// $chi_arr[(string)$w_score][(string)$mz]["w_chi"] = $willi_chi;

					$max_peaks--;
				}
			}
		}
	}

	// if there's multiple best chi scores for alternate mono's, compare what the other scores have to say.

	// krsort($chi_arr);

	// $chi_scores = array_keys($chi_arr);
	// // print "CHi scores\n";
	// // print_r(array_values($chi_scores));

	// $best = -1;
	// $mono = -1;

	// foreach($chi_scores as $chi_single){
	// 	if($best == -1){
	// 		$best = $chi_single;
	// 	}
	// 	foreach($chi_arr[$chi_single] as $mz => $chi_scores_arr){
	// 		if($mono == -1){
	// 			$mono = $mz;
	// 		}

	// 		if ($chi_single >= $best){
	// 			if($mz >= $mono){
	// 				print "New vals $chi_single $mz\n";
	// 				$best = $chi_single;
	// 				$mono = $mz;
	// 			}
	// 		}
	// 		//print "Current score $chi_single and mz $mz\n";
	// 	}

	// }

	krsort($mono_arr);
	$score_arr = array_keys($mono_arr);

	$best = -1;
	$mono = -1;
	foreach ($score_arr as $score) {
		if ($best == -1) {
			$best = $score;
		}
		foreach ($mono_arr[$score] as $mz => $scores_arr) {
			$tally = 0;
			if ($mono == -1) {
				$mono = $mz;
			}
			/*
			 print "best is $best; mono is $mono; corresponding scores array:\n";
			 var_dump($scores_arr);

			 */
			if ($best != -1 && $best == $score && $mono == $mz) {
				$mono_arr[$best][$mono]["tally"] = floor(count($mono_arr[$best][$mono])/2)+1; // tally score to beat is half the number of scores + 1; for the best score to be upstaged, other monos must have better values for more than half of the scores; this procedure shouldn't count "tally" as one of the scores
			}else
			if ($best != -1 && ((abs($best - $score)/$best < 0.2)) && $mz != $mono) { // evaluate multiple possible monos if they have similar correlation scores (20% similarity)
				if ($debug_mode) {
					print "close one...comparing $mz  with score $score to $mono with score $best\n";
				}

				foreach ($scores_arr as $scorename => $val) {
					if ($scorename == "tally") {
						// do nothing; don't use it!
					}else
					if ($scorename == "w_chi" || $scorename == "chi") { // less is better
						if ($val < $mono_arr[$best][$mono][$scorename]) {
							if ($debug_mode) {
								print "new $scorename ($val) is better than old value (".$mono_arr[$best][$mono][$scorename].")\n";
							}
							$tally++;
						}else
						if ($debug_mode) {
							print "old $scorename (".$mono_arr[$best][$mono][$scorename].") is better than new value ($val)\n";
						}

					}else{ // greater is better
						if ($val > $mono_arr[$best][$mono][$scorename]) {

							if ($debug_mode) {
								print "new $scorename ($val) is bettern than old value (".$mono_arr[$best][$mono][$scorename].")\n";
							}
							$tally++;
						}else
						if ($debug_mode) {
							print "old $scorename (".$mono_arr[$best][$mono][$scorename].") is better than new value ($val)\n";
						}

					}
				}

				if ($tally > $mono_arr[$best][$mono]["tally"]) {

					if ($debug_mode) {
						print "new mono $mz wins: tally $tally beats ".$mono_arr[$best][$mono]["tally"]."\n";
					}

					$best = $score;
					$mono = $mz;
					$mono_arr[$best][$mono]["tally"] = $tally;


				}else if ($debug_mode) {
					print "old mono $mono remains: tally ".$mono_arr[$best][$mono]["tally"]." beats $tally\n";
				}
			}
		}
	}



	/*
	 //both scores agree
	 if ($correl > $best_score && $chi < $best_chi || $best_chi == 0) {
	 $mono = $mz;
	 $best_score = $correl;
	 $best_chi = $chi;
	 if ($debug_mode) {
	 print "setting mono to $mono because of a better chi and correl\n";
	 }

	 //chi is better but score isn't -- favor score
	 else if ($correl > $best_score && $chi > $best_chi || $best_chi == 0) {
	 $mono = $mz;
	 $best_chi = $chi;
	 $best_score = $correl;
	 if ($debug_mode) {
	 print "setting mono to $mono because of a better correl\n";
	 }
	 }
	 }

	 */

	//check to make sure mono was set --


	if (isset($mono)) {
		return array($mono,$best);
	}
	else {
		return -1;
	}
}

// // Correl
// function correl($obs_arr, $exp_arr) {

// 	$avg_obs = array_sum($obs_arr)/count($obs_arr);
// 	$avg_exp = array_sum($exp_arr)/count($exp_arr);

// 	$numerator = $exp_denominator = $obs_denominator  = 0;
// 	foreach ($exp_arr as $idx => $exp_int) {
// 		if (isset($obs_arr[$idx])) {
// 			$obs_int = $obs_arr[$idx];
// 		}else{
// 			$obs_int = 0;
// 		}
// 		$numerator+= ($obs_int - $avg_obs)*($exp_int - $avg_exp);
// 		$obs_denominator += ($obs_int - $avg_obs)*($obs_int - $avg_obs);
// 		$exp_denominator += ($exp_int - $avg_exp)*($exp_int - $avg_exp);

// 	}
// 	if ($obs_denominator*$exp_denominator == 0) {
// 		return 0;
// 	}else{

// 		return($numerator/sqrt($obs_denominator*$exp_denominator));
// 	}
// }

/*
 Function which calculates the slope for obs and exp spectra
 Returns slope

 Assumes pre-normalized, pre-sized matching input arrays
 */

function calc_slope($obs_arr, $exp_arr) {

	/*
	 $count_obs = count($obs_arr);
	 if ($count_obs == 1) {
	 $obs_arrp[] = 0; // have at least two data points to calculate slope
	 $count_obs++;
	 }
	 //normalize first
	 */

	$max_exp = array_sum($exp_arr);
	$max_obs = array_sum($obs_arr);

	$sum_exp = $m_top = $m_bot = $r2_top = $r2_bot1 = $r2_bot2 = 0;
	foreach ($obs_arr as $idx => $int) {

		/*
		 if ($max_exp != 0 && $max_obs != 0) {
		 */
		if (isset($exp_arr[$idx])) {
			$m_top += ($int * $exp_arr[$idx] );
			$m_bot += ($int) * ($int);
			$sum_exp += ($exp_arr[$idx]/$max_exp);
		}
	}
	if ($m_bot == 0) {
		$m = 0;
	}else{
		$m = $m_top/$m_bot;
	}
	foreach ($obs_arr as $idx => $int) {
		if ($max_obs != 0 && $max_exp != 0 && isset($exp_arr[$idx])) {
			$r2_top += ($int * $exp_arr[$idx]);
			$r2_bot1 += ($int * $int);
			$r2_bot2 += ($exp_arr[$idx]) * ($exp_arr[$idx]);
		}
	}
	if ($r2_bot1 == 0 || $r2_bot2 == 0) {
		$r2 = 0;
	}else{
		$r2 = ($r2_top * $r2_top)/($r2_bot1 * $r2_bot2);
	}
	$score = (1-abs($m-1))*$r2; //slope component represents deviation from unity

	return array($score, $m, $r2);
}

/*
 Function which calculates the chi square value for obs and exp spectra
 Returns chi value

 Assumes pre-normalized arrays (JE)
 */

function calc_chi($obs_arr, $exp_arr) {

	/*
	 print "calculating chi value for:\n";
	 print "obs_arr:\n";
	 var_dump($obs_arr);
	 print "exp_arr:\n";
	 var_dump($exp_arr);
	 */

	//get obs sum
	$chi = -1;

	/*
	 print "observed sum: $sum_obs\n";
	 */

	foreach ($obs_arr as $idx => $val) {
		if (isset($exp_arr[$idx]) && $exp_arr[$idx]>0) {
			$chi += (($val - ($exp_arr[$idx])) * ($val - ($exp_arr[$idx])) )/ ($exp_arr[$idx]);
		}
	}
	if ($chi != -1) {
		$chi += 1; // correct starting off with -1
	}

	return $chi; //return -1 if error
}

/**
 *	calc_diff: calculates difference between two numbers given a specific ppm
 *
 *	given an expected value and ppm, calculates diff to 6 decimal places
 *  
 *      @author	Cocksman
 *      @param  FLOAT   $exp            float, expected value
 *      @param  INT     $ppm            int, ppm value
 *	    @return	FLOAT                   diff
 */
function calc_diff ($exp, $ppm) {
  return (sprintf("%.6f", ($exp*$ppm)/1000000));
}

// Deisotope
function deisotope($mzInt, $prec, $prec_charge){
	global $leftIsoWindow, $rightIsoWindow, $proton, $neutron, $ppmWindow, $multWindow, $noise_flag, $debug_mode, $major_debug;
	global $stopIntThres, $stopCountThres;

	$stopThres = $stopIntThres*median(array_values($mzInt));
	
	$stopCount = 0;
	$newMzInt = array();
	$indMzInt = array();
	while (true){
		#Find highest peak
		$maxHeight = -1;
		$maxMz = 0;
		foreach ($mzInt as $mz => $int){
			if ($int > $maxHeight){
				$maxMz = $mz;
				$maxHeight = $int;
			}
		}
		#If highest peak below threshold, break
		if ($maxHeight < $stopThres){
			if ($debug_mode){ print "Reached height threshold\n"; }
			break;
		}
		if ($stopCount > $stopCountThres){
			if ($debug_mode){ print "Reached count threshold\n"; }
			break;
		}
		
		#Get peaks within window
		$windowPeaks = array();
		if ($debug_mode){ print "Found:\n"; }
		foreach ($mzInt as $mz => $int){
			if (($mz < $maxMz + $rightIsoWindow) && ($mz > $maxMz - $leftIsoWindow)){
				$windowPeaks[$mz] = $int;
				if ($debug_mode){ print $mz . "\t" . $int . "\n"; }
			}
		}
		if ($debug_mode){ print "Done Found\n\n"; }
		
		$mono = array();
		$charge = array();
		#Find charge state
		list($mono, $charge) = get_charge_state($windowPeaks, $maxMz, null, $prec_charge, $major_debug);
		if ($debug_mode){
			print "make_mono found:\n";
			foreach($mono as $monos){
				print $monos . "\t";
			}print "\n";
			foreach($charge as $charges){
				print $charges . "\t";
			}print "\n";
		}
		$theCharge = max($charge);
		#(m+zh)/z to m+h 
		$oldMass = end($mono);
		$newMass = end($mono)*floatval($theCharge) - floatval($theCharge-1)*$proton;
		#Make sure new mass is less than precursor (cannot be larger)
		if ($mono[0] != 0 && $newMass <= $prec){
			if ($stopCount > 0){
				$stopCount--;
			}
			if ($debug_mode){ echo "Removing: " . $oldMass . "\t Charge: " . $theCharge . "\t" . $windowPeaks[$oldMass] . "\n"; }
			#Check if there's a peak already where we're moving
			$newPeak = $newMass;
			foreach ($newMzInt as $mz => $int){
				$window = calc_diff($newMass, $multWindow);
				if (($mz < $newMass + $window) && ($mz > $newMass - $window)){
					$newPeak = $mz;
				}
			}
			#convert to string, php doesn't like float keys
			$newKey = strval($newPeak);
			if ($debug_mode){ print "Int for " . $newPeak . " before was: " . $newMzInt[$newKey] . "\n"; }
			if (isset($indMzInt[$newKey])){
				$indMzInt[$newKey]++;
			}else{
				$indMzInt[$newKey] = 1;
			}
			for ($i=0; $i<=4; $i++){
				#Find all peaks involved
				$testMz = $oldMass+$i*$neutron/$theCharge;
				$window = calc_diff($testMz, $ppmWindow);
				if ($debug_mode){ print "Looking for: " . $testMz . " with window " . $window . "\n"; }
				foreach ($windowPeaks as $mz => $int){
					#If it's there, add it too
					if (($mz < $testMz + $window) && ($mz > $testMz - $window)){
						if ($debug_mode){ print "Removing: " . $mz . "\t" . $int . "\n"; }
						if (isset($newMzInt[$newKey])){
							$newMzInt[$newKey] += $int;
						}else{
							$newMzInt[$newKey] = $int;
						}
						unset($mzInt[$mz]);
						if ($debug_mode){ print "Int is now: " . $newMzInt[$newKey] . "\n"; }
					}
				}
			}
			if ($debug_mode){ 
				print "Moved to: mz=" . $newPeak . "\tint=" . $newMzInt[$newKey] . "\n";
				print "Done Removed\n\n";
			}
		}else{
			if ($debug_mode){ print "make_mono failed (dropping: " . $maxMz . ")\n\n"; }
			$stopCount++;
			#If not found, get rid of it (and possibly add to noise)
			if ($noise_flag){
				if (isset($indMzInt[$maxMz])){
					$indMzInt[$maxMz]++;
					$newMzInt[$maxMz] += $maxHeight;
				}else{
					$indMzInt[$maxMz] = 1;
					$newMzInt[$maxMz] = $maxHeight;
				}
			}
			unset($mzInt[$maxMz]);
		}
	}

	if ($noise_flag){
		#Transfer all other peaks
		foreach ($mzInt as $mz => $int){
			if (isset($indMzInt[$mz])){
				$indMzInt[$mz]++;
				$newMzInt[$mz] += $int;
			}else{
				$indMzInt[$mz] = 1;
				$newMzInt[$mz] = $int;
			}
		}
	}

	foreach ($newMzInt as $mz => $int){
		$outMzInt[$mz] = $newMzInt[$mz]/$indMzInt[$mz];
	}
	return $outMzInt;

	// Xdebug call - remove for publication!
	xdebug_pring_function_stack('End of deisotope function');
}

function extractMS($precursor, $ms2_peaks, $tic, $charge) {
	$debug_mode = FALSE;	
	if ($debug_mode) {
		print "Running with these vars parent scan: precursor: $precursor, tic: $tic, z: $charge\n";
	}
	$mono_arr = array();
	$charge_arr = array();
	//decode ms2 peaks
	//    if (!$debug_mode) {
	$mz_int_array = mzxml_decode_peaks($ms2_peaks);
	// }

	//check for one plus by tic
	//if (!$debug_mode) {
	if ( check_one_plus($mz_int_array, $tic, $precursor)  && $charge == NULL) {
		//one plus
		$charge_arr[] = 1;
		$parent_arr[] = $precursor;
		return array($parent_arr, $charge_arr);
	}

	else {
		//not one plus
		if ($charge == NULL) {
			//print "charge was null\n";
			$mono_arr[] = $precursor;
			$mono_arr[] = $precursor;
			$charge_arr[] = 2;
			$charge_arr[] = 3;			
			return array($mono_arr, $charge_arr);
		}
		else {
			$mono_arr[] = $precursor;
			$charge_arr[] = $charge;			
			return array($mono_arr, $charge_arr);
		}		 
	}
}

function check_one_plus($mz_int_array, $tic, $precursor) {
	//tic below precursor
	$low_tic = 0;
	if ($tic <=0) {
		return TRUE;
	}

	ksort($mz_int_array);
	foreach ($mz_int_array as $mz => $int) {
		if ($mz <= $precursor) {
			$low_tic += $int;
		}
		else {
			break;
		}
	}

	$ratio = $low_tic/$tic;
	if ($ratio >= 0.95) {
		return TRUE;
	}
	else {
		return FALSE;
	}
}

//xdebug_stop_trace();
?>
