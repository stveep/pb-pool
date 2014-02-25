#!/usr/bin/ruby
# Pipeline for analysis of barcodes from transposon libraries.
# Steve Pettitt, ICR
# spettitt@icr.ac.uk

# Location of configuration file to use:
# Config file format is (tab-separated):
#<fastq file name, full path>	<sequencing run ID>	<seq. index>	<library>	<drug>
config = ARGV.last

# Default output file name - better to pass this direct to summary task from cmd line?
outf = "summary.tsv"

# The fastq.rb file needs to be available.
require '/home/breakthr/spettitt/data/scripts/pool/fastq.rb'

# Barcode match/capture regexps can be changed if required:
sbc = /TGAATTC(.{20,30})(TACATC)/
wbc = /CAGACTG(.{20,30})(GGTCGA)/

# Separators for output file and column names respectively:
sep = "\t"
fieldsep = "_"


# Initialise the hash that will contain data about each file
fileh = Hash.new{|h,k| h[k] = Hash.new}
# Sort out the config file:
IO.foreach(config) do |line|
	fields = line.chomp.split("\t")
	fileh[fields[0]] = {'run' => fields[1], 'index' => fields[2], 'library' => fields[3], 'drug' => fields[4]}
end
# Fileh is now a hash keyed by file name, pointing to another hash with all the sample information.
# The rest of the script will use the list of keys from the files hash to decide which files to process.
files = fileh.keys
puts files

desc "Count barcodes in each file"
task :barcodes => files.each.map{|a| a.sub(/\.fastq$/,".sims.bc")} + files.each.map{|a| a.sub(/\.fastq$/,".west.bc")}

desc "Associate barcodes with inverse PCR mappings"
task :map => files.each.map{|a| a.sub(/\.fastq$/,".sims.map")} + files.each.map{|a| a.sub(/\.fastq$/,".west.map")}

desc "Do barcode processing in parallel on HPC (recommended for lots of large files)"
task :parallel => files do |t|
	files.each do |f|
	puts f
	#Make a temporary config file for just this one file
	cfname = f + ".configtmp"
	cf = File.open(cfname,'w')
	cf.puts f + "\t" + fileh[f]["run"] + "\t" + fileh[f]["index"] + "\t" + fileh[f]["library"] + "\t" + fileh[f]["drug"]
	sh "bsub -P BCAVAS -J parallelbc rake -f ~/data/scripts/pool/barcode-pipe.rake barcodes #{cfname}"
	end
end

# Sims and West barcodes from the same sample will have the same index, but can be identified from the surrounding sequence.
# This rule identifies barcode reads from the fastq file and outputs them, with counts, to the appropriate barcode (.bc) file:
# The rule says "a foo.*.bc file [* = sims or west in this case] can be made from a foo.fastq file by doing the following...." 
rule(/\.\w+\.bc$/ => proc{ |task_name| task_name.sub(/\.\w+\.bc$/, ".fastq") } ) do |a|
	f = Fastq.new(a.prerequisites[0])
	
	# Set up hashes for counting barcodes. This code sets the default value of any new key to be zero.
	simsh = Hash.new{|h,k| h[k] = 0}
	westh = Hash.new{|h,k| h[k] = 0}
	
	# Go through each read (referred to as r) in the fastq file
	f.each{|r|
		# Is the read a sims barcode?
		if r.seq =~ sbc
			simsh[$1] += 1 #Increment hash entry for that barcode, which has been captured into $1 by the regexp match	
		# ... or if not, maybe a "west" barcode?
		elsif r.seq =~ wbc
			westh[$1] += 1
		end #if-else
	} # fastq each
	# After looping through all the reads in the fastq file, we have two hashes containing all the counts for each barcode. Now to write these to new files:
	osb = File.open(a.prerequisites[0].sub(".fastq",".sims.bc"),'w')
	simsh.each{|k,v| osb.puts k+"\t"+v.to_s}
	osb.close

	owb = File.open(a.prerequisites[0].sub(".fastq",".west.bc"),'w')
	westh.each{|k,v| owb.puts k+"\t"+v.to_s}
	owb.close
end # bc rule


rule(".map" => '.bc' ) do |a|
	sh "cp #{a.name.sub('.map','.bc')} #{a.name}" 
#	m = File.open(a.name,'w')
#	$1 here is sims or west:
#	m.puts "chr\tstrand\tpos\tgene\tbc\t" + fileh[a.name.sub(/\.(\w+)\.map$/,'.fastq')]['library'] + "_" + $1 + "_" + fileh[a.name.sub(/\.(\w+)\.map$/,'.fastq')]['drug']
	# This is to be modified to lookup barcodes in the inverse PCR data.
#	IO.foreach(a.prerequisites[0]){|l| m.puts "chr\tstrand\tpos\tgene\t" + l.chomp } 
end # map rule

task :summary => files.map{|a| a.sub('.fastq','.sims.map')} + files.map{|a| a.sub('.fastq','.west.map')} do |t|
# Counter for file number
i = 1
hash = {}
titles = []
t.prerequisites.each do |f|
	input = File.open(f,'r')
	input.each do |l|
# Split into fields
		(barcode, count) = l.split(" ") 

# Use hash to associate coverage values array (PB5+3 ends) with site
		if (hash[barcode])
		hash[barcode] << count
		else
			hash[barcode] = []
# If it's a new site, we can add if first file...
			if (i == 1)
				hash[barcode] << count
# ...if not, need to add an appropriate no. of zeroes first to keep table format
			else
				j = i-1
				j.times {hash[barcode] << 0}
				hash[barcode] << count
			end
		end
	end
# Add zero array to any barcode that wasn't seen in this file
	hash.each {|k,v| hash[k] << 0 if (v.length < i)}
	# Now to make the corresponding title for the column in the final table
	# Use the config file data to construct the column names:
	# Column name will be Lib_drug_sims/west_run_index
	# Need to reconstruct the fastq filename for the hash:
	entry = fileh[f.sub(/\.(\w+).map/,'.fastq')]
	# reminder of format: fileh[fields[0]] = {'run' => fields[1], 'index' => fields[2], 'library' => fields[3], 'drug' => fields[4]}
	puts f.sub(/\.(\w+).map/,'.fastq') 
	titles << entry["library"] + fieldsep + entry["drug"] + fieldsep + $1 + fieldsep + entry["run"] + fieldsep + entry["index"]
# That's enough looping through files
	i += 1
end

# Output a table from the hash

o = File.open(outf, 'w')

# [Need to add column for gene/mapping too]
o.puts "barcode" + sep + titles.join(sep)

# Order in which barcodes are accessed isn't predictable, but the the counts are in an array and therefore in correct order:
hash.each{|k,v| o.puts k + sep + v.join(sep)}

end # summary task
