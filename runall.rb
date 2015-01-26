#!/usr/bin/env ruby
#require 'pathname'
#require 'fileutils'
#require 'ftools'
#require 'find'
#require 'securerandom'
## require 'LoggerDecorator'
#include FileUtils

script="./execute.sh"
mode="loadrun"
code="NTUPmaker"
type="analysis"
outdir="$HOME/data/bout/"
source="muons"
version="CUTS:ALL"
split="0:0"

build="#{script} build #{code} #{type} #{outdir} Wtaunu_3mu #{source} #{version} #{split}"

procs = Hash.new
procs = {
	'Wtaunu_3mu'    => "#{script} #{mode} #{code} #{type} #{outdir} Wtaunu_3mu    #{source} #{version} #{split}",
	'bbTotau10_3mu' => "#{script} #{mode} #{code} #{type} #{outdir} bbTotau10_3mu #{source} #{version} #{split}",
	'ccTotau10_3mu' => "#{script} #{mode} #{code} #{type} #{outdir} ccTotau10_3mu #{source} #{version} #{split}",

	'bbTomu15'      => "#{script} #{mode} #{code} #{type} #{outdir} bbTomu15      #{source} #{version} #{split}",
	'ccTomu15'      => "#{script} #{mode} #{code} #{type} #{outdir} ccTomu15      #{source} #{version} #{split}",
	'bb_mu4mu4'     => "#{script} #{mode} #{code} #{type} #{outdir} bb_mu4mu4     #{source} #{version} #{split}",
	'bb_Jpsimu4mu4' => "#{script} #{mode} #{code} #{type} #{outdir} bb_Jpsimu4mu4 #{source} #{version} #{split}",
	'JZ2W'          => "#{script} #{mode} #{code} #{type} #{outdir} JZ2W          #{source} #{version} #{split}",
	#'JZ3W'          => "#{script} #{mode} #{code} #{type} #{outdir} JZ3W          #{source} #{version} #{split}",

	#'periodA'       => "#{script} #{mode} #{code} #{type} #{outdir} periodA       #{source} #{version} #{split}",
	'periodB'       => "#{script} #{mode} #{code} #{type} #{outdir} periodB       #{source} #{version} #{split}",
	'periodC'       => "#{script} #{mode} #{code} #{type} #{outdir} periodC       #{source} #{version} #{split}",
	'periodD'       => "#{script} #{mode} #{code} #{type} #{outdir} periodD       #{source} #{version} #{split}",
	'periodE'       => "#{script} #{mode} #{code} #{type} #{outdir} periodE       #{source} #{version} #{split}",
	'periodG'       => "#{script} #{mode} #{code} #{type} #{outdir} periodG       #{source} #{version} #{split}",
	'periodH'       => "#{script} #{mode} #{code} #{type} #{outdir} periodH       #{source} #{version} #{split}",
	#'periodI'       => "#{script} #{mode} #{code} #{type} #{outdir} periodI       #{source} #{version} #{split}",
	#'periodJ'       => "#{script} #{mode} #{code} #{type} #{outdir} periodJ       #{source} #{version} #{split}",
	#'periodL'       => "#{script} #{mode} #{code} #{type} #{outdir} periodL       #{source} #{version} #{split}",
}

def runproc(command)
	system(command)
	sleep(rand(10))
	puts "Done... #{command}"
end

# end of definition - start of run
#system(build) # compile first

# Set the threads going
tasklist = []
procs.each{|key,val|
	procname = key
	command  = val
	task = Thread.new(command) { |comm| runproc(comm) }
	tasklist << task
}

# Wait for the threads to finish
tasklist.each { |task|
	task.join
}

# join all the outputs
sleep(rand(20))
hadd="hadd -f /tmp/hod/flatout.periodall+MC.muons.cuts.analysis.n0.j0.root "
procs.each{|key,val|
	procname = key
	hadd = hadd + " flatout.#{procname}.muons.cuts.analysis.n0.j0.root"
}
puts "\n\n#{hadd}"
#system(hadd)
