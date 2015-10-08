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
	'Wtaunu_3mu'    => "#{script} #{mode} #{code} #{type} #{outdir} Wtaunu_3mu    #{source} #{version} #{split}  >  #{outdir}Wtaunu_3mu.log 2>&1",
	'bbTotau10_3mu' => "#{script} #{mode} #{code} #{type} #{outdir} bbTotau10_3mu #{source} #{version} #{split}  >  #{outdir}bbTotau10_3mu.log 2>&1",
	'ccTotau10_3mu' => "#{script} #{mode} #{code} #{type} #{outdir} ccTotau10_3mu #{source} #{version} #{split}  >  #{outdir}ccTotau10_3mu.log 2>&1",

	'bbTomu15'      => "#{script} #{mode} #{code} #{type} #{outdir} bbTomu15      #{source} #{version} #{split}  >  #{outdir}bbTomu15.log 2>&1",
	'ccTomu15'      => "#{script} #{mode} #{code} #{type} #{outdir} ccTomu15      #{source} #{version} #{split}  >  #{outdir}ccTomu15.log 2>&1",
	'bb_mu4mu4'     => "#{script} #{mode} #{code} #{type} #{outdir} bb_mu4mu4     #{source} #{version} #{split}  >  #{outdir}bb_mu4mu4.log 2>&1",
	'bb_Jpsimu4mu4' => "#{script} #{mode} #{code} #{type} #{outdir} bb_Jpsimu4mu4 #{source} #{version} #{split}  >  #{outdir}bb_Jpsimu4mu4.log 2>&1",
	'JZ2W'          => "#{script} #{mode} #{code} #{type} #{outdir} JZ2W          #{source} #{version} #{split}  >  #{outdir}JZ2W.log 2>&1",
	#'JZ3W'          => "#{script} #{mode} #{code} #{type} #{outdir} JZ3W          #{source} #{version} #{split}  >  #{outdir}JZ3W.log 2>&1",

	#'periodA'       => "#{script} #{mode} #{code} #{type} #{outdir} periodA       #{source} #{version} #{split}  >  #{outdir}periodA.log 2>&1",
	'periodB'       => "#{script} #{mode} #{code} #{type} #{outdir} periodB       #{source} #{version} #{split}  >  #{outdir}periodB.log 2>&1",
	'periodC'       => "#{script} #{mode} #{code} #{type} #{outdir} periodC       #{source} #{version} #{split}  >  #{outdir}periodC.log 2>&1",
	'periodD'       => "#{script} #{mode} #{code} #{type} #{outdir} periodD       #{source} #{version} #{split}  >  #{outdir}periodD.log 2>&1",
	'periodE'       => "#{script} #{mode} #{code} #{type} #{outdir} periodE       #{source} #{version} #{split}  >  #{outdir}periodE.log 2>&1",
	'periodG'       => "#{script} #{mode} #{code} #{type} #{outdir} periodG       #{source} #{version} #{split}  >  #{outdir}periodG.log 2>&1",
	'periodH'       => "#{script} #{mode} #{code} #{type} #{outdir} periodH       #{source} #{version} #{split}  >  #{outdir}periodH.log 2>&1",
	#'periodI'       => "#{script} #{mode} #{code} #{type} #{outdir} periodI       #{source} #{version} #{split}  >  #{outdir}periodI.log 2>&1",
	#'periodJ'       => "#{script} #{mode} #{code} #{type} #{outdir} periodJ       #{source} #{version} #{split}  >  #{outdir}periodJ.log 2>&1",
	#'periodL'       => "#{script} #{mode} #{code} #{type} #{outdir} periodL       #{source} #{version} #{split}  >  #{outdir}periodL.log 2>&1",
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
%x("cd #{outdir}")
system(hadd)
%x("cd -")
