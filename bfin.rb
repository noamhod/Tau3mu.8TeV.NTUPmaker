#!/usr/bin/env ruby
require 'pathname'
require 'fileutils'
require 'ftools'
require 'find'
# require 'LoggerDecorator'
include FileUtils

$thisdir = Dir.pwd
$homedir = "/afs/cern.ch/user/h/hod"
$outdir  = "#{$homedir}/data/bout"
$logdir  = "#{$homedir}/data/logs"
$jobdir  = "#{$homedir}/data/jobs"
$figdir  = "#{$homedir}/data/figs"

masters = Array.new
# masters = ['muons', 'muid']
masters = ['muons']


insets = Hash.new
insets = {
  # 'periodA' =>       22, # nfiles=44
  # 'periodB' =>      17, # nfiles=167
  # 'periodC' =>       9, # nfiles=83
  # 'periodD' =>      12, # nfiles=112
  # 'periodE' =>       8, # nfiles=80
  # 'periodG' =>       5, # nfiles=42
  # 'periodH' =>       6, # nfiles=52
  # 'periodI' =>      50, # nfiles=99
  # 'periodJ' =>      86, # nfiles=257
  # 'periodL' =>        46 # nfiles=91
  # 
  # 'bb_mu4mu4'     => 3, # nfiles=27
  # 'bb_Jpsimu4mu4' => 2, # nfiles=19
  # 'bbTomu15'      => 1, # nfiles=6
  # 'ccTomu15'      => 1, # nfiles=6
  # 
  # 'JZ2W' => 1,  # nfiles=4
  # 'JZ3W' => 3,  # nfiles=11
 
    #'WtaunuNp0' => 3, #nfiles=23
    #'WtaunuNp1' => 2, #nfiles=16
    'WtaunuNp2' => 3, #nfiles=25
    #'WtaunuNp3' => 1, #nfiles=7
    #'WtaunuNp4' => 1, #nfiles=4
    #'WtaunuNp5' => 1, #nfiles=3

    #'WmunuNp0' => 3, #nfiles=21
    #'WmunuNp1' => 2, #nfiles=17
    #'WmunuNp2' => 3, #nfiles=23
    #'WmunuNp3' => 1, #nfiles=5
    #'WmunuNp4' => 1, #nfiles=6
    'WmunuNp5' => 1, #nfiles=2

    #'ZtautauNp0' => 4, #nfiles=40
    ##'ZtautauNp1' => , #nfiles=
    #'ZtautauNp2' => 1, #nfiles=4
    #'ZtautauNp3' => 1, #nfiles=2
    #'ZtautauNp4' => 1, #nfiles=2
    #'ZtautauNp5' => 1, #nfiles=1

    #'ZmumuNp0' => 4, #nfiles=39
    #'ZmumuNp1' => 1, #nfiles=10
    #'ZmumuNp2' => 1, #nfiles=5
    ##'ZmumuNp3' => 1, #nfiles=2
    #'ZmumuNp4' => 1, #nfiles=2
    #'ZmumuNp5' => 1, #nfiles=1
  
  # 'bbTotau10_3mu' => 1, # nfiles=2
  # 'ccTotau10_3mu' => 1, # nfiles=2
  # 'Wtaunu_3mu'    => 1  # nfiles=2
}

fall=0
fmiss=0
fexist=0
fwarn=0
fgood=0

logfiles = Array.new
insets.each {|key,val|
  inset = key
  njobs = val
  masters.each do |master|
    for j in 1..njobs
      job = (j<10) ? "0#{j}" : "#{j}"
      logname = "#{$logdir}/#{inset}.#{master}.n#{njobs}.j#{job}.log";
      if(!File.file?(logname)) then
        puts "File does not (yet) exist: #{logname}"
        fmiss += 1
        next
      end
      logfiles << logname
    end
    fall += njobs
  end
}
logfiles.sort
fexist = logfiles.length


lastline=-1
firstline=-30
phrase1="EVENTS PROCESSED"
phrase2="ALL DONE CORRECTLY"
phrase3="otrees[MUONS_TRIPLET] is NULL"
phrase4="otrees[MUID_TRIPLET] is NULL"

logfiles.each do |file|
  
  nlines = %x{wc -l #{file}}.split[0].to_i
  if(nlines<firstline.abs) then abort end
  
  lines = IO.readlines(file)[firstline..lastline]
  ok1 = false
  ok2 = false
  ok3 = false
  ok4 = false
  lines.each do |line|
    ok1 = (!ok1) ? (line.include? phrase1) : true
    ok2 = (!ok2) ? (line.include? phrase2) : true
    ok3 = (!ok3) ? (line.include? phrase3) : true
    ok4 = (!ok4) ? (line.include? phrase4) : true
  end
  ok = ((ok1 and ok2) or (ok3 or ok4))
  if(!ok) then
    puts "Bad file: #{file}, missing patterns:"
    puts "`#{phrase1}`  and  `#{phrase2}`"
    puts "or  `#{phrase3}`"
  end
  if(ok3 or ok4) then
    puts "Note! unprocessed sample: #{file}"
  end
  fwarn = (ok3) ? fwarn+1 : fwarn 
  fwarn = (ok4) ? fwarn+1 : fwarn 
  fgood = (ok)  ? fgood+1 : fgood
  # fgood = (ok1 and ok2) ? fgood+1 : fgood
end

puts "<<<<<<<<< validation summary >>>>>>>>>"
puts "All files:      #{fall}"
puts "Missing files:  #{fmiss}"
puts "Existing files: #{fexist}"
puts "Warning files:  #{fwarn}"
puts "Good files:     #{fgood}"
puts "<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>"


if(fgood==fall) then
  Dir.chdir $outdir
  %x(mkdir data)
  %x(mkdir data/bphys)
  %x(mkdir data/bphys/data)
  %x(mkdir data/bphys/mc)
  masters.each do |master|
    %x(mkdir data/bphys/data/#{master})
    %x(mkdir data/bphys/mc/#{master})
    datapattern = ""
    mcspattern  = ""
    insets.each {|key,val|
      inset = key
      pattern = "#{$figdir}/fig.#{inset}*#{master}*.root"
      if(inset.include? "period") then
        datapattern = datapattern + " #{pattern}"
      else
        mcspattern  = mcspattern  + " #{pattern}"
      end
    }
    # # datapattern = "#{$figdir}/fig.*#{master}*.root"
    # # mcspattern  = "#{$figdir}/fig.bb*#{master}*.root #{$figdir}/fig.cc*#{master}*.root"
    # mcspattern = mcspattern + " #{$figdir}/fig.J*#{master}*.root"
    # mcspattern = mcspattern + " #{$figdir}/fig.W*#{master}*.root"
    
    %x(hadd -f #{$figdir}/fig.data.#{master}.root #{datapattern})
    %x(hadd -f #{$figdir}/fig.mcs.#{master}.root  #{mcspattern})
    %x(hadd -f #{$figdir}/fig.all.#{master}.root  #{datapattern} #{mcspattern})
    insets.each {|key,val|
      inset = key
      target = ""
      if(inset.include? "period") then
        target = "data/bphys/data/#{master}/#{inset}"
      else
        target = "data/bphys/mc/#{master}/#{inset}"
      end
      %x(mkdir #{target})
      %x(mv out.#{inset}*#{master}* #{target}/)
    }
  end
  Dir.chdir $thisdir
else
  puts "Cannot merge since not all jobs are done"
  puts "Call `bjobs` to monitor the running jobs"
end
