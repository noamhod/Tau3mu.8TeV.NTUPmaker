#!/usr/bin/env ruby
require 'pathname'
require 'fileutils'
require 'ftools'
require 'find'
require 'securerandom'
# require 'LoggerDecorator'
include FileUtils

$thisdir = Dir.pwd
$homedir = "/afs/cern.ch/user/h/hod"
$outdir  = "#{$homedir}/data/bout"
$logdir  = "#{$homedir}/data/logs"
$jobdir  = "#{$homedir}/data/jobs"
$figdir  = "#{$homedir}/data/figs"

def areyousure()
  puts "Are you sure you want to run bsub ??? [yes/no]";
  option = gets();
  option = option.gsub(" ","");
  option = option.gsub("\n","");
  if(option=="yes") then
    puts "OK !"
    puts "to list all, do: `bjobs -w`"
    puts "to kill all, do: `bkill -u hod 0`"
  else
    abort("Cancelling");
  end
end

def iscompiled()
  puts "did you remember to compile before submission ??? [yes/no]";
  option = gets();
  option = option.gsub(" ","");
  option = option.gsub("\n","");
  if(option=="yes") then
    puts "OK !";
  else
    puts "please run ./execute.sh and follow the instructions";
    abort("Cancelling");
  end
end

def cleandir(path)
  puts "erase previous output from: #{path} ??? [yes/no]";
  option = gets();
  option = option.gsub(" ","");
  option = option.gsub("\n","");
  if(option=="yes") then
    puts "your choice is to clean: #{option}"
    r = SecureRandom.random_number(100);
    puts "to erase previous output from: #{path},";
    puts "please enter this number: #{r}";
    randnum = gets();
    randnum = randnum.gsub(" ", "");
    randnum = randnum.gsub("\n", "");
    if(randnum.to_i==r) then
      if(File.directory?(path)) then
        %x(rm -rf #{path}/*)
        puts "OK !";
      else
        puts "dir does not exist: #{path}";
        abort("Cancelling");
      end
    else
      puts "Wrong number: #{randnum} -> not erasing !";
    end
  elsif(option=="no") then
    puts "your choice is NOT to clean: #{option}"
  else
    puts "unknown option: #{option}"
    abort("Cancelling");
  end
end

def makejobfile(jobname,inset,master,method,split)
  jobfile = File.open(jobname,'w') { |f|
    f.puts "#!/bin/bash"
    f.puts "echo \"host = $HOSTNAME\""
    f.puts "cd #{$thisdir}"
    f.puts "source setup.sh"
    f.puts "./execute.sh  loadrun  NTUPmaker  skim  #{$outdir}/  #{inset}  #{master}  #{method}  #{split}"
    f.puts "echo \"host = $HOSTNAME\""
    f.puts "cd -"
  }
end

masters = Array.new
# masters = ["muons", "muid"]
masters = ["muons"]

insets = Hash.new
insets = {
  # 'periodA' =>      22, # nfiles=44
  # 'periodB' =>      17, # nfiles=167
  # 'periodC' =>       9, # nfiles=83
  # 'periodD' =>      12, # nfiles=112
  # 'periodE' =>       8, # nfiles=80
  # 'periodG' =>       5, # nfiles=42
  # 'periodH' =>       6, # nfiles=52
  # 'periodI' =>      50, # nfiles=99
  # 'periodJ' =>     86, # nfiles=257
  # 'periodL' =>       46 # nfiles=91
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
    #'ZmumuNp3' => 1, #nfiles=2
    #'ZmumuNp4' => 1, #nfiles=2
    #'ZmumuNp5' => 1, #nfiles=1
 
  # 'bbTotau10_3mu' => 1, # nfiles=2
  # 'ccTotau10_3mu' => 1, # nfiles=2
  # 'Wtaunu_3mu'    => 1  # nfiles=2
}


#################################################
#################################################
#################################################

areyousure();
iscompiled();
cleandir($logdir);
cleandir($jobdir);
cleandir($outdir);
cleandir($figdir);

Dir.chdir "#{$jobdir}"
insets.each{|key,val|
  inset = key
  njobs = val
  for j in 1..njobs
    job = (j<10) ? "0#{j}" : "#{j}"
    masters.each do |master|
      jobname = "#{$jobdir}/#{inset}.#{master}.n#{njobs}.j#{job}.sh";
      logname = "#{$logdir}/#{inset}.#{master}.n#{njobs}.j#{job}.log";
      split = "#{njobs}:#{j}"
      method = "CUTS:ALL" # MVA is irrelevant here
      makejobfile(jobname,inset,master,method,split)
      %x(chmod 755 #{jobname})
      # command = "bsub -q1nh -W 240 -o #{logname} #{jobname}"
      # command = "bsub -q 8nm -o #{logname} #{jobname}"
      command = "bsub -q 1nd -o #{logname} #{jobname}"
      puts "running: #{command}"
      %x(#{command})
    end
  end
}
puts "run `bjobs -w` to list all running jobs"
Dir.chdir $thisdir
