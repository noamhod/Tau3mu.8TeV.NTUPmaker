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
  'periodA' =>     5, # nfiles=41
  'periodB' =>    14, # nfiles=132 (133 with one corrupted, available in $EOSHOME/data/backup/user.hod.4143407.d3pdstream._000078.root)
  'periodC' =>     6, # nfiles=52
  'periodD' =>    11, # nfiles=104
  #'periodE' =>     , # nfiles=33
  'periodG' =>     4, # nfiles=31
  'periodH' =>     4, # nfiles=39
  'periodI' =>     3, # nfiles=26
  #'periodJ' =>     , # nfiles=
  'periodL' =>     2, # nfiles=18
  
  'bb_mu4mu4'     => 5, # nfiles=46
  'pp_Jpsimu4mu4' => 2, # nfiles=15
  'bb_Jpsimu4mu4' => 3, # nfiles=21
  'bbTomu15'      => 1, # nfiles=5
  'ccTomu15'      => 1, # nfiles=4
  
  'JZ2W' => 1,  # nfiles=5
  'JZ3W' => 1,  # nfiles=4
  
  'bbTotau10_3mu' => 1, # nfiles=2
  'ccTotau10_3mu' => 1, # nfiles=2
  'Wtaunu_3mu'    => 1  # nfiles=3
  
  # 'ZmumuNp0'  => 1,
  # 'ZmumuNp1'  => 1,
  # 'ZmumuNp2'  => 1,
  # 'ZmumuNp3'  => 1,
  # 'ZmumuNp4'  => 1,
  # 'ZmumuNp5'  => 1,
  # 'WmunuNp0'  => 1,
  # 'WmunuNp1'  => 1,
  # 'WmunuNp2'  => 1,
  # 'WmunuNp3'  => 1,
  # 'WmunuNp4'  => 1,
  # 'WmunuNp5'  => 1,
  # 'WtaunuNp0' => 1,
  # 'WtaunuNp1' => 1,
  # 'WtaunuNp2' => 1,
  # 'WtaunuNp3' => 1,
  # 'WtaunuNp4' => 1,
  # 'WtaunuNp5' => 1
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
      method = "CUTS:vertex,MET" # MVA is irrelevant here
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
