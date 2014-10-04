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
    f.puts "source #{$homedir}/setROOT.sh"
    f.puts "cd #{$thisdir}"
    f.puts "./execute.sh  loadrun  vtxing  skim  #{$outdir}/  #{inset}  #{master}  #{method}  #{split}"
    f.puts "echo \"host = $HOSTNAME\""
    f.puts "cd -"
  }
end

masters = Array.new
# masters = ["muons", "muid"]
masters = ["muons"]

insets = Hash.new
insets = {
  'periodA' =>     8, # nfiles=29
  'periodB' =>    22, # nfiles=87
  'periodC' =>     8, # nfiles=30
  'periodD' =>    14, # nfiles=55
  'periodE' =>     9, # nfiles=33
  'periodG' =>     5, # nfiles=20
  'periodH' =>     6, # nfiles=22
  'periodI' =>     4, # nfiles=16
  'periodJ' =>     9, # nfiles=33 (is 34 but one file is corrupted)
  'periodL' =>     3, # nfiles=11
  
  'bb_mu4mu4'     => 6, # nfiles=21
  'bb_Jpsimu4mu4' => 4, # nfiles=14
  'bbTomu15'      => 1, # nfiles=2
  'ccTomu15'      => 1, # nfiles=1
  
  'JZ2W' => 1,  # nfiles=4
  'JZ3W' => 1,  # nfiles=3
  
  'bbTotau10_3mu' => 1, # nfiles=2
  'ccTotau10_3mu' => 1, # nfiles=2
  'Wtaunu_3mu'    => 1  # nfiles=2
  
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
