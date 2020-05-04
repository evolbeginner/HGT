#! /usr/bin/env ruby


##############################################################################
DIR = File.dirname($0)
$: << File.join(DIR, 'lib')


##############################################################################
require 'getoptlong'
require 'parallel'
require 'bio'

require 'Dir'
require 'process_species_name'


##############################################################################
infile = nil
cpu = 1
outdir = nil
is_force = false
is_tolerate = false


##############################################################################
def checkFas(infile, orgn0, level)
  seq_objs = Hash.new
  in_fh = Bio::FlatFile.open(infile)
  in_fh.each_entry do |f|
    f.definition =~ /\[ ([^\[\]]+) \]$/x
    orgn2 = $1
    if level == 'species'
      next if orgn0 !~ /#{orgn2}/ and orgn2 !~ /#{orgn0}/
    elsif level == 'genus'
      next if orgn0.nil? or orgn2.nil?
      genus0 = orgn0.split(' ')[0]
      genus2 = orgn2.split(' ')[0]
      next if genus0 != genus2
    elsif level == 'any'
      ;
    else
      break
    end
    seq_objs[f.definition] = f
  end
  in_fh.close
  return(seq_objs)
end


def outputSeqObjs(outfile, seq_objs)
  out_fh = File.open(outfile, 'w')
  seq_objs.each_pair do |seq_title, f|
    out_fh.puts '>' + seq_title
    out_fh.puts f.seq
  end
  out_fh.close
end


##############################################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--tolerate', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--cpu'
      cpu = value.to_i
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
    when '--tolerate'
      is_tolerate = true
  end
end


##############################################################################
mkdir_with_force(outdir, is_force, is_tolerate)


##############################################################################
in_fh = File.open(infile, 'r')
lines = in_fh.readlines()
in_fh.close


##############################################################################
Parallel.map(lines, in_threads: cpu) do |line|
  line.chomp!
  line_arr = line.split("\t")
  orgn0 = line_arr[0]
  orgn = Marshal.load(Marshal.dump(orgn0))
  orgn.process_species_name!

  outfile = File.join(outdir, orgn+'.fas')
  #if not File.exists?(outfile)
  #  puts orgn0
  #end

  `esearch -db protein -query "#{orgn0}[orgn] AND refseq[filter]" 2>/dev/null | efetch -format fasta 2>/dev/null >#{outfile}`
  seq_objs = checkFas(outfile, orgn0, 'any')

  if seq_objs.size < 100
    `esearch -db protein -query "#{orgn0}[orgn] AND refseq[filter]" 2>/dev/null | efetch -format fasta 2>/dev/null >#{outfile}`
    seq_objs = checkFas(outfile, orgn0, 'genus')
  end

  puts [orgn0, orgn, seq_objs.size].join("\t")

  if seq_objs.size < 100
    `rm #{outfile}`
  else
    outputSeqObjs(outfile, seq_objs)
  end
end


