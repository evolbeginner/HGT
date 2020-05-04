#! /usr/bin/env ruby


##############################################################################
DIR = File.dirname($0)
$: << File.join(DIR, 'lib')


##############################################################################
require 'getoptlong'

require 'process_species_name'


##############################################################################
infile = nil


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
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
  end
end


##############################################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  orgn0 = line_arr[0]
  orgn = Marshal.load(Marshal.dump(orgn0)).process_species_name!
  puts [line_arr, orgn].flatten.join("\t")
end
in_fh.close


