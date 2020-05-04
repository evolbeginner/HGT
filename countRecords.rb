#! /usr/bin/env ruby


################################################
require 'getoptlong'
require 'parallel'


################################################
infile = nil
cpu = 1

orgn_info = Hash.new


################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--cpu'
      cpu = value.to_i
  end
end


################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  orgn = line_arr[0]
  orgn_info[orgn] = line_arr[1, line_arr.size-1]
end
in_fh.close


################################################
results = Parallel.map(orgn_info, in_processes: cpu) do |orgn, line_arr|
  #orgn_query = ['"', '"'].join(orgn)
  a = `esearch -db protein -query "#{orgn}[orgn] AND refseq[filter]" 2>/dev/null  | grep Count`.chomp
  a =~ /^\s+ \<Count\> (\d+) \<\/Count\>/x
  [orgn, $1.to_i]
end

results.sort_by{|orgn, count| count}.each do |orgn, count|
  puts [orgn, orgn_info[orgn], count].flatten.join("\t")
end


