#! /usr/bin/env ruby


#############################################
require 'getoptlong'
require 'parallel'

require 'Hash'


#############################################
infile = nil
ref_infiles = Array.new
cpu = 1
window_size = 10
$is_alpha = false
$is_fungi = false

taxonomic_info = Hash.new


#############################################
class RANK
  attr_accessor :domain, :kingdom, :phylum, :class, :order, :family, :genus, :species
  def initialize(arr)
    @species = arr[0]
    @domain = arr[1]
    @kingdom = arr[2]
    @phylum = arr[3]
    @class = arr[4]
    @order = arr[5]
    @family = arr[6]
    @genus = arr[7]
  end
end


#############################################
def read_ref_infile(ref_infile, taxonomic_info)
  in_fh = File.open(ref_infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    orgn = line_arr[-1]
    rank = RANK.new(line_arr[0, line_arr.size-2])
    #p rank
    rank.phylum = rank.class if $is_alpha and rank.phylum == 'Proteobacteria'
    rank.phylum = rank.kingdom if $is_fungi and rank.kingdom == 'Fungi'
    taxonomic_info[orgn] = rank
  end
  return(taxonomic_info)
end


def getDominantRank(taxa_tmp, taxonomic_info)
  rv = nil
  #p taxa_tmp;exit

  counter = getCounts(taxa_tmp.map{|i|taxonomic_info[i].phylum})

  largestPair = counter.sort_by{|orgn, count| count}[-1]
  rank = largestPair[0]
  if largestPair[1].to_f/counter.values.sum > 0.5
    rv = rank
  else
    rv = 'NA'
  end
  return(rv)
end


#############################################
if $0 == __FILE__
  opts = GetoptLong.new(
    ['-i', GetoptLong::REQUIRED_ARGUMENT],
    ['--ref', GetoptLong::REQUIRED_ARGUMENT],
    ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
    ['--window', GetoptLong::REQUIRED_ARGUMENT],
    ['--alpha', GetoptLong::NO_ARGUMENT],
    ['--fungi', GetoptLong::NO_ARGUMENT],
  )

  opts.each do |opt, value|
    case opt
      when '-i'
        infile = value
      when '--ref'
        ref_infiles << value.split(',')
      when '--cpu'
        cpu = value.to_i
      when '--window'
        window_size = value.to_i
      when '--alpha'
        $is_alpha = true
      when '--fungi'
        $is_fungi = true
    end
  end


  #############################################
  ref_infiles.flatten!.each do |ref_infile|
    taxonomic_info = read_ref_infile(ref_infile, taxonomic_info)
  end


  #############################################
  in_fh = File.open(infile, 'r')
  lines = in_fh.readlines()
  results = Parallel.map(lines, in_processes: cpu) do |line|
    line.chomp!
    line_arr = line.split("\t")
    gene = line_arr[0]
    taxa = line_arr[1, line_arr.size-1]
    dominant_ranks = Array.new
    taxa.each_with_index do |taxon, index|
      taxa_tmp = taxa[index, window_size]
      next unless taxonomic_info.include?(taxon)
      next if taxa_tmp.size < window_size
      dominant_rank = getDominantRank(taxa_tmp, taxonomic_info)
      dominant_ranks << dominant_rank
    end
    [gene, dominant_ranks].flatten.join("\t")
  end
  in_fh.close


  results.each do |str|
    puts str
  end

end


