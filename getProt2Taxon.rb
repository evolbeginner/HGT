#! /usr/bin/env ruby


###########################################################################
DIR = File.dirname(__FILE__)
$: << File.join(DIR, 'lib')


###########################################################################
require 'getoptlong'
require 'parallel'
require 'bio'

require 'Dir'
require 'prot2taxonFromIndir'


###########################################################################
indirs = Array.new
blast_file = nil
evalue_cutoff = 1
cpu = 1
seq_outdir = nil
is_seqObjs = false
is_force = false
is_tolerate = false

infiles = Array.new


###########################################################################
class BLAST
  attr_accessor :subject, :evalue, :bit_score, :q_start, :q_end, :s_start, :s_end, :identity
  def initialize(subject, evalue, bit_score)
    @subject = subject
    @evalue = evalue
    @bit_score = bit_score
  end
end


###########################################################################
def output_result(blast_file, seq_outdir, prot2taxon, seqObjs, evalue_cutoff, cpu, is_getQuery2Blast=false)
  query2blast = Hash.new{|h,k|h[k]=[]}
  in_fh = File.open(blast_file, 'r')
  lines = in_fh.readlines
  in_fh.close

  results = Parallel.map(lines, in_processes: cpu) do |line|
    query2blast = Hash.new
    line.chomp!
    line_arr = line.split("\t")
    query, subject = line_arr[0, 2]
    identity = line_arr[2].to_i
    q_start, q_end, s_start, s_end = line_arr[6, 4].map{|i|i.to_f}
    evalue, bit_score = line_arr[-2, 2].map{|i|i.to_f}
    next if evalue > evalue_cutoff
    query2blast[query] = Array.new if not query2blast.include?(query)
    blast = BLAST.new(subject, evalue, bit_score)
    blast.q_start = q_start; blast.q_end = q_end; blast.s_start = s_start; blast.s_end = s_end; blast.identity = identity
    query2blast[query] << blast
    query2blast
  end

  results.each do |h|
    h.each_pair do |query, blasts|
      query2blast[query] << blasts
      query2blast[query].flatten!
    end
  end

  ###########################################################################
  if not seq_outdir.nil?
    mkdir_with_force(seq_outdir, is_force, is_tolerate)
    query2blast.each_pair do |prot, blasts|
      outputtedSeqs = Hash.new
      outfile = File.join(seq_outdir, prot+'.fas')
      out_fh = File.open(outfile, 'w')
      blasts.each do |blast|
        outputtedSeqs.include?(blast.subject) ? next : outputtedSeqs[blast.subject] = ''
        out_fh.puts '>' + [prot2taxon[blast.subject], blast.subject].join('|')
        out_fh.puts seqObjs[blast.subject].seq
      end
      out_fh.close
    end
  elsif not is_getQuery2Blast
    query2blast.each_pair do |prot, blasts|
      puts [prot, blasts.sort_by{|i|i.bit_score}.reverse.map{|i|prot2taxon[i.subject]}.compact].flatten.join("\t")
    end
  else
    return(query2blast)
  end
end


###########################################################################
if $0 == __FILE__
  opts = GetoptLong.new(
    ['--indir', GetoptLong::REQUIRED_ARGUMENT],
    ['-b', GetoptLong::REQUIRED_ARGUMENT],
    ['-e', GetoptLong::REQUIRED_ARGUMENT],
    ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
    ['--seq_outdir', GetoptLong::REQUIRED_ARGUMENT],
    ['--force', GetoptLong::NO_ARGUMENT],
    ['--tolerate', GetoptLong::NO_ARGUMENT],
  )


  opts.each do |opt, value|
    case opt
      when '--indir'
        indirs << value.split(',')
      when '-b'
        blast_file = value
      when '-e'
        evalue_cutoff = value.to_f
      when '--cpu'
        cpu = value.to_i
      when '--seq_outdir'
        seq_outdir = value
        is_seqObjs = true
      when '--force'
        is_force = true
      when '--tolerate'
        is_tolerate = true
    end
  end


  indirs.flatten!


  ###########################################################################
  prot2taxon, seqObjs = getProt2Taxon(indirs, cpu, is_seqObjs)

  output_result(blast_file, seq_outdir, prot2taxon, seqObjs, evalue_cutoff, cpu)

end

