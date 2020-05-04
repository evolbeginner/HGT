#! /usr/bin/env ruby


###################################################################
require 'parallel'
require 'bio'

require 'Dir'


###################################################################
def getProt2Taxon(indirs, cpu, is_seqObj=false)
  infiles = Array.new
  prot2taxon = Hash.new
  seqObjs = Hash.new

  indirs.each do |indir|
    infiles << read_infiles(indir)
  end
  infiles.flatten!
  infiles.select!{|i| File.ftype(i) != 'directory'}

  results = Parallel.map(infiles, in_processes: cpu) do |infile|
    prot2taxon = Hash.new
    b = File.basename(infile)
    b =~ /^[^.]+/
    taxon = $&
    in_fh = Bio::FlatFile.open(infile, 'r') 
    in_fh.each_entry do |f|
      f.definition =~ /^[^ ]+/
      prot = $&
      prot2taxon[prot] = taxon
      seqObjs[prot] = f if is_seqObj
    end
    in_fh.close
    [prot2taxon, seqObjs]
  end

  results.each do |h1, h2|
    prot2taxon.merge!(h1)
    seqObjs.merge!(h2)
  end

  return([prot2taxon, seqObjs])
end


