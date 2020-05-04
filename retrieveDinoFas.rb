#! /usr/bin/env ruby


####################################################################
DIR = File.dirname($0)
$: << File.join(DIR, 'lib')


####################################################################
require 'getoptlong'
require 'bio'

require 'process_species_name'
require 'Dir'
require 'processbar'


####################################################################
infile = nil
seqfile = nil
outdir = nil
is_force = false
is_tolerate = false


project_id_2_orgn = Hash.new
orgn_2_project_id = Hash.new{|k,v|k[v]=[]}
outdirs = Hash.new
out_fhs = Hash.new


####################################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--tolerate', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--seq'
      seqfile = value
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
    when '--tolerate'
      is_tolerate = true
  end
end


####################################################################
outdirs[:seq] = File.join(outdir, 'seq')
mkdir_with_force(outdirs[:seq], is_force, is_tolerate)


####################################################################
out_fh = File.open(File.join(outdir, 'ref.dino'), 'w')
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  next if $. == 1
  line.chomp!
  #MMETSP0328	Jakobsen	Kjetill	NA	Dinophyta	Dinophyceae	Goniodomataceae	Alexandrium	minutum	CCMP113	39455	1	0
  line_arr = line.split("\t")
  project_id = line_arr[0]
  rank_tmp = Hash.new
  [:phylum, :class, :family, :genus, :species, :strain].zip(line_arr[4,6]).each do |i, j|
    j = $1 if j =~ /(.+) \(=/
    rank_tmp[i] = j
  end
  rank_tmp[:order] = 'None'
  orgn0 = [rank_tmp[:genus], rank_tmp[:species], rank_tmp[:strain]].join(' ')
  orgn = Marshal.load(Marshal.dump(orgn0))
  orgn.process_species_name!
  project_id_2_orgn[project_id] = orgn
  orgn_2_project_id[orgn] << project_id
  #Mus musculus	Eukaryota	Metazoa	Chordata	Mammalia	Rodentia	Muridae	Mus	77425	Mus_musculus
  out_fh.puts [orgn0, 'Eukaryota', 'None', [:phylum, :class, :order, :family, :genus].map{|i|rank_tmp[i]}, '', orgn].flatten.join("\t")
end
in_fh.close
out_fh.close


####################################################################
# output project orgn rela
out_fh = File.open(File.join(outdir, 'project-orgn.rela'), 'w')
project_id_2_orgn.each_pair do |project_id, orgn|
  out_fh.puts [project_id, orgn].join("\t")
end
out_fh.close


####################################################################
in_fh = Bio::FlatFile.open(seqfile)
in_fh.each_entry do |f|
  #>CAMPEP_0113835384 /NCGR_PEP_ID=MMETSP0328-20130328|8914_1 /TAXON_ID=39455 /ORGANISM="Alexandrium minutum" /LENGTH=90 /DNA_ID=CAMNT_0000803721 /DNA_START=1 /DNA_END=269 /DNA_ORIENTATION=+ /assembly_acc=CAM_ASM_000350
  f.definition =~ /^(CAMPEP_\d+) \/NCGR_PEP_ID=([^-]+)/
  pep_id = $1
  project_id = $2
  next if not project_id_2_orgn.include?(project_id)

  orgn = project_id_2_orgn[project_id]
  if not out_fhs.include?(project_id)
    outfile = File.join(outdirs[:seq], project_id + '.fas')
    out_fh = File.open(outfile, 'w')
    out_fhs[project_id] = out_fh
  end
  out_fh = out_fhs[project_id]
  out_fh.puts '>' + pep_id
  out_fh.puts f.seq
  processbar(out_fhs.size, project_id_2_orgn.size)
end
in_fh.close


####################################################################
out_fhs.each_value do |out_fh|
  out_fh.close
end


