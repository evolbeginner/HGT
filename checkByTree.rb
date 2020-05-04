#! /usr/bin/env ruby


##################################################
DIR = File.dirname($0)
$: << File.join(DIR, 'lib')


##################################################
require 'getoptlong'
require 'parallel'
require 'tmpdir'
require "sqlite3"
require 'bio'

require 'prot2taxonFromIndir'
require 'Dir'
require 'util'
require 'processbar'


##################################################
$FASTTREE = 'FastTree'
$TRIMAL=File.expand_path("/home-user/software/trimAl/trimal-1.4.1/source/trimal")


#REGEXP = Regexp.new('^[^\t]+ \t (Fungi\t){0,30} (NA\t){0,10} (Alphaproteobacteria\t){20,}', Regexp::EXTENDED)
#REGEXP = Regexp.new('^[^\t]+ \t (Streptophyta\t){0,6} (NA\t){0,3} (Alphaproteobacteria\t){20,}', Regexp::EXTENDED)
#REGEXP = Regexp.new('^[^\t]+ \t (Dinophyta\t){0,4} (NA\t){0,2} (Proteobacteria\t){30,}', Regexp::EXTENDED)
REGEXP = Regexp.new('^[^\t]+ \t (Cyanobacteria\t){0,50} (NA\t){0,5} (Planctomycetes\t){50,}', Regexp::EXTENDED)
#REGEXP = Regexp.new('^[^\t]+ \t (Fungi\t){0,30} (NA\t){0,10} (Planctomycetes\t){20,}', Regexp::EXTENDED)


##################################################
taxa_res_file = nil
db_files = Array.new
indirs = Array.new
blast_file = nil
cpu = 1
n = 500
is_iqtree = false
outdir = nil
is_force = false
is_tolerate = false


query2subjects = Hash.new{|h,k|h[k]=[]}
prot2taxon = Hash.new
seqObjs = Hash.new


##################################################
def parse_taxa_res_file(infile, cpu, regexp)
  genes = Array.new

  in_fh = File.open(infile, 'r')
  lines = in_fh.readlines
  in_fh.close

  results = Parallel.map(lines, in_processes: cpu) do |line|
    line.chomp!
    line_arr = line.split("\t")
    gene = line_arr[0]
    #next unless line =~ /^[^\t]+ \t (Streptophyta\t){3,8} (Proteobacteria\t){30,}/x
    next unless line =~ regexp 
    genes << gene
    genes
  end

  results.compact.each do |arr|
    genes << arr
  end

  if genes.flatten.size == 0
    puts "No genes found!"
    exit 0
  end

  genes.flatten!.uniq!

  return(genes)
end


def getProt2TaxonFromDb(db_file, query2subjects)
  seqObjs = Hash.new
  prot2taxon = Hash.new
  db = SQLite3::Database.new(db_file)
  #stm = db.prepare ("select * from seq where prot = \"#{subject}\"")
  stm = db.prepare ("select * from seq")
  rs = stm.execute

  while (row = rs.next) do
    #puts prot2taxon.values.uniq.size if count % 100000 == 0
    prot, taxon, seq = row
    f = Bio::Sequence::AA.new(seq)
    seqObjs[prot] = f
    prot2taxon[prot] = taxon
  end

  return([prot2taxon, seqObjs])
end


##################################################
opts = GetoptLong.new(
  ['-t', GetoptLong::REQUIRED_ARGUMENT],
  ['--db', GetoptLong::REQUIRED_ARGUMENT],
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['-b', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['-n', GetoptLong::REQUIRED_ARGUMENT],
  ['--iqtree', GetoptLong::NO_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--tolerate', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-t'
      taxa_res_file = value
    when '--db'
      db_files << value.split(',')
    when '--indir'
      indirs << value.split(',')
    when '-b'
      blast_file = value
    when '--cpu'
      cpu = value.to_i
    when '--outdir'
      outdir = value
    when '-n'
      n = value.to_i
    when '--iqtree'
      is_iqtree = true
    when '--force'
      is_force = true
    when '--tolerate'
      is_tolerate = true
  end
end


db_files.flatten!


##################################################
mkdir_with_force(outdir, is_force, is_tolerate)
param_outfile = File.join(outdir, 'param')
param_out_fh = File.open(param_outfile, 'w')
param_out_fh.puts REGEXP
param_out_fh.close


##################################################
genes = parse_taxa_res_file(taxa_res_file, cpu, REGEXP)
puts "#{genes.size} genes detected!"

Dir.mktmpdir do |tmpdir|
  outfile1 = File.join(tmpdir, 'genes.list')
  out_fh = File.open(outfile1, 'w')
  genes.each do |gene|
    out_fh.puts "^#{gene}"
  end
  out_fh.close  
  outfile2 = File.join(tmpdir, 'tmp.blast8')
  `grep -f #{outfile1} #{blast_file} > #{outfile2}`
  #`cp #{outfile1} ./`
  #`cp #{outfile2} ./`

  in_fh = File.open(outfile2, 'r')
  in_fh.each_line do |line|
    line.chomp!
    next if line =~ /^$/
    line_arr = line.split("\t")
    query, subject = line_arr[0, 2]
    query2subjects[query] << subject
  end
  in_fh.close

  query2subjects.each_pair do |query, subjects|
    subjects = subjects[0, n]
    subjects.uniq!
  end
end


##################################################
if not indirs.empty? or not db_files.empty?
  prot2taxon_1, seqObjs_1, prot2taxon_2, seqObjs_2 = Hash.new, Hash.new, Hash.new, Hash.new
  if not indirs.empty?
    indirs.flatten!
    prot2taxon_1, seqObjs_1 = getProt2Taxon(indirs, cpu, true)
  end
  if not db_files.empty?
    db_files.each do |db_file|
      a, b = getProt2TaxonFromDb(db_file, query2subjects)
      prot2taxon_2.merge!(a)
      seqObjs_2.merge!(b)
    end
  end
  prot2taxon = prot2taxon_1.merge(prot2taxon_2)
  seqObjs = seqObjs_1.merge(seqObjs_2)
else
  STDERR.puts "Wrong! Either indir or db_file has to be provided! Exiting ......"
  exit 1
end


##################################################
outdir_fas = File.join(outdir, 'fas')
outdir_aln = File.join(outdir, 'aln')
outdir_iqtree_aln = File.join(outdir, 'iqtree_aln')
outdir_trimal = File.join(outdir, 'trimal')
outdir_iqtree_trimal = File.join(outdir, 'iqtree_trimal')
outdir_fasttree = File.join(outdir, 'FastTree')
outdir_iqtree = File.join(outdir, 'iqtree')
mkdir_with_force(outdir_fas, is_force, is_tolerate)
mkdir_with_force(outdir_aln, is_force, is_tolerate)
mkdir_with_force(outdir_iqtree_aln, is_force, is_tolerate)
mkdir_with_force(outdir_trimal, is_force, is_tolerate)
mkdir_with_force(outdir_iqtree_trimal, is_force, is_tolerate)
mkdir_with_force(outdir_fasttree, is_force, is_tolerate)
mkdir_with_force(outdir_iqtree, is_force, is_tolerate)


# generate fasta sequences
query2subjects.each_pair do |query, subjects|
  outfile = File.join(outdir_fas, query.split('|')[0] + '.fas')
  out_fh = File.open(outfile, 'w')
  subjects.uniq.each do |subject|
    next unless seqObjs.include?(subject)
    seqObj = seqObjs[subject]
    taxon = prot2taxon[subject] 
    seqTitle = [taxon, subject].join('|')
    out_fh.puts '>' + seqTitle
    out_fh.puts seqObjs[subject].seq
  end
  out_fh.close
end


##################################################
# do alignment
puts "Alignment and tree construction is starting ......"
fas_files = read_infiles(outdir_fas)
Parallel.map(fas_files, in_processes: cpu) do |fas_file|
  c = getCorename(fas_file)
  puts c
  aln_file = File.join(outdir_aln, c+'.aln')

  `mafft --thread 8 --quiet #{fas_file} > #{aln_file}`

  iqtree_aln_file = File.join(outdir_iqtree_aln, c+'.aln')
  `cp #{aln_file} #{iqtree_aln_file}`
  `sed -i '/^>/s/[|]/!/g' #{iqtree_aln_file}`
  
  # do trimal
  trimal_file = File.join(outdir_trimal, c+'.aln')
  iqtree_trimal_file = File.join(outdir_iqtree_trimal, c+'.aln')
  system("#{$TRIMAL} -in #{aln_file} -out #{trimal_file} -st 0.001")
  system("#{$TRIMAL} -in #{iqtree_aln_file} -out #{iqtree_trimal_file} -st 0.001")

  fasttree_file = File.join(outdir_fasttree, c+'.FastTree.tre')
  #`#{$FASTTREE} -quiet < #{trimal_file} > #{fasttree_file}`
  `#{$FASTTREE} -quiet < #{aln_file} > #{fasttree_file}`

  if is_iqtree
    outdir_iqtree2 = File.join(outdir_iqtree, c)
    mkdir_with_force(outdir_iqtree2, true, true)
    iqtree_pre = File.join(outdir_iqtree2, c)
    #`iqtree -s #{iqtree_trimal_file} -m MFP -mset LG,JTT,WAG -mrate E,G,I,G+I -mfreq F,FU -bb 1000 -nt AUTO -redo -pre #{iqtree_pre}`
    `iqtree -s #{iqtree_aln_file} -m MFP -mset LG,JTT,WAG -mrate E,G,I,G+I -mfreq F,FU -bb 1000 -nt AUTO -redo -pre #{iqtree_pre}`
  end
end


