#! /usr/bin/env ruby


##################################################
DIR = File.dirname($0)
$: << File.join(DIR, 'lib')


##################################################
require 'getoptlong'
require 'parallel'
require "sqlite3"

require 'util'
require 'prot2taxonFromIndir'
require 'processbar'


##################################################
db_file = nil
indirs = Array.new
cpu = 1

query2subjects = Hash.new{|h,k|h[k]=[]}


##################################################
def create_db(seqObjs, prot2taxon, db_file)
  db = SQLite3::Database.new(db_file)
  rows = db.execute <<-SQL
  create table if not exists seq (
    prot varchar(100),
    orgn varchar(150),
    seq varchar(5000)
  );
  SQL

  count = 0
  seqObjs.each_pair do |prot, f|
    orgn = prot2taxon[prot]
    db.execute("INSERT OR IGNORE INTO seq (prot, orgn, seq) 
            VALUES (?, ?, ?)", [prot, orgn, f.seq])
    count += 1
    processbar(count, seqObjs.size)
  end
end


##################################################
opts = GetoptLong.new(
  ['--db', GetoptLong::REQUIRED_ARGUMENT],
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--db'
      db_file = value
    when '--indir'
      indirs << value.split(',')
    when '--cpu'
      cpu = value.to_i
  end
end


##################################################
indirs.flatten!


##################################################
prot2taxon, seqObjs = getProt2Taxon(indirs, cpu, true)

create_db(seqObjs, prot2taxon, db_file)


