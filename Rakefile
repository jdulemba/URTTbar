$project_dir = ENV['URA_PROJECT']
$fwk_dir = ENV['URA']+'/AnalysisFW'

tools = "#{$fwk_dir}/rake/tools.rb"
require tools

analysis_rules = "#{$fwk_dir}/rake/analysis.rake"
import analysis_rules

meta_rules = "#{$fwk_dir}/rake/meta.rake"
import meta_rules

task :local_build => make_libs($project_dir, $fwk_dir)
task :build => make_libs($fwk_dir)

def psub(target, sub)
  return proc {|name| name.sub(target, sub)}
end

def psubs(target, *subs)
  return subs.map{|x| psub(target, x)}
end

rule /\.model\.root$/ => psub(/\.model\.root$/, '.txt') do |t|
  dir = File.dirname(t.name)
  chdir(dir) do
    sh "text2workspace.py #{File.basename(t.source)} -P URAnalysis.AnalysisTools.statistics.TTBarXsecFittter:ttxsecfitter -o #{File.basename(t.name)}"
  end
end

rule /\.mlfit\.root$/ => psub(/\.mlfit\.root$/, '.model.root') do |t|
  dir = File.dirname(t.name)
  chdir(dir) do
    sh "combine #{File.basename(t.source)} -M MaxLikelihoodFit --saveNormalizations --saveWithUncertainties --saveNLL"
    sh "cp mlfit.root #{File.basename(t.name)}"
  end
end

rule /\.harvested\.root$/ => psub(/\.harvested\.root$/, '.mlfit.root') do |t|
  sh "./harvest_fit.py #{t.source} #{t.source.sub(/\.mlfit\.root$/, '.binning.json')} -o #{t.name}"
end

