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

rule  /\.model\.root$/ => [
    proc {|name| name.sub(/\.model\.root$/, '.txt') }] do |t|
  dir = File.dirname(t.name)
  chdir(dir) do
    sh "text2workspace.py #{File.basename(t.source)} -P HiggsAnalysis.CombinedLimit.TTBarXsecFittter:ttxsecfitter -o #{File.basename(t.name)}"
  end
end

rule /\.mlfit\.root$/ => [ 
     proc {|name| name.sub(/\.mlfit\.root$/, '.model.root') }] do |t|
  dir = File.dirname(t.name)
  chdir(dir) do
    sh "combine #{File.basename(t.source)} -M MaxLikelihoodFit"
    sh "cp mlfit.root #{File.basename(t.name)}"
  end
end

