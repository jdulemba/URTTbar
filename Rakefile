$project_dir = ENV['URA_PROJECT']
$fwk_dir = ENV['URA']+'/AnalysisFW'
$jobid = ENV['jobid']

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

task :new_ttbar_trial, [:info] do |t, args|
  res_dir = "results/#{$jobid}/ttbarxsec"
  if File.exist? res_dir
    if not File.symlink? res_dir
      throw "#{res_dir} MUST be a symlink in the first place to work!"
    end
  sh "rm -f #{res_dir}"
  end

  plot_dir = "plots/#{$jobid}/ttxsec"
  if File.exist? plot_dir 
    if not File.symlink? plot_dir
      throw "#{plot_dir} MUST be a symlink in the first place to work!"
    end
    sh "rm -f #{plot_dir}"    
  end

  #add new dir and new link
  timestamp = Time.now.strftime("%Y%b%d")
  chdir("results/#{$jobid}") do
    new_res = "ttbarxsec_#{timestamp}_#{args.info}"
    sh "mkdir -p #{new_res}" 
    sh "ln -s #{new_res} ttbarxsec"
  end

  chdir("plots/#{$jobid}") do
    new_plot = "ttxsec_#{timestamp}_#{args.info}"
    sh "mkdir -p #{new_plot}" 
    sh "ln -s #{new_plot} ttxsec"
  end
end

rule /\.model\.root$/ => psub(/\.model\.root$/, '.txt') do |t|
  dir = File.dirname(t.name)
  chdir(dir) do
    sh "text2workspace.py #{File.basename(t.source)} -P URAnalysis.AnalysisTools.statistics.TTBarXsecFitter:ttxsecfitter -o #{File.basename(t.name)}"
  end
end

rule /(:?\.toy)?\.mlfit\.root$/ => psub(/(:?\.toy)?\.mlfit\.root$/, '.model.root') do |t|
  dir = File.dirname(t.name)
  toy_cmd = ''
  if t.name.include? '.toy.'
    toy_cmd = '--saveToys --expectSignal 1 -t 200'
  end
  chdir(dir) do
    sh "combine #{File.basename(t.source)} -M MaxLikelihoodFit --saveNormalizations --saveWithUncertainties --saveNLL #{toy_cmd}"
    sh "cp mlfit.root #{File.basename(t.name)}"
  end
end

rule /\.harvested\.root$/ => psub(/\.harvested\.root$/, '.mlfit.root') do |t|
  toy_cmd = ''
  if t.name.include? '.toy.'
    toy_cmd = '-t'
  end
  sh "./harvest_fit.py #{t.source} #{t.source.sub(/(:?\.toy)?\.mlfit\.root$/, '.binning.json')} -o #{t.name} #{toy_cmd}"
end

###########################
# Charm-tagging tasks
###########################

rule /fitModel.root$/ => psub(/fitModel.root$/, 'datacard.txt') do |t|
  dir = File.dirname(t.name)
  chdir(dir) do
    sh "sed -i 's|$MASS||g' datacard.txt"
    sh "sed -i 's|\x1b\[?1034h||g' datacard.txt"
    puts 'creating workspace'
    sh "text2workspace.py datacard.txt -P HiggsAnalysis.CombinedLimit.CTagEfficiencies:ctagEfficiency -o #{File.basename(t.name)} --PO inclusive=True"
  end
end

rule /MultiDimFit(:?Toy)?.root$/ => psub(/MultiDimFit(:?Toy)?.root$/, 'fitModel.root') do |t|
  toy_cmd = ''
  if t.name.include? '.toy.'
    toy_cmd = '--saveToys --expectSignal 1 -t 200'
  end
  dir = File.dirname(t.name)
  chdir(dir) do
    puts 'running Multi-dimensional with Profile-Likelyhood errors'
    sh "combine fitModel.root -M MultiDimFit --algo=singles --setPhysicsModelParameterRanges charmSF=0,2:lightSF=0,2 #{toy_cmd} > #{File.basename(t.name).sub(/\.root$/,'.log')}"
    sh "mv higgsCombineTest.MultiDimFit.mH120.root #{File.basename(t.name)}"
  end
end

rule /MultiDimScan(:?Toy)?.root$/ => psub(/MultiDimScan(:?Toy)?.root$/, 'fitModel.root') do |t|
  toy_cmd = ''
  if t.name.include? '.toy.'
    toy_cmd = '--saveToys --expectSignal 1 -t 200'
  end
  dir = File.dirname(t.name)
  chdir(dir) do
    puts 'running Multi-dimensional with Profile-Likelyhood errors'
    sh "combine -M MultiDimFit fitModel.root --algo=grid --points=400 --setPhysicsModelParameterRanges charmSF=0.05,2.05:lightSF=0.05,2.05 #{toy_cmd} > #{File.basename(t.name).sub(/\.root$/,'.log')}"
    sh "mv higgsCombineTest.MultiDimFit.mH120.root #{File.basename(t.name)}"
  end
end

rule /MaxLikeFit(:?Toy)?.root$/ => psub(/MaxLikeFit(:?Toy)?.root$/, 'fitModel.root') do |t|
  toy_cmd = ''
  if t.name.include? '.toy.'
    toy_cmd = '--saveToys --expectSignal 1 -t 200'
  end
  dir = File.dirname(t.name)
  chdir(dir) do
    puts 'running MaxLikelihood fit with Profile-Likelyhood errors'
    sh "combine fitModel.root -M MaxLikelihoodFit --saveNormalizations --saveWithUncertainties --saveNLL #{toy_cmd}"
    sh "mv mlfit.root #{File.basename(t.name)}"
    #sh "mv higgsCombineTest.MultiDimFit.mH120.root MultiDimFit.root"
  end
end
