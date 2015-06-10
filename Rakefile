$project_dir = ENV['URA_PROJECT']
$fwk_dir = ENV['URA']+'/AnalysisFW'
$jobid = ENV['jobid']

tools = "#{$fwk_dir}/rake/tools.rb"
require tools

analysis_rules = "#{$fwk_dir}/rake/analysis.rake"
import analysis_rules

meta_rules = "#{$fwk_dir}/rake/meta.rake"
import meta_rules

require 'json'

task :local_build => make_libs($project_dir, $fwk_dir)
task :build => make_libs($fwk_dir)

def psub(target, sub)
  return proc {|name| name.sub(target, sub)}
end

def psubs(target, *subs)
  return subs.map{|x| psub(target, x)}
end

task :new_ttbar_trial, [:info] do |t, args|
  new_trial('ttbarxsec', 'ttxsec', args.info)
end

task :new_ttbar_plots, [:info] do |t, args|
  new_trial('', 'ttxsec', args.info)
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
    sh "combine #{File.basename(t.source)} -M MaxLikelihoodFit --saveNormalizations --saveWithUncertainties --saveNLL --skipBOnlyFit #{toy_cmd}"
    sh "cp mlfit.root #{File.basename(t.name)}"
  end
end

rule /\.harvested\.root$/ => psub(/\.harvested\.root$/, '.mlfit.root') do |t|
  toy_cmd = ''
  if t.name.include? '.toy.'
    toy_cmd = '-t'
  end
  sh "./harvest_fit.py #{t.source} #{t.source.sub(/(:?\.toy)?\.mlfit\.root$/, '.binning.json')} #{t.source.sub(/(:?\.toy)?\.mlfit\.root$/, '.root')} -o #{t.name} #{toy_cmd}"
end

def get_toy_harvest(fname)
  vardir = fname.sub("plots/#{$jobid}/ttxsec/", "").split('/')[0]
  var = vardir.sub(/_\d+\.\d_\d+\.\d/, '')
  return "plots/#{$jobid}/ttxsec/#{vardir}/#{var}.toy.harvested.root"
end

rule /toys\/shapes\/tt_right.json$/ => proc {|name| get_toy_harvest(name)} do |t|
  vardir = File.dirname(t.source)
  var = File.basename(t.source).sub('.toy.harvested.root','')
  sh "./toy_diagnostics.py #{var} #{t.source} #{vardir}/#{var}.toy.mlfit.root -o #{vardir}/toys"
end

task :toy_diagnostics, [:vardir] do |t, args|
  Rake::Task["plots/#{$jobid}/ttxsec/#{args.vardir}/toys/shapes/tt_right.json"].invoke
end

task :fit, [:var] do |t, args|
  Rake::Task["plots/#{$jobid}/ttxsec/#{args.var}/#{args.var}.harvested.root"].invoke
end

task :fit_toys, [:var] do |t, args|
  Rake::Task["plots/#{$jobid}/ttxsec/#{args.var}/#{args.var}.toy.harvested.root"].invoke
end

task :optimize_binning, [:var, :varmin, :vmin, :vmax, :vstep, :terminate] do |t, args|
  terminate = args.terminate
  if not terminate
    terminate = 0.1
  end
  varmin = Float(args.varmin)

  (Float(args.vmin)..Float(args.vmax)).step(Float(args.vstep)) do |varmax|
    sh "python runTTXSecPlotter.py --optimize_binning '#{args.var}:[#{varmin},#{varmax}]'"
  end

  jsons = []
  (Float(args.vmin)..Float(args.vmax)).step(Float(args.vstep)) do |varmax|
    vardir = "#{args.var}_#{format("%.1f", varmin)}_#{format("%.1f", varmax)}"
    jfile = "plots/#{$jobid}/ttxsec/#{vardir}/toys/shapes/tt_right.json"
    Rake::Task[jfile].invoke
    jsons << jfile
  end
  sh "./plot_optimization.py #{args.var} plots/#{$jobid}/ttxsec #{jsons.join(' ')}"
end

###########################
# Charm-tagging tasks
###########################

task :new_ctag_trial, [:info] do |t, args|
  new_trial('btag_efficiency', 'btageff', args.info)
end

rule /fitModel.root$/ => psub(/fitModel.root$/, 'datacard.txt') do |t|
  dir = File.dirname(t.name)
  chdir(dir) do
    #sh "sed -i 's|$MASS||g' datacard.txt"
    #sh "sed -i 's|\x1b\[?1034h||g' datacard.txt"
    puts 'creating workspace'
    sh "text2workspace.py datacard.txt -P HiggsAnalysis.CombinedLimit.CTagEfficiencies:ctagEfficiency -o #{File.basename(t.name)}"
  end
end

rule /MultiDimFit(:?Toy)?.root$/ => psub(/MultiDimFit(:?Toy)?.root$/, 'fitModel.root') do |t|
  toy_cmd = ''
  if t.name.include? 'Toy'
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
  if t.name.include? 'Toy'
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
  dir = File.dirname(t.name)
  chdir(dir) do
    if t.name.include? 'Toy'
      toy_cmd = '--saveToys --expectSignal 1 -t 200'
      sh "cp #{ENV['URA']}/AnalysisTools/scripts/batch_job.sh ."
      File.open('condor.mltoys.jdl', 'w') do |file|
        file << "universe = vanilla\n"
        file << "Executable = batch_job.sh\n"
        file << "Should_Transfer_Files = YES\n"
        file << "WhenToTransferOutput = ON_EXIT\n"
        file << "Transfer_Input_Files = #{Dir.pwd}/fitModel.root\n"
        
        (2345678...2545678).step(10000) do |seed|
          file << "\n"
          file << "Output = con_#{seed}.stdout\n"
          file << "Error = con_#{seed}.stderr\n"
          file << "Arguments = combine fitModel.root -M MaxLikelihoodFit --saveNormalizations --saveWithUncertainties --saveNLL --saveToys --expectSignal 1 -t 10 -n #{seed} -s #{seed}\n"
          file << "Queue\n"
        end
      end
      sh 'condor_submit condor.mltoys.jdl'
    else
      puts 'running MaxLikelihood fit with Profile-Likelyhood errors'
      sh "combine fitModel.root -M MaxLikelihoodFit --saveNormalizations --saveWithUncertainties --minos=all --saveNLL #{toy_cmd}"
      sh "mv mlfit.root #{File.basename(t.name)}"      
    end
    #sh "mv higgsCombineTest.MultiDimFit.mH120.root MultiDimFit.root"
  end
end
